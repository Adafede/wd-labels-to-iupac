"""
Wikidata-Labels-to-IUPAC-Converter (WLIC)

Converts chemical labels from Wikidata to IUPAC names using OPSIN.
Validates results by comparing InChIKeys.

Author: Adriano Rutz
License: MIT (code) and CC0 (data)
Version: 0.0.2
"""

import argparse
import csv
import json
import logging
import os
import requests
import subprocess
import sys
import tarfile
import time

from dataclasses import dataclass, asdict
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import inchi
from SPARQLWrapper import SPARQLWrapper, JSON
from tqdm import tqdm
from typing import Dict, List, Tuple

__version__ = "0.0.2"
__author__ = "Adriano Rutz"
__email__ = "adafede@gmail.com"
__license__ = "MIT (code) and CC0 (data)"


@dataclass
class Config:
    opsin_url: str = "https://github.com/dan2097/opsin/releases/download/2.8.0/opsin-cli-2.8.0-jar-with-dependencies.jar"
    opsin_jar: str = "opsin-cli-2.8.0-jar-with-dependencies.jar"
    sparql_endpoint: str = "https://qlever.dev/api/wikidata"
    max_workers: int = 4
    timeout: int = 300
    retry_attempts: int = 3
    output_dir: str = "output"
    export_sdf: bool = False


@dataclass
class ChemicalMatch:
    label: str
    expected_inchikey: str
    actual_inchikey: str
    smiles: str
    confidence: float = 1.0
    source: str = "wikidata"
    processing_time: float = 0.0


class WLICError(Exception):
    pass


class Logger:
    def __init__(self, name: str = "WLIC", level: int = logging.INFO):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        if not self.logger.handlers:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(level)

            file_handler = logging.FileHandler("WLIC_converter.log")
            file_handler.setLevel(logging.DEBUG)

            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            console_handler.setFormatter(formatter)
            file_handler.setFormatter(formatter)

            self.logger.addHandler(console_handler)
            self.logger.addHandler(file_handler)

    def info(self, msg: str) -> None:
        self.logger.info(msg)

    def error(self, msg: str) -> None:
        self.logger.error(msg)

    def warning(self, msg: str) -> None:
        self.logger.warning(msg)

    def debug(self, msg: str) -> None:
        self.logger.debug(msg)


class DependencyManager:
    def __init__(self, config: Config, logger: Logger):
        self.config = config
        self.logger = logger

    def download_opsin(self) -> bool:
        if os.path.exists(self.config.opsin_jar):
            self.logger.info(f"OPSIN JAR found: {self.config.opsin_jar}")
            return True

        try:
            self.logger.info(f"Downloading OPSIN from {self.config.opsin_url}")
            response = requests.get(
                self.config.opsin_url,
                stream=True,
                timeout=self.config.timeout,
                headers={"User-Agent": "WLIC-Converter/0.0.2"},
            )
            response.raise_for_status()

            total_size = int(response.headers.get("content-length", 0))
            with open(self.config.opsin_jar, "wb") as f:
                with tqdm(
                    total=total_size, unit="B", unit_scale=True, desc="Downloading"
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

            self.logger.info("Download complete")
            return True

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Download failed: {e}")
            return False

    def verify_java(self) -> bool:
        try:
            result = subprocess.run(
                ["java", "-version"], capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                version = result.stderr.split("\n")[0] if result.stderr else "Unknown"
                self.logger.info(f"Java verified: {version}")
                return True
            return False
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.logger.error(f"Java verification failed: {e}")
            return False


class WikidataQuerier:
    def __init__(self, config: Config, logger: Logger):
        self.config = config
        self.logger = logger
        self.sparql = SPARQLWrapper(config.sparql_endpoint)
        self.sparql.setReturnFormat(JSON)
        self.sparql.addCustomHttpHeader("User-Agent", "WLIC-Converter/0.0.2")

    def query_chemical_data(self) -> Dict[str, str] | None:
        query = """
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX wdt: <http://www.wikidata.org/prop/direct/>
        
        SELECT DISTINCT ?item ?inchikey ?label WHERE { 
            ?item wdt:P235 ?inchikey .
            ?item rdfs:label ?label .
            FILTER(lang(?label) = "en")
        }
        ORDER BY ?label
        """

        for attempt in range(self.config.retry_attempts):
            try:
                self.logger.info(f"Querying Wikidata (attempt {attempt + 1})")
                start = time.time()
                self.sparql.setQuery(query)
                results = self.sparql.query().convert()
                self.logger.info(f"Query completed in {time.time() - start:.2f}s")
                return self._process_results(results)

            except Exception as e:
                self.logger.warning(f"Query failed: {e}")
                if attempt < self.config.retry_attempts - 1:
                    time.sleep(2**attempt)
                else:
                    raise WLICError(
                        f"Query failed after {self.config.retry_attempts} attempts"
                    )
        return None

    def _process_results(self, results: dict) -> Dict[str, str]:
        label_dict = {}
        duplicates = set()
        invalid = 0

        for result in results["results"]["bindings"]:
            try:
                label = result["label"]["value"].lower().strip()
                inchikey = result["inchikey"]["value"].strip()

                if not self._validate_inchikey(inchikey):
                    invalid += 1
                    continue

                if label in label_dict and label_dict[label] != inchikey:
                    duplicates.add(label)
                    continue

                label_dict[label] = inchikey

            except KeyError:
                continue

        if duplicates:
            self.logger.warning(f"{len(duplicates)} duplicate labels (kept first)")
        if invalid:
            self.logger.warning(f"{invalid} invalid InChIKeys filtered")

        self.logger.info(f"Processed {len(label_dict)} compounds")
        return label_dict

    @staticmethod
    def _validate_inchikey(key: str) -> bool:
        if not key or len(key) != 27:
            return False
        parts = key.split("-")
        return (
            len(parts) == 3
            and len(parts[0]) == 14
            and len(parts[1]) == 10
            and len(parts[2]) == 1
            and all(c.isalnum() for p in parts for c in p)
        )


class OpsinProcessor:
    def __init__(self, config: Config, logger: Logger):
        self.config = config
        self.logger = logger

    def process_labels(self, labels: List[str]) -> List[str]:
        if not labels:
            return []

        input_file = Path("opsin_input.txt")
        output_file = Path("opsin_output.txt")

        try:
            with open(input_file, "w", encoding="utf-8") as f:
                f.write("\n".join(labels))

            self.logger.info(f"Processing {len(labels)} labels through OPSIN")
            start = time.time()

            subprocess.run(
                [
                    "java",
                    "-jar",
                    self.config.opsin_jar,
                    "-osmi",
                    str(input_file),
                    str(output_file),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            self.logger.info(f"Processing completed in {time.time() - start:.2f}s")

            if not output_file.exists():
                return [""] * len(labels)

            with open(output_file, "r", encoding="utf-8") as f:
                smiles = [line.strip() for line in f]

            while len(smiles) < len(labels):
                smiles.append("")

            success = sum(1 for s in smiles if s)
            self.logger.info(f"Success rate: {(success/len(labels))*100:.1f}%")
            return smiles

        except subprocess.CalledProcessError as e:
            self.logger.error(f"OPSIN failed: {e}")
            return [""] * len(labels)
        finally:
            for f in [input_file, output_file]:
                if f.exists():
                    f.unlink()


class MoleculeValidator:
    def __init__(self, logger: Logger):
        self.logger = logger

    def validate_and_compare(
        self, labels: List[str], expected: List[str], smiles: List[str]
    ) -> Tuple[List[ChemicalMatch], Dict[str, int]]:
        matches = []
        stats = {
            "total_processed": 0,
            "valid_names": 0,
            "inchikey_matches": 0,
            "invalid_names": 0,
            "conversion_errors": 0,
        }

        self.logger.info("Validating structures")

        for label, exp_key, smi in tqdm(
            zip(labels, expected, smiles), total=len(labels), desc="Validating"
        ):
            start = time.time()
            stats["total_processed"] += 1

            if not smi:
                stats["invalid_names"] += 1
                continue

            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    stats["invalid_names"] += 1
                    continue

                stats["valid_names"] += 1
                actual_key = inchi.MolToInchiKey(mol)

                if actual_key == exp_key:
                    stats["inchikey_matches"] += 1
                    matches.append(
                        ChemicalMatch(
                            label=label,
                            expected_inchikey=exp_key,
                            actual_inchikey=actual_key,
                            smiles=smi,
                            processing_time=time.time() - start,
                        )
                    )
            except Exception as e:
                stats["conversion_errors"] += 1
                self.logger.debug(f"Error validating {label}: {e}")

        if stats["total_processed"] > 0:
            rate = (stats["inchikey_matches"] / stats["total_processed"]) * 100
            self.logger.info(f"Success rate: {rate:.2f}%")

        return matches, stats


class ResultExporter:
    def __init__(self, config: Config, logger: Logger):
        self.config = config
        self.logger = logger
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(exist_ok=True)

    def export_results(
        self,
        matches: List[ChemicalMatch],
        stats: Dict[str, int],
        original: Dict[str, str],
        export_sdf: bool,
    ) -> None:
        self.logger.info("Exporting results")

        self._export_csv(matches, "chemical_matches.csv")
        self._export_json(matches, "chemical_matches.json")
        if export_sdf:
            self._export_sdf(matches, "chemical_matches.sdf")

        self._export_statistics(stats, len(original))
        self._export_metadata()
        self._export_original(original)
        self._create_readme(stats, original)
        self._compress_large_files()

        self.logger.info(f"Exported to {self.output_dir}")

    def _export_csv(self, matches: List[ChemicalMatch], filename: str) -> None:
        with open(self.output_dir / filename, "w", newline="", encoding="utf-8") as f:
            if matches:
                writer = csv.DictWriter(f, fieldnames=list(asdict(matches[0]).keys()))
                writer.writeheader()
                for m in matches:
                    writer.writerow(asdict(m))
            else:
                csv.writer(f).writerow(
                    [
                        "label",
                        "expected_inchikey",
                        "actual_inchikey",
                        "smiles",
                        "confidence",
                        "source",
                        "processing_time",
                    ]
                )

    def _export_json(self, matches: List[ChemicalMatch], filename: str) -> None:
        data = {
            "metadata": {
                "version": __version__,
                "total_matches": len(matches),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            },
            "data": [asdict(m) for m in matches],
        }
        with open(self.output_dir / filename, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

    def _export_sdf(self, matches: List[ChemicalMatch], filename: str) -> None:
        try:
            with open(self.output_dir / filename, "w", encoding="utf-8") as f:
                for m in matches:
                    try:
                        mol = Chem.MolFromSmiles(m.smiles)
                        if mol:
                            mol.SetProp("_Name", m.label)
                            mol.SetProp("InChIKey", m.actual_inchikey)
                            f.write(Chem.MolToMolBlock(mol))
                            f.write("$$$$\n")
                    except Exception as e:
                        self.logger.debug(f"SDF write failed for {m.label}: {e}")
        except Exception as e:
            self.logger.warning(f"SDF export failed: {e}")

    def _export_statistics(self, stats: Dict[str, int], total: int) -> None:
        rate = (
            (stats["inchikey_matches"] / stats["total_processed"]) * 100
            if stats["total_processed"] > 0
            else 0
        )
        data = {
            "summary": {
                "total_input": total,
                "processed": stats["total_processed"],
                "matches": stats["inchikey_matches"],
                "success_rate": round(rate, 2),
            },
            "details": stats,
        }
        with open(self.output_dir / "statistics.json", "w") as f:
            json.dump(data, f, indent=2)

    def _export_metadata(self) -> None:
        metadata = {
            "title": "Wikidata-Labels-to-IUPAC Conversion Dataset",
            "description": "Conversion of Wikidata labels to IUPAC names using OPSIN with InChIKey validation",
            "version": __version__,
            "creators": [
                {
                    "name": __author__,
                    "orcid": "0000-0003-0443-9902",
                    "affiliation": {
                        "name": "ETH Zurich",
                        "department": "Institute for Molecular Systems Biology",
                        "address": "Otto-Stern-Weg, 3",
                        "city": "Zurich",
                        "country": "Switzerland",
                        "postal_code": "8093",
                        "ror": "https://ror.org/03j5gm982",
                    },
                }
            ],
            "license": __license__,
            "keywords": ["cheminformatics", "IUPAC", "OPSIN", "wikidata"],
            "related_identifiers": [
                {
                    "identifier": "https://github.com/dan2097/opsin",
                    "relation": "requires",
                    "resource_type": "software",
                    "scheme": "url",
                },
                {
                    "identifier": "https://github.com/rdkit/rdkit",
                    "relation": "requires",
                    "resource_type": "software",
                    "scheme": "url",
                },
                {
                    "identifier": "https://www.wikidata.org",
                    "relation": "isSourceOf",
                    "resource_type": "dataset",
                    "scheme": "url",
                },
                {
                    "identifier": "10.59350/dycsw-qeq51",
                    "relation": "cites",
                    "resource_type": "other",
                    "scheme": "doi",
                },
            ],
        }
        with open(self.output_dir / ".zenodo.json", "w") as f:
            json.dump(metadata, f, indent=2)

    def _export_original(self, data: Dict[str, str]) -> None:
        with open(
            self.output_dir / "wikidata_reference.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(["label", "inchikey"])
            for label, key in data.items():
                writer.writerow([label, key])

    def _create_readme(
        self,
        stats: Dict[str, int],
        original: Dict[str, str],
    ) -> None:
        rate = (
            (stats["inchikey_matches"] / stats["total_processed"]) * 100
            if stats["total_processed"] > 0
            else 0
        )
        content = f"""# Wikidata-Labels-to-IUPAC Dataset

Conversion of Wikidata chemical labels to structures using OPSIN with InChIKey validation.

## Statistics

- Total input: {len(original):,}
- Processed: {stats['total_processed']:,}
- Matches: {stats['inchikey_matches']:,}
- Success rate: {rate:.2f}%
- Date: {time.strftime('%Y-%m-%d', time.gmtime())}

## Files

- `chemical_matches.csv/json`: Successful conversions
- `wikidata_reference.csv`: Original Wikidata data
- `statistics.json`: Processing statistics
- `.zenodo.json`: Dataset metadata

## Usage

```bash
docker build -t wd-labels-to-iupac .
docker run -v $(pwd):/app/output wd-labels-to-iupac
```

## Citation

Adriano Rutz. ({time.strftime('%Y', time.gmtime())}). Wikidata-Labels-to-IUPAC Conversion Dataset ({__version__}). Zenodo. <https://doi.org/10.5281/zenodo.17378415>

## Acknowledgments

- Egon Willighagen ([0000-0001-7542-0286](https://orcid.org/0000-0001-7542-0286)) for the original idea
- Wikidata contributors
- OPSIN and RDKit developers
"""
        with open(self.output_dir / "README.md", "w") as f:
            f.write(content)

    def _compress_large_files(self, min_size_mb: float = 50.0) -> None:
        """Compress large CSV/JSON files individually as .tar.gz archives."""
        for file in self.output_dir.glob("*.csv"):
            self._maybe_tar(file, min_size_mb)
        for file in self.output_dir.glob("*.json"):
            self._maybe_tar(file, min_size_mb)

    def _maybe_tar(self, file_path: Path, min_size_mb: float) -> None:
        size_mb = file_path.stat().st_size / (1024 * 1024)
        if size_mb >= min_size_mb:
            tar_path = file_path.with_suffix(file_path.suffix + ".tar.gz")
            self.logger.info(f"Compressing {file_path.name} ({size_mb:.1f} MB)")
            with tarfile.open(tar_path, "w:gz") as tar:
                tar.add(file_path, arcname=file_path.name)
            file_path.unlink()
            self.logger.info(f"→ {tar_path.name}")


class WLIC:
    def __init__(self, config: Config = None):
        self.config = config or Config()
        self.logger = Logger()
        self.deps = DependencyManager(self.config, self.logger)
        self.querier = WikidataQuerier(self.config, self.logger)
        self.processor = OpsinProcessor(self.config, self.logger)
        self.validator = MoleculeValidator(self.logger)
        self.exporter = ResultExporter(self.config, self.logger)

    def run(self) -> bool:
        try:
            self.logger.info(f"WLIC v{__version__}")

            if not self.deps.verify_java() or not self.deps.download_opsin():
                return False

            data = self.querier.query_chemical_data()
            if not data:
                return False

            labels = list(data.keys())
            expected = list(data.values())
            smiles = self.processor.process_labels(labels)
            matches, stats = self.validator.validate_and_compare(
                labels, expected, smiles
            )

            self.exporter.export_results(matches, stats, data, self.config.export_sdf)

            self.logger.info(
                f"Complete. {len(matches)} matches exported to {self.config.output_dir}"
            )
            return True

        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            return False


def main():
    parser = argparse.ArgumentParser(description="Wikidata-Labels-to-IUPAC-Converter")
    parser.add_argument("--max-workers", type=int, default=4, help="Parallel threads")
    parser.add_argument(
        "--output-dir", type=str, default="output", help="Output directory"
    )
    parser.add_argument("--timeout", type=int, default=300, help="Timeout (seconds)")
    parser.add_argument("--retry-attempts", type=int, default=3, help="Retry attempts")
    parser.add_argument("--export-sdf", action="store_true", help="Export SDF format")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument("--quiet", "-q", action="store_true", help="Quiet mode")
    parser.add_argument("--version", action="version", version=f"WLIC {__version__}")

    args = parser.parse_args()

    if args.verbose and args.quiet:
        parser.error("--verbose and --quiet are mutually exclusive")

    config = Config(
        max_workers=args.max_workers,
        output_dir=args.output_dir,
        timeout=args.timeout,
        retry_attempts=args.retry_attempts,
        export_sdf=args.export_sdf,
    )

    level = (
        logging.DEBUG
        if args.verbose
        else (logging.ERROR if args.quiet else logging.INFO)
    )
    logging.getLogger().setLevel(level)

    try:
        converter = WLIC(config)
        sys.exit(0 if converter.run() else 1)
    except KeyboardInterrupt:
        print("\nCancelled")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
