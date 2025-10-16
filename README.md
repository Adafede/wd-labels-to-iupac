# Wikidata-Labels-to-IUPAC Dataset

Conversion of Wikidata chemical labels to structures using OPSIN with InChIKey validation.

## Statistics

- Total input: 1,341,787
- Processed: 1,341,787
- Matches: 838,633
- Success rate: 62.50%
- Date: 2025-10-16

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

Adriano Rutz (2025). Wikidata-Labels-to-IUPAC Dataset. Version 0.0.1.

## Acknowledgments

- Egon Willighagen ([0000-0001-7542-0286](https://orcid.org/0000-0001-7542-0286)) for the original idea
- Wikidata contributors
- OPSIN and RDKit developers
