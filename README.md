# Wikidata-Labels-to-IUPAC Conversion Dataset

## Overview

This dataset contains the results of converting labels from Wikidata to molecular structures using OPSIN (Open Parser for Systematic IUPAC Nomenclature). The conversion process validates results by comparing generated InChIKeys with values from Wikidata.

## Dataset Statistics

- **Total input compounds**: 1,341,287
- **Successfully processed**: 1,341,287
- **Successful matches**: 838,452
- **Success rate**: 62.51%
- **Processing date**: 2025-06-10

## Files Description

### Primary Data Files

- **`chemical_matches.csv`**: Successful matches in CSV format
  - Contains: chemical names, InChIKeys, SMILES structures, processing metadata
  - Suitable for: Spreadsheet analysis, statistical processing

- **`chemical_matches.json`**: Successful matches in JSON format
  - Contains: Same data as CSV with additional metadata and schema
  - Suitable for: Programmatic access, web applications

- **`chemical_matches.sdf`** (optional): Molecular structures in SDF format
  - Contains: 3D molecular structures with properties
  - Suitable for: Chemical visualization, molecular modeling

### Reference and Metadata Files

- **`wikidata_reference.csv`**: Original Wikidata compounds
  - Contains: All chemical names and InChIKeys from Wikidata query
  - Purpose: Complete reference dataset for reproducibility

- **`processing_statistics.json`**: Detailed processing statistics
  - Contains: Success rates, error counts, quality metrics
  - Purpose: Quality assessment and method validation

- **`metadata.json`**: Dataset metadata
  - Contains: Complete dataset description, methodology, provenance
  - Purpose: Zenodo compliance and data citation

## Methodology

### Data Source
Chemical compound data was retrieved from Wikidata using SPARQL queries targeting:
- Compounds with English language labels (`rdfs:label`)
- Compounds with InChIKey identifiers (`wdt:P235`)

### Processing Pipeline
1. **Data Retrieval**: SPARQL query to Wikidata endpoint
2. **Name Processing**: Batch conversion using OPSIN v2.8.0
3. **Structure Validation**: SMILES validation using RDKit
4. **Quality Control**: InChIKey comparison for accuracy verification
5. **Result Export**: Multi-format output generation

### Quality Assurance
- InChIKey format validation
- SMILES structure validation using RDKit
- Exact InChIKey matching for success determination
- Comprehensive error logging and statistics

## Software Requirements

### Runtime Dependencies
- Python 3.10+
- Java Runtime Environment (JRE) 8+
- OPSIN v2.8.0 (automatically downloaded)

### Python Packages
- RDKit (chemistry toolkit)
- SPARQLWrapper (SPARQL queries)
- requests (HTTP requests)
- tqdm (progress bars)

## Limitations and Considerations

1. **Coverage**: Results represent only compounds with both Wikidata names and InChIKeys
2. **Accuracy**: Success determined by exact InChIKey matching
3. **Scope**: Limited to systematic chemical nomenclature parseable by OPSIN
4. **Language**: Only English language chemical names processed

## Use

```bash
docker build -t wd-labels-to-iupac .

# Run (assuming script outputs to current directory)
docker run -v $(pwd):/app/output wd-labels-to-iupac
```

## Reproducibility

To reproduce this dataset:
1. Install required dependencies
2. Run the conversion script with default settings
3. Compare results with provided reference data

## Citation

If you use this dataset in your research, please cite:

```
Adriano Rutz (2025). Wikidata-Labels-to-IUPAC Conversion Dataset. 
Version 0.0.1. [Dataset]. Zenodo. https://doi.org/[DOI]
```

## License

This dataset is released under the MIT (code) and CC0 (data) License. The original chemical data is from Wikidata (CC0 1.0 Universal).

## Contact

For questions or issues regarding this dataset:
- Email: adafede@gmail.com
- GitHub: https://github.com/Adafede/wd-labels-to-iupac

## Acknowledgments

- Egon Willighagen (0000-0001-7542-0286) for the original idea (see https://doi.org/10.59350/dycsw-qeq51)
- Wikidata contributors for chemical compound data
- OPSIN developers for the nomenclature parsing tool
- RDKit developers for chemical informatics capabilities
