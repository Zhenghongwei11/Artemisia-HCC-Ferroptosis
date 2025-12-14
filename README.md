# CHB Ferroptosis Analysis - Review Package

## Overview
This is a reproducible analysis package for the manuscript "Network Pharmacology Reveals Artemisia capillaris Targets Ferroptosis-Related Genes in Hepatocellular Carcinoma".

## Directory Structure
```
fx_review_package/
├── run_complete_pipeline.R     # Main analysis pipeline script
├── README.md                  # This file
├── scripts/                   # R analysis scripts
├── data/                      # Reference data
│   └── references/            # Small reference files (ferroptosis genes, TCM targets)
└── plots/                     # Output directory for figures
```

## Data Sources
The analysis uses data from public repositories:
- GSE14520: GEO database (https://www.ncbi.nlm.nih.gov/geo/)
- GSE83148: GEO database (https://www.ncbi.nlm.nih.gov/geo/)
- Ferroptosis genes: FerrDb database (http://www.zhounan.org/ferrdb/)
- TCM targets: TCMSP database (https://tcmsp-e.com/)

## Running the Analysis

1. Make sure you have R (>=4.0) installed with internet access for downloading packages
2. Navigate to the project directory
3. Run the complete analysis pipeline:
   ```bash
   Rscript run_complete_pipeline.R
   ```

## Note
- This package contains analysis scripts and reference data only
- Raw data is downloaded from public databases during execution
- Large intermediate and result files are generated during execution
- Final figures are saved to the `plots/publication/` directory