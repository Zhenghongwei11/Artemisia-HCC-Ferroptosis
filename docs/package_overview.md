# CHB Ferroptosis Analysis - Review Package

## Overview
This is a reproducible analysis package for the manuscript "Network Pharmacology Reveals Artemisia capillaris Targets Ferroptosis-Related Genes in Hepatocellular Carcinoma".

## Directory Structure
```
fx_review_package/
├── run_complete_pipeline.R     # Main analysis pipeline script
├── README.md                  # This file
├── .gitignore                 # Git ignore rules
├── scripts/                   # R analysis scripts
├── data/                      # Data directories
│   └── references/            # Small reference files (ferroptosis genes, TCM targets)
├── docs/                      # Documentation
├── plots/                     # Output directory for figures
└── results/                   # Output directory for results
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

## Important Notes
- This package contains analysis scripts and reference data only
- Figures 2–6 are regenerated from curated result tables included under `results/`.
- Raw datasets are not redistributed here; some upstream steps require manual downloads and file placement (see `docs/DATA_MANIFEST.md`).
- Final figures are saved to `plots/publication/`.

## Contents
- All analysis scripts (WGCNA, ML, survival analysis, immune infiltration, network pharmacology, molecular docking)
- Reference gene lists (ferroptosis genes, TCM targets)
- Pipeline script to execute the complete analysis
- Documentation for installation and usage

This package regenerates Figures 2–6 from curated result tables and documents the manual inputs needed for an optional full re-run from public sources.
