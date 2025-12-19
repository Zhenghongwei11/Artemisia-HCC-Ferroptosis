# Installation and Usage Guide

## Prerequisites

- R (version >= 4.0)
- Internet access is recommended for installing packages (data downloads are not fully automated in this package)

## Installation

1. Clone or download this repository to your local machine
2. Navigate to the project directory:
   ```bash
   cd fx_review_package
   ```

## Environment Setup

Install required R packages by running:
```R
# Start R/RStudio
R
# Then run:
source("scripts/00_setup_env.R")
```

## Running the Complete Analysis

Execute the main pipeline script:
```bash
Rscript run_complete_pipeline.R
```

## What the Pipeline Does

This entry script regenerates Figures 2–6 from the curated result tables in `results/`.

If you want to re-run upstream analyses from public raw data, follow:
- `docs/DATA_MANIFEST.md`
- `docs/REPRODUCIBILITY.md`

## Expected Runtime

Regenerating Figures 2–6 from included results typically completes in minutes.

A full re-run from public raw data may take several hours and requires manual downloads (and, for some steps, external services).

## Output Files

- Results (inputs to figure regeneration): `results/`
- Figures: `plots/publication/`
- Optional raw/processed inputs (if re-running upstream): `data/raw/`, `data/processed/`

## Troubleshooting

- If package installation fails, ensure you have internet access and try installing packages individually.
- If you attempt a full re-run and a required file is missing, check `docs/DATA_MANIFEST.md` for the exact filename and placement path.
