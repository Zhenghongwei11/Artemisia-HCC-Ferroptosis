# Installation and Usage Guide

## Prerequisites

- R (version >= 4.0)
- Internet access for downloading packages and data

## Installation

1. Clone or download this repository to your local machine
2. Navigate to the project directory:
   ```bash
   cd fx_review_package
   ```

## Environment Setup

The analysis pipeline will automatically install required R packages. You can also install them manually by running:
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

1. **Environment Setup**: Installs required R packages
2. **Data Download**: Downloads GSE14520 and GSE83148 data from GEO database
3. **Batch Correction**: Combines and corrects for batch effects
4. **Differential Expression**: Identifies DEGs between conditions
5. **Ferroptosis Analysis**: Identifies ferroptosis-related genes
6. **Prognostic Model**: Builds a 15-gene ferroptosis-based prognostic model
7. **Immune Analysis**: Analyzes immune infiltration patterns
8. **Network Pharmacology**: Identifies intersection with TCM targets
9. **Molecular Docking**: Performs computational docking validation
10. **Figure Generation**: Creates publication-ready figures

## Expected Runtime

The complete analysis can take several hours to complete depending on your system specifications and internet connection speed, as it involves:
- Downloading large datasets from GEO
- Running computationally intensive analyses (WGCNA, ML algorithms)
- Performing molecular docking calculations

## Output Files

- Processed data: `data/processed/`
- Analysis results: `results/`
- Figures: `plots/publication/`

## Troubleshooting

- If R package installation fails, ensure you have internet access and try installing packages individually
- Large data downloads may take time; be patient during the initial steps
- If memory errors occur, ensure your system has sufficient RAM (16GB+ recommended)