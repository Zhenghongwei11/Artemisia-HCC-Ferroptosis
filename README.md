# Artemisia capillaris Targets Ferroptosis in Hepatocellular Carcinoma

**Network Pharmacology and Molecular Mechanism Study**

> 📄 **Reproducible analysis package for peer review** - Contains all scripts, data, and figures for independent verification of the manuscript analysis.

---

## 📋 Overview

A systematic network pharmacology investigation revealing how Artemisia capillaris (绵茵陈) regulates ferroptosis-related genes in hepatocellular carcinoma (HCC), with integrated multi-omics analysis and molecular validation.

**Key Findings:**
- 3 Hub genes identified: **ACSL4, TFRC, NQO1**
- 15-gene prognostic risk model (C-index: 0.717)
- Validated immune infiltration and checkpoint mechanisms
- Dual molecular docking validation (CB-Dock2 + NVIDIA Boltz-2)

---

## 📊 Study Workflow & Key Figures

**For detailed figure visualization**, please download the repository and view:
- High-resolution figures in `plots/publication/` (PNG, PDF, TIFF formats)
- Figure legends and specifications in `plots/publication/README_figures.md`

**Key Results Summary:**
- **Figure 1:** Comprehensive research workflow (5 analytical modules)
- **Figure 2:** DEG analysis (Volcano plots + Venn intersection) → 3 Hub genes
- **Figure 3:** Prognostic model (KM curves, ROC, Risk scores)
- **Figure 4:** Clinical utility (Cox regression, DCA, C-index)
- **Figure 5:** Immune microenvironment (MDSC r=0.48, checkpoint analysis)
- **Figure 6:** Drug sensitivity & molecular docking (18 protein-ligand pairs)

---

## 📁 Directory Structure
```
fx_review_package/
├── run_complete_pipeline.R          # Main analysis pipeline
├── README.md                        # This file
├── FIGURE_GENERATION_METADATA.md    # Figure 1 generation prompt
├── scripts/                         # R analysis scripts (15 modules)
│   ├── 00_setup_env.R
│   ├── 01b_download_GSE14520.R
│   ├── 02_batch_correction.R
│   ├── 02_bulk_WGCNA_ML_CHB.R
│   ├── 02b_multi_cohort_DEG.R
│   ├── 02c_ml_feature_selection_*.R (3 variants)
│   ├── 02c_prognostic_model_v2.R
│   ├── 02d_nomogram_calibration_dca.R
│   ├── 02e_external_validation.R
│   ├── 03a_immune_infiltration_v2.R
│   ├── 03b_immune_checkpoint.R
│   ├── 04_scrna_CHB.R
│   ├── 05_network_MianYinChen_CHB.R
│   ├── 05c_drug_sensitivity_v2.R
│   ├── 06_molecular_docking_heatmap.R
│   ├── 06_publication_figures_v2.R
│   └── 07_final_publication_figures.R
├── data/
│   └── references/
│       ├── ferroptosis_genes_expanded.csv    (86 genes)
│       └── tcm_targets_CHB.csv              (42 targets)
└── plots/
    └── publication/
        ├── Figure1_flowchart.jpeg
        ├── Figure2-6_*.png (publication quality)
        ├── Figure2-6_*.tiff (300 dpi, journal submission)
        └── README_figures.md (detailed specifications)
```

## 📊 Data Sources

| Source | Type | Size | Reference |
|--------|------|------|-----------|
| **GSE14520** | HCC cohort (Tumor vs Non-tumor) | n=445 (225 tumor) | GEO Database |
| **TCGA-LIHC** | External validation cohort | n=369 | TCGA |
| **FerrDb V2** | Ferroptosis gene database | 108 genes | Zhou et al. |
| **TCMSP** | Traditional Chinese medicine targets | 42 targets | Artemisia capillaris |

## 🚀 Quick Start

### System Requirements
- **R version:** ≥ 4.0.0
- **Internet connection** required for data download
- **RAM:** ≥ 8GB recommended
- **Disk space:** ≥ 10GB for intermediate files

### Installation & Execution

1. Clone this repository:
   ```bash
   git clone https://github.com/Zhenghongwei11/Artemisia-HCC-Ferroptosis.git
   cd fx_review_package
   ```

2. Install dependencies (automated):
   ```bash
   Rscript run_complete_pipeline.R
   ```
   The script automatically installs required packages:
   - tidyverse, ggplot2, ggpubr
   - survival, survminer, timeROC
   - ComplexHeatmap, VennDiagram
   - GEOquery, limma, WGCNA, glmnet

3. Final figures generated in: `plots/publication/`

### Advanced: Run Individual Modules
```bash
cd scripts/
Rscript 02c_prognostic_model_v2.R    # Risk model only
Rscript 03a_immune_infiltration_v2.R # Immune analysis only
Rscript 06_molecular_docking_heatmap.R # Docking results
```

---

## 📈 Key Statistical Results

| Metric | Value | 95% CI |
|--------|-------|--------|
| **C-index** | 0.717 | 0.681-0.753 |
| **Log-rank P** | 3.1×10⁻⁸ | - |
| **1-year AUC** | 0.742 | 0.665-0.819 |
| **3-year AUC** | 0.765 | 0.695-0.835 |
| **5-year AUC** | 0.760 | 0.691-0.829 |
| **MDSC correlation** | r=0.48 | p<1×10⁻¹³ |
| **M2 Macrophage** | r=0.42 | p<1×10⁻¹⁰ |

---

## � About This Repository

This repository contains the **complete computational pipeline** for the manuscript manuscript (currently under review). It is designed to ensure:
- **Full reproducibility** of all analyses
- **Transparency** in methodology and results
- **Accessibility** for peer reviewers and researchers

### What's Included
- ✅ All 15 R analysis scripts (complete source code)
- ✅ Reference data (ferroptosis genes, TCM targets)
- ✅ Publication-quality figures (PNG, PDF, TIFF 300dpi)
- ✅ Figure specifications and generation metadata
- ✅ Detailed statistical results and computational logs

### What's NOT Included
- ⚠️ Raw datasets are downloaded automatically from public repositories (GEO, TCGA)
- ⚠️ Large intermediate RDS files are generated during pipeline execution

---

## 🔄 For Reviewers & Reproducibility

To verify all analyses:
1. Clone this repository
2. Run `Rscript run_complete_pipeline.R` 
3. All figures and statistical results will be regenerated

Estimated runtime: 2-4 hours (depending on internet speed and system)

---

## �🔬 Molecular Docking Methods

### NVIDIA Boltz-2 (AI-based)
- **Model:** MIT Boltz-2
- **Best result:** ACSL4 + Scoparone (ipTM=0.833)
- **API:** https://build.nvidia.com/mit/boltz2

### CB-Dock2 (Traditional)
- **Method:** AutoDock Vina
- **Best result:** NQO1 + Isorhamnetin (-9.1 kcal/mol)
- **Coverage:** 18 combinations (3 proteins × 6 ligands)

---

## � References & Acknowledgments

### Key Databases Used
- **GEO Database:** GSE14520, GSE83148
- **TCGA:** TCGA-LIHC (external validation)
- **FerrDb V2:** Ferroptosis gene database
- **TCMSP:** Traditional Chinese medicine targets

### Computational Tools
- R packages: tidyverse, ggplot2, survival, WGCNA, ComplexHeatmap
- Molecular docking: CB-Dock2, NVIDIA Boltz-2
- Statistical analysis: Univariate/multivariate Cox regression, LASSO-Cox, SSGSEA

---

## ⚠️ Important Notes

- **Manuscript Status:** Under review (not yet published)
- **Data Access:** Raw data automatically downloaded from public repositories during execution
- **Reproducibility:** All analyses are fully version-controlled and reproducible
- **Code Quality:** Production-ready R scripts with error handling and logging
- **Figure 1:** AI-assisted visualization (workflow diagram); scientific design and analysis entirely author-designed

---

## 📧 Support

For technical questions about reproduction:
- Check individual script headers for specific dependencies
- Review `plots/publication/README_figures.md` for detailed figure specifications
- Verify internet connection for public database downloads
