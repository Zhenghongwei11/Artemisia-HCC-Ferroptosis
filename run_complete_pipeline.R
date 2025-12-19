#!/usr/bin/env Rscript

# Entry point for peer-review reproducibility.
#
# This script regenerates manuscript Figures 2–6 from the curated result tables
# shipped in this package (see `results/`).
#
# Usage:
#   Rscript run_complete_pipeline.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

{
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    setwd(normalizePath(dirname(script_path)))
  }
}

message("============================================================")
message("fx_review_package: regenerate Figures 2–6")
message("Working directory: ", getwd())
message("============================================================")

required <- c(
  "results/deg_GSE14520_all.csv",
  "results/deg_GSE83148_all.csv",
  "results/risk_score_data_v2.csv",
  "results/prognostic_model_coef_v2.csv",
  "results/prognostic_model_stats_v2.csv",
  "results/univariate_cox_clinical.csv",
  "results/multivariate_cox_clinical.csv",
  "results/immune_risk_correlation.csv",
  "results/ssgsea_immune_scores.csv",
  "results/checkpoint_risk_diff.csv",
  "results/drug_risk_correlation.csv",
  "results/hcc_drug_risk_diff.csv",
  "results/cbdock2_docking_results.csv",
  "results/molecular_docking_results.csv",
  "data/processed/GSE14520_tumor_clinical.rds",
  "data/references/ferroptosis_genes_expanded.csv",
  "data/references/tcm_targets_CHB.csv"
)

missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  message("Missing required files:")
  for (m in missing) message("  - ", m)
  message("")
  message("See `docs/DATA_MANIFEST.md` for manual download / placement details.")
  quit(save = "no", status = 1)
}

source("scripts/07_final_publication_figures.R")

message("============================================================")
message("Done. Outputs are in `plots/publication/`.")
message("============================================================")
