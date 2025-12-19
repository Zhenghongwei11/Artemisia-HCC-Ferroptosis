# Reproducibility guide (for peer review)

## 1. What this package is meant to reproduce

This package focuses on reproducing the manuscript’s computational outputs in a reviewer-friendly way:

- Figures 2–6 can be regenerated from the included result tables.
- Figure 1 is an AI-assisted workflow diagram; the exact prompt is provided in `FIGURE_GENERATION_METADATA.md`.

## 2. Quick reproduction (recommended)

From the `fx_review_package/` directory:

```bash
Rscript run_complete_pipeline.R
```

Outputs:
- Regenerated figures: `plots/publication/Figure2_DEG_analysis.*`, `plots/publication/Figure3_prognostic_model.*`, `plots/publication/Figure4_clinical_utility.*`, `plots/publication/Figure5_immune_analysis.*`, `plots/publication/Figure6_drug_docking.*`

## 3. File-to-figure map

See:
- `plots/publication/figure_index.csv`
- `results/README_RESULTS.md`

## 4. Full re-run (optional, requires manual downloads)

If you want to re-run upstream steps from public raw data, follow `docs/DATA_MANIFEST.md` and then run scripts in this order (as needed):

1. `scripts/00_setup_env.R` (install packages)
2. `scripts/01b_download_GSE14520.R` (requires GEO series matrix files placed under `data/raw/GSE14520/`)
3. `scripts/02b_multi_cohort_DEG.R` (requires processed GSE83148 + GSE14520 files)
4. `scripts/02c_prognostic_model_v2.R`
5. `scripts/03a_immune_infiltration_v2.R` and `scripts/03b_immune_checkpoint.R`
6. `scripts/05c_drug_sensitivity_v2.R` (requires GDSC2 training data)
7. `scripts/07_final_publication_figures.R`
