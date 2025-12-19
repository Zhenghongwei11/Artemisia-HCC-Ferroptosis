# Artemisia capillaris targets ferroptosis in hepatocellular carcinoma (HCC)

This is the peer-review reproduction package for the manuscript’s computational results.

## What you can reproduce with this package

- Figures 2–6 can be regenerated from curated result tables included under `results/`.
- Figure 1 is an AI-assisted workflow diagram; the exact prompt is provided in `FIGURE_GENERATION_METADATA.md`.

The regenerated figures are written to `plots/publication/`.

## Quick start (recommended)

From the `fx_review_package/` directory:

```bash
Rscript run_complete_pipeline.R
```

This runs `scripts/07_final_publication_figures.R` and regenerates Figures 2–6.

## Manual downloads (only needed for a full re-run)

Some upstream steps in the original project require manual downloads (GEO series matrix files, GPL annotation, GDSC training data, TCGA access). This package does not claim “automatic data download”.

See:
- `docs/DATA_MANIFEST.md`
- `docs/REPRODUCIBILITY.md`

## Directory layout

```
fx_review_package/
├── run_complete_pipeline.R
├── FIGURE_GENERATION_METADATA.md
├── scripts/
├── data/
│   ├── processed/              # minimal processed inputs (see DATA_MANIFEST)
│   ├── raw/                    # raw files (manual placement if re-running upstream)
│   └── references/             # ferroptosis genes, TCM targets
├── results/                    # curated result tables (used by Figures 2–6)
└── plots/publication/          # publication-ready figures and index
```

## Notes

- `data/references/ferroptosis_genes_expanded.csv` contains 108 genes.
- `data/references/tcm_targets_CHB.csv` contains 42 targets.
