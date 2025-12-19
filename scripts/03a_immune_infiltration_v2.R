#!/usr/bin/env Rscript

# 03a_immune_infiltration_v2.R
# Immune infiltration analysis (ssGSEA) for GSE14520.
#
# Usage:
#   cd fx_review_package
#   Rscript scripts/03a_immune_infiltration_v2.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(pheatmap)
  library(ggpubr)
})

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"
set.seed(123)

message("[Immune V2] Start immune infiltration analysis...")

# ============================================
# 1. Load data and map probe IDs to gene symbols
# ============================================
message("[Immune V2] Load GPL571 annotation and map probe IDs...")

anno_file <- "data/raw/GSE14520/GPL571.annot.gz"
if (!file.exists(anno_file)) {
  stop(
    "Missing file: data/raw/GSE14520/GPL571.annot.gz\n",
    "See docs/DATA_MANIFEST.md for where to download and where to place it.\n"
  )
}

all_lines <- readLines(gzfile(anno_file))
header_line <- grep("^ID\t", all_lines)[1]
if (is.na(header_line)) stop("Cannot find GPL571 annotation table header (line starting with 'ID\\t').")

anno <- read.delim(
  text = all_lines[header_line:length(all_lines)],
  stringsAsFactors = FALSE,
  sep = "\t"
)

gene_col <- grep("Gene.symbol|Gene symbol", colnames(anno), value = TRUE, ignore.case = TRUE)[1]
if (is.na(gene_col)) stop("Cannot find gene symbol column in GPL571 annotation.")

probe2gene <- anno %>%
  dplyr::select(ID, all_of(gene_col)) %>%
  rename(Gene.Symbol = all_of(gene_col)) %>%
  filter(Gene.Symbol != "" & !is.na(Gene.Symbol)) %>%
  mutate(Gene.Symbol = sapply(strsplit(Gene.Symbol, "///"), function(x) trimws(x[1])))

message("[Immune V2] Valid probe-to-gene mappings: ", nrow(probe2gene))

expr_file <- file.path(proc_dir, "GSE14520_expr.rds")
if (!file.exists(expr_file)) {
  stop(
    "Missing file: data/processed/GSE14520_expr.rds\n",
    "If you want to re-generate it, run scripts/01b_download_GSE14520.R after placing GEO series matrix files.\n",
    "See docs/DATA_MANIFEST.md.\n"
  )
}

expr_14 <- readRDS(expr_file)
message("[Immune V2] Expression matrix: ", nrow(expr_14), " x ", ncol(expr_14))

expr_df <- as.data.frame(expr_14)
expr_df$ProbeID <- rownames(expr_df)

expr_gene <- expr_df %>%
  inner_join(probe2gene, by = c("ProbeID" = "ID")) %>%
  dplyr::select(-ProbeID) %>%
  group_by(Gene.Symbol) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames("Gene.Symbol")

expr_gene <- as.matrix(expr_gene)
message("[Immune V2] Gene-level matrix: ", nrow(expr_gene), " x ", ncol(expr_gene))

saveRDS(expr_gene, file.path(proc_dir, "GSE14520_expr_symbol.rds"))

# Risk Score (manuscript uses v2)
risk_file <- file.path(res_dir, "risk_score_data_v2.csv")
if (!file.exists(risk_file)) {
  stop(
    "Missing file: results/risk_score_data_v2.csv\n",
    "If you only need to regenerate Figures 2–6, run: Rscript run_complete_pipeline.R\n",
    "If you want to re-build the prognostic model, run: Rscript scripts/02c_prognostic_model_v2.R\n"
  )
}
risk_data <- read.csv(risk_file)
message("[Immune V2] Risk Score (v2) samples: ", nrow(risk_data))

# ============================================
# 2. Immune signatures (ssGSEA)
# ============================================
message("[Immune V2] Prepare immune signatures...")

immune_signatures <- list(
  "Activated_B_cell" = c("CD19", "CD79A", "CD79B", "MS4A1", "CD22"),
  "Activated_CD4_T_cell" = c("CD4", "IL2RA", "CD40LG", "ICOS", "CTLA4"),
  "Activated_CD8_T_cell" = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG"),
  "Central_memory_CD4_T_cell" = c("CD4", "CCR7", "SELL", "IL7R"),
  "Central_memory_CD8_T_cell" = c("CD8A", "CCR7", "SELL", "IL7R"),
  "Effector_memory_CD4_T_cell" = c("CD4", "GZMK", "CCL5"),
  "Effector_memory_CD8_T_cell" = c("CD8A", "GZMK", "CCL5", "NKG7"),
  "Gamma_delta_T_cell" = c("TRGC1", "TRGC2", "TRDC"),
  "Regulatory_T_cell" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"),
  "T_follicular_helper_cell" = c("CXCR5", "ICOS", "BCL6", "CD40LG"),
  "Type_1_T_helper_cell" = c("TBX21", "IFNG", "IL12RB2"),
  "Type_2_T_helper_cell" = c("GATA3", "IL4", "IL5", "IL13"),
  "Type_17_T_helper_cell" = c("RORC", "IL17A", "IL17F", "IL22"),
  "Activated_dendritic_cell" = c("CD80", "CD86", "CD83", "CCR7"),
  "Immature_dendritic_cell" = c("CD1A", "CD1C", "CLEC10A"),
  "Plasmacytoid_dendritic_cell" = c("CLEC4C", "IL3RA", "NRP1"),
  "Macrophage" = c("CD68", "CD163", "MSR1", "MRC1"),
  "M1_Macrophage" = c("NOS2", "IL1B", "IL6", "TNF", "CD80"),
  "M2_Macrophage" = c("CD163", "MRC1", "ARG1", "IL10"),
  "Monocyte" = c("CD14", "FCGR3A", "CSF1R"),
  "Natural_killer_cell" = c("NCAM1", "NKG7", "KLRD1", "KLRB1"),
  "Neutrophil" = c("FCGR3B", "CXCR1", "CXCR2", "CSF3R"),
  "Mast_cell" = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  "Eosinophil" = c("CCR3", "SIGLEC8", "PRG2"),
  "MDSC" = c("CD33", "ITGAM", "ARG1", "NOS2"),
  "Natural_killer_T_cell" = c("CD3D", "NCAM1", "KLRB1")
)

all_immune_genes <- unique(unlist(immune_signatures))
matched_genes <- intersect(all_immune_genes, rownames(expr_gene))
message("[Immune V2] Immune gene matched: ", length(matched_genes), "/", length(all_immune_genes))

message("[Immune V2] Run ssGSEA...")
ssgsea_scores <- gsva(expr_gene, immune_signatures, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
ssgsea_scores <- as.matrix(ssgsea_scores)

if (nrow(ssgsea_scores) == 0) stop("ssGSEA returned empty scores matrix.")

write.csv(ssgsea_scores, file.path(res_dir, "ssgsea_immune_scores.csv"))

# ============================================
# 3. Correlation with risk score
# ============================================
message("[Immune V2] Correlate immune cells with risk score...")

common_samples <- intersect(risk_data$sample, colnames(ssgsea_scores))
message("[Immune V2] Matched samples: ", length(common_samples))

ssgsea_matched <- ssgsea_scores[, common_samples, drop = FALSE]
risk_matched <- risk_data[match(common_samples, risk_data$sample), ]

cor_results <- data.frame(
  Cell_Type = rownames(ssgsea_matched),
  Correlation = NA_real_,
  P_value = NA_real_
)

for (i in seq_len(nrow(ssgsea_matched))) {
  ct <- cor.test(as.numeric(ssgsea_matched[i, ]), risk_matched$risk_score, method = "spearman")
  cor_results$Correlation[i] <- unname(ct$estimate)
  cor_results$P_value[i] <- ct$p.value
}

cor_results <- cor_results %>% arrange(P_value)
write.csv(cor_results, file.path(res_dir, "immune_risk_correlation.csv"), row.names = FALSE)

# ============================================
# 4. Figure output (optional)
# ============================================
pdf(file.path(plot_dir, "Figure8_immune_infiltration.pdf"), width = 16, height = 14)
par(mfrow = c(2, 2))

heatmap_data <- ssgsea_scores
if (ncol(heatmap_data) > 50) {
  set.seed(123)
  heatmap_data <- heatmap_data[, sample(seq_len(ncol(heatmap_data)), 50), drop = FALSE]
}
heatmap_data <- heatmap_data[complete.cases(heatmap_data), , drop = FALSE]
if (nrow(heatmap_data) > 0) {
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  heatmap_data_scaled[is.na(heatmap_data_scaled)] <- 0
  heatmap_data_scaled[heatmap_data_scaled > 2] <- 2
  heatmap_data_scaled[heatmap_data_scaled < -2] <- -2
  heatmap(
    heatmap_data_scaled,
    main = "A. Immune Cell Infiltration (ssGSEA)",
    col = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    cexRow = 0.6,
    cexCol = 0.5
  )
}

top_cor <- head(cor_results, 15)
barplot(
  top_cor$Correlation,
  names.arg = top_cor$Cell_Type,
  col = ifelse(top_cor$Correlation > 0, "red", "blue"),
  las = 2,
  cex.names = 0.6,
  main = "B. Immune Cell - Risk Score Correlation",
  ylab = "Spearman Correlation"
)
abline(h = 0, lty = 2)

dev.off()

message("[Immune V2] Done.")
