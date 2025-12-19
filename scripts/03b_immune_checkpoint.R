#!/usr/bin/env Rscript

# 03b_immune_checkpoint.R
# Immune checkpoint analysis + simplified TIDE-like scores (for manuscript figures).
#
# Usage:
#   cd fx_review_package
#   Rscript scripts/03b_immune_checkpoint.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(corrplot)
})

if (!dir.exists("data/processed")) stop("Please run this script from the fx_review_package/ directory.")

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"
set.seed(123)

message("[Checkpoint] Start immune checkpoint analysis...")

expr_symbol_file <- file.path(proc_dir, "GSE14520_expr_symbol.rds")
if (!file.exists(expr_symbol_file)) {
  stop("Missing file: data/processed/GSE14520_expr_symbol.rds (run scripts/03a_immune_infiltration_v2.R first).")
}
expr_14 <- readRDS(expr_symbol_file)

risk_file <- file.path(res_dir, "risk_score_data_v2.csv")
if (!file.exists(risk_file)) {
  stop("Missing file: results/risk_score_data_v2.csv (run scripts/02c_prognostic_model_v2.R first).")
}
risk_data <- read.csv(risk_file)

checkpoint_genes <- list(
  Inhibitory = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA", "VSIR", "IDO1", "SIGLEC15"),
  Stimulatory = c("CD28", "ICOS", "TNFRSF4", "TNFRSF9", "CD40", "CD40LG", "TNFRSF18", "CD27", "CD70")
)
all_checkpoints <- unlist(checkpoint_genes)

available_checkpoints <- intersect(all_checkpoints, rownames(expr_14))
checkpoint_expr <- expr_14[available_checkpoints, , drop = FALSE]

common_samples <- intersect(risk_data$sample, colnames(checkpoint_expr))
checkpoint_matched <- checkpoint_expr[, common_samples, drop = FALSE]
risk_matched <- risk_data[match(common_samples, risk_data$sample), ]

diff_results <- data.frame(
  Gene = rownames(checkpoint_matched),
  High_Mean = NA_real_,
  Low_Mean = NA_real_,
  Log2FC = NA_real_,
  P_value = NA_real_
)

for (i in seq_len(nrow(checkpoint_matched))) {
  high_vals <- as.numeric(checkpoint_matched[i, risk_matched$risk_group == "High"])
  low_vals <- as.numeric(checkpoint_matched[i, risk_matched$risk_group == "Low"])
  diff_results$High_Mean[i] <- mean(high_vals, na.rm = TRUE)
  diff_results$Low_Mean[i] <- mean(low_vals, na.rm = TRUE)
  diff_results$Log2FC[i] <- log2((diff_results$High_Mean[i] + 0.01) / (diff_results$Low_Mean[i] + 0.01))
  if (sum(!is.na(high_vals)) > 3 && sum(!is.na(low_vals)) > 3) {
    diff_results$P_value[i] <- wilcox.test(high_vals, low_vals)$p.value
  }
}

diff_results <- diff_results %>% arrange(P_value)
write.csv(diff_results, file.path(res_dir, "checkpoint_risk_diff.csv"), row.names = FALSE)

cor_results <- data.frame(
  Gene = rownames(checkpoint_matched),
  Correlation = NA_real_,
  P_value = NA_real_
)

for (i in seq_len(nrow(checkpoint_matched))) {
  ct <- cor.test(as.numeric(checkpoint_matched[i, ]), risk_matched$risk_score, method = "spearman")
  cor_results$Correlation[i] <- unname(ct$estimate)
  cor_results$P_value[i] <- ct$p.value
}

cor_results <- cor_results %>% arrange(P_value)
write.csv(cor_results, file.path(res_dir, "checkpoint_risk_correlation.csv"), row.names = FALSE)

# simplified TIDE-like score (for internal comparison only)
dysfunction_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
available_dys <- intersect(dysfunction_genes, rownames(expr_14))
dysfunction_score <- if (length(available_dys) > 0) colMeans(expr_14[available_dys, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(expr_14))

exclusion_genes <- c("FAP", "ACTA2", "COL1A1", "COL1A2", "TGFB1")
available_exc <- intersect(exclusion_genes, rownames(expr_14))
exclusion_score <- if (length(available_exc) > 0) colMeans(expr_14[available_exc, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(expr_14))

tide_score <- scale(dysfunction_score) + scale(exclusion_score)
names(tide_score) <- colnames(expr_14)

tide_data <- data.frame(
  sample = names(tide_score),
  dysfunction_score = dysfunction_score,
  exclusion_score = exclusion_score,
  tide_score = as.numeric(tide_score)
) %>%
  left_join(risk_data[, c("sample", "risk_score", "risk_group")], by = "sample")

write.csv(tide_data, file.path(res_dir, "tide_prediction.csv"), row.names = FALSE)

pdf(file.path(plot_dir, "Figure9_immune_checkpoint.pdf"), width = 14, height = 12)
par(mfrow = c(2, 2))

checkpoint_scaled <- t(scale(t(as.matrix(checkpoint_matched))))
checkpoint_scaled[checkpoint_scaled > 2] <- 2
checkpoint_scaled[checkpoint_scaled < -2] <- -2
order_idx <- order(risk_matched$risk_group, risk_matched$risk_score)
checkpoint_ordered <- checkpoint_scaled[, order_idx, drop = FALSE]
group_colors <- ifelse(risk_matched$risk_group[order_idx] == "High", "red", "blue")
heatmap(checkpoint_ordered, main = "A. Immune Checkpoint Expression",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none", Colv = NA, cexRow = 0.7, ColSideColors = group_colors)

sig_checkpoints <- diff_results %>% filter(!is.na(P_value) & P_value < 0.1) %>% head(10)
if (nrow(sig_checkpoints) > 0) {
  diff_mat <- as.matrix(sig_checkpoints[, c("High_Mean", "Low_Mean")])
  rownames(diff_mat) <- sig_checkpoints$Gene
  barplot(t(diff_mat), beside = TRUE, col = c("red", "blue"), las = 2,
          main = "B. Checkpoint Expression: High vs Low Risk", ylab = "Expression Level")
  legend("topright", legend = c("High Risk", "Low Risk"), fill = c("red", "blue"))
}

top_cor <- head(cor_results, 10)
barplot(top_cor$Correlation, names.arg = top_cor$Gene,
        col = ifelse(top_cor$Correlation > 0, "red", "blue"),
        las = 2, main = "C. Checkpoint - Risk Score Correlation", ylab = "Spearman Correlation")
abline(h = 0, lty = 2)

tide_complete <- tide_data[complete.cases(tide_data[, c("tide_score", "risk_group")]), ]
if (nrow(tide_complete) > 0) {
  boxplot(tide_score ~ risk_group, data = tide_complete, col = c("blue", "red"),
          main = "D. TIDE-like Score by Risk Group", xlab = "Risk Group", ylab = "Score")
}

dev.off()

message("[Checkpoint] Done.")
