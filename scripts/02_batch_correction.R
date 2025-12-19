#!/usr/bin/env Rscript

# 02_batch_correction.R
# 合并GSE83148和GSE14520表达矩阵，去除批次效应
# 输入: GSE83148_exprs_log2.rds, GSE14520_expr.rds
# 输出: combined_expr_batch_corrected.rds, batch_correction_check.pdf

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(sva)
  library(limma)
  library(ggplot2)
  library(tidyverse)
})

# 设置工作目录（支持命令行和RStudio）
# 使用方法: cd fx_review_package && Rscript scripts/02_batch_correction.R
if (!dir.exists("data/processed")) {
  stop("错误: 请在项目根目录(fx_review_package)运行此脚本\n使用方法: cd fx_review_package && Rscript scripts/02_batch_correction.R")
}

proc_dir <- "data/processed"
plot_dir <- "plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("[批次校正] 开始加载数据...")

# 1. 加载表达矩阵
expr_83 <- readRDS(file.path(proc_dir, "GSE83148_exprs_log2.rds"))
expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))

message("[批次校正] GSE83148: ", dim(expr_83)[1], " x ", dim(expr_83)[2])
message("[批次校正] GSE14520: ", dim(expr_14)[1], " x ", dim(expr_14)[2])

# 2. 找到共同基因
common_genes <- intersect(rownames(expr_83), rownames(expr_14))
message("[批次校正] 共同基因数: ", length(common_genes))

if (length(common_genes) < 1000) {
  warning("[批次校正] 共同基因数过少，可能存在基因符号不一致的问题")
}

# 3. 提取共同基因的表达矩阵
expr_83_sub <- expr_83[common_genes, ]
expr_14_sub <- expr_14[common_genes, ]

# 4. 合并表达矩阵
combined_expr <- cbind(expr_83_sub, expr_14_sub)
message("[批次校正] 合并后表达矩阵: ", dim(combined_expr)[1], " x ", dim(combined_expr)[2])

# 5. 创建批次标签
batch <- c(rep("GSE83148", ncol(expr_83_sub)), rep("GSE14520", ncol(expr_14_sub)))
message("[批次校正] 批次分布: GSE83148=", sum(batch=="GSE83148"), ", GSE14520=", sum(batch=="GSE14520"))

# 6. 校正前PCA分析
message("[批次校正] 进行校正前PCA分析...")
pca_before <- prcomp(t(combined_expr), scale. = TRUE)
pca_before_df <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  Batch = batch
)

# 7. 应用ComBat批次校正
message("[批次校正] 应用ComBat算法...")
combined_expr_corrected <- ComBat(combined_expr, batch = batch, par.prior = TRUE, prior.plots = FALSE)
message("[批次校正] ComBat校正完成")

# 8. 校正后PCA分析
message("[批次校正] 进行校正后PCA分析...")
pca_after <- prcomp(t(combined_expr_corrected), scale. = TRUE)
pca_after_df <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Batch = batch
)

# 9. 生成对比图
message("[批次校正] 生成对比图...")
pdf(file.path(plot_dir, "batch_correction_check.pdf"), width = 12, height = 5)

# 校正前
p1 <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.6) +
  theme_bw() +
  ggtitle("Before Batch Correction") +
  xlab(paste0("PC1 (", round(summary(pca_before)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_before)$importance[2, 2] * 100, 1), "%)")) +
  theme(legend.position = "right")

# 校正后
p2 <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.6) +
  theme_bw() +
  ggtitle("After Batch Correction") +
  xlab(paste0("PC1 (", round(summary(pca_after)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_after)$importance[2, 2] * 100, 1), "%)")) +
  theme(legend.position = "right")

print(gridExtra::grid.arrange(p1, p2, ncol = 2))

dev.off()

message("[批次校正] 对比图已保存: ", file.path(plot_dir, "batch_correction_check.pdf"))

# 10. 保存校正后的表达矩阵
saveRDS(combined_expr_corrected, file.path(proc_dir, "combined_expr_batch_corrected.rds"))
message("[批次校正] 校正后表达矩阵已保存: ", file.path(proc_dir, "combined_expr_batch_corrected.rds"))

# 11. 生成校正效果报告
report <- data.frame(
  Metric = c("Total Genes", "Total Samples", "GSE83148 Samples", "GSE14520 Samples", 
             "PC1 Variance Before", "PC1 Variance After",
             "PC2 Variance Before", "PC2 Variance After"),
  Value = c(
    nrow(combined_expr_corrected),
    ncol(combined_expr_corrected),
    sum(batch == "GSE83148"),
    sum(batch == "GSE14520"),
    round(summary(pca_before)$importance[2, 1] * 100, 2),
    round(summary(pca_after)$importance[2, 1] * 100, 2),
    round(summary(pca_before)$importance[2, 2] * 100, 2),
    round(summary(pca_after)$importance[2, 2] * 100, 2)
  )
)

write.csv(report, file.path("results", "batch_correction_report.csv"), row.names = FALSE)
message("[批次校正] 校正报告已保存: results/batch_correction_report.csv")

message("[批次校正] ✅ 批次校正完成！")
message("[批次校正] 输出文件:")
message("  ✅ ", file.path(proc_dir, "combined_expr_batch_corrected.rds"))
message("  ✅ ", file.path(plot_dir, "batch_correction_check.pdf"))
message("  ✅ results/batch_correction_report.csv")
