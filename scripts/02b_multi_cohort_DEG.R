#!/usr/bin/env Rscript

# 02b_multi_cohort_DEG.R
# 多队列DEG分析 + Meta-analysis + 与铁死亡基因交集
# 输入: combined_expr_batch_corrected.rds, GSE83148_pdata.rds, GSE14520_clinical.rds
# 输出: DEG结果、Meta-analysis结果、Venn图、火山图

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(VennDiagram)
  library(grid)  # 用于grid.newpage()
  library(ggplot2)
  library(AnnotationDbi)
  library(hgu133plus2.db)  # Affymetrix HG-U133 Plus 2.0 注释包
})

# 禁用VennDiagram的日志输出
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir <- "results"
plot_dir <- "plots"

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("[多队列DEG] 开始加载数据...")

# 探针ID到基因符号的转换函数
probe_to_symbol <- function(probe_ids) {
  # 尝试使用hgu133plus2.db注释包
  tryCatch({
    symbols <- mapIds(hgu133plus2.db, 
                      keys = probe_ids,
                      column = "SYMBOL",
                      keytype = "PROBEID",
                      multiVals = "first")
    return(symbols)
  }, error = function(e) {
    message("  注释包转换失败，使用备用方法...")
    # 备用方法：从探针ID提取基因符号（如果有的话）
    return(rep(NA, length(probe_ids)))
  })
}

# 1. 加载表达矩阵和临床信息
expr_83 <- readRDS(file.path(proc_dir, "GSE83148_exprs_log2.rds"))
expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
pdata_83 <- readRDS(file.path(proc_dir, "GSE83148_pdata.rds"))
pdata_14 <- readRDS(file.path(proc_dir, "GSE14520_clinical.rds"))

# 2. 加载铁死亡基因集 (扩展版，108个基因)
ferro_genes <- read.csv(file.path("data/references", "ferroptosis_genes_expanded.csv"))$Gene

message("[多队列DEG] GSE83148: ", dim(expr_83)[1], " x ", dim(expr_83)[2])
message("[多队列DEG] GSE14520: ", dim(expr_14)[1], " x ", dim(expr_14)[2])
message("[多队列DEG] 铁死亡基因: ", length(ferro_genes))

# ============================================
# GSE83148 DEG分析
# ============================================
message("[多队列DEG] 开始GSE83148 DEG分析...")

# 构造分组信息
group_83 <- ifelse(grepl("control|normal|healthy", pdata_83$title, ignore.case = TRUE), "Control", "CHB")
group_83 <- factor(group_83, levels = c("Control", "CHB"))

message("[多队列DEG] GSE83148分组: Control=", sum(group_83=="Control"), ", CHB=", sum(group_83=="CHB"))

# Limma DEG分析
design_83 <- model.matrix(~0 + group_83)
colnames(design_83) <- levels(group_83)
fit_83 <- lmFit(expr_83, design_83)
contrast_83 <- makeContrasts(CHB - Control, levels = design_83)
fit_83 <- eBayes(contrasts.fit(fit_83, contrast_83))

deg_83 <- topTable(fit_83, adjust.method = "fdr", number = Inf)

# 添加基因符号列
message("[多队列DEG] 转换探针ID到基因符号...")
deg_83$ProbeID <- rownames(deg_83)
deg_83$Gene <- probe_to_symbol(rownames(deg_83))
deg_83 <- deg_83 %>% filter(!is.na(Gene) & Gene != "")

deg_83_sig <- deg_83 %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)

message("[多队列DEG] GSE83148显著DEG: ", nrow(deg_83_sig))

write.csv(deg_83, file.path(res_dir, "deg_GSE83148_all.csv"), row.names = FALSE)
write.csv(deg_83_sig, file.path(res_dir, "deg_GSE83148_sig.csv"), row.names = FALSE)

# ============================================
# GSE14520 DEG分析
# ============================================
message("[多队列DEG] 开始GSE14520 DEG分析...")

# 构造分组信息（HCC vs Normal）
# GSE14520样本标记为 "Liver Tumor Tissue" 和 "Liver Non-Tumor Tissue"
group_14 <- ifelse(grepl("Non-Tumor|non-tumor|Normal|normal|control|healthy", pdata_14$title, ignore.case = TRUE), "Normal", "HCC")
group_14 <- factor(group_14, levels = c("Normal", "HCC"))

message("[多队列DEG] GSE14520分组: Normal=", sum(group_14=="Normal"), ", HCC=", sum(group_14=="HCC"))

# Limma DEG分析
design_14 <- model.matrix(~0 + group_14)
colnames(design_14) <- levels(group_14)
fit_14 <- lmFit(expr_14, design_14)
contrast_14 <- makeContrasts(HCC - Normal, levels = design_14)
fit_14 <- eBayes(contrasts.fit(fit_14, contrast_14))

deg_14 <- topTable(fit_14, adjust.method = "fdr", number = Inf)

# 添加基因符号列
message("[多队列DEG] 转换GSE14520探针ID到基因符号...")
deg_14$ProbeID <- rownames(deg_14)
deg_14$Gene <- probe_to_symbol(rownames(deg_14))
deg_14 <- deg_14 %>% filter(!is.na(Gene) & Gene != "")

deg_14_sig <- deg_14 %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)

message("[多队列DEG] GSE14520显著DEG: ", nrow(deg_14_sig))

write.csv(deg_14, file.path(res_dir, "deg_GSE14520_all.csv"), row.names = FALSE)
write.csv(deg_14_sig, file.path(res_dir, "deg_GSE14520_sig.csv"), row.names = FALSE)

# ============================================
# 多队列交集分析
# ============================================
message("[多队列DEG] 进行多队列交集分析...")

# 两个队列都显著的基因 (使用基因符号)
common_deg <- intersect(deg_83_sig$Gene, deg_14_sig$Gene)
message("[多队列DEG] 两队列共同DEG: ", length(common_deg))

# 检查方向一致性
deg_83_common <- deg_83_sig %>% filter(Gene %in% common_deg)
deg_14_common <- deg_14_sig %>% filter(Gene %in% common_deg)

# 按基因符号排序后比较
deg_83_common <- deg_83_common %>% arrange(Gene)
deg_14_common <- deg_14_common %>% filter(Gene %in% deg_83_common$Gene) %>% arrange(Gene)

direction_consistent <- sum(sign(deg_83_common$logFC) == sign(deg_14_common$logFC))
message("[多队列DEG] 方向一致的基因: ", direction_consistent, "/", length(common_deg))

# ============================================
# 与铁死亡基因交集
# ============================================
message("[多队列DEG] 进行与铁死亡基因的交集...")

# 两队列都显著且与铁死亡基因相关 (使用Gene列而非rownames)
ferro_deg_83 <- intersect(deg_83_sig$Gene, ferro_genes)
ferro_deg_14 <- intersect(deg_14_sig$Gene, ferro_genes)
ferro_deg_common <- intersect(ferro_deg_83, ferro_deg_14)

message("[多队列DEG] GSE83148铁死亡DEG: ", length(ferro_deg_83))
message("[多队列DEG] GSE14520铁死亡DEG: ", length(ferro_deg_14))
message("[多队列DEG] 两队列共同铁死亡DEG: ", length(ferro_deg_common))

# 保存结果
ferro_deg_result <- data.frame(
  Gene = ferro_deg_common,
  logFC_83 = deg_83_sig$logFC[match(ferro_deg_common, deg_83_sig$Gene)],
  logFC_14 = deg_14_sig$logFC[match(ferro_deg_common, deg_14_sig$Gene)],
  adj.P.Val_83 = deg_83_sig$adj.P.Val[match(ferro_deg_common, deg_83_sig$Gene)],
  adj.P.Val_14 = deg_14_sig$adj.P.Val[match(ferro_deg_common, deg_14_sig$Gene)]
)

write.csv(ferro_deg_result, file.path(res_dir, "ferroptosis_DEG_intersection.csv"), row.names = FALSE)

# ============================================
# 生成Venn图
# ============================================
message("[多队列DEG] 生成Venn图...")

pdf(file.path(plot_dir, "multi_cohort_venn.pdf"), width = 8, height = 6)

venn_list <- list(
  GSE83148_DEG = rownames(deg_83_sig),
  GSE14520_DEG = rownames(deg_14_sig),
  Ferroptosis = ferro_genes
)

draw.triple.venn(
  area1 = length(venn_list[[1]]),
  area2 = length(venn_list[[2]]),
  area3 = length(venn_list[[3]]),
  n12 = length(intersect(venn_list[[1]], venn_list[[2]])),
  n23 = length(intersect(venn_list[[2]], venn_list[[3]])),
  n13 = length(intersect(venn_list[[1]], venn_list[[3]])),
  n123 = length(Reduce(intersect, venn_list)),
  category = c("GSE83148 DEG", "GSE14520 DEG", "Ferroptosis"),
  fill = c("skyblue", "pink", "lightgreen"),
  main = "Multi-cohort DEG and Ferroptosis Intersection"
)

dev.off()

message("[多队列DEG] Venn图已保存: ", file.path(plot_dir, "multi_cohort_venn.pdf"))

# ============================================
# 生成火山图
# ============================================
message("[多队列DEG] 生成火山图...")

# GSE83148火山图
pdf(file.path(plot_dir, "volcano_GSE83148.pdf"), width = 8, height = 6)

deg_83$color <- ifelse(deg_83$adj.P.Val < 0.05 & abs(deg_83$logFC) > 1, 
                       ifelse(rownames(deg_83) %in% ferro_genes, "Ferroptosis", "DEG"),
                       "Not Sig")

p1 <- ggplot(deg_83, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("DEG" = "blue", "Ferroptosis" = "red", "Not Sig" = "grey")) +
  theme_bw() +
  ggtitle("GSE83148: CHB vs Control") +
  xlab("log2 Fold Change") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

print(p1)
dev.off()

# GSE14520火山图
pdf(file.path(plot_dir, "volcano_GSE14520.pdf"), width = 8, height = 6)

deg_14$color <- ifelse(deg_14$adj.P.Val < 0.05 & abs(deg_14$logFC) > 1,
                       ifelse(rownames(deg_14) %in% ferro_genes, "Ferroptosis", "DEG"),
                       "Not Sig")

p2 <- ggplot(deg_14, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("DEG" = "blue", "Ferroptosis" = "red", "Not Sig" = "grey")) +
  theme_bw() +
  ggtitle("GSE14520: HCC vs Normal") +
  xlab("log2 Fold Change") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

print(p2)
dev.off()

message("[多队列DEG] 火山图已保存")

# ============================================
# 总结
# ============================================
message("[多队列DEG] ✅ 多队列DEG分析完成！")
message("[多队列DEG] 输出文件:")
message("  ✅ ", file.path(res_dir, "deg_GSE83148_sig.csv"))
message("  ✅ ", file.path(res_dir, "deg_GSE14520_sig.csv"))
message("  ✅ ", file.path(res_dir, "ferroptosis_DEG_intersection.csv"))
message("  ✅ ", file.path(plot_dir, "multi_cohort_venn.pdf"))
message("  ✅ ", file.path(plot_dir, "volcano_GSE83148.pdf"))
message("  ✅ ", file.path(plot_dir, "volcano_GSE14520.pdf"))

