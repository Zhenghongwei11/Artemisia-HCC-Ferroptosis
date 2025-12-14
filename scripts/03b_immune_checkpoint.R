#!/usr/bin/env Rscript

# 03b_immune_checkpoint.R
# 免疫检查点分析 + TIDE免疫治疗预测
# 使用方法: Rscript scripts_final/03b_immune_checkpoint.R
# 
# 输入: 表达矩阵, Risk Score数据
# 输出: 免疫检查点分析结果, Figure 9

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(corrplot)
})

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"

set.seed(123)

message("[免疫检查点] 开始免疫检查点分析...")

# ============================================
# 1. 加载数据
# ============================================
message("[免疫检查点] 加载数据...")

# 使用基因符号版本的表达矩阵
expr_symbol_file <- file.path(proc_dir, "GSE14520_expr_symbol.rds")
if (file.exists(expr_symbol_file)) {
  expr_14 <- readRDS(expr_symbol_file)
  message("[免疫检查点] 使用基因符号版本表达矩阵")
} else {
  stop("请先运行03a_immune_infiltration_v2.R生成基因符号版本表达矩阵")
}

risk_file <- file.path(res_dir, "risk_score_data.csv")
if (file.exists(risk_file)) {
  risk_data <- read.csv(risk_file)
} else {
  stop("找不到Risk Score数据")
}

message("[免疫检查点] 样本数: ", ncol(expr_14))

# ============================================
# 2. 定义免疫检查点基因
# ============================================
message("[免疫检查点] 定义免疫检查点基因...")

# 免疫检查点基因列表
checkpoint_genes <- list(
  # 抑制性检查点
  "Inhibitory" = c(
    "PDCD1",      # PD-1
    "CD274",      # PD-L1
    "PDCD1LG2",   # PD-L2
    "CTLA4",      # CTLA-4
    "LAG3",       # LAG-3
    "HAVCR2",     # TIM-3
    "TIGIT",      # TIGIT
    "BTLA",       # BTLA
    "VSIR",       # VISTA
    "IDO1",       # IDO1
    "SIGLEC15"    # Siglec-15
  ),
  # 刺激性检查点
  "Stimulatory" = c(
    "CD28",       # CD28
    "ICOS",       # ICOS
    "TNFRSF4",    # OX40
    "TNFRSF9",    # 4-1BB
    "CD40",       # CD40
    "CD40LG",     # CD40L
    "TNFRSF18",   # GITR
    "CD27",       # CD27
    "CD70"        # CD70
  )
)

all_checkpoints <- unlist(checkpoint_genes)

# ============================================
# 3. 提取免疫检查点表达
# ============================================
message("[免疫检查点] 提取免疫检查点表达...")

available_checkpoints <- intersect(all_checkpoints, rownames(expr_14))
message("[免疫检查点] 可用检查点基因: ", length(available_checkpoints), "/", length(all_checkpoints))

checkpoint_expr <- expr_14[available_checkpoints, , drop = FALSE]

# ============================================
# 4. 高低风险组检查点表达差异
# ============================================
message("[免疫检查点] 分析高低风险组差异...")

common_samples <- intersect(risk_data$sample, colnames(checkpoint_expr))
checkpoint_matched <- checkpoint_expr[, common_samples]
risk_matched <- risk_data[match(common_samples, risk_data$sample), ]

diff_results <- data.frame(
  Gene = rownames(checkpoint_matched),
  High_Mean = numeric(nrow(checkpoint_matched)),
  Low_Mean = numeric(nrow(checkpoint_matched)),
  Log2FC = numeric(nrow(checkpoint_matched)),
  P_value = numeric(nrow(checkpoint_matched)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(checkpoint_matched)) {
  high_vals <- as.numeric(checkpoint_matched[i, risk_matched$risk_group == "High"])
  low_vals <- as.numeric(checkpoint_matched[i, risk_matched$risk_group == "Low"])
  
  diff_results$High_Mean[i] <- mean(high_vals, na.rm = TRUE)
  diff_results$Low_Mean[i] <- mean(low_vals, na.rm = TRUE)
  diff_results$Log2FC[i] <- log2((diff_results$High_Mean[i] + 0.01) / (diff_results$Low_Mean[i] + 0.01))
  
  if (length(high_vals) > 3 && length(low_vals) > 3) {
    diff_results$P_value[i] <- wilcox.test(high_vals, low_vals)$p.value
  }
}

diff_results <- diff_results %>% arrange(P_value)
write.csv(diff_results, file.path(res_dir, "checkpoint_risk_diff.csv"), row.names = FALSE)

message("[免疫检查点] 显著差异的检查点: ", sum(diff_results$P_value < 0.05))

# ============================================
# 5. 检查点与Risk Score相关性
# ============================================
message("[免疫检查点] 计算与Risk Score相关性...")

cor_results <- data.frame(
  Gene = rownames(checkpoint_matched),
  Correlation = numeric(nrow(checkpoint_matched)),
  P_value = numeric(nrow(checkpoint_matched))
)

for (i in 1:nrow(checkpoint_matched)) {
  cor_test <- cor.test(as.numeric(checkpoint_matched[i, ]), risk_matched$risk_score,
                       method = "spearman")
  cor_results$Correlation[i] <- cor_test$estimate
  cor_results$P_value[i] <- cor_test$p.value
}

cor_results <- cor_results %>% arrange(P_value)
write.csv(cor_results, file.path(res_dir, "checkpoint_risk_correlation.csv"), row.names = FALSE)

# ============================================
# 6. TIDE评分计算 (简化版)
# ============================================
message("[免疫检查点] 计算TIDE评分...")

# TIDE评分简化计算
# 实际应使用TIDE官网 (http://tide.dfci.harvard.edu/)
# 这里使用简化的替代方法

# 计算免疫排斥评分 (基于T细胞功能障碍相关基因)
dysfunction_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
available_dys <- intersect(dysfunction_genes, rownames(expr_14))

if (length(available_dys) > 0) {
  dysfunction_score <- colMeans(expr_14[available_dys, , drop = FALSE], na.rm = TRUE)
} else {
  dysfunction_score <- rep(NA, ncol(expr_14))
}

# 计算免疫排斥评分 (基于CAF和MDSC相关基因)
exclusion_genes <- c("FAP", "ACTA2", "COL1A1", "COL1A2", "TGFB1")
available_exc <- intersect(exclusion_genes, rownames(expr_14))

if (length(available_exc) > 0) {
  exclusion_score <- colMeans(expr_14[available_exc, , drop = FALSE], na.rm = TRUE)
} else {
  exclusion_score <- rep(NA, ncol(expr_14))
}

# 简化TIDE评分 = 功能障碍 + 排斥
tide_score <- scale(dysfunction_score) + scale(exclusion_score)
names(tide_score) <- colnames(expr_14)

# 保存TIDE评分
tide_data <- data.frame(
  sample = names(tide_score),
  dysfunction_score = dysfunction_score,
  exclusion_score = exclusion_score,
  tide_score = as.numeric(tide_score)
)

# 匹配Risk数据
tide_data <- tide_data %>%
  left_join(risk_data[, c("sample", "risk_score", "risk_group")], by = "sample")

write.csv(tide_data, file.path(res_dir, "tide_prediction.csv"), row.names = FALSE)

# ============================================
# 7. 生成Figure 9: 免疫检查点多面板图
# ============================================
message("[免疫检查点] 生成Figure 9...")

pdf(file.path(plot_dir, "Figure9_immune_checkpoint.pdf"), width = 14, height = 12)

par(mfrow = c(2, 2))

# A: 免疫检查点表达热图
checkpoint_scaled <- t(scale(t(as.matrix(checkpoint_matched))))
checkpoint_scaled[checkpoint_scaled > 2] <- 2
checkpoint_scaled[checkpoint_scaled < -2] <- -2

# 按Risk分组排序
order_idx <- order(risk_matched$risk_group, risk_matched$risk_score)
checkpoint_ordered <- checkpoint_scaled[, order_idx]

# 颜色条
group_colors <- ifelse(risk_matched$risk_group[order_idx] == "High", "red", "blue")

heatmap(checkpoint_ordered,
        main = "A. Immune Checkpoint Expression",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none", Colv = NA, cexRow = 0.7,
        ColSideColors = group_colors)

# B: 高低风险组检查点差异
sig_checkpoints <- diff_results %>% filter(P_value < 0.1) %>% head(10)

if (nrow(sig_checkpoints) > 0) {
  diff_mat <- as.matrix(sig_checkpoints[, c("High_Mean", "Low_Mean")])
  rownames(diff_mat) <- sig_checkpoints$Gene
  
  barplot(t(diff_mat), beside = TRUE,
          col = c("red", "blue"), las = 2,
          main = "B. Checkpoint Expression: High vs Low Risk",
          ylab = "Expression Level")
  legend("topright", legend = c("High Risk", "Low Risk"),
         fill = c("red", "blue"))
}

# C: 检查点与Risk Score相关性
top_cor <- head(cor_results, 10)
barplot(top_cor$Correlation, names.arg = top_cor$Gene,
        col = ifelse(top_cor$Correlation > 0, "red", "blue"),
        las = 2, main = "C. Checkpoint - Risk Score Correlation",
        ylab = "Spearman Correlation")
abline(h = 0, lty = 2)

# D: TIDE评分与Risk分组
if (!all(is.na(tide_data$tide_score))) {
  tide_complete <- tide_data[complete.cases(tide_data[, c("tide_score", "risk_group")]), ]
  
  boxplot(tide_score ~ risk_group, data = tide_complete,
          col = c("blue", "red"),
          main = "D. TIDE Score by Risk Group",
          xlab = "Risk Group", ylab = "TIDE Score")
  
  # 添加p值
  if (nrow(tide_complete) > 10) {
    p_val <- wilcox.test(tide_score ~ risk_group, data = tide_complete)$p.value
    mtext(paste("p =", format(p_val, digits = 3)), side = 3, line = -2)
  }
}

dev.off()

message("[免疫检查点] Figure 9已保存")

# ============================================
# 8. 免疫治疗响应预测
# ============================================
message("[免疫检查点] 预测免疫治疗响应...")

# 基于TIDE评分预测响应
tide_complete <- tide_data[complete.cases(tide_data$tide_score), ]
tide_complete$predicted_response <- ifelse(tide_complete$tide_score < median(tide_complete$tide_score, na.rm = TRUE),
                                           "Responder", "Non-responder")

# 统计
response_by_risk <- table(tide_complete$risk_group, tide_complete$predicted_response)
message("[免疫检查点] 免疫治疗响应预测:")
print(response_by_risk)

write.csv(response_by_risk, file.path(res_dir, "immunotherapy_response_prediction.csv"))

# ============================================
# 完成
# ============================================
message("[免疫检查点] ✅ 免疫检查点分析完成！")
message("[免疫检查点] 输出文件:")
message("  ✅ ", file.path(res_dir, "checkpoint_risk_diff.csv"))
message("  ✅ ", file.path(res_dir, "checkpoint_risk_correlation.csv"))
message("  ✅ ", file.path(res_dir, "tide_prediction.csv"))
message("  ✅ ", file.path(res_dir, "immunotherapy_response_prediction.csv"))
message("  ✅ ", file.path(plot_dir, "Figure9_immune_checkpoint.pdf"))
