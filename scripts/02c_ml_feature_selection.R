#!/usr/bin/env Rscript

# 02c_ml_feature_selection.R
# 机器学习特征筛选 - 主控脚本
# 协调预筛选和主分析两个阶段
# 使用方法: Rscript scripts/02c_ml_feature_selection.R

message("========================================")
message("机器学习特征筛选 - 完整流程")
message("========================================")

# 第一阶段：预筛选
message("\n[阶段1] 预筛选 - 单变量分析...")
source("scripts/02c_ml_feature_selection_prefilter.R")

# 检查预筛选结果
if (!file.exists("results/ml_candidate_genes_top50.csv")) {
  stop("预筛选失败：输出文件不存在")
}

message("\n[检查] 预筛选结果文件已生成 ✓")

# 第二阶段：主分析
message("\n[阶段2] 主分析 - SVM-RFE/RF/LASSO...")
source("scripts/02c_ml_feature_selection_main.R")

# 检查主分析结果
if (!file.exists("results/diagnostic_core_genes_CHB.csv")) {
  stop("主分析失败：输出文件不存在")
}

message("\n[检查] 主分析结果文件已生成 ✓")

# 最终验证
core_genes <- read.csv("results/diagnostic_core_genes_CHB.csv")$Gene
message("\n========================================")
message("✅ 机器学习特征筛选完成！")
message("核心基因数: ", length(core_genes))
message("输出文件: results/diagnostic_core_genes_CHB.csv")
message("========================================")
