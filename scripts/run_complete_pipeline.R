#!/usr/bin/env Rscript

# run_complete_pipeline.R (optional full re-run driver)
# 用途：按“原始工程”脚本顺序串联执行（可选）。
#
# 重要：
# - 本脚本不等同于“零数据一键全自动”。多处输入需要手动下载/放置（见 docs/DATA_MANIFEST.md）。
# - 审稿复现 Figures 2–6 建议只运行仓库根目录的 `run_complete_pipeline.R`（从 results/ 复现出图）。
# 
# 使用方法:
#   cd fx_review_package
#   Rscript scripts/run_complete_pipeline.R

# ============================================================================
# 禁用交互式图形设备
# ============================================================================
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

message("================================================================================")
message("CHB 绵茵陈-铁死亡-免疫研究 - 可选全流程串联脚本")
message("================================================================================")
message("该脚本需要手动准备部分输入（见 docs/DATA_MANIFEST.md）。")
message("开始时间: ", Sys.time())

# 确保在包根目录运行（允许从 scripts/ 目录启动）
if (!dir.exists("data/processed")) {
  if (dir.exists("../data/processed")) {
    setwd("..")
  } else {
    stop("错误: 请在 fx_review_package/ 根目录运行此脚本")
  }
}

# ============================================================================
# Stage 1: 多队列整合与验证
# ============================================================================
message("\n[Stage 1] 数据准备与基础分析...")
source("scripts/00_setup_env.R")
# source("scripts/01b_download_GSE14520.R")
source("scripts/02_batch_correction.R")
source("scripts/02b_multi_cohort_DEG.R")
source("scripts/02_bulk_WGCNA_ML_CHB.R")

# ============================================================================
# Stage 2: 临床关联分析 (核心更新部分)
# ============================================================================
message("\n[Stage 2] 预后模型构建 (V2)...")
# 使用 V2 版本的预后模型脚本
source("scripts/02c_prognostic_model_v2.R")
source("scripts/02d_nomogram_calibration_dca.R")
source("scripts/02e_external_validation.R")

# ============================================================================
# Stage 3: 免疫微环境分析
# ============================================================================
message("\n[Stage 3] 免疫微环境分析...")
# 使用 V2 或正确的脚本名称
if (file.exists("scripts/03a_immune_infiltration_v2.R")) {
  source("scripts/03a_immune_infiltration_v2.R")
} else {
  source("scripts/03a_immune_infiltration_enhanced.R")
}
source("scripts/03b_immune_checkpoint.R")

# ============================================================================
# Stage 4: 绵茵陈网络药理学与分子对接
# ============================================================================
message("\n[Stage 4] 网络药理与分子对接...")
source("scripts/05_network_MianYinChen_CHB.R")

# 药物敏感性 (V2)
if (file.exists("scripts/05c_drug_sensitivity_v2.R")) {
  source("scripts/05c_drug_sensitivity_v2.R")
} else {
  source("scripts/05c_drug_sensitivity.R")
}

# 分子对接热图绘制
if (file.exists("scripts/06_molecular_docking_heatmap.R")) {
  source("scripts/06_molecular_docking_heatmap.R")
}

# ============================================================================
# Stage 5: 最终出版级图表生成
# ============================================================================
message("\n[Stage 5] 生成最终出版级图表 (Figures 2-6)...")
source("scripts/07_final_publication_figures.R")

message("\n================================================================================")
message("✅ 分析管道全部完成！")
message("请查看 'plots/publication/' 目录获取最终图表。")
message("================================================================================")
