#!/usr/bin/env Rscript

# run_complete_pipeline.R (Fixed Version)
# 完整的CHB分析管道 - 从头到尾的可复现分析
# 修正：引用了正确的 _v2 版本脚本
# 
# 使用方法: Rscript scripts_final/run_complete_pipeline.R

# ============================================================================
# 禁用交互式图形设备
# ============================================================================
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

message("================================================================================")
message("CHB 绵茵陈-铁死亡-免疫研究 - 完整分析管道 (Fixed)")
message("================================================================================")
message("开始时间: ", Sys.time())

# 检查工作目录
if (!dir.exists("data/processed")) {
  stop("错误: 请在项目根目录运行此脚本")
}

# ============================================================================
# Stage 1: 多队列整合与验证
# ============================================================================
message("\n[Stage 1] 数据准备与基础分析...")
source("scripts_final/00_setup_env.R")
# source("scripts_final/01b_download_GSE14520.R") # 数据已下载，跳过以节省时间
source("scripts_final/02_batch_correction.R")
source("scripts_final/02b_multi_cohort_DEG.R")
source("scripts_final/02_bulk_WGCNA_ML_CHB.R")

# ============================================================================
# Stage 2: 临床关联分析 (核心更新部分)
# ============================================================================
message("\n[Stage 2] 预后模型构建 (V2)...")
# 使用 V2 版本的预后模型脚本
source("scripts_final/02c_prognostic_model_v2.R") 
source("scripts_final/02d_nomogram_calibration_dca.R")
source("scripts_final/02e_external_validation.R")

# ============================================================================
# Stage 3: 免疫微环境分析
# ============================================================================
message("\n[Stage 3] 免疫微环境分析...")
# 使用 V2 或正确的脚本名称
if(file.exists("scripts_final/03a_immune_infiltration_v2.R")) {
  source("scripts_final/03a_immune_infiltration_v2.R")
} else {
  source("scripts_final/03a_immune_infiltration_enhanced.R")
}
source("scripts_final/03b_immune_checkpoint.R")

# ============================================================================
# Stage 4: 绵茵陈网络药理学与分子对接
# ============================================================================
message("\n[Stage 4] 网络药理与分子对接...")
source("scripts_final/05_network_MianYinChen_CHB.R")

# 药物敏感性 (V2)
if(file.exists("scripts_final/05c_drug_sensitivity_v2.R")) {
  source("scripts_final/05c_drug_sensitivity_v2.R")
} else {
  source("scripts_final/05c_drug_sensitivity.R")
}

# 分子对接热图绘制
if(file.exists("scripts_final/06_molecular_docking_heatmap.R")) {
  source("scripts_final/06_molecular_docking_heatmap.R")
}

# ============================================================================
# Stage 5: 最终出版级图表生成
# ============================================================================
message("\n[Stage 5] 生成最终出版级图表 (Figures 2-6)...")
source("scripts_final/07_final_publication_figures.R")

message("\n================================================================================")
message("✅ 分析管道全部完成！")
message("请查看 'plots/publication/' 目录获取最终图表。")
message("================================================================================")