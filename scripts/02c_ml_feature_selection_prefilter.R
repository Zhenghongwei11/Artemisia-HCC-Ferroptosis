#!/usr/bin/env Rscript

# 02c_ml_feature_selection_prefilter.R
# 机器学习特征筛选 - 预筛选阶段
# 使用单变量分析将704个候选基因筛选到50个
# 输入: wgcna_module_genes_CHB.csv, deg_GSE83148_sig.csv
# 输出: ml_candidate_genes_top50.csv

suppressPackageStartupMessages({
  library(tidyverse)
})

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir <- "results"

message("[ML预筛选] 开始...")

# 1. 加载数据
message("[ML预筛选] 加载数据...")
expr_83 <- readRDS(file.path(proc_dir, "GSE83148_exprs_log2.rds"))
pdata_83 <- readRDS(file.path(proc_dir, "GSE83148_pdata.rds"))

# 读取WGCNA模块基因 (探针ID格式)
wgcna_df <- read.csv(file.path(res_dir, "wgcna_module_genes_CHB.csv"))
wgcna_probes <- wgcna_df[, 1]  # 第一列是探针ID

# 读取DEG结果 (包含ProbeID和Gene列)
deg_83_sig <- read.csv(file.path(res_dir, "deg_GSE83148_sig.csv"))

# 如果WGCNA为空或交集为空，直接使用DEG的探针ID
if (length(wgcna_probes) == 0 || all(is.na(wgcna_probes))) {
  message("[ML预筛选] WGCNA模块基因为空，使用全部DEG基因")
  wgcna_probes <- deg_83_sig$ProbeID
}

# 2. 构造分组
group_83 <- ifelse(grepl("control|normal|healthy", pdata_83$title, ignore.case = TRUE), "Control", "CHB")
group_83 <- factor(group_83, levels = c("Control", "CHB"))

message("[ML预筛选] 样本分组: Control=", sum(group_83=="Control"), ", CHB=", sum(group_83=="CHB"))

# 3. 获取候选基因 (使用探针ID匹配)
candidate_probes <- intersect(deg_83_sig$ProbeID, wgcna_probes)

# 如果交集为空，使用DEG的全部探针
if (length(candidate_probes) == 0) {
  message("[ML预筛选] WGCNA与DEG交集为空，使用全部DEG探针")
  candidate_probes <- deg_83_sig$ProbeID
}

# 过滤只保留表达矩阵中存在的探针
candidate_probes <- intersect(candidate_probes, rownames(expr_83))
message("[ML预筛选] 候选基因数: ", length(candidate_probes))

# 4. 单变量分析：对每个探针进行t检验
message("[ML预筛选] 进行单变量t检验...")
expr_sub <- expr_83[candidate_probes, ]

pvalues <- apply(expr_sub, 1, function(x) {
  tryCatch({
    t.test(x[group_83=="CHB"], x[group_83=="Control"])$p.value
  }, error = function(e) {
    return(1)
  })
})

# 5. 按p值排序，选择top 50
message("[ML预筛选] 选择p值最小的50个基因...")
n_select <- min(50, length(pvalues))
top_probes <- names(sort(pvalues))[1:n_select]

message("[ML预筛选] 选中基因数: ", length(top_probes))
if (length(top_probes) > 0) {
  message("[ML预筛选] p值范围: ", min(pvalues[top_probes]), " - ", max(pvalues[top_probes]))
}

# 6. 保存结果 (使用探针ID，因为主分析脚本需要探针ID)
result_df <- data.frame(
  Gene = top_probes,  # 这里实际是探针ID
  pvalue = pvalues[top_probes],
  row.names = NULL
)

write.csv(result_df, file.path(res_dir, "ml_candidate_genes_top50.csv"), row.names = FALSE)

message("[ML预筛选] ✅ 完成！")
message("[ML预筛选] 输出文件: ", file.path(res_dir, "ml_candidate_genes_top50.csv"))

