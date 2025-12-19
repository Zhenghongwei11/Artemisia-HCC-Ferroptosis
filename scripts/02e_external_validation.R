#!/usr/bin/env Rscript

# 02e_external_validation.R
# 外部验证 - 使用TCGA-LIHC队列验证预后模型
# 使用方法: Rscript scripts/02e_external_validation.R
# 
# 输入: results/prognostic_model_coef_v2.csv
# 输出: TCGA验证结果

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"

set.seed(123)

message("[外部验证] 开始TCGA-LIHC外部验证...")

# ============================================
# 1. 加载模型系数
# ============================================
message("[外部验证] 加载模型系数...")

coef_file <- file.path(res_dir, "prognostic_model_coef_v2.csv")
if (!file.exists(coef_file)) {
  message("[外部验证] 未找到模型系数文件: ", coef_file)
  message("[外部验证] 请先运行 scripts/02c_prognostic_model_v2.R 生成该文件，或使用本包提供的结果表。")
  quit(save = "no", status = 1)
}

model_coef <- read.csv(coef_file)
message("[外部验证] 模型基因数: ", nrow(model_coef))

# ============================================
# 2. 下载TCGA-LIHC数据
# ============================================
message("[外部验证] 检查TCGA-LIHC数据...")

tcga_expr_file <- file.path(proc_dir, "TCGA_LIHC_expr.rds")
tcga_clin_file <- file.path(proc_dir, "TCGA_LIHC_clinical.rds")

if (!file.exists(tcga_expr_file) || !file.exists(tcga_clin_file)) {
  message("[外部验证] 下载TCGA-LIHC数据...")
  
  tryCatch({
    # 查询TCGA-LIHC
    query <- GDCquery(
      project = "TCGA-LIHC",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    # 下载数据
    GDCdownload(query, method = "api", files.per.chunk = 50)
    
    # 准备数据
    data <- GDCprepare(query)
    
    # 提取表达矩阵（不同版本的TCGAbiolinks可能使用不同assay名称）
    assays <- assayNames(data)
    preferred <- c("tpm_unstrand", "tpm_unstranded", "HTSeq - Counts")
    use_assay <- intersect(preferred, assays)
    if (length(use_assay) == 0) use_assay <- assays[1]
    expr_mat <- assay(data, use_assay[1])
    
    # 转换基因ID为Symbol
    gene_info <- rowData(data)
    if ("gene_name" %in% colnames(gene_info)) {
      rownames(expr_mat) <- gene_info$gene_name
    } else if ("external_gene_name" %in% colnames(gene_info)) {
      rownames(expr_mat) <- gene_info$external_gene_name
    } else if ("gene_id" %in% colnames(gene_info)) {
      rownames(expr_mat) <- gene_info$gene_id
    }
    
    # 去重
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
    
    # 保存
    saveRDS(expr_mat, tcga_expr_file)
    
    # 获取临床数据
    clinical <- GDCquery_clinic("TCGA-LIHC", type = "clinical")
    saveRDS(clinical, tcga_clin_file)
    
    message("[外部验证] TCGA数据下载完成")
    
  }, error = function(e) {
    message("[外部验证] TCGA下载失败: ", e$message)
    message("[外部验证] 如果希望执行该步骤，请手动准备并放置以下文件：")
    message("  - ", tcga_expr_file)
    message("  - ", tcga_clin_file)
    message("[外部验证] 详见 docs/DATA_MANIFEST.md (TCGA-LIHC 外部验证部分)。")
    quit(save = "no", status = 1)
  })
}

# 加载数据
tcga_expr <- readRDS(tcga_expr_file)
tcga_clin <- readRDS(tcga_clin_file)

message("[外部验证] TCGA-LIHC样本数: ", ncol(tcga_expr))
message("[外部验证] TCGA-LIHC基因数: ", nrow(tcga_expr))

# ============================================
# 3. 准备验证数据
# ============================================
message("[外部验证] 准备验证数据...")

# 筛选肿瘤样本
tumor_samples <- colnames(tcga_expr)[grepl("-01A-", colnames(tcga_expr))]
tcga_expr_tumor <- tcga_expr[, tumor_samples]

message("[外部验证] 肿瘤样本数: ", length(tumor_samples))

# 匹配模型基因
available_genes <- intersect(model_coef$Gene, rownames(tcga_expr_tumor))
message("[外部验证] 可用模型基因: ", length(available_genes), "/", nrow(model_coef))

if (length(available_genes) < 2) {
  message("[外部验证] 可用基因过少，无法进行验证")
  quit(save = "no", status = 1)
}

# 提取表达矩阵
expr_sub <- t(tcga_expr_tumor[available_genes, , drop = FALSE])
expr_sub <- log2(expr_sub + 1)  # log2转换

# ============================================
# 4. 计算Risk Score
# ============================================
message("[外部验证] 计算Risk Score...")

# 使用模型系数计算Risk Score
coef_matched <- model_coef %>% filter(Gene %in% available_genes)

risk_scores <- as.matrix(expr_sub[, coef_matched$Gene]) %*% coef_matched$Coefficient
risk_scores <- as.numeric(risk_scores)

# 创建验证数据框
valid_data <- data.frame(
  sample = rownames(expr_sub),
  risk_score = risk_scores,
  stringsAsFactors = FALSE
)

# 提取样本ID用于匹配临床数据
valid_data$patient_id <- substr(valid_data$sample, 1, 12)

# ============================================
# 5. 匹配生存数据
# ============================================
message("[外部验证] 匹配生存数据...")

# TCGA临床数据字段
if ("days_to_death" %in% colnames(tcga_clin) && "days_to_last_follow_up" %in% colnames(tcga_clin)) {
  tcga_clin$time <- ifelse(!is.na(tcga_clin$days_to_death), 
                           tcga_clin$days_to_death,
                           tcga_clin$days_to_last_follow_up)
  tcga_clin$status <- ifelse(!is.na(tcga_clin$days_to_death), 1, 0)
}

# 匹配
valid_data <- valid_data %>%
  left_join(tcga_clin[, c("submitter_id", "time", "status")], 
            by = c("patient_id" = "submitter_id"))

# 移除缺失值
valid_data <- valid_data[complete.cases(valid_data[, c("time", "status", "risk_score")]), ]
valid_data <- valid_data[valid_data$time > 0, ]

message("[外部验证] 有效样本数: ", nrow(valid_data))

# Risk分群
median_risk <- median(valid_data$risk_score)
valid_data$risk_group <- ifelse(valid_data$risk_score > median_risk, "High", "Low")
valid_data$risk_group <- factor(valid_data$risk_group, levels = c("Low", "High"))

# ============================================
# 6. 验证集生存分析
# ============================================
message("[外部验证] 进行验证集生存分析...")

# KM曲线
surv_obj <- Surv(valid_data$time, valid_data$status)
km_fit <- survfit(surv_obj ~ risk_group, data = valid_data)

# Log-rank检验
logrank_test <- survdiff(surv_obj ~ risk_group, data = valid_data)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = 1)

message("[外部验证] 验证集Log-rank p值: ", format(logrank_p, digits = 4))

# C-index
cox_valid <- coxph(Surv(time, status) ~ risk_score, data = valid_data)
c_index_valid <- summary(cox_valid)$concordance[1]

message("[外部验证] 验证集C-index: ", round(c_index_valid, 3))

# ============================================
# 7. 生成验证图
# ============================================
message("[外部验证] 生成验证图...")

pdf(file.path(plot_dir, "Figure6_external_validation.pdf"), width = 12, height = 10)

par(mfrow = c(2, 2))

# A: 验证集KM曲线
plot(km_fit, col = c("#2E9FDF", "#E7B800"), lwd = 2,
     xlab = "Time (days)", ylab = "Survival Probability",
     main = paste0("A. TCGA-LIHC Validation (p = ", format(logrank_p, digits = 3), ")"))
legend("topright", legend = c("Low Risk", "High Risk"), 
       col = c("#2E9FDF", "#E7B800"), lwd = 2)

# B: Risk Score分布比较
# 加载训练集数据
train_risk_file <- file.path(res_dir, "risk_score_data_v2.csv")
if (!file.exists(train_risk_file)) {
  message("[外部验证] 缺少训练集 risk score 文件: ", train_risk_file)
  quit(save = "no", status = 1)
}
train_risk <- read.csv(train_risk_file)

boxplot(list(Training = train_risk$risk_score, Validation = valid_data$risk_score),
        col = c("steelblue", "coral"),
        main = "B. Risk Score Distribution Comparison",
        ylab = "Risk Score")

# C: 验证集Risk Score分布
hist(valid_data$risk_score, breaks = 30, col = "coral", border = "white",
     main = "C. TCGA-LIHC Risk Score Distribution",
     xlab = "Risk Score", ylab = "Frequency")
abline(v = median_risk, col = "red", lwd = 2, lty = 2)

# D: C-index比较
c_index_train <- summary(coxph(Surv(time, status) ~ risk_score, data = train_risk))$concordance[1]
barplot(c(Training = c_index_train, Validation = c_index_valid),
        col = c("steelblue", "coral"), ylim = c(0, 1),
        main = "D. C-index Comparison",
        ylab = "C-index")
abline(h = 0.7, lty = 2, col = "red")

dev.off()

message("[外部验证] 验证图已保存")

# ============================================
# 8. 保存验证结果
# ============================================
message("[外部验证] 保存验证结果...")

valid_stats <- data.frame(
  Dataset = c("Training (GSE14520)", "Validation (TCGA-LIHC)"),
  Sample_N = c(nrow(train_risk), nrow(valid_data)),
  High_Risk_N = c(sum(train_risk$risk_group == "High"), sum(valid_data$risk_group == "High")),
  Low_Risk_N = c(sum(train_risk$risk_group == "Low"), sum(valid_data$risk_group == "Low")),
  Logrank_P = c(NA, logrank_p),
  C_index = c(c_index_train, c_index_valid)
)

write.csv(valid_stats, file.path(res_dir, "external_validation_stats.csv"), row.names = FALSE)
write.csv(valid_data, file.path(res_dir, "TCGA_LIHC_risk_score.csv"), row.names = FALSE)

# ============================================
# 完成
# ============================================
message("[外部验证] ✅ 外部验证完成！")
message("[外部验证] 输出文件:")
message("  ✅ ", file.path(plot_dir, "Figure6_external_validation.pdf"))
message("  ✅ ", file.path(res_dir, "external_validation_stats.csv"))
message("  ✅ ", file.path(res_dir, "TCGA_LIHC_risk_score.csv"))
