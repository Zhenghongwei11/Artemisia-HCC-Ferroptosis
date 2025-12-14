#!/usr/bin/env Rscript

# 02c_prognostic_model_v2.R
# 预后模型构建 V2 - 扩大候选基因范围
# 策略: WGCNA关键模块基因 ∩ 铁死亡基因 + 单因素Cox筛选
# 
# 使用方法: Rscript scripts_final/02c_prognostic_model_v2.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(glmnet)
  library(timeROC)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

proc_dir <- "data/processed"
ref_dir  <- "data/references"
res_dir  <- "results"
plot_dir <- "plots"

dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

set.seed(123)

message("[预后模型V2] 开始构建预后模型（扩大候选基因范围）...")

# ============================================
# 1. 加载数据
# ============================================
message("[预后模型V2] 加载数据...")

expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
clinical_14 <- readRDS(file.path(proc_dir, "GSE14520_tumor_clinical.rds"))
ferro_genes <- read.csv(file.path(ref_dir, "ferroptosis_genes_expanded.csv"))$Gene

message("[预后模型V2] 表达矩阵: ", nrow(expr_14), " x ", ncol(expr_14))
message("[预后模型V2] 铁死亡基因: ", length(ferro_genes))

# ============================================
# 2. 探针ID转基因符号
# ============================================
message("[预后模型V2] 转换探针ID到基因符号...")

probe_ids <- rownames(expr_14)
gene_symbols <- mapIds(hgu133plus2.db, keys = probe_ids, 
                       column = "SYMBOL", keytype = "PROBEID", 
                       multiVals = "first")

expr_gene <- expr_14
rownames(expr_gene) <- gene_symbols
expr_gene <- expr_gene[!is.na(rownames(expr_gene)), ]

# 合并重复基因
expr_df <- as.data.frame(expr_gene)
expr_df$gene <- rownames(expr_df)
expr_agg <- aggregate(. ~ gene, data = expr_df, FUN = mean)
rownames(expr_agg) <- expr_agg$gene
expr_agg$gene <- NULL
expr_gene <- as.matrix(expr_agg)

message("[预后模型V2] 转换后基因数: ", nrow(expr_gene))


# ============================================
# 3. 匹配表达数据和临床数据
# ============================================
message("[预后模型V2] 匹配表达数据和临床数据...")

matched_gsm <- intersect(colnames(expr_gene), clinical_14$Affy_GSM)
message("[预后模型V2] 匹配样本数: ", length(matched_gsm))

expr_matched <- expr_gene[, matched_gsm]
clinical_matched <- clinical_14[match(matched_gsm, clinical_14$Affy_GSM), ]

# 构建生存数据
surv_data <- data.frame(
  sample = matched_gsm,
  time = as.numeric(clinical_matched$Survival.months),
  status = as.numeric(clinical_matched$Survival.status),
  stringsAsFactors = FALSE
)
surv_data <- surv_data[complete.cases(surv_data), ]

message("[预后模型V2] 有效样本数: ", nrow(surv_data))
message("[预后模型V2] 事件数: ", sum(surv_data$status))

# ============================================
# 4. 扩大候选基因范围
# ============================================
message("[预后模型V2] 筛选候选基因...")

# 策略1: 所有铁死亡基因（在表达矩阵中存在的）
ferro_in_expr <- intersect(ferro_genes, rownames(expr_gene))
message("[预后模型V2] 表达矩阵中的铁死亡基因: ", length(ferro_in_expr))

# 策略2: 对所有铁死亡基因进行单因素Cox筛选
message("[预后模型V2] 对铁死亡基因进行单因素Cox筛选...")

cox_results <- data.frame(
  Gene = character(),
  HR = numeric(),
  HR_lower = numeric(),
  HR_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (gene in ferro_in_expr) {
  gene_expr <- as.numeric(expr_matched[gene, surv_data$sample])
  
  if (sd(gene_expr, na.rm = TRUE) > 0.1) {
    tryCatch({
      cox_fit <- coxph(Surv(time, status) ~ gene_expr, 
                       data = data.frame(time = surv_data$time, 
                                        status = surv_data$status,
                                        gene_expr = gene_expr))
      
      summary_cox <- summary(cox_fit)
      cox_results <- rbind(cox_results, data.frame(
        Gene = gene,
        HR = summary_cox$conf.int[1, 1],
        HR_lower = summary_cox$conf.int[1, 3],
        HR_upper = summary_cox$conf.int[1, 4],
        P_value = summary_cox$coefficients[1, 5]
      ))
    }, error = function(e) {})
  }
}

cox_results <- cox_results %>% arrange(P_value)
write.csv(cox_results, file.path(res_dir, "cox_univariate_ferroptosis_v2.csv"), row.names = FALSE)

sig_genes <- cox_results %>% filter(P_value < 0.1) %>% pull(Gene)
message("[预后模型V2] 单因素Cox显著基因(p<0.1): ", length(sig_genes))

if (length(sig_genes) < 3) {
  # 如果显著基因太少，取top 10
  sig_genes <- head(cox_results$Gene, 10)
  message("[预后模型V2] 使用top 10基因")
}

print(head(cox_results, 15))


# ============================================
# 5. LASSO-Cox回归构建模型
# ============================================
message("[预后模型V2] 进行LASSO-Cox回归...")

# 准备数据
expr_sub <- t(expr_matched[sig_genes, surv_data$sample])
expr_sub <- as.data.frame(expr_sub)

model_data <- cbind(surv_data, expr_sub)
model_data <- model_data[complete.cases(model_data), ]

message("[预后模型V2] 模型数据样本数: ", nrow(model_data))

# LASSO-Cox
x <- as.matrix(model_data[, sig_genes])
y <- Surv(model_data$time, model_data$status)

cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
best_lambda <- cv_fit$lambda.min

lasso_coef <- coef(cv_fit, s = best_lambda)
selected_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]

message("[预后模型V2] LASSO选择的基因: ", length(selected_genes))
print(selected_genes)

# 如果LASSO没选到基因，使用单因素Cox top基因
if (length(selected_genes) == 0) {
  selected_genes <- head(sig_genes, 5)
  message("[预后模型V2] LASSO未选到基因，使用top 5")
}

# ============================================
# 6. 计算Risk Score
# ============================================
message("[预后模型V2] 计算Risk Score...")

# 多因素Cox回归获取系数
if (length(selected_genes) > 1) {
  formula_str <- paste("Surv(time, status) ~", paste(paste0("`", selected_genes, "`"), collapse = " + "))
  multi_cox <- coxph(as.formula(formula_str), data = model_data)
  coef_df <- data.frame(
    Gene = names(coef(multi_cox)),
    Coefficient = as.numeric(coef(multi_cox))
  )
} else {
  # 单基因
  coef_df <- data.frame(
    Gene = selected_genes,
    Coefficient = as.numeric(lasso_coef[selected_genes, ])
  )
}

write.csv(coef_df, file.path(res_dir, "prognostic_model_coef_v2.csv"), row.names = FALSE)

# 计算Risk Score
risk_score <- rep(0, nrow(model_data))
for (i in 1:nrow(coef_df)) {
  gene <- coef_df$Gene[i]
  coef <- coef_df$Coefficient[i]
  gene_clean <- gsub("`", "", gene)
  if (gene_clean %in% colnames(model_data)) {
    risk_score <- risk_score + coef * model_data[[gene_clean]]
  }
}

model_data$risk_score <- risk_score
model_data$risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
model_data$risk_group <- factor(model_data$risk_group, levels = c("Low", "High"))

message("[预后模型V2] High Risk: ", sum(model_data$risk_group == "High"))
message("[预后模型V2] Low Risk: ", sum(model_data$risk_group == "Low"))

write.csv(model_data[, c("sample", "time", "status", "risk_score", "risk_group")],
          file.path(res_dir, "risk_score_data_v2.csv"), row.names = FALSE)


# ============================================
# 7. 生存分析
# ============================================
message("[预后模型V2] 进行生存分析...")

# KM曲线
surv_obj <- Surv(model_data$time, model_data$status)
km_fit <- survfit(surv_obj ~ risk_group, data = model_data)

# Log-rank检验
logrank_test <- survdiff(surv_obj ~ risk_group, data = model_data)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = 1)

message("[预后模型V2] Log-rank p值: ", format(logrank_p, digits = 4))

# 时间依赖ROC
roc_res <- timeROC(
  T = model_data$time,
  delta = model_data$status,
  marker = model_data$risk_score,
  cause = 1,
  times = c(12, 36, 60),
  iid = TRUE
)

auc_1y <- roc_res$AUC[1]
auc_3y <- roc_res$AUC[2]
auc_5y <- roc_res$AUC[3]

message("[预后模型V2] 1年AUC: ", round(auc_1y, 3))
message("[预后模型V2] 3年AUC: ", round(auc_3y, 3))
message("[预后模型V2] 5年AUC: ", round(auc_5y, 3))

# C-index
cox_model <- coxph(Surv(time, status) ~ risk_score, data = model_data)
c_index <- summary(cox_model)$concordance[1]
message("[预后模型V2] C-index: ", round(c_index, 3))

# ============================================
# 8. 生成图表
# ============================================
message("[预后模型V2] 生成图表...")

pdf(file.path(plot_dir, "Figure5_prognostic_model_v2.pdf"), width = 14, height = 10)
par(mfrow = c(2, 2))

# A: KM曲线
km_plot <- ggsurvplot(
  km_fit,
  data = model_data,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#2E9FDF", "#E7B800"),
  title = paste0("A. Kaplan-Meier Survival Curve\n(p = ", format(logrank_p, digits = 3), ")"),
  xlab = "Time (months)",
  ylab = "Survival Probability"
)
print(km_plot)

# B: ROC曲线
plot(roc_res, time = 12, col = "red", lwd = 2, 
     main = "B. Time-dependent ROC Curves")
plot(roc_res, time = 36, col = "blue", lwd = 2, add = TRUE)
plot(roc_res, time = 60, col = "green", lwd = 2, add = TRUE)
legend("bottomright", 
       legend = c(paste0("1-year AUC = ", round(auc_1y, 3)),
                  paste0("3-year AUC = ", round(auc_3y, 3)),
                  paste0("5-year AUC = ", round(auc_5y, 3))),
       col = c("red", "blue", "green"), lwd = 2)

# C: Risk Score分布
hist(model_data$risk_score, breaks = 30, col = "steelblue",
     main = "C. Risk Score Distribution",
     xlab = "Risk Score", ylab = "Frequency")
abline(v = median(model_data$risk_score), col = "red", lwd = 2, lty = 2)

# D: 森林图
if (nrow(coef_df) > 0) {
  barplot(coef_df$Coefficient, names.arg = coef_df$Gene,
          col = ifelse(coef_df$Coefficient > 0, "red", "blue"),
          main = "D. Model Coefficients",
          ylab = "Coefficient", las = 2, cex.names = 0.8)
  abline(h = 0, lty = 2)
}

dev.off()

# 保存统计结果
stats_df <- data.frame(
  Metric = c("Log-rank p-value", "1-year AUC", "3-year AUC", "5-year AUC", "C-index",
             "High Risk N", "Low Risk N", "Total Events", "Selected Genes"),
  Value = c(format(logrank_p, digits = 4), round(auc_1y, 3), round(auc_3y, 3), 
            round(auc_5y, 3), round(c_index, 3),
            sum(model_data$risk_group == "High"), sum(model_data$risk_group == "Low"),
            sum(model_data$status), paste(selected_genes, collapse = ", "))
)
write.csv(stats_df, file.path(res_dir, "prognostic_model_stats_v2.csv"), row.names = FALSE)

message("[预后模型V2] ✅ 预后模型V2构建完成！")
message("[预后模型V2] 输出文件:")
message("  ✅ results/cox_univariate_ferroptosis_v2.csv")
message("  ✅ results/prognostic_model_coef_v2.csv")
message("  ✅ results/risk_score_data_v2.csv")
message("  ✅ results/prognostic_model_stats_v2.csv")
message("  ✅ plots/Figure5_prognostic_model_v2.pdf")
