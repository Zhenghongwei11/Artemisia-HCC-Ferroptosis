#!/usr/bin/env Rscript

# 05c_drug_sensitivity_v2.R
# 药物敏感性分析 - 基于Ridge回归的自定义实现
# 参考oncoPredict/pRRophetic方法
# 使用方法: Rscript scripts_final/05c_drug_sensitivity_v2.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(sva)
  library(glmnet)
  library(pheatmap)
})

proc_dir <- "data/processed"
ref_dir  <- "data/references"
res_dir  <- "results"
plot_dir <- "plots"

set.seed(123)

message("[药物敏感性V2] 开始药物敏感性分析...")

# ============================================
# 1. 加载数据
# ============================================
message("[药物敏感性V2] 加载数据...")

# 表达矩阵
expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr_symbol.rds"))

# GDSC数据
GDSC2_Expr <- readRDS(file.path(ref_dir, "GDSC/GDSC2_Expr.rds"))
GDSC2_Res <- readRDS(file.path(ref_dir, "GDSC/GDSC2_Res.rds"))

# Risk Score数据
risk_data <- read.csv(file.path(res_dir, "risk_score_data.csv"))

message("[药物敏感性V2] GSE14520: ", nrow(expr_14), " x ", ncol(expr_14))
message("[药物敏感性V2] GDSC2表达: ", nrow(GDSC2_Expr), " x ", ncol(GDSC2_Expr))
message("[药物敏感性V2] GDSC2药物: ", ncol(GDSC2_Res))

# ============================================
# 2. 数据预处理和批次校正
# ============================================
message("[药物敏感性V2] 数据预处理...")

# 统一基因名为大写
test_expr <- as.matrix(expr_14)
rownames(test_expr) <- toupper(rownames(test_expr))
rownames(GDSC2_Expr) <- toupper(rownames(GDSC2_Expr))

# 找到共同基因
commonGenes <- intersect(rownames(GDSC2_Expr), rownames(test_expr))
message("[药物敏感性V2] 共同基因: ", length(commonGenes))

# 子集化
trainExpr <- GDSC2_Expr[commonGenes, ]
testExpr <- test_expr[commonGenes, ]

# 合并数据进行批次校正
combined <- cbind(trainExpr, testExpr)
batch <- c(rep(1, ncol(trainExpr)), rep(2, ncol(testExpr)))

# 移除低变异基因
gene_var <- apply(combined, 1, var, na.rm = TRUE)
keep_genes <- gene_var > quantile(gene_var, 0.2, na.rm = TRUE)
combined_filtered <- combined[keep_genes, ]

message("[药物敏感性V2] 过滤后基因数: ", nrow(combined_filtered))

# ComBat批次校正
message("[药物敏感性V2] 运行ComBat批次校正...")
corrected <- ComBat(dat = combined_filtered, batch = batch, par.prior = TRUE)

# 分离训练集和测试集
trainExpr_corrected <- corrected[, 1:ncol(trainExpr)]
testExpr_corrected <- corrected[, (ncol(trainExpr)+1):ncol(corrected)]

message("[药物敏感性V2] 批次校正完成")


# ============================================
# 3. 使用Ridge回归预测IC50
# ============================================
message("[药物敏感性V2] 使用Ridge回归预测IC50...")

# 匹配训练集样本
train_samples <- intersect(colnames(trainExpr_corrected), rownames(GDSC2_Res))
message("[药物敏感性V2] 训练样本数: ", length(train_samples))

trainExpr_final <- t(trainExpr_corrected[, train_samples])
trainIC50 <- GDSC2_Res[train_samples, ]

# 选择关键药物进行预测（肝癌相关+常用化疗药）
key_drug_keywords <- c("Sorafenib", "Regorafenib", "Cabozantinib", "Lenvatinib",
                       "Doxorubicin", "Cisplatin", "Fluorouracil", "Gemcitabine",
                       "Oxaliplatin", "Erlotinib", "Gefitinib", "Lapatinib",
                       "Sunitinib", "Pazopanib", "Axitinib", "Imatinib",
                       "Paclitaxel", "Docetaxel", "Etoposide", "Vincristine",
                       "Methotrexate", "Temozolomide", "Rapamycin", "Everolimus",
                       "Bortezomib", "Vorinostat", "Panobinostat", "Olaparib",
                       "Vemurafenib", "Dabrafenib", "Trametinib", "Selumetinib",
                       "Crizotinib", "Ceritinib", "Afatinib", "Osimertinib",
                       "Ibrutinib", "Venetoclax", "Navitoclax", "Nutlin",
                       "JQ1", "I-BET", "AZD", "GSK", "PD-", "MK-", "BMS-")

key_drugs_idx <- grep(paste(key_drug_keywords, collapse = "|"), 
                      colnames(trainIC50), ignore.case = TRUE)

# 如果匹配太少，选择数据最完整的前50个药物
if (length(key_drugs_idx) < 30) {
  drug_completeness <- colSums(!is.na(trainIC50))
  key_drugs_idx <- order(drug_completeness, decreasing = TRUE)[1:min(80, ncol(trainIC50))]
}

message("[药物敏感性V2] 选择 ", length(key_drugs_idx), " 个关键药物进行预测")

# 预测选定药物的IC50
n_drugs <- length(key_drugs_idx)
test_samples <- colnames(testExpr_corrected)
predicted_IC50 <- matrix(NA, nrow = length(test_samples), ncol = n_drugs)
rownames(predicted_IC50) <- test_samples
colnames(predicted_IC50) <- colnames(trainIC50)[key_drugs_idx]

testExpr_final <- t(testExpr_corrected)

message("[药物敏感性V2] 开始预测 ", n_drugs, " 个药物的IC50...")

pb <- txtProgressBar(min = 0, max = n_drugs, style = 3)

for (i in 1:n_drugs) {
  drug_idx <- key_drugs_idx[i]
  drug_name <- colnames(trainIC50)[drug_idx]
  y <- trainIC50[, drug_idx]
  
  # 移除NA样本
  valid_idx <- !is.na(y)
  if (sum(valid_idx) < 30) {
    setTxtProgressBar(pb, i)
    next
  }
  
  y_valid <- y[valid_idx]
  X_train <- trainExpr_final[valid_idx, ]
  
  # Ridge回归 (alpha = 0)，减少交叉验证折数加速
  tryCatch({
    cv_fit <- cv.glmnet(X_train, y_valid, alpha = 0, nfolds = 3)
    
    # 预测测试集
    predicted_IC50[, i] <- predict(cv_fit, newx = testExpr_final, s = "lambda.min")
  }, error = function(e) {
    # 跳过失败的药物
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

# 移除全NA的药物
valid_drugs <- colSums(!is.na(predicted_IC50)) > 0
predicted_IC50 <- predicted_IC50[, valid_drugs]

message("[药物敏感性V2] 成功预测 ", ncol(predicted_IC50), " 个药物")

# 保存预测结果
write.csv(predicted_IC50, file.path(res_dir, "predicted_IC50_GDSC2.csv"))
saveRDS(predicted_IC50, file.path(res_dir, "predicted_IC50_GDSC2.rds"))

# ============================================
# 4. 分析IC50与Risk Score的关系
# ============================================
message("[药物敏感性V2] 分析IC50与Risk Score的关系...")

common_samples <- intersect(risk_data$sample, rownames(predicted_IC50))
message("[药物敏感性V2] 匹配样本数: ", length(common_samples))

ic50_matched <- predicted_IC50[common_samples, ]
risk_matched <- risk_data[match(common_samples, risk_data$sample), ]

# 计算每个药物IC50与Risk Score的相关性
cor_results <- data.frame(
  Drug = colnames(ic50_matched),
  Correlation = numeric(ncol(ic50_matched)),
  P_value = numeric(ncol(ic50_matched)),
  stringsAsFactors = FALSE
)

for (i in 1:ncol(ic50_matched)) {
  ic50_vals <- ic50_matched[, i]
  valid_idx <- !is.na(ic50_vals)
  
  if (sum(valid_idx) > 20) {
    cor_test <- cor.test(ic50_vals[valid_idx], risk_matched$risk_score[valid_idx], 
                         method = "spearman")
    cor_results$Correlation[i] <- cor_test$estimate
    cor_results$P_value[i] <- cor_test$p.value
  }
}

cor_results <- cor_results %>% 
  filter(!is.na(Correlation)) %>%
  arrange(P_value) %>%
  mutate(FDR = p.adjust(P_value, method = "BH"))

write.csv(cor_results, file.path(res_dir, "drug_risk_correlation.csv"), row.names = FALSE)

sig_drugs <- sum(cor_results$P_value < 0.05)
message("[药物敏感性V2] 与Risk Score显著相关的药物: ", sig_drugs)


# ============================================
# 5. 高低风险组IC50差异
# ============================================
message("[药物敏感性V2] 分析高低风险组IC50差异...")

diff_results <- data.frame(
  Drug = colnames(ic50_matched),
  High_Mean = numeric(ncol(ic50_matched)),
  Low_Mean = numeric(ncol(ic50_matched)),
  Diff = numeric(ncol(ic50_matched)),
  P_value = numeric(ncol(ic50_matched)),
  stringsAsFactors = FALSE
)

for (i in 1:ncol(ic50_matched)) {
  high_vals <- ic50_matched[risk_matched$risk_group == "High", i]
  low_vals <- ic50_matched[risk_matched$risk_group == "Low", i]
  
  high_vals <- high_vals[!is.na(high_vals)]
  low_vals <- low_vals[!is.na(low_vals)]
  
  diff_results$High_Mean[i] <- mean(high_vals, na.rm = TRUE)
  diff_results$Low_Mean[i] <- mean(low_vals, na.rm = TRUE)
  diff_results$Diff[i] <- diff_results$High_Mean[i] - diff_results$Low_Mean[i]
  
  if (length(high_vals) > 10 && length(low_vals) > 10) {
    diff_results$P_value[i] <- wilcox.test(high_vals, low_vals)$p.value
  }
}

diff_results <- diff_results %>% 
  filter(!is.na(P_value)) %>%
  arrange(P_value) %>%
  mutate(FDR = p.adjust(P_value, method = "BH"))

write.csv(diff_results, file.path(res_dir, "drug_risk_group_diff.csv"), row.names = FALSE)

# 高风险组更敏感的药物 (IC50更低)
high_sensitive <- diff_results %>% 
  filter(P_value < 0.05 & Diff < 0) %>%
  arrange(Diff)

# 低风险组更敏感的药物
low_sensitive <- diff_results %>% 
  filter(P_value < 0.05 & Diff > 0) %>%
  arrange(desc(Diff))

message("[药物敏感性V2] 高风险组更敏感的药物: ", nrow(high_sensitive))
message("[药物敏感性V2] 低风险组更敏感的药物: ", nrow(low_sensitive))

# 保存推荐药物
if (nrow(high_sensitive) > 0) {
  write.csv(high_sensitive, file.path(res_dir, "recommended_drugs_high_risk.csv"), row.names = FALSE)
}
if (nrow(low_sensitive) > 0) {
  write.csv(low_sensitive, file.path(res_dir, "recommended_drugs_low_risk.csv"), row.names = FALSE)
}

# ============================================
# 6. 肝癌特异性药物分析
# ============================================
message("[药物敏感性V2] 分析肝癌特异性药物...")

hcc_drug_keywords <- c("Sorafenib", "Regorafenib", "Cabozantinib", "Doxorubicin", "Cisplatin", "Fluorouracil", "Gemcitabine", "Oxaliplatin", "Erlotinib", "Gefitinib", "Lapatinib", "Sunitinib", "Pazopanib", "Axitinib")

hcc_drugs <- grep(paste(hcc_drug_keywords, collapse = "|"), 
                  colnames(predicted_IC50), ignore.case = TRUE, value = TRUE)

if (length(hcc_drugs) > 0) {
  message("[药物敏感性V2] 找到肝癌相关药物: ", length(hcc_drugs))
  
  hcc_ic50 <- predicted_IC50[, hcc_drugs, drop = FALSE]
  write.csv(hcc_ic50, file.path(res_dir, "hcc_drug_ic50.csv"))
  
  # 肝癌药物在高低风险组的差异
  hcc_diff <- diff_results %>% filter(Drug %in% hcc_drugs)
  if (nrow(hcc_diff) > 0) {
    message("[药物敏感性V2] 肝癌药物高低风险组差异:")
    print(hcc_diff[, c("Drug", "High_Mean", "Low_Mean", "Diff", "P_value")])
    write.csv(hcc_diff, file.path(res_dir, "hcc_drug_risk_diff.csv"), row.names = FALSE)
  }
}

# ============================================
# 7. 生成Figure 10: 药物敏感性多面板图
# ============================================
message("[药物敏感性V2] 生成Figure 10...")

pdf(file.path(plot_dir, "Figure10_drug_sensitivity.pdf"), width = 16, height = 14)
par(mfrow = c(2, 2))

# A: 药物-Risk Score相关性Top20
top_cor <- head(cor_results, 20)
barplot(top_cor$Correlation, names.arg = gsub("_.*", "", top_cor$Drug),
        col = ifelse(top_cor$Correlation > 0, "red", "blue"),
        las = 2, cex.names = 0.6,
        main = "A. Drug IC50 - Risk Score Correlation (Top 20)",
        ylab = "Spearman Correlation")
abline(h = 0, lty = 2)
legend("topright", legend = c("Positive (resistant)", "Negative (sensitive)"),
       fill = c("red", "blue"), cex = 0.8)

# B: 高低风险组IC50差异Top20
top_diff <- head(diff_results, 20)
diff_mat <- as.matrix(top_diff[, c("High_Mean", "Low_Mean")])
rownames(diff_mat) <- gsub("_.*", "", top_diff$Drug)
barplot(t(diff_mat), beside = TRUE,
        col = c("red", "blue"), las = 2, cex.names = 0.6,
        main = "B. IC50 Difference: High vs Low Risk (Top 20)",
        ylab = "Predicted IC50")
legend("topright", legend = c("High Risk", "Low Risk"),
       fill = c("red", "blue"), cex = 0.8)

# C: 肝癌药物IC50分布
if (length(hcc_drugs) > 0 && exists("hcc_diff") && nrow(hcc_diff) > 0) {
  hcc_mat <- as.matrix(hcc_diff[, c("High_Mean", "Low_Mean")])
  rownames(hcc_mat) <- gsub("_.*", "", hcc_diff$Drug)
  barplot(t(hcc_mat), beside = TRUE,
          col = c("red", "blue"), las = 2, cex.names = 0.7,
          main = "C. HCC Drugs: High vs Low Risk",
          ylab = "Predicted IC50")
  legend("topright", legend = c("High Risk", "Low Risk"),
         fill = c("red", "blue"), cex = 0.8)
} else {
  plot.new()
  text(0.5, 0.5, "No HCC drugs found", cex = 1.5)
}

# D: 高风险患者推荐药物
if (nrow(high_sensitive) > 0) {
  top_rec <- head(high_sensitive, 15)
  barplot(-top_rec$Diff, names.arg = gsub("_.*", "", top_rec$Drug),
          col = "darkgreen", las = 2, cex.names = 0.6,
          main = "D. Recommended Drugs for High Risk Patients",
          ylab = "IC50 Difference (Low - High)")
} else {
  plot.new()
  text(0.5, 0.5, "No significant drugs for high risk", cex = 1.2)
}

dev.off()

message("[药物敏感性V2] ✅ 药物敏感性分析完成！")
message("[药物敏感性V2] 输出文件:")
message("  ✅ ", file.path(res_dir, "predicted_IC50_GDSC2.csv"))
message("  ✅ ", file.path(res_dir, "drug_risk_correlation.csv"))
message("  ✅ ", file.path(res_dir, "drug_risk_group_diff.csv"))
if (nrow(high_sensitive) > 0) message("  ✅ ", file.path(res_dir, "recommended_drugs_high_risk.csv"))
if (length(hcc_drugs) > 0) message("  ✅ ", file.path(res_dir, "hcc_drug_ic50.csv"))
message("  ✅ ", file.path(plot_dir, "Figure10_drug_sensitivity.pdf"))
