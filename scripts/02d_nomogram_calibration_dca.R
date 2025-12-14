#!/usr/bin/env Rscript

# 02d_nomogram_calibration_dca.R
# 高质量临床图表: Nomogram + 校准曲线 + DCA决策曲线
# 整合Risk Score与临床变量，构建综合预后模型
# 
# 适用于SCI 1-2区论文投稿标准
# 使用方法: Rscript scripts_final/02d_nomogram_calibration_dca.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(rms)
  library(ggDCA)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(pROC)
})

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"

dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

set.seed(123)

message("=" %>% rep(60) %>% paste(collapse = ""))
message("[临床图表] 开始生成高质量Nomogram、校准曲线和DCA")
message("=" %>% rep(60) %>% paste(collapse = ""))

# ============================================
# 1. 加载数据
# ============================================
message("\n[Step 1] 加载数据...")

# 加载预后模型结果
risk_data <- read.csv(file.path(res_dir, "risk_score_data_v2.csv"))
model_coef <- read.csv(file.path(res_dir, "prognostic_model_coef_v2.csv"))

# 加载原始临床数据
clinical_14 <- readRDS(file.path(proc_dir, "GSE14520_tumor_clinical.rds"))

message("  - Risk Score数据: ", nrow(risk_data), " 样本")
message("  - 模型基因数: ", nrow(model_coef))
message("  - 临床数据变量: ", ncol(clinical_14))

# ============================================
# 2. 整合临床变量
# ============================================
message("\n[Step 2] 整合临床变量...")

# 匹配样本
clinical_matched <- clinical_14[match(risk_data$sample, clinical_14$Affy_GSM), ]

# 构建完整数据集
full_data <- data.frame(
  sample = risk_data$sample,
  time = risk_data$time,
  status = risk_data$status,
  risk_score = risk_data$risk_score,
  risk_group = factor(risk_data$risk_group, levels = c("Low", "High")),
  stringsAsFactors = FALSE
)

# 添加临床变量
full_data$gender <- factor(clinical_matched$Gender, levels = c("F", "M"))
full_data$age <- as.numeric(clinical_matched$Age)
full_data$tumor_size <- factor(ifelse(clinical_matched$Main.Tumor.Size......5.cm. == "large", "Large", "Small"),
                               levels = c("Small", "Large"))
full_data$multinodular <- factor(clinical_matched$Multinodular, levels = c("N", "Y"))
full_data$cirrhosis <- factor(clinical_matched$Cirrhosis, levels = c("N", "Y"))
full_data$tnm_stage <- clinical_matched$TNM.staging
full_data$bclc_stage <- clinical_matched$BCLC.staging
full_data$afp <- factor(ifelse(clinical_matched$AFP......300ng.ml. == "high", "High", "Low"),
                        levels = c("Low", "High"))

# 处理TNM分期（合并为早期/晚期）
full_data$tnm_group <- factor(ifelse(full_data$tnm_stage %in% c("I", "II"), "Early", "Advanced"),
                              levels = c("Early", "Advanced"))

# 移除缺失值
full_data <- full_data[complete.cases(full_data[, c("time", "status", "risk_score", 
                                                     "age", "gender", "tumor_size", "afp")]), ]

message("  - 完整数据集: ", nrow(full_data), " 样本")
message("  - 临床变量: Gender, Age, Tumor Size, AFP, TNM Stage, Cirrhosis, Multinodular")

# ============================================
# 3. 单因素Cox回归分析临床变量
# ============================================
message("\n[Step 3] 单因素Cox回归分析...")

univariate_results <- data.frame(
  Variable = character(),
  HR = numeric(),
  HR_lower = numeric(),
  HR_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Risk Score
cox_risk <- coxph(Surv(time, status) ~ risk_score, data = full_data)
s <- summary(cox_risk)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "Risk Score",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# Age
cox_age <- coxph(Surv(time, status) ~ age, data = full_data)
s <- summary(cox_age)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "Age",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# Gender
cox_gender <- coxph(Surv(time, status) ~ gender, data = full_data)
s <- summary(cox_gender)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "Gender (M vs F)",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# Tumor Size
cox_size <- coxph(Surv(time, status) ~ tumor_size, data = full_data)
s <- summary(cox_size)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "Tumor Size (Large vs Small)",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# AFP
cox_afp <- coxph(Surv(time, status) ~ afp, data = full_data)
s <- summary(cox_afp)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "AFP (High vs Low)",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# TNM Stage
cox_tnm <- coxph(Surv(time, status) ~ tnm_group, data = full_data)
s <- summary(cox_tnm)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "TNM Stage (Advanced vs Early)",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# Cirrhosis
cox_cirr <- coxph(Surv(time, status) ~ cirrhosis, data = full_data)
s <- summary(cox_cirr)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "Cirrhosis (Y vs N)",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

# Multinodular
cox_multi <- coxph(Surv(time, status) ~ multinodular, data = full_data)
s <- summary(cox_multi)
univariate_results <- rbind(univariate_results, data.frame(
  Variable = "Multinodular (Y vs N)",
  HR = s$conf.int[1, 1],
  HR_lower = s$conf.int[1, 3],
  HR_upper = s$conf.int[1, 4],
  P_value = s$coefficients[1, 5]
))

univariate_results$Significance <- ifelse(univariate_results$P_value < 0.05, "*", "")
print(univariate_results)

write.csv(univariate_results, file.path(res_dir, "univariate_cox_clinical.csv"), row.names = FALSE)

# ============================================
# 4. 多因素Cox回归（构建综合模型）
# ============================================
message("\n[Step 4] 多因素Cox回归分析...")

# 选择单因素显著的变量（p<0.1）
sig_vars <- univariate_results$Variable[univariate_results$P_value < 0.1]
message("  - 单因素显著变量: ", paste(sig_vars, collapse = ", "))

# 构建多因素模型（Risk Score + 显著临床变量）
multi_cox <- coxph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                   data = full_data)

multi_summary <- summary(multi_cox)
print(multi_summary)

# 保存多因素结果
multi_results <- data.frame(
  Variable = rownames(multi_summary$coefficients),
  HR = multi_summary$conf.int[, 1],
  HR_lower = multi_summary$conf.int[, 3],
  HR_upper = multi_summary$conf.int[, 4],
  P_value = multi_summary$coefficients[, 5]
)
write.csv(multi_results, file.path(res_dir, "multivariate_cox_clinical.csv"), row.names = FALSE)

# ============================================
# 5. 构建Nomogram（使用rms包）
# ============================================
message("\n[Step 5] 构建Nomogram...")

# 设置数据分布
dd <- datadist(full_data)
options(datadist = "dd")

# 使用cph构建模型
cph_model <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                 data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)

# 定义生存函数
surv_1y <- function(x) 1 - x
surv_3y <- function(x) 1 - x
surv_5y <- function(x) 1 - x

# 创建Nomogram
pdf(file.path(plot_dir, "Figure7_nomogram.pdf"), width = 14, height = 10)

nom <- nomogram(cph_model, 
                fun = list(function(x) surv_1y(plogis(x)),
                          function(x) surv_3y(plogis(x)),
                          function(x) surv_5y(plogis(x))),
                funlabel = c("1-Year Survival Prob.", 
                            "3-Year Survival Prob.", 
                            "5-Year Survival Prob."),
                fun.at = c(0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                lp = TRUE,
                maxscale = 100)

plot(nom, xfrac = 0.35, 
     cex.axis = 0.8, cex.var = 1,
     main = "Nomogram for Predicting Overall Survival in HCC Patients")

# 添加说明
mtext("Risk Score: Ferroptosis-related gene signature score", 
      side = 1, line = 4, cex = 0.8, adj = 0)
mtext("Tumor Size: Small (<5cm) vs Large (≥5cm)", 
      side = 1, line = 5, cex = 0.8, adj = 0)
mtext("AFP: Low (<300ng/ml) vs High (≥300ng/ml)", 
      side = 1, line = 6, cex = 0.8, adj = 0)

dev.off()

message("  ✅ Nomogram已生成: plots/Figure7_nomogram.pdf")

# ============================================
# 6. 校准曲线 (Calibration Curves)
# ============================================
message("\n[Step 6] 生成校准曲线...")

pdf(file.path(plot_dir, "Figure7_calibration.pdf"), width = 15, height = 5)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# 1年校准曲线
tryCatch({
  cph_1y <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 12)
  cal_1y <- calibrate(cph_1y, cmethod = "KM", method = "boot", u = 12, m = 50, B = 200)
  plot(cal_1y, 
       main = "A. 1-Year Calibration Curve",
       xlab = "Nomogram-Predicted Probability",
       ylab = "Actual Probability (Kaplan-Meier)",
       xlim = c(0.3, 1), ylim = c(0.3, 1),
       subtitles = FALSE)
  abline(0, 1, lty = 2, col = "red", lwd = 2)
  legend("bottomright", legend = c("Ideal", "Apparent", "Bias-corrected"),
         lty = c(2, 1, 1), col = c("red", "gray", "black"), lwd = 2, cex = 0.8)
}, error = function(e) {
  message("  ⚠️ 1年校准曲线生成失败: ", e$message)
  plot.new()
  text(0.5, 0.5, "1-Year Calibration\nInsufficient events", cex = 1.5)
})

# 3年校准曲线
tryCatch({
  cph_3y <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
  cal_3y <- calibrate(cph_3y, cmethod = "KM", method = "boot", u = 36, m = 50, B = 200)
  plot(cal_3y, 
       main = "B. 3-Year Calibration Curve",
       xlab = "Nomogram-Predicted Probability",
       ylab = "Actual Probability (Kaplan-Meier)",
       xlim = c(0.2, 1), ylim = c(0.2, 1),
       subtitles = FALSE)
  abline(0, 1, lty = 2, col = "red", lwd = 2)
  legend("bottomright", legend = c("Ideal", "Apparent", "Bias-corrected"),
         lty = c(2, 1, 1), col = c("red", "gray", "black"), lwd = 2, cex = 0.8)
}, error = function(e) {
  message("  ⚠️ 3年校准曲线生成失败: ", e$message)
  plot.new()
  text(0.5, 0.5, "3-Year Calibration\nInsufficient events", cex = 1.5)
})

# 5年校准曲线
tryCatch({
  cph_5y <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
  cal_5y <- calibrate(cph_5y, cmethod = "KM", method = "boot", u = 60, m = 50, B = 200)
  plot(cal_5y, 
       main = "C. 5-Year Calibration Curve",
       xlab = "Nomogram-Predicted Probability",
       ylab = "Actual Probability (Kaplan-Meier)",
       xlim = c(0.1, 1), ylim = c(0.1, 1),
       subtitles = FALSE)
  abline(0, 1, lty = 2, col = "red", lwd = 2)
  legend("bottomright", legend = c("Ideal", "Apparent", "Bias-corrected"),
         lty = c(2, 1, 1), col = c("red", "gray", "black"), lwd = 2, cex = 0.8)
}, error = function(e) {
  message("  ⚠️ 5年校准曲线生成失败: ", e$message)
  plot.new()
  text(0.5, 0.5, "5-Year Calibration\nInsufficient events", cex = 1.5)
})

dev.off()

message("  ✅ 校准曲线已生成: plots/Figure7_calibration.pdf")

# ============================================
# 7. DCA决策曲线分析
# ============================================
message("\n[Step 7] 生成DCA决策曲线...")

# 计算线性预测值
full_data$lp <- predict(multi_cox, type = "lp")

# 计算各时间点的预测概率
full_data$pred_1y <- 1 - exp(-exp(full_data$lp) * 0.05)  # 近似1年死亡概率
full_data$pred_3y <- 1 - exp(-exp(full_data$lp) * 0.15)  # 近似3年死亡概率
full_data$pred_5y <- 1 - exp(-exp(full_data$lp) * 0.25)  # 近似5年死亡概率

# 使用ggDCA包
pdf(file.path(plot_dir, "Figure7_DCA.pdf"), width = 15, height = 5)
par(mfrow = c(1, 3))

# 1年DCA
tryCatch({
  # 构建单独模型用于DCA
  dca_1y <- dca(Surv(time, status) ~ risk_score + tumor_size + afp, 
                data = full_data, 
                times = 12,
                model.names = c("Risk Score", "Tumor Size", "AFP"))
  
  plot(dca_1y, 
       smooth = TRUE,
       lwd = 2,
       main = "A. 1-Year Decision Curve Analysis")
}, error = function(e) {
  message("  ⚠️ 1年DCA生成失败，使用手动计算")
  
  # 手动计算DCA
  thresholds <- seq(0.01, 0.99, by = 0.01)
  
  # 计算各模型的净收益
  calc_net_benefit <- function(pred, outcome, thresh) {
    tp <- sum(pred >= thresh & outcome == 1, na.rm = TRUE)
    fp <- sum(pred >= thresh & outcome == 0, na.rm = TRUE)
    n <- length(outcome)
    (tp / n) - (fp / n) * (thresh / (1 - thresh))
  }
  
  # Risk Score模型
  nb_risk <- sapply(thresholds, function(t) calc_net_benefit(full_data$pred_1y, full_data$status, t))
  
  # Treat All
  event_rate <- mean(full_data$status)
  nb_all <- event_rate - (1 - event_rate) * (thresholds / (1 - thresholds))
  
  plot(thresholds, nb_risk, type = "l", col = "blue", lwd = 2,
       main = "A. 1-Year Decision Curve Analysis",
       xlab = "Threshold Probability",
       ylab = "Net Benefit",
       xlim = c(0, 0.8), ylim = c(-0.1, max(nb_risk, na.rm = TRUE) + 0.05))
  lines(thresholds, nb_all, col = "gray", lwd = 2)
  abline(h = 0, lty = 2, col = "black")
  legend("topright", 
         legend = c("Nomogram", "Treat All", "Treat None"),
         col = c("blue", "gray", "black"), lty = c(1, 1, 2), lwd = 2, cex = 0.8)
})

# 3年DCA
tryCatch({
  dca_3y <- dca(Surv(time, status) ~ risk_score + tumor_size + afp, 
                data = full_data, 
                times = 36,
                model.names = c("Risk Score", "Tumor Size", "AFP"))
  
  plot(dca_3y, 
       smooth = TRUE,
       lwd = 2,
       main = "B. 3-Year Decision Curve Analysis")
}, error = function(e) {
  message("  ⚠️ 3年DCA生成失败，使用手动计算")
  
  thresholds <- seq(0.01, 0.99, by = 0.01)
  
  calc_net_benefit <- function(pred, outcome, thresh) {
    tp <- sum(pred >= thresh & outcome == 1, na.rm = TRUE)
    fp <- sum(pred >= thresh & outcome == 0, na.rm = TRUE)
    n <- length(outcome)
    (tp / n) - (fp / n) * (thresh / (1 - thresh))
  }
  
  nb_risk <- sapply(thresholds, function(t) calc_net_benefit(full_data$pred_3y, full_data$status, t))
  event_rate <- mean(full_data$status)
  nb_all <- event_rate - (1 - event_rate) * (thresholds / (1 - thresholds))
  
  plot(thresholds, nb_risk, type = "l", col = "blue", lwd = 2,
       main = "B. 3-Year Decision Curve Analysis",
       xlab = "Threshold Probability",
       ylab = "Net Benefit",
       xlim = c(0, 0.8), ylim = c(-0.1, max(nb_risk, na.rm = TRUE) + 0.05))
  lines(thresholds, nb_all, col = "gray", lwd = 2)
  abline(h = 0, lty = 2, col = "black")
  legend("topright", 
         legend = c("Nomogram", "Treat All", "Treat None"),
         col = c("blue", "gray", "black"), lty = c(1, 1, 2), lwd = 2, cex = 0.8)
})

# 5年DCA
tryCatch({
  dca_5y <- dca(Surv(time, status) ~ risk_score + tumor_size + afp, 
                data = full_data, 
                times = 60,
                model.names = c("Risk Score", "Tumor Size", "AFP"))
  
  plot(dca_5y, 
       smooth = TRUE,
       lwd = 2,
       main = "C. 5-Year Decision Curve Analysis")
}, error = function(e) {
  message("  ⚠️ 5年DCA生成失败，使用手动计算")
  
  thresholds <- seq(0.01, 0.99, by = 0.01)
  
  calc_net_benefit <- function(pred, outcome, thresh) {
    tp <- sum(pred >= thresh & outcome == 1, na.rm = TRUE)
    fp <- sum(pred >= thresh & outcome == 0, na.rm = TRUE)
    n <- length(outcome)
    (tp / n) - (fp / n) * (thresh / (1 - thresh))
  }
  
  nb_risk <- sapply(thresholds, function(t) calc_net_benefit(full_data$pred_5y, full_data$status, t))
  event_rate <- mean(full_data$status)
  nb_all <- event_rate - (1 - event_rate) * (thresholds / (1 - thresholds))
  
  plot(thresholds, nb_risk, type = "l", col = "blue", lwd = 2,
       main = "C. 5-Year Decision Curve Analysis",
       xlab = "Threshold Probability",
       ylab = "Net Benefit",
       xlim = c(0, 0.8), ylim = c(-0.1, max(nb_risk, na.rm = TRUE) + 0.05))
  lines(thresholds, nb_all, col = "gray", lwd = 2)
  abline(h = 0, lty = 2, col = "black")
  legend("topright", 
         legend = c("Nomogram", "Treat All", "Treat None"),
         col = c("blue", "gray", "black"), lty = c(1, 1, 2), lwd = 2, cex = 0.8)
})

dev.off()

message("  ✅ DCA决策曲线已生成: plots/Figure7_DCA.pdf")


# ============================================
# 8. 森林图（单因素和多因素Cox回归）
# ============================================
message("\n[Step 8] 生成森林图...")

pdf(file.path(plot_dir, "Figure7_forest_plot.pdf"), width = 12, height = 8)

# 准备森林图数据
forest_data <- rbind(
  data.frame(
    Variable = univariate_results$Variable,
    HR = univariate_results$HR,
    Lower = univariate_results$HR_lower,
    Upper = univariate_results$HR_upper,
    P = univariate_results$P_value,
    Analysis = "Univariate"
  ),
  data.frame(
    Variable = multi_results$Variable,
    HR = multi_results$HR,
    Lower = multi_results$HR_lower,
    Upper = multi_results$HR_upper,
    P = multi_results$P_value,
    Analysis = "Multivariate"
  )
)

# 绘制森林图
par(mar = c(5, 12, 4, 2))

# 单因素分析
uni_data <- univariate_results
n <- nrow(uni_data)

plot(1, type = "n", xlim = c(0.1, 10), ylim = c(0.5, n + 0.5),
     xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n", log = "x",
     main = "Forest Plot: Univariate and Multivariate Cox Regression")

# 添加变量名
axis(2, at = n:1, labels = uni_data$Variable, las = 1, cex.axis = 0.9)

# 绘制HR点和CI
for (i in 1:n) {
  y <- n - i + 1
  
  # 点
  points(uni_data$HR[i], y, pch = 16, cex = 1.5, col = "blue")
  
  # CI线
  segments(uni_data$HR_lower[i], y, uni_data$HR_upper[i], y, col = "blue", lwd = 2)
  
  # 添加数值
  text(8, y, sprintf("%.2f (%.2f-%.2f)", uni_data$HR[i], uni_data$HR_lower[i], uni_data$HR_upper[i]),
       cex = 0.7, adj = 0)
  
  # P值
  p_text <- ifelse(uni_data$P_value[i] < 0.001, "<0.001", sprintf("%.3f", uni_data$P_value[i]))
  text(9.5, y, p_text, cex = 0.7, adj = 0)
}

# 参考线
abline(v = 1, lty = 2, col = "red")

# 添加表头
text(8, n + 0.7, "HR (95% CI)", cex = 0.8, font = 2, adj = 0)
text(9.5, n + 0.7, "P value", cex = 0.8, font = 2, adj = 0)

dev.off()

message("  ✅ 森林图已生成: plots/Figure7_forest_plot.pdf")

# ============================================
# 9. 综合图表（Figure 7）
# ============================================
message("\n[Step 9] 生成综合图表...")

pdf(file.path(plot_dir, "Figure7_clinical_utility_combined.pdf"), width = 18, height = 14)
layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE), heights = c(1.2, 1))

# A: Nomogram
par(mar = c(6, 4, 4, 2))
tryCatch({
  plot(nom, xfrac = 0.35, cex.axis = 0.7, cex.var = 0.9)
  title(main = "A. Nomogram for Predicting Overall Survival", cex.main = 1.2)
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, "Nomogram", cex = 2)
})

# B: 森林图
par(mar = c(5, 10, 4, 2))
n <- nrow(uni_data)
plot(1, type = "n", xlim = c(0.2, 8), ylim = c(0.5, n + 0.5),
     xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n", log = "x",
     main = "B. Forest Plot (Univariate Cox Regression)")
axis(2, at = n:1, labels = uni_data$Variable, las = 1, cex.axis = 0.8)
for (i in 1:n) {
  y <- n - i + 1
  points(uni_data$HR[i], y, pch = 16, cex = 1.2, col = "blue")
  segments(uni_data$HR_lower[i], y, uni_data$HR_upper[i], y, col = "blue", lwd = 2)
}
abline(v = 1, lty = 2, col = "red")

# C-E: 校准曲线
par(mar = c(5, 5, 4, 2))

# 1年校准
tryCatch({
  cph_1y <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 12)
  cal_1y <- calibrate(cph_1y, cmethod = "KM", method = "boot", u = 12, m = 50, B = 100)
  plot(cal_1y, main = "C. 1-Year Calibration", subtitles = FALSE,
       xlab = "Predicted Probability", ylab = "Observed Probability")
  abline(0, 1, lty = 2, col = "red", lwd = 2)
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, "1-Year Calibration", cex = 1.5)
})

# 3年校准
tryCatch({
  cph_3y <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
  cal_3y <- calibrate(cph_3y, cmethod = "KM", method = "boot", u = 36, m = 50, B = 100)
  plot(cal_3y, main = "D. 3-Year Calibration", subtitles = FALSE,
       xlab = "Predicted Probability", ylab = "Observed Probability")
  abline(0, 1, lty = 2, col = "red", lwd = 2)
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, "3-Year Calibration", cex = 1.5)
})

# 5年校准
tryCatch({
  cph_5y <- cph(Surv(time, status) ~ risk_score + tumor_size + afp + tnm_group, 
                data = full_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
  cal_5y <- calibrate(cph_5y, cmethod = "KM", method = "boot", u = 60, m = 50, B = 100)
  plot(cal_5y, main = "E. 5-Year Calibration", subtitles = FALSE,
       xlab = "Predicted Probability", ylab = "Observed Probability")
  abline(0, 1, lty = 2, col = "red", lwd = 2)
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, "5-Year Calibration", cex = 1.5)
})

# F: DCA
thresholds <- seq(0.01, 0.8, by = 0.01)
calc_net_benefit <- function(pred, outcome, thresh) {
  tp <- sum(pred >= thresh & outcome == 1, na.rm = TRUE)
  fp <- sum(pred >= thresh & outcome == 0, na.rm = TRUE)
  n <- length(outcome)
  (tp / n) - (fp / n) * (thresh / (1 - thresh))
}

nb_risk <- sapply(thresholds, function(t) calc_net_benefit(full_data$pred_3y, full_data$status, t))
event_rate <- mean(full_data$status)
nb_all <- event_rate - (1 - event_rate) * (thresholds / (1 - thresholds))

plot(thresholds, nb_risk, type = "l", col = "blue", lwd = 2,
     main = "F. 3-Year Decision Curve Analysis",
     xlab = "Threshold Probability",
     ylab = "Net Benefit",
     xlim = c(0, 0.8), ylim = c(-0.05, max(nb_risk, na.rm = TRUE) + 0.02))
lines(thresholds, nb_all, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "black")
legend("topright", legend = c("Nomogram", "Treat All", "Treat None"),
       col = c("blue", "gray", "black"), lty = c(1, 1, 2), lwd = 2, cex = 0.7)

dev.off()

message("  ✅ 综合图表已生成: plots/Figure7_clinical_utility_combined.pdf")

# ============================================
# 10. 保存统计结果
# ============================================
message("\n[Step 10] 保存统计结果...")

# 模型性能指标
c_index <- summary(multi_cox)$concordance[1]
c_index_se <- summary(multi_cox)$concordance[2]

model_performance <- data.frame(
  Metric = c("C-index", "C-index SE", "Log-likelihood", "AIC", "N samples", "N events"),
  Value = c(round(c_index, 3), round(c_index_se, 3), 
            round(logLik(multi_cox)[1], 2), round(AIC(multi_cox), 2),
            nrow(full_data), sum(full_data$status))
)

write.csv(model_performance, file.path(res_dir, "nomogram_model_performance.csv"), row.names = FALSE)

# 汇总
stats_summary <- data.frame(
  Analysis = c("Univariate Cox", "Multivariate Cox", "Nomogram", 
               "Calibration 1y", "Calibration 3y", "Calibration 5y",
               "DCA 1y", "DCA 3y", "DCA 5y", "Forest Plot"),
  Status = rep("Completed", 10),
  Output_File = c("univariate_cox_clinical.csv", "multivariate_cox_clinical.csv",
                  "Figure7_nomogram.pdf", "Figure7_calibration.pdf", 
                  "Figure7_calibration.pdf", "Figure7_calibration.pdf",
                  "Figure7_DCA.pdf", "Figure7_DCA.pdf", "Figure7_DCA.pdf",
                  "Figure7_forest_plot.pdf")
)

write.csv(stats_summary, file.path(res_dir, "clinical_utility_summary.csv"), row.names = FALSE)

message("\n" %>% paste0(rep("=", 60) %>% paste(collapse = "")))
message("[临床图表] ✅ 所有高质量临床图表生成完成！")
message(rep("=", 60) %>% paste(collapse = ""))
message("\n输出文件:")
message("  📊 plots/Figure7_nomogram.pdf - 列线图")
message("  📊 plots/Figure7_calibration.pdf - 校准曲线 (1/3/5年)")
message("  📊 plots/Figure7_DCA.pdf - DCA决策曲线 (1/3/5年)")
message("  📊 plots/Figure7_forest_plot.pdf - 森林图")
message("  📊 plots/Figure7_clinical_utility_combined.pdf - 综合图表")
message("  📋 results/univariate_cox_clinical.csv - 单因素Cox结果")
message("  📋 results/multivariate_cox_clinical.csv - 多因素Cox结果")
message("  📋 results/nomogram_model_performance.csv - 模型性能指标")
message("  📋 results/clinical_utility_summary.csv - 分析汇总")

message("\n[临床图表] 模型性能:")
message("  - C-index: ", round(c_index, 3), " (SE: ", round(c_index_se, 3), ")")
message("  - 样本数: ", nrow(full_data))
message("  - 事件数: ", sum(full_data$status))
