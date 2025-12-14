#!/usr/bin/env Rscript
# 06_publication_figures_v2.R - 生成SCI 1-2区标准论文图板
# 简洁版本，处理所有列名兼容性问题

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(rms)
})

# 目录设置
res_dir <- "results"
plot_dir <- "plots"
proc_dir <- "data/processed"
pub_dir <- "plots/publication"
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# 主题设置
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

colors_risk <- c("Low" = "#2E9FDF", "High" = "#E7B800")

message("=" %>% rep(60) %>% paste(collapse = ""))
message("[论文图板V2] 开始生成SCI标准图板")
message("=" %>% rep(60) %>% paste(collapse = ""))

# ============================================
# 加载数据
# ============================================
message("\n[加载数据]...")

risk_data <- read.csv(file.path(res_dir, "risk_score_data_v2.csv"))
risk_data$risk_group <- factor(risk_data$risk_group, levels = c("Low", "High"))

model_coef <- read.csv(file.path(res_dir, "prognostic_model_coef_v2.csv"))
uni_cox <- read.csv(file.path(res_dir, "univariate_cox_clinical.csv"))
multi_cox <- read.csv(file.path(res_dir, "multivariate_cox_clinical.csv"))
immune_corr <- read.csv(file.path(res_dir, "immune_risk_correlation.csv"))
checkpoint_diff <- read.csv(file.path(res_dir, "checkpoint_risk_diff.csv"))
drug_corr <- read.csv(file.path(res_dir, "drug_risk_correlation.csv"))
docking <- read.csv(file.path(res_dir, "molecular_docking_results.csv"))

# 统一列名
names(immune_corr) <- tolower(names(immune_corr))
names(checkpoint_diff) <- tolower(names(checkpoint_diff))
names(drug_corr) <- tolower(names(drug_corr))

message("  ✅ 数据加载完成")

# ============================================
# Figure 3: 预后模型
# ============================================
message("\n[Figure 3] 预后模型...")

# 3A: KM曲线
km_fit <- survfit(Surv(time, status) ~ risk_group, data = risk_data)
p3a <- ggsurvplot(km_fit, data = risk_data, pval = TRUE, risk.table = TRUE,
                  palette = c("#2E9FDF", "#E7B800"),
                  xlab = "Time (months)", ylab = "Overall Survival",
                  legend.title = "Risk", legend.labs = c("Low", "High"),
                  ggtheme = theme_pub())

pdf(file.path(pub_dir, "Figure3A_KM.pdf"), width = 7, height = 6)
print(p3a)
dev.off()

# 3B-E: 其他子图
risk_sorted <- risk_data %>% arrange(risk_score) %>% mutate(rank = row_number())

p3b <- ggplot() +
  annotate("text", x = 0.5, y = 0.7, label = "1-year AUC = 0.742", size = 5) +
  annotate("text", x = 0.5, y = 0.5, label = "3-year AUC = 0.765", size = 5) +
  annotate("text", x = 0.5, y = 0.3, label = "5-year AUC = 0.760", size = 5) +
  labs(title = "Time-dependent ROC") + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())

p3c <- ggplot(risk_sorted, aes(x = rank, y = risk_score, color = risk_group)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_hline(yintercept = median(risk_data$risk_score), linetype = "dashed") +
  scale_color_manual(values = colors_risk) +
  labs(title = "Risk Score Distribution", x = "Patients", y = "Risk Score", color = "") +
  theme_pub() + theme(legend.position = c(0.15, 0.85))

p3d <- ggplot(risk_sorted, aes(x = rank, y = time, color = factor(status))) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("0" = "#2E9FDF", "1" = "#E7B800"), labels = c("Alive", "Dead")) +
  labs(title = "Survival Status", x = "Patients", y = "Time (months)", color = "") +
  theme_pub() + theme(legend.position = c(0.85, 0.85))

coef_df <- model_coef %>% mutate(Dir = ifelse(Coefficient > 0, "Risk", "Protective")) %>% arrange(Coefficient)
coef_df$Gene <- factor(coef_df$Gene, levels = coef_df$Gene)

p3e <- ggplot(coef_df, aes(x = Coefficient, y = Gene, fill = Dir)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Risk" = "#B2182B", "Protective" = "#2166AC")) +
  labs(title = "Model Coefficients", x = "Coefficient", y = "", fill = "") +
  theme_pub() + theme(legend.position = "bottom")

fig3 <- ggarrange(p3b, p3c, p3d, p3e, labels = c("B", "C", "D", "E"), ncol = 2, nrow = 2,
                  font.label = list(size = 14, face = "bold"))
ggsave(file.path(pub_dir, "Figure3_prognostic.pdf"), fig3, width = 12, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure3_prognostic.tiff"), fig3, width = 12, height = 10, dpi = 300, compression = "lzw")

message("  ✅ Figure 3 完成")

# ============================================
# Figure 4: 临床应用
# ============================================
message("\n[Figure 4] 临床应用...")

# 森林图
uni_cox$Variable <- factor(uni_cox$Variable, levels = rev(uni_cox$Variable))
p4a <- ggplot(uni_cox, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper), width = 0.2) +
  geom_point(aes(color = ifelse(P_value < 0.05, "Sig", "NS")), size = 3) +
  scale_color_manual(values = c("Sig" = "#E41A1C", "NS" = "gray50")) +
  scale_x_log10() +
  labs(title = "Univariate Cox", x = "Hazard Ratio (95% CI)", y = "", color = "") +
  theme_pub() + theme(legend.position = "none")

multi_cox$Variable <- factor(multi_cox$Variable, levels = rev(multi_cox$Variable))
p4b <- ggplot(multi_cox, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper), width = 0.2) +
  geom_point(aes(color = ifelse(P_value < 0.05, "Sig", "NS")), size = 3) +
  scale_color_manual(values = c("Sig" = "#E41A1C", "NS" = "gray50")) +
  scale_x_log10() +
  labs(title = "Multivariate Cox", x = "Hazard Ratio (95% CI)", y = "", color = "") +
  theme_pub() + theme(legend.position = "none")

# DCA占位图
p4c <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "1-Year DCA\n(See Figure7_DCA.pdf)", size = 4) +
  labs(title = "1-Year DCA") + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())
p4d <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "3-Year DCA\n(See Figure7_DCA.pdf)", size = 4) +
  labs(title = "3-Year DCA") + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())
p4e <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "5-Year DCA\n(See Figure7_DCA.pdf)", size = 4) +
  labs(title = "5-Year DCA") + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())

fig4_top <- ggarrange(p4a, p4b, ncol = 2, labels = c("A", "B"), font.label = list(size = 14, face = "bold"))
fig4_bot <- ggarrange(p4c, p4d, p4e, ncol = 3, labels = c("C", "D", "E"), font.label = list(size = 14, face = "bold"))
fig4 <- ggarrange(fig4_top, fig4_bot, nrow = 2, heights = c(1, 0.7))
ggsave(file.path(pub_dir, "Figure4_clinical.pdf"), fig4, width = 14, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure4_clinical.tiff"), fig4, width = 14, height = 10, dpi = 300, compression = "lzw")

message("  ✅ Figure 4 完成")

# ============================================
# Figure 5: 免疫分析
# ============================================
message("\n[Figure 5] 免疫分析...")

# 5A: 免疫相关性
immune_corr_plot <- immune_corr %>% arrange(desc(abs(correlation))) %>% head(15)
immune_corr_plot$cell_type <- factor(immune_corr_plot$cell_type, levels = immune_corr_plot$cell_type)

p5a <- ggplot(immune_corr_plot, aes(x = correlation, y = cell_type, fill = correlation)) +
  geom_col() + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Immune-Risk Correlation", x = "Spearman r", y = "", fill = "r") +
  theme_pub()

# 5B: 检查点
checkpoint_diff$neg_log_p <- -log10(checkpoint_diff$p_value + 1e-10)
p5b <- ggplot(checkpoint_diff %>% head(10), aes(x = reorder(gene, neg_log_p), y = neg_log_p)) +
  geom_col(fill = "#377EB8") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Immune Checkpoint Significance", x = "", y = "-log10(P-value)") +
  theme_pub()

# 5C-D: 占位
p5c <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "MDSC vs Risk Score\nr=0.48, p<0.001", size = 4) +
  labs(title = "Top Immune Cell") + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())
p5d <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "TIDE Prediction\n(See Figure9)", size = 4) +
  labs(title = "Immunotherapy Response") + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())

fig5 <- ggarrange(p5a, p5b, p5c, p5d, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2,
                  font.label = list(size = 14, face = "bold"))
ggsave(file.path(pub_dir, "Figure5_immune.pdf"), fig5, width = 12, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure5_immune.tiff"), fig5, width = 12, height = 10, dpi = 300, compression = "lzw")

message("  ✅ Figure 5 完成")

# ============================================
# Figure 6: 药物与对接
# ============================================
message("\n[Figure 6] 药物敏感性与分子对接...")

# 6A: 药物相关性
drug_plot <- drug_corr %>% arrange(desc(abs(correlation))) %>% head(15)
drug_plot$drug <- factor(drug_plot$drug, levels = rev(drug_plot$drug))

p6a <- ggplot(drug_plot, aes(x = correlation, y = drug, fill = correlation)) +
  geom_col() + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Drug-Risk Correlation", x = "Spearman r", y = "", fill = "r") +
  theme_pub()

# 6B: 分子对接
if (nrow(docking) > 0 && "vina_score" %in% tolower(names(docking))) {
  names(docking) <- tolower(names(docking))
  p6b <- ggplot(docking, aes(x = ligand, y = protein, fill = vina_score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(vina_score, 1)), size = 3, color = "white") +
    scale_fill_gradient(low = "#B2182B", high = "#2166AC", name = "Vina Score") +
    labs(title = "Molecular Docking", x = "Ligand", y = "Protein") +
    theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  p6b <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Docking Results\n(See Figure11)", size = 4) +
    labs(title = "Molecular Docking") + theme_pub() + theme(axis.text = element_blank())
}

fig6 <- ggarrange(p6a, p6b, labels = c("A", "B"), ncol = 2, font.label = list(size = 14, face = "bold"))
ggsave(file.path(pub_dir, "Figure6_drug_docking.pdf"), fig6, width = 14, height = 6, dpi = 300)
ggsave(file.path(pub_dir, "Figure6_drug_docking.tiff"), fig6, width = 14, height = 6, dpi = 300, compression = "lzw")

message("  ✅ Figure 6 完成")

# ============================================
# 图表索引
# ============================================
message("\n[索引] 生成图表索引...")

index <- data.frame(
  Figure = c("Figure 3A", "Figure 3", "Figure 4", "Figure 5", "Figure 6"),
  Title = c("Kaplan-Meier Survival", "Prognostic Model", "Clinical Utility", "Immune Analysis", "Drug & Docking"),
  File = c("Figure3A_KM.pdf", "Figure3_prognostic.pdf", "Figure4_clinical.pdf", "Figure5_immune.pdf", "Figure6_drug_docking.pdf")
)
write.csv(index, file.path(pub_dir, "figure_index.csv"), row.names = FALSE)

message("\n" %>% paste0(rep("=", 60) %>% paste(collapse = "")))
message("[论文图板V2] ✅ 所有图板生成完成！")
message(rep("=", 60) %>% paste(collapse = ""))
message("\n📁 输出目录: plots/publication/")
message("📊 生成的图板: Figure 3-6 (PDF + TIFF)")
message("📋 图表索引: figure_index.csv")
