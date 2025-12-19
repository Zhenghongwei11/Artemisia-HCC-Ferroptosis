#!/usr/bin/env Rscript
# 07_final_publication_figures.R - 最终版论文图板 (修正版 V3 - 消除硬编码数据)

# Ensure non-interactive execution
if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(rms)
  library(timeROC)
  library(gridExtra)
  library(grid)
  library(cowplot)
})

detect_base_dir <- function() {
  candidates <- c(".", "fx_review_package")
  for (d in candidates) {
    if (dir.exists(file.path(d, "results")) &&
        dir.exists(file.path(d, "plots")) &&
        dir.exists(file.path(d, "data/processed")) &&
        dir.exists(file.path(d, "data/references"))) {
      return(d)
    }
  }
  stop(
    "Cannot locate fx_review_package base directory.\n",
    "Run either:\n",
    "  1) from fx_review_package/ directory, or\n",
    "  2) from its parent directory (so fx_review_package/ is a subfolder).\n"
  )
}

# 目录设置
base_dir <- detect_base_dir()
res_dir <- file.path(base_dir, "results")
proc_dir <- file.path(base_dir, "data/processed")
ref_dir <- file.path(base_dir, "data/references") # For gene lists
pub_dir <- file.path(base_dir, "plots/publication")
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

# 主题设置
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold")
    )
}

colors_risk <- c("Low" = "#2E9FDF", "High" = "#E7B800")
colors_reg <- c("Up" = "#B2182B", "Down" = "#2166AC", "NS" = "gray70")

message("============================================================")
message("[最终图板 V3] 开始生成SCI标准论文图板 (消除硬编码数据)")
message("============================================================")

# ============================================ 
# 加载数据
# ============================================ 
message("\n[1/7] 加载数据...")

risk_data <- read.csv(file.path(res_dir, "risk_score_data_v2.csv"))
risk_data$risk_group <- factor(risk_data$risk_group, levels = c("Low", "High"))

model_coef <- read.csv(file.path(res_dir, "prognostic_model_coef_v2.csv"))
model_stats <- read.csv(file.path(res_dir, "prognostic_model_stats_v2.csv")) # NEW: For C-index
uni_cox <- read.csv(file.path(res_dir, "univariate_cox_clinical.csv"))
multi_cox <- read.csv(file.path(res_dir, "multivariate_cox_clinical.csv"))

immune_corr <- read.csv(file.path(res_dir, "immune_risk_correlation.csv") )
names(immune_corr) <- tolower(names(immune_corr)) 
if( !"cell_type" %in% names(immune_corr)) names(immune_corr)[1] <- "cell_type"

immune_scores <- read.csv(file.path(res_dir, "ssgsea_immune_scores.csv"), row.names = 1)
checkpoint_diff <- read.csv(file.path(res_dir, "checkpoint_risk_diff.csv") )
names(checkpoint_diff) <- tolower(names(checkpoint_diff))

drug_corr <- read.csv(file.path(res_dir, "drug_risk_correlation.csv") )
names(drug_corr) <- tolower(names(drug_corr))

docking <- read.csv(file.path(res_dir, "molecular_docking_results.csv") )
names(docking) <- tolower(names(docking))

clinical_14 <- readRDS(file.path(proc_dir, "GSE14520_tumor_clinical.rds") )

message("  ✅ 数据加载完成")

# ============================================ 
# Figure 2: DEG (A-D) - Fix Dynamic Venn Count
# ============================================ 
message("\n[2/7] Figure 2: DEG分析...")

deg_14520 <- read.csv(file.path(res_dir, "deg_GSE14520_all.csv") )
deg_83148 <- read.csv(file.path(res_dir, "deg_GSE83148_all.csv") )

# Helper for volcano
plot_volcano <- function(df, title) {
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < 0.05 & df$logFC > 1] <- "Up"
  df$Sig[df$adj.P.Val < 0.05 & df$logFC < -1] <- "Down"
  df$nlogp <- -log10(df$adj.P.Val + 1e-300)
  
  ggplot(df, aes(x = logFC, y = nlogp, color = Sig)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = colors_reg) +
    geom_vline(xintercept = c(-1, 1), linetype="dashed", color="grey") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(Adj. P)") +
    theme_pub() + theme(legend.position = "top")
}

p2a <- plot_volcano(deg_14520, "GSE14520 (HCC)")
p2b <- plot_volcano(deg_83148, "GSE83148 (CHB)")

# C: Venn Diagram (Manual Construction with White Background & Dynamic Counts)
circle_fun <- function(center = c(0,0), diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
c1 <- circle_fun(c(0.35, 0.6), 0.6); c1$group <- "DEGs"
c2 <- circle_fun(c(0.65, 0.6), 0.6); c2$group <- "Ferroptosis"
c3 <- circle_fun(c(0.5, 0.35), 0.6); c3$group <- "TCM Targets"
venn_dat <- rbind(c1, c2, c3)

# Dynamic calculation of intersection counts
ferroptosis_genes <- read.csv(file.path(ref_dir, "ferroptosis_genes_expanded.csv") )$Gene
tcm_targets <- read.csv(file.path(ref_dir, "tcm_targets_CHB.csv") )$Target
deg_sig <- deg_14520 %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
deg_genes <- deg_sig$Gene

hub_genes_actual <- Reduce(intersect, list(deg_genes, ferroptosis_genes, tcm_targets))
n_hub_actual <- length(hub_genes_actual)

p2c <- ggplot(venn_dat, aes(x, y)) +
  geom_polygon(aes(fill = group, group = group), alpha = 0.3) +
  geom_path(aes(color = group, group = group), linewidth = 1) +
  scale_fill_manual(values = c("DEGs"="#E41A1C", "Ferroptosis"="#377EB8", "TCM Targets"="#4DAF4A")) +
  scale_color_manual(values = c("DEGs"="#E41A1C", "Ferroptosis"="#377EB8", "TCM Targets"="#4DAF4A")) +
  annotate("text", x=0.2, y=0.8, label="DEGs", color="#E41A1C", fontface="bold") +
  annotate("text", x=0.8, y=0.8, label="Ferroptosis", color="#377EB8", fontface="bold") +
  annotate("text", x=0.5, y=0.15, label="TCM Targets", color="#4DAF4A", fontface="bold") +
  annotate("text", x=0.5, y=0.5, label=paste0(n_hub_actual, "\n(Hub)"), fontface="bold", size=5) + # Dynamic count
  coord_equal() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none") +
  labs(title = "Intersection Analysis")

# D: Hub Genes (使用动态计算的hub_genes_actual)
hub_dat <- deg_14520 %>% filter(Gene %in% hub_genes_actual)
hub_dat$Direction <- ifelse(hub_dat$logFC > 0, "Up", "Down")
p2d <- ggplot(hub_dat, aes(x = logFC, y = reorder(Gene, logFC), fill = Direction)) +
  geom_col(width=0.6) +
  scale_fill_manual(values = c("Up"="#B2182B", "Down"="#2166AC")) +
  geom_vline(xintercept = 0) +
  labs(title = "Hub Gene Expression", x = "log2 Fold Change", y = "", fill = "") +
  theme_pub() + theme(legend.position = "bottom")

fig2 <- plot_grid(p2a, p2b, p2c, p2d, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 16)
ggsave(file.path(pub_dir, "Figure2_DEG_analysis.pdf"), fig2, width = 12, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure2_DEG_analysis.tiff"), fig2, width = 12, height = 10, dpi = 300, compression = "lzw", bg="white")
message("  ✅ Figure 2 完成")

# ============================================ 
# Figure 3: 预后模型 (A-D)
# ============================================ 
message("\n[3/7] Figure 3: 预后模型...")

# A: KM Curve
km_fit <- survfit(Surv(time, status) ~ risk_group, data = risk_data)
p3a_obj <- ggsurvplot(km_fit, data = risk_data, 
                  pval = TRUE, 
                  palette = c("#2E9FDF", "#E7B800"),
                  xlab = "Time (months)", ylab = "Overall Survival",
                  legend.title = "", legend.labs = c("Low Risk", "High Risk"),
                  ggtheme = theme_pub(),
                  risk.table = FALSE) # Disable table for panel integration to avoid complexity
p3a <- p3a_obj$plot + labs(title = "Kaplan-Meier Survival")

# B: Time-dependent ROC
roc_res <- timeROC(T = risk_data$time, delta = risk_data$status, marker = risk_data$risk_score,
                   cause = 1, weighting = "marginal", times = c(12, 36, 60), iid = FALSE)

roc_df <- data.frame()
for (t in c(12, 36, 60)) {
  idx <- which(roc_res$times == t)
  tmp <- data.frame(FPR = roc_res$FP[, idx], TPR = roc_res$TP[, idx], 
                    Time = paste0(t/12, "-year (AUC=", round(roc_res$AUC[idx], 3), ")"))
  roc_df <- rbind(roc_df, tmp)
}
p3b <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Time)) +
  geom_path(linewidth = 1.2) + # Thicker lines
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) + 
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  labs(title = "Time-dependent ROC", x = "1 - Specificity", y = "Sensitivity", color = "") + 
  theme_pub() + 
  theme(legend.position = c(0.65, 0.25), 
        legend.background = element_rect(fill="white", color="black", size=0.2)) # Box around legend

# C: Risk Score
risk_sorted <- risk_data %>% arrange(risk_score) %>% mutate(rank = row_number())
p3c <- ggplot(risk_sorted, aes(x = rank, y = risk_score, color = risk_group)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = median(risk_data$risk_score), linetype = "dashed", color = "black") +
  scale_color_manual(values = colors_risk) +
  labs(title = "Risk Score Distribution", x = "Patients (ranked)", y = "Risk Score", color = "") +
  theme_pub() + theme(legend.position = c(0.15, 0.85))

# D: Coefficients
cof_df <- model_coef %>%
  mutate(Direction = ifelse(Coefficient > 0, "Risk", "Protective")) %>%
  arrange(Coefficient)
cof_df$Gene <- factor(cof_df$Gene, levels = cof_df$Gene)
p3d <- ggplot(cof_df, aes(x = Coefficient, y = Gene, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Risk" = "#B2182B", "Protective" = "#2166AC")) +
  labs(title = "Model Coefficients", x = "Coefficient", y = "", fill = "") +
  theme_pub() + theme(legend.position = "bottom", axis.text.y = element_text(size = 8))

fig3 <- plot_grid(p3a, p3b, p3c, p3d, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 16)
ggsave(file.path(pub_dir, "Figure3_prognostic_model.pdf"), fig3, width = 12, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure3_prognostic_model.tiff"), fig3, width = 12, height = 10, dpi = 300, compression = "lzw", bg="white")
message("  ✅ Figure 3 完成")

# ============================================ 
# Figure 4: 临床应用 (A-D) - Fix C-index hardcoding
# ============================================ 
message("\n[4/7] Figure 4: 临床应用价值...")

# A & B Cox
uni_cox$Variable <- factor(uni_cox$Variable, levels = rev(uni_cox$Variable))
p4a <- ggplot(uni_cox, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper), width = 0.25) +
  geom_point(aes(color = ifelse(P_value < 0.05, "p<0.05", "p>=0.05")), size = 3) +
  scale_color_manual(values = c("p<0.05" = "#E41A1C", "p>=0.05" = "gray50")) +
  scale_x_log10() +
  labs(title = "Univariate Cox", x = "Hazard Ratio", y = "", color = "") +
  theme_pub() + theme(legend.position = "bottom")

multi_cox$Variable <- factor(multi_cox$Variable, levels = rev(multi_cox$Variable))
p4b <- ggplot(multi_cox, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper), width = 0.25) +
  geom_point(aes(color = ifelse(P_value < 0.05, "p<0.05", "p>=0.05")), size = 3) +
  scale_color_manual(values = c("p<0.05" = "#E41A1C", "p>=0.05" = "gray50")) +
  scale_x_log10() +
  labs(title = "Multivariate Cox", x = "Hazard Ratio", y = "", color = "") +
  theme_pub() + theme(legend.position = "none")

# C: DCA - Robust Data Prep
clinical_matched <- clinical_14[match(risk_data$sample, clinical_14$Affy_GSM), ]
full_data <- data.frame(
  time = risk_data$time,
  status = risk_data$status,
  risk_score = risk_data$risk_score,
  tumor_size = ifelse(clinical_matched$Main.Tumor.Size......5.cm. == "large", 1, 0),
  afp = ifelse(clinical_matched$AFP......300ng.ml. == "high", 1, 0)
)
full_data <- full_data[complete.cases(full_data), ]

cox_model <- coxph(Surv(time, status) ~ risk_score + tumor_size + afp, data = full_data)
full_data$pred <- 1 - exp(-predict(cox_model, type = "expected"))
thresholds <- seq(0.01, 0.9, by = 0.02)
dca_res <- data.frame()
for(pt in thresholds) {
  tp <- sum(full_data$pred >= pt & full_data$status == 1)
  fp <- sum(full_data$pred >= pt & full_data$status == 0)
  n <- nrow(full_data)
  nb <- (tp/n) - (fp/n) * (pt / (1-pt))
  all_nb <- (sum(full_data$status==1)/n) - (sum(full_data$status==0)/n) * (pt / (1 - pt))
  dca_res <- rbind(dca_res, data.frame(Threshold = pt, Benefit = nb, Type = "Nomogram"))
  dca_res <- rbind(dca_res, data.frame(Threshold = pt, Benefit = all_nb, Type = "Treat All"))
  dca_res <- rbind(dca_res, data.frame(Threshold = pt, Benefit = 0, Type = "Treat None"))
}
ymax <- max(dca_res$Benefit[dca_res$Type=="Nomogram"], na.rm=T) + 0.05
p4c <- ggplot(dca_res, aes(x = Threshold, y = Benefit, color = Type, linetype = Type)) +
  geom_line(linewidth = 1.2) +
  coord_cartesian(ylim = c(-0.05, ymax), xlim = c(0, 0.8)) +
  scale_color_manual(values = c("Nomogram"="#E41A1C", "Treat All"="gray50", "Treat None"="black")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  labs(title = "Decision Curve Analysis", x = "Threshold Probability", y = "Net Benefit") +
  theme_pub() + theme(legend.position = c(0.75, 0.8))

# D: C-index (Dynamic Value - Calculate Combined Model C-index)
risk_cindex <- model_stats %>% filter(Metric == "C-index") %>% pull(Value) %>% as.numeric()

# Calculate Combined Model C-index from multivariate Cox model
combined_cox <- coxph(Surv(time, status) ~ risk_score + tumor_size + afp, data = full_data)
combined_cindex <- summary(combined_cox)$concordance[1]

cindex_df <- data.frame(
  Model = c("Risk Score", "Combined Model"),
  C_index = c(risk_cindex, combined_cindex)
)
cindex_df$Model <- factor(cindex_df$Model, levels = rev(cindex_df$Model))

p4d <- ggplot(cindex_df, aes(x = C_index, y = Model, fill = Model)) +
  geom_col(width = 0.6, alpha=0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%.3f", C_index)), hjust = -0.3, fontface = "bold") +
  scale_fill_manual(values = c("Risk Score" = "#377EB8", "Combined Model" = "#E41A1C")) +
  coord_cartesian(xlim = c(0.4, 0.9)) +
  labs(title = "C-index Comparison", x = "C-index", y = "") +
  theme_pub() + theme(legend.position = "none")

fig4 <- plot_grid(p4a, p4b, p4c, p4d, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 16)
ggsave(file.path(pub_dir, "Figure4_clinical_utility.pdf"), fig4, width = 12, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure4_clinical_utility.tiff"), fig4, width = 12, height = 10, dpi = 300, compression = "lzw", bg="white")
message("  ✅ Figure 4 完成")

# ============================================ 
# Figure 5: 免疫 (A-D)
# ============================================ 
message("\n[5/7] Figure 5: 免疫微环境...")

# A: Correlation
immune_corr_plot <- immune_corr %>% arrange(desc(abs(correlation))) %>% head(15)
immune_corr_plot$cell_type <- factor(immune_corr_plot$cell_type, levels = immune_corr_plot$cell_type)
p5a <- ggplot(immune_corr_plot, aes(x = correlation, y = cell_type, fill = correlation)) +
  geom_col() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits=c(-0.6, 0.6)) +
  labs(title = "Immune Correlation", x = "Correlation", y = "") +
  theme_pub()

# B: Boxplot
sig_cells <- immune_corr %>% filter(p_value < 0.05) %>% arrange(p_value) %>% head(4) %>% pull(cell_type)
if(length(sig_cells) == 0) sig_cells <- immune_corr %>% arrange(desc(abs(correlation))) %>% head(4) %>% pull(cell_type)
plot_data <- immune_scores %>%
  rownames_to_column("sample") %>%
  inner_join(risk_data %>% select(sample, risk_group), by = "sample") %>%
  select(sample, risk_group, all_of(sig_cells)) %>%
  pivot_longer(cols = -c(sample, risk_group), names_to = "Cell", values_to = "Score")
p5b <- ggplot(plot_data, aes(x = Cell, y = Score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(aes(group = risk_group), label = "p.signif") +
  scale_fill_manual(values = colors_risk) +
  labs(title = "Immune Differences", x = "", y = "Score") +
  theme_pub() + theme(axis.text.x = element_text(angle=30, hjust=1))

# C: Checkpoint
check_plot <- checkpoint_diff %>% arrange(p_value) %>% head(10)
check_plot$neg_log_p <- -log10(check_plot$p_value + 1e-300)
check_plot$gene <- factor(check_plot$gene, levels = rev(check_plot$gene))
p5c <- ggplot(check_plot, aes(x = neg_log_p, y = gene, fill = neg_log_p)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="red") +
  scale_fill_gradient(low = "#FEE0D2", high = "#CB181D") +
  labs(title = "Checkpoint Expression", x = "-log10(P-value)", y = "") +
  theme_pub() + theme(legend.position = "none")

# D: Scatter
top_cell <- sig_cells[1]
merged <- inner_join(risk_data, immune_scores %>% rownames_to_column("sample"), by="sample")
p5d <- ggplot(merged, aes_string(x = "risk_score", y = top_cell)) +
  geom_point(aes(color = risk_group), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = 0.5) +
  scale_color_manual(values = colors_risk) +
  labs(title = paste0("Risk Score vs ", top_cell), y = "Immune Score") +
  theme_pub() + theme(legend.position = "none")

fig5 <- plot_grid(p5a, p5b, p5c, p5d, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 16)
ggsave(file.path(pub_dir, "Figure5_immune_analysis.pdf"), fig5, width = 12, height = 10, dpi = 300)
ggsave(file.path(pub_dir, "Figure5_immune_analysis.tiff"), fig5, width = 12, height = 10, dpi = 300, compression = "lzw", bg="white")
message("  ✅ Figure 5 完成")

# ============================================ 
# Figure 6: 药物 (A-D)
# ============================================ 
message("\n[6/7] Figure 6: 药物与分子对接...")
drug_sub <- drug_corr %>% arrange(desc(abs(correlation))) %>% head(15)
drug_sub$drug <- factor(drug_sub$drug, levels = drug_sub$drug)
p6a <- ggplot(drug_sub, aes(x = correlation, y = drug, fill = correlation)) +
  geom_col() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(title = "Drug Sensitivity Correlation", x = "Correlation", y = "") +
  theme_pub() + theme(axis.text.y = element_text(size=8))

hcc_drug_file <- file.path(res_dir, "hcc_drug_risk_diff.csv")
if(file.exists(hcc_drug_file)) {
  hcc_dat <- read.csv(hcc_drug_file)
  if("High_Mean" %in% names(hcc_dat)) hcc_dat <- hcc_dat %>% rename(mean_high = High_Mean, mean_low = Low_Mean)
  names(hcc_dat) <- tolower(names(hcc_dat))
  hcc_plot <- hcc_dat %>% head(6) %>%
    select(drug, mean_high, mean_low) %>%
    pivot_longer(cols = c(mean_high, mean_low), names_to = "Group", values_to = "IC50") %>%
    mutate(Group = ifelse(Group == "mean_high", "High", "Low"))
  p6b <- ggplot(hcc_plot, aes(x = drug, y = IC50, fill = Group)) +
    geom_bar(stat="identity", position="dodge", width=0.7) +
    scale_fill_manual(values = colors_risk) +
    labs(title = "HCC Drug Sensitivity", x = "", y = "IC50") +
    theme_pub() + theme(axis.text.x = element_text(angle=30, hjust=1), legend.position = "top")
} else { p6b <- ggplot() + theme_void() }

# ===== CB-Dock2 传统分子对接结果 =====
cbdock_file <- file.path(res_dir, "cbdock2_docking_results.csv")
if(file.exists(cbdock_file)) {
  cbdock <- read.csv(cbdock_file)
  
  # CB-Dock2 热图
  cbdock$Protein <- factor(cbdock$Protein, levels = c("ACSL4", "TFRC", "NQO1"))
  cbdock$Ligand <- factor(cbdock$Ligand, levels = c("Capillarisin", "Chlorogenic_acid", 
                                                     "Isorhamnetin", "Kaempferol", 
                                                     "Quercetin", "Scoparone"))
  p6c <- ggplot(cbdock, aes(x = Ligand, y = Protein, fill = Vina_Score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", Vina_Score)), color = "white", size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "#67A9CF", high = "#D1E5F0",
                         midpoint = -8, name = "Vina Score\n(kcal/mol)",
                         limits = c(-10, -6)) +
    labs(title = "Molecular Docking (CB-Dock2)", x = "", y = "") +
    theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # CB-Dock2 最佳结合柱状图
  best_cb <- cbdock %>% group_by(Protein) %>% slice_min(Vina_Score, n = 1)
  best_cb$Label <- paste0(best_cb$Protein, "\n+ ", best_cb$Ligand)
  p6d <- ggplot(best_cb, aes(x = reorder(Label, Vina_Score), y = -Vina_Score, fill = Protein)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 7, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%.1f", Vina_Score)), vjust = 1.5, color = "white", fontface = "bold") +
    labs(title = "Top Binding Pairs (CB-Dock2)", x = "", y = "-Vina Score (kcal/mol)") +
    scale_fill_brewer(palette = "Set2") +
    coord_cartesian(ylim = c(0, 10)) +
    theme_pub() + theme(legend.position = "none")
} else { 
  p6c <- ggplot() + annotate("text", x=0.5, y=0.5, label="CB-Dock2 data\nnot found") + theme_void()
  p6d <- ggplot() + theme_void() 
}

# ===== Boltz-2 AI分子对接结果 (补充) =====
if(nrow(docking) > 0) {
  p6e <- ggplot(docking, aes(x = ligand, y = protein, fill = iptm)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = sprintf("%.2f", iptm)), color = "white", fontface="bold", size=3.5) +
    scale_fill_gradient(low = "#FEE0D2", high = "#B2182B", name = "ipTM\n(Boltz-2)") +
    labs(title = "Molecular Docking (Boltz-2)", x = "", y = "") + 
    theme_pub() + theme(axis.text.x = element_text(angle=45, hjust=1))
  
  best <- docking %>% group_by(protein) %>% slice_max(iptm, n=1)
  best$Label <- paste0(best$protein, "\n+ ", best$ligand)
  p6f <- ggplot(best, aes(x = reorder(Label, iptm), y = iptm, fill = protein)) +
    geom_col(width=0.5) + 
    geom_hline(yintercept = 0.5, linetype="dashed", color="red") +
    geom_text(aes(label=sprintf("%.2f", iptm)), vjust=-0.5, fontface="bold") +
    labs(title = "Top Binding Pairs (Boltz-2)", x = "", y = "ipTM Score") + 
    scale_fill_brewer(palette = "Set2") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_pub() + theme(legend.position = "none")
} else { p6e <- ggplot() + theme_void(); p6f <- ggplot() + theme_void() }

# 组合6个子图: A-B药物敏感性, C-D CB-Dock2, E-F Boltz-2
fig6 <- plot_grid(p6a, p6b, p6c, p6d, p6e, p6f, 
                  labels = c("A", "B", "C", "D", "E", "F"), 
                  ncol = 2, label_size = 16)
ggsave(file.path(pub_dir, "Figure6_drug_docking.pdf"), fig6, width = 14, height = 15, dpi = 300)
ggsave(file.path(pub_dir, "Figure6_drug_docking.tiff"), fig6, width = 14, height = 15, dpi = 300, compression = "lzw", bg="white")
message("  ✅ Figure 6 完成")

message("\n============================================================")
message("所有图板生成完毕! 输出目录: plots/publication/")
message("============================================================")
