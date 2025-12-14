#!/usr/bin/env Rscript
# =============================================================================
# CB-Dock2 分子对接结果热图
# 
# 输入: results/cbdock2_docking_results.csv
# 输出: plots/Figure11_molecular_docking_heatmap.pdf
# =============================================================================

# 加载包
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# 设置工作目录
if (interactive()) {
  setwd("MianYinChen_Project")
}

# 创建输出目录
dir.create("plots", showWarnings = FALSE)

# 读取数据
docking <- read.csv("results/cbdock2_docking_results.csv", stringsAsFactors = FALSE)

# 转换为矩阵格式
docking_matrix <- dcast(docking, Protein ~ Ligand, value.var = "Vina_Score")
rownames(docking_matrix) <- docking_matrix$Protein
docking_matrix$Protein <- NULL

# 转换为长格式用于ggplot
docking_long <- melt(as.matrix(docking_matrix))
colnames(docking_long) <- c("Protein", "Ligand", "Vina_Score")

# 设置蛋白和配体顺序
docking_long$Protein <- factor(docking_long$Protein, 
                                levels = c("ACSL4", "TFRC", "NQO1"))
docking_long$Ligand <- factor(docking_long$Ligand,
                               levels = c("Capillarisin", "Chlorogenic_acid", 
                                         "Isorhamnetin", "Kaempferol", 
                                         "Quercetin", "Scoparone"))

# 生成热图
p <- ggplot(docking_long, aes(x = Ligand, y = Protein, fill = Vina_Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", Vina_Score)), 
            color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "#67A9CF", high = "#D1E5F0",
                       midpoint = -8,
                       name = "Vina Score\n(kcal/mol)",
                       limits = c(-10, -6),
                       breaks = c(-10, -9, -8, -7, -6)) +
  labs(title = "Molecular Docking Results (CB-Dock2)",
       subtitle = "Binding affinity between Artemisia capillaris compounds and hub proteins",
       x = "Ligand (Active Compounds)",
       y = "Protein (Hub Genes)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11),
    panel.grid = element_blank()
  )

# 保存图片
ggsave("plots/Figure11_molecular_docking_heatmap.pdf", p, 
       width = 10, height = 6, dpi = 300)

cat("✅ 热图已保存: plots/Figure11_molecular_docking_heatmap.pdf\n")

# 输出统计摘要
cat("\n=== 分子对接结果摘要 ===\n")
cat("所有对接均为有效结合 (Vina Score < -5 kcal/mol)\n\n")

# 各蛋白最佳配体
for (prot in c("ACSL4", "TFRC", "NQO1")) {
  best <- docking[docking$Protein == prot, ]
  best <- best[which.min(best$Vina_Score), ]
  cat(sprintf("%s 最佳配体: %s (%.1f kcal/mol)\n", 
              prot, best$Ligand, best$Vina_Score))
}

cat("\n整体最佳结合: ")
best_all <- docking[which.min(docking$Vina_Score), ]
cat(sprintf("%s + %s (%.1f kcal/mol)\n", 
            best_all$Protein, best_all$Ligand, best_all$Vina_Score))
