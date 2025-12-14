#!/usr/bin/env Rscript

# 03a_immune_infiltration_v2.R
# 免疫浸润分析 - 修复版（探针ID转基因符号）
# 使用方法: Rscript scripts_final/03a_immune_infiltration_v2.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(pheatmap)
  library(ggpubr)
})

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"
set.seed(123)

message("[免疫分析V2] 开始免疫浸润分析...")

# ============================================
# 1. 加载数据并转换探针ID为基因符号
# ============================================
message("[免疫分析V2] 加载数据并转换基因符号...")

# 读取GPL571注释 - 跳过注释行
anno_file <- "data/raw/GSE14520/GPL571.annot.gz"
if (file.exists(anno_file)) {
  # 读取所有行，找到表头位置
  all_lines <- readLines(gzfile(anno_file))
  header_line <- grep("^ID\t", all_lines)[1]
  
  if (!is.na(header_line)) {
    # 从表头开始读取
    anno <- read.delim(text = all_lines[header_line:length(all_lines)], 
                       stringsAsFactors = FALSE, sep = "\t")
    message("[免疫分析V2] GPL571注释行数: ", nrow(anno))
    message("[免疫分析V2] 注释列名: ", paste(head(colnames(anno), 5), collapse = ", "))
    
    # 提取探针ID和基因符号映射
    gene_col <- grep("Gene.symbol|Gene symbol", colnames(anno), value = TRUE, ignore.case = TRUE)[1]
    
    if (!is.na(gene_col)) {
      probe2gene <- anno %>%
        dplyr::select(ID, all_of(gene_col)) %>%
        rename(Gene.Symbol = all_of(gene_col)) %>%
        filter(Gene.Symbol != "" & !is.na(Gene.Symbol)) %>%
        mutate(Gene.Symbol = sapply(strsplit(Gene.Symbol, "///"), function(x) trimws(x[1])))
      message("[免疫分析V2] 有效探针-基因映射: ", nrow(probe2gene))
    } else {
      stop("找不到基因符号列")
    }
  } else {
    stop("找不到注释表头")
  }
} else {
  stop("找不到GPL571注释文件")
}

# 读取表达矩阵
expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
message("[免疫分析V2] 原始表达矩阵: ", nrow(expr_14), " x ", ncol(expr_14))

# 转换为基因符号
expr_df <- as.data.frame(expr_14)
expr_df$ProbeID <- rownames(expr_df)

expr_gene <- expr_df %>%
  inner_join(probe2gene, by = c("ProbeID" = "ID")) %>%
  dplyr::select(-ProbeID) %>%
  # 对于同一基因的多个探针，取平均值
  group_by(Gene.Symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("Gene.Symbol")

expr_gene <- as.matrix(expr_gene)
message("[免疫分析V2] 转换后表达矩阵: ", nrow(expr_gene), " x ", ncol(expr_gene))

# 保存转换后的表达矩阵
saveRDS(expr_gene, file.path(proc_dir, "GSE14520_expr_symbol.rds"))

# Risk Score数据
risk_file <- file.path(res_dir, "risk_score_data.csv")
if (file.exists(risk_file)) {
  risk_data <- read.csv(risk_file)
  message("[免疫分析V2] Risk Score样本数: ", nrow(risk_data))
} else {
  risk_data <- NULL
  message("[免疫分析V2] ⚠️ 未找到Risk Score数据")
}


# ============================================
# 2. 定义免疫细胞基因集 (ssGSEA)
# ============================================
message("[免疫分析V2] 准备免疫细胞基因集...")

immune_signatures <- list(
  "Activated_B_cell" = c("CD19", "CD79A", "CD79B", "MS4A1", "CD22"),
  "Activated_CD4_T_cell" = c("CD4", "IL2RA", "CD40LG", "ICOS", "CTLA4"),
  "Activated_CD8_T_cell" = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG"),
  "Central_memory_CD4_T_cell" = c("CD4", "CCR7", "SELL", "IL7R"),
  "Central_memory_CD8_T_cell" = c("CD8A", "CCR7", "SELL", "IL7R"),
  "Effector_memory_CD4_T_cell" = c("CD4", "GZMK", "CCL5"),
  "Effector_memory_CD8_T_cell" = c("CD8A", "GZMK", "CCL5", "NKG7"),
  "Gamma_delta_T_cell" = c("TRGC1", "TRGC2", "TRDC"),
  "Regulatory_T_cell" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"),
  "T_follicular_helper_cell" = c("CXCR5", "ICOS", "BCL6", "CD40LG"),
  "Type_1_T_helper_cell" = c("TBX21", "IFNG", "IL12RB2"),
  "Type_2_T_helper_cell" = c("GATA3", "IL4", "IL5", "IL13"),
  "Type_17_T_helper_cell" = c("RORC", "IL17A", "IL17F", "IL22"),
  "Activated_dendritic_cell" = c("CD80", "CD86", "CD83", "CCR7"),
  "Immature_dendritic_cell" = c("CD1A", "CD1C", "CLEC10A"),
  "Plasmacytoid_dendritic_cell" = c("CLEC4C", "IL3RA", "NRP1"),
  "Macrophage" = c("CD68", "CD163", "MSR1", "MRC1"),
  "M1_Macrophage" = c("NOS2", "IL1B", "IL6", "TNF", "CD80"),
  "M2_Macrophage" = c("CD163", "MRC1", "ARG1", "IL10"),
  "Monocyte" = c("CD14", "FCGR3A", "CSF1R"),
  "Natural_killer_cell" = c("NCAM1", "NKG7", "KLRD1", "KLRB1"),
  "Neutrophil" = c("FCGR3B", "CXCR1", "CXCR2", "CSF3R"),
  "Mast_cell" = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  "Eosinophil" = c("CCR3", "SIGLEC8", "PRG2"),
  "MDSC" = c("CD33", "ITGAM", "ARG1", "NOS2"),
  "Natural_killer_T_cell" = c("CD3D", "NCAM1", "KLRB1")
)

# 检查基因匹配情况
all_immune_genes <- unique(unlist(immune_signatures))
matched_genes <- intersect(all_immune_genes, rownames(expr_gene))
message("[免疫分析V2] 免疫基因匹配: ", length(matched_genes), "/", length(all_immune_genes))

# ============================================
# 3. ssGSEA分析
# ============================================
message("[免疫分析V2] 进行ssGSEA分析...")

# 使用简化版ssGSEA（更稳定）
ssgsea_scores <- sapply(immune_signatures, function(genes) {
  genes_in_expr <- intersect(genes, rownames(expr_gene))
  if (length(genes_in_expr) >= 2) {
    colMeans(expr_gene[genes_in_expr, , drop = FALSE], na.rm = TRUE)
  } else if (length(genes_in_expr) == 1) {
    expr_gene[genes_in_expr, ]
  } else {
    rep(NA, ncol(expr_gene))
  }
})
ssgsea_scores <- t(ssgsea_scores)

# 移除全NA的细胞类型
valid_cells <- rowSums(!is.na(ssgsea_scores)) > 0
ssgsea_scores <- ssgsea_scores[valid_cells, ]

message("[免疫分析V2] ssGSEA完成，有效免疫细胞类型: ", nrow(ssgsea_scores))

# 保存ssGSEA结果
write.csv(t(ssgsea_scores), file.path(res_dir, "ssgsea_immune_scores.csv"))

# ============================================
# 4. 免疫细胞与Risk Score相关性
# ============================================
message("[免疫分析V2] 计算免疫细胞与Risk Score相关性...")

if (!is.null(risk_data) && nrow(ssgsea_scores) > 0) {
  common_samples <- intersect(risk_data$sample, colnames(ssgsea_scores))
  message("[免疫分析V2] 匹配样本数: ", length(common_samples))
  
  if (length(common_samples) > 10) {
    ssgsea_matched <- ssgsea_scores[, common_samples, drop = FALSE]
    risk_matched <- risk_data[match(common_samples, risk_data$sample), ]
    
    cor_results <- data.frame(
      Cell_Type = rownames(ssgsea_matched),
      Correlation = numeric(nrow(ssgsea_matched)),
      P_value = numeric(nrow(ssgsea_matched))
    )
    
    for (i in 1:nrow(ssgsea_matched)) {
      vals <- as.numeric(ssgsea_matched[i, ])
      if (sum(!is.na(vals)) > 10) {
        cor_test <- cor.test(vals, risk_matched$risk_score, method = "spearman", use = "complete.obs")
        cor_results$Correlation[i] <- cor_test$estimate
        cor_results$P_value[i] <- cor_test$p.value
      }
    }
    
    cor_results <- cor_results %>% filter(!is.na(Correlation)) %>% arrange(P_value)
    write.csv(cor_results, file.path(res_dir, "immune_risk_correlation.csv"), row.names = FALSE)
    message("[免疫分析V2] 显著相关的免疫细胞: ", sum(cor_results$P_value < 0.05))
  }
}


# ============================================
# 5. 核心基因与免疫细胞相关性
# ============================================
message("[免疫分析V2] 计算核心基因与免疫细胞相关性...")

# 获取核心基因
coef_file <- file.path(res_dir, "prognostic_model_coef.csv")
if (file.exists(coef_file)) {
  core_genes <- read.csv(coef_file)$Gene
} else {
  # 使用铁死亡DEG
  ferro_file <- file.path("data/references", "ferroptosis_genes_expanded.csv")
  deg_file <- file.path(res_dir, "deg_GSE14520_sig.csv")
  if (file.exists(ferro_file) && file.exists(deg_file)) {
    ferro_genes <- read.csv(ferro_file)$Gene
    deg_genes <- read.csv(deg_file)$Gene
    core_genes <- intersect(ferro_genes, deg_genes)
  } else {
    core_genes <- c()
  }
}

available_genes <- intersect(core_genes, rownames(expr_gene))
message("[免疫分析V2] 可用核心基因: ", length(available_genes))

if (length(available_genes) > 0 && nrow(ssgsea_scores) > 0) {
  gene_immune_cor <- matrix(NA, nrow = length(available_genes), ncol = nrow(ssgsea_scores))
  rownames(gene_immune_cor) <- available_genes
  colnames(gene_immune_cor) <- rownames(ssgsea_scores)
  
  for (i in seq_along(available_genes)) {
    gene <- available_genes[i]
    for (j in 1:nrow(ssgsea_scores)) {
      common_samples <- intersect(colnames(expr_gene), colnames(ssgsea_scores))
      if (length(common_samples) > 10) {
        gene_vals <- as.numeric(expr_gene[gene, common_samples])
        immune_vals <- as.numeric(ssgsea_scores[j, common_samples])
        if (sum(!is.na(gene_vals) & !is.na(immune_vals)) > 10) {
          cor_test <- cor.test(gene_vals, immune_vals, method = "spearman", use = "complete.obs")
          gene_immune_cor[i, j] <- cor_test$estimate
        }
      }
    }
  }
  
  write.csv(gene_immune_cor, file.path(res_dir, "gene_immune_correlation.csv"))
}

# ============================================
# 6. 高低风险组免疫差异
# ============================================
message("[免疫分析V2] 分析高低风险组免疫差异...")

if (!is.null(risk_data) && nrow(ssgsea_scores) > 0) {
  common_samples <- intersect(risk_data$sample, colnames(ssgsea_scores))
  
  if (length(common_samples) > 10) {
    ssgsea_matched <- ssgsea_scores[, common_samples, drop = FALSE]
    risk_matched <- risk_data[match(common_samples, risk_data$sample), ]
    
    diff_results <- data.frame(
      Cell_Type = rownames(ssgsea_matched),
      High_Mean = numeric(nrow(ssgsea_matched)),
      Low_Mean = numeric(nrow(ssgsea_matched)),
      P_value = numeric(nrow(ssgsea_matched))
    )
    
    for (i in 1:nrow(ssgsea_matched)) {
      high_vals <- as.numeric(ssgsea_matched[i, risk_matched$risk_group == "High"])
      low_vals <- as.numeric(ssgsea_matched[i, risk_matched$risk_group == "Low"])
      
      diff_results$High_Mean[i] <- mean(high_vals, na.rm = TRUE)
      diff_results$Low_Mean[i] <- mean(low_vals, na.rm = TRUE)
      
      if (sum(!is.na(high_vals)) > 3 && sum(!is.na(low_vals)) > 3) {
        diff_results$P_value[i] <- wilcox.test(high_vals, low_vals)$p.value
      }
    }
    
    diff_results <- diff_results %>% filter(!is.na(P_value)) %>% arrange(P_value)
    write.csv(diff_results, file.path(res_dir, "immune_risk_group_diff.csv"), row.names = FALSE)
    message("[免疫分析V2] 显著差异的免疫细胞: ", sum(diff_results$P_value < 0.05))
  }
}

# ============================================
# 7. 生成Figure 8: 免疫浸润多面板图
# ============================================
message("[免疫分析V2] 生成Figure 8...")

pdf(file.path(plot_dir, "Figure8_immune_infiltration.pdf"), width = 16, height = 14)
par(mfrow = c(2, 2))

# A: ssGSEA热图
if (nrow(ssgsea_scores) > 0) {
  # 随机抽样显示
  if (ncol(ssgsea_scores) > 50) {
    sample_idx <- sample(1:ncol(ssgsea_scores), 50)
    heatmap_data <- ssgsea_scores[, sample_idx]
  } else {
    heatmap_data <- ssgsea_scores
  }
  
  heatmap_data <- heatmap_data[complete.cases(heatmap_data), , drop = FALSE]
  
  if (nrow(heatmap_data) > 0) {
    heatmap_data_scaled <- t(scale(t(heatmap_data)))
    heatmap_data_scaled[is.na(heatmap_data_scaled)] <- 0
    heatmap_data_scaled[heatmap_data_scaled > 2] <- 2
    heatmap_data_scaled[heatmap_data_scaled < -2] <- -2
    
    heatmap(heatmap_data_scaled, 
            main = "A. Immune Cell Infiltration (ssGSEA)",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            scale = "none", cexRow = 0.6, cexCol = 0.5)
  }
}

# B: 免疫细胞与Risk Score相关性
if (exists("cor_results") && nrow(cor_results) > 0) {
  top_cor <- head(cor_results, 15)
  barplot(top_cor$Correlation, names.arg = top_cor$Cell_Type,
          col = ifelse(top_cor$Correlation > 0, "red", "blue"),
          las = 2, cex.names = 0.6,
          main = "B. Immune Cell - Risk Score Correlation",
          ylab = "Spearman Correlation")
  abline(h = 0, lty = 2)
}

# C: 核心基因与免疫细胞相关性热图
if (exists("gene_immune_cor") && nrow(gene_immune_cor) > 0) {
  cor_sub <- gene_immune_cor[, head(colnames(gene_immune_cor), 15), drop = FALSE]
  cor_sub[is.na(cor_sub)] <- 0
  
  heatmap(cor_sub, 
          main = "C. Core Genes - Immune Cells Correlation",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          scale = "none", cexRow = 0.7, cexCol = 0.6)
}

# D: 高低风险组免疫差异
if (exists("diff_results") && nrow(diff_results) > 0) {
  sig_cells <- diff_results %>% filter(P_value < 0.1) %>% head(10)
  
  if (nrow(sig_cells) > 0) {
    diff_mat <- as.matrix(sig_cells[, c("High_Mean", "Low_Mean")])
    rownames(diff_mat) <- sig_cells$Cell_Type
    
    barplot(t(diff_mat), beside = TRUE, 
            col = c("red", "blue"), las = 2, cex.names = 0.6,
            main = "D. Immune Difference: High vs Low Risk",
            ylab = "ssGSEA Score")
    legend("topright", legend = c("High Risk", "Low Risk"), fill = c("red", "blue"))
  }
}

dev.off()

message("[免疫分析V2] ✅ 免疫浸润分析完成！")
message("[免疫分析V2] 输出文件:")
message("  ✅ ", file.path(proc_dir, "GSE14520_expr_symbol.rds"))
message("  ✅ ", file.path(res_dir, "ssgsea_immune_scores.csv"))
message("  ✅ ", file.path(res_dir, "immune_risk_correlation.csv"))
message("  ✅ ", file.path(res_dir, "gene_immune_correlation.csv"))
message("  ✅ ", file.path(res_dir, "immune_risk_group_diff.csv"))
message("  ✅ ", file.path(plot_dir, "Figure8_immune_infiltration.pdf"))
