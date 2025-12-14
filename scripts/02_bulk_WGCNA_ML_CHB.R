#!/usr/bin/env Rscript

# 02_bulk_WGCNA_ML_CHB.R
# CHB Bulk DEG + ferroptosis 交集 + WGCNA + 机器学习筛选核心基因
# 使用方法: Rscript scripts/02_bulk_WGCNA_ML_CHB.R

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(WGCNA)
  library(caret)
  library(randomForest)
  library(glmnet)
  library(pROC)
  library(VennDiagram)
  library(grid)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# 禁用VennDiagram的日志输出
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

options(stringsAsFactors = FALSE)

# 设置工作目录（支持从项目根目录运行）
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
ref_dir  <- "data/references"
res_dir  <- "results"
plot_dir <- "plots"

dir.create(res_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1. 选择主要训练队列（CHB 肝组织）
# ---------------------------

train_gse <- "GSE83148"  # 与 01_download_data_CHB.R 中保持一致（可按需修改）
exprs_log <- readRDS(file.path(proc_dir, paste0(train_gse, "_exprs_log2.rds")))
pdata     <- readRDS(file.path(proc_dir, paste0(train_gse, "_pdata.rds")))

message("pdata 中的列名如下，请在第一次运行时注意查看：")
message(paste(colnames(pdata), collapse = ", "))

# 从 pdata 自动推断 Group（Control vs CHB），逻辑：
# 1）优先使用字符型列中的 "title"；
# 2）若无 title，则在所有字符型列中选择最可能包含疾病/对照信息的一列；
# 3）在该列中，包含 control/normal/healthy 归为 Control，其余归为 CHB。

# 仅保留字符型列用于模式匹配
char_cols <- sapply(pdata, is.character)
pdata_char <- pdata[, char_cols, drop = FALSE]

if (ncol(pdata_char) == 0) {
  stop("pdata 中没有字符型列，无法自动推断 Group，请手动在脚本中构造 Group 向量。")
}

if ("title" %in% colnames(pdata_char)) {
  src_col <- "title"
  src_vec <- pdata_char[["title"]]
} else {
  # 在所有字符列中，按是否包含 HBV/hepatitis/CHB/control/healthy 等关键词打分，选分最高的一列
  score_col <- function(x) {
    sum(grepl("HBV|hepatitis|chronic hepatitis b|CHB|control|healthy|normal", x, ignore.case = TRUE), na.rm = TRUE)
  }
  scores <- sapply(pdata_char, score_col)
  src_col <- names(which.max(scores))[1]
  src_vec <- pdata_char[[src_col]]
}

message("自动选择列用于构造 Group：", src_col)

grp <- ifelse(grepl("control|normal|healthy", src_vec, ignore.case = TRUE),
              "Control", "CHB")
Group <- factor(grp, levels = c("Control", "CHB"))

message("自动推断的分组频数如下（请核对，如有不符请在脚本中手动修改 Group 构造部分）：")
print(table(Group))

# 保存样本分组信息，供免疫浸润等后续脚本使用
sample_group <- data.frame(Sample = rownames(pdata), Group = Group)
saveRDS(sample_group, file.path(res_dir, "sample_group_CHB.rds"))

# ---------------------------
# 2. 差异分析 (Limma) + 铁死亡交集
# ---------------------------

# 使用扩展版铁死亡基因集 (108个基因，来自FerrDb+KEGG+MSigDB)
ferro_genes <- read.csv(file.path(ref_dir, "ferroptosis_genes_expanded.csv"))$Gene

# 保证表达矩阵列名与 pdata 行名一致顺序
exprs_log <- exprs_log[, sample_group$Sample]

design <- model.matrix(~0 + Group)
colnames(design) <- levels(Group)
fit <- lmFit(exprs_log, design)
contrast.matrix <- makeContrasts(CHB - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))

deg_all <- topTable(fit2, adjust.method = "fdr", number = Inf)

# 添加基因符号列 (探针ID到基因符号转换)
message("[DEG] 转换探针ID到基因符号...")
tryCatch({
  library(hgu133plus2.db)
  deg_all$ProbeID <- rownames(deg_all)
  deg_all$Gene <- mapIds(hgu133plus2.db, 
                         keys = rownames(deg_all),
                         column = "SYMBOL",
                         keytype = "PROBEID",
                         multiVals = "first")
  # 过滤掉没有基因符号的探针
  deg_all <- deg_all %>% filter(!is.na(Gene) & Gene != "")
  message("[DEG] 成功转换 ", nrow(deg_all), " 个探针到基因符号")
}, error = function(e) {
  message("[DEG] 注释包转换失败: ", e$message)
  deg_all$ProbeID <- rownames(deg_all)
  deg_all$Gene <- rownames(deg_all)
})

write.csv(deg_all, file.path(res_dir, "bulk_deg_all_CHB.csv"), row.names = FALSE)

# 筛选显著 DEGs
sig_deg <- deg_all %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(sig_deg, file.path(res_dir, "bulk_deg_sig_CHB.csv"), row.names = FALSE)

message("[DEG] 显著DEG数: ", nrow(sig_deg))

# 与铁死亡基因交集 (使用基因符号)
ferro_deg <- sig_deg %>% filter(Gene %in% ferro_genes)
write.csv(ferro_deg, file.path(res_dir, "ferroptosis_deg_CHB.csv"), row.names = FALSE)

message("[DEG] 铁死亡相关DEG数: ", nrow(ferro_deg))

# 简单火山图（铁死亡基因高亮）
ferro_flag <- rownames(deg_all) %in% ferro_genes

pdf(file.path(plot_dir, "bulk_volcano_ferroptosis_CHB.pdf"), width = 6, height = 5)
plot(deg_all$logFC, -log10(deg_all$adj.P.Val), pch = 20,
     col = ifelse(ferro_flag, "red", "grey"),
     xlab = "log2FC (CHB vs Control)",
     ylab = "-log10(FDR)",
     main = "CHB vs Control (Ferroptosis genes in red)")
abline(h = -log10(0.05), v = c(-1, 1), lty = 2, col = "blue")
dev.off()

# ---------------------------
# 3. WGCNA 寻找与 CHB 状态相关模块
# ---------------------------

enableWGCNAThreads()

# 重要：使用WGCNA的cor函数，避免与stats包冲突
cor <- WGCNA::cor

# WGCNA 要求：行为样本，列为基因
expr_t <- t(exprs_log)

# 取方差最高的前 5000 基因简化计算
vars <- apply(expr_t, 2, var)
expr_t_filt <- expr_t[, order(vars, decreasing = TRUE)[1:min(5000, ncol(expr_t))]]

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(expr_t_filt, powerVector = powers, verbose = 5, corFnc = "cor")
softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

net <- blockwiseModules(expr_t_filt, power = softPower,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 3,
                        corType = "pearson")

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# 恢复stats包的cor函数
cor <- stats::cor

traits <- data.frame(CHB_status = as.numeric(Group == "CHB"))
rownames(traits) <- rownames(expr_t_filt)

moduleTraitCor <- WGCNA::cor(MEs, traits, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(expr_t_filt))

saveRDS(net, file.path(proc_dir, "wgcna_net_CHB.rds"))

pdf(file.path(plot_dir, "wgcna_module_trait_heatmap_CHB.pdf"), width = 8, height = 6)
par(mar = c(6, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleTraitP, 2),
               main = "Module–trait relationships (CHB status)")
dev.off()

# 选取与 CHB 状态相关性最高的模块
module_of_interest <- names(which.max(abs(moduleTraitCor[, "CHB_status"])))
module_genes <- names(which(net$colors == gsub("ME", "", module_of_interest)))

# 保存模块基因，使用标准列名 "Gene"
write.csv(data.frame(Gene = module_genes), file.path(res_dir, "wgcna_module_genes_CHB.csv"), row.names = FALSE)

# 与铁死亡 DEGs 进一步交集，得到候选核心基因池
candidate_genes <- intersect(module_genes, rownames(ferro_deg))
write.csv(candidate_genes, file.path(res_dir, "candidate_genes_for_ML_CHB.csv"), row.names = FALSE)

# ---------------------------
# 4. 机器学习筛选 Diagnostic Core Genes
# ---------------------------

if (length(candidate_genes) < 5) {
  warning("候选基因太少，无法稳定运行 ML。请适当放宽 DEG/WGCNA 筛选阈值。")
} else {
  expr_sub <- t(exprs_log[candidate_genes, ])  # 样本 x 基因
  expr_sub <- scale(expr_sub)
  expr_sub <- as.data.frame(expr_sub)
  expr_sub$Group <- Group

  set.seed(123)
  train_index <- createDataPartition(expr_sub$Group, p = 0.7, list = FALSE)
  train_dat <- expr_sub[train_index, ]
  test_dat  <- expr_sub[-train_index, ]

  # SVM-RFE
  ctrl_rfe <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
  svm_rfe <- rfe(x = train_dat[, candidate_genes, drop = FALSE],
                 y = train_dat$Group,
                 sizes = 1:min(20, length(candidate_genes)),
                 rfeControl = ctrl_rfe,
                 method = "svmRadial")
  genes_svm <- predictors(svm_rfe)

  # Random Forest 特征重要性
  rf_model <- randomForest(Group ~ ., data = train_dat[, c(candidate_genes, "Group")])
  imp <- importance(rf_model)
  genes_rf <- rownames(imp)[order(imp[, 1], decreasing = TRUE)[1:min(15, nrow(imp))]]

  # LASSO
  x_mat <- as.matrix(train_dat[, candidate_genes])
  y_bin <- ifelse(train_dat$Group == "CHB", 1, 0)
  cv_lasso <- cv.glmnet(x_mat, y_bin, family = "binomial", alpha = 1)
  coef_lasso <- coef(cv_lasso, s = cv_lasso$lambda.min)
  genes_lasso <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
  genes_lasso <- setdiff(genes_lasso, "(Intercept)")

  venn_list <- list(SVM = genes_svm, RF = genes_rf, LASSO = genes_lasso)
  core_genes <- Reduce(intersect, venn_list)

  pdf(file.path(plot_dir, "ML_feature_venn_CHB.pdf"))
  draw.triple.venn(area1 = length(genes_svm), area2 = length(genes_rf), area3 = length(genes_lasso),
                   n12 = length(intersect(genes_svm, genes_rf)),
                   n23 = length(intersect(genes_rf, genes_lasso)),
                   n13 = length(intersect(genes_svm, genes_lasso)),
                   n123 = length(core_genes),
                   category = c("SVM-RFE", "RF", "LASSO"),
                   fill = c("pink", "lightblue", "lightgreen"))
  dev.off()

  if (length(core_genes) > 0) {
    write.csv(core_genes, file.path(res_dir, "diagnostic_core_genes_CHB.csv"), row.names = FALSE)

    # 逻辑回归模型 + ROC
    train_glm <- glm(Group ~ ., data = train_dat[, c(core_genes, "Group")], family = binomial)
    prob_train <- predict(train_glm, type = "response")
    roc_train <- roc(train_dat$Group, prob_train)

    prob_test <- predict(train_glm, newdata = test_dat, type = "response")
    roc_test <- roc(test_dat$Group, prob_test)

    pdf(file.path(plot_dir, "ML_core_genes_ROC_train_test_CHB.pdf"), width = 6, height = 6)
    plot(roc_train, col = "red", main = "Diagnostic model ROC (CHB)")
    plot(roc_test, col = "blue", add = TRUE)
    legend("bottomright", legend = c(paste0("Train AUC=", round(auc(roc_train), 3)),
                                      paste0("Test AUC=", round(auc(roc_test), 3))),
           col = c("red", "blue"), lwd = 2)
    dev.off()
  } else {
    warning("三种算法交集为空，可考虑放宽筛选标准或用并集 + 手动挑选。")
  }
}

message("CHB Bulk + WGCNA + ML 分析完成。")
