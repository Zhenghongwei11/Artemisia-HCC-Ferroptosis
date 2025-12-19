#!/usr/bin/env Rscript

# 04_scrna_CHB.R
# 肝脏单细胞转录组分析：定位 CHB 核心铁死亡基因的细胞类型与拟时序变化
# 
# ⚠️ 注意：此脚本为可选分析，需要大量内存（建议16GB+）
# 如果内存不足，可以跳过此脚本
# 使用方法: Rscript scripts/04_scrna_CHB.R

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(monocle)
  library(Matrix)
})

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

raw_sc_dir   <- "data/raw/scRNA_CHB"   # 需手动放置肝脏单细胞数据
proc_dir     <- "data/processed"
plot_dir     <- "plots"
ref_dir      <- "data/references"

dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1. 读取单细胞数据
# ---------------------------

if (!dir.exists(raw_sc_dir)) {
  message("⚠️ 未找到 data/raw/scRNA_CHB 目录")
  message("单细胞分析为可选步骤，跳过...")
  message("如需运行，请下载肝脏单细胞数据并放入该目录")
  quit(save = "no", status = 0)
}

files_h5 <- list.files(raw_sc_dir, pattern = "\\.h5$", full.names = TRUE)

if (length(files_h5) > 0) {
  message("检测到 .h5 文件，使用 Read10X_h5 读取...")
  counts <- Read10X_h5(files_h5[1])
} else {
  message("未检测到 .h5 文件，尝试按 10X matrix 目录结构读取...")
  counts <- Read10X(data.dir = raw_sc_dir)
}

seu <- CreateSeuratObject(counts = counts, project = "Liver_CHB_like", min.cells = 3, min.features = 200)

# 线粒体比例（人类基因一般以 "MT-" 开头）
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# 简单 QC 过滤，可按数据量微调
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# ---------------------------
# 2. 归一化、降维与聚类
# ---------------------------

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(seu))
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

p_umap_cluster <- DimPlot(seu, reduction = "umap", label = TRUE) +
  ggtitle("Liver single-cell UMAP (clusters)")

ggsave(file.path(plot_dir, "scCHB_umap_clusters.pdf"), p_umap_cluster, width = 7, height = 6)

# ---------------------------
# 3. （占位）细胞类型注释
# ---------------------------
# 不同数据集的 marker 不同，这里只提供一个示意：
# 实际使用时建议：
#   markers <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25)
# 根据已知肝细胞 / Kupffer / stellate / 内皮等 marker 手动注释。

# 为了保证脚本可运行，这里仅按 cluster ID 命名：
new_ids <- paste0("Cluster_", levels(seu))
names(new_ids) <- levels(seu)
seu <- RenameIdents(seu, new_ids)
seu$cell_type <- Idents(seu)

p_umap_celltype <- DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE) +
  ggtitle("Liver single-cell UMAP (cell types placeholder)")

ggsave(file.path(plot_dir, "scCHB_umap_celltypes_placeholder.pdf"), p_umap_celltype, width = 7, height = 6)

# ---------------------------
# 4. 铁死亡模块评分（含 TFRC 等基因）
# ---------------------------

# 使用扩展版铁死亡基因集 (108个基因)
ferro_file <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")
if (!file.exists(ferro_file)) {
  stop("找不到 ferroptosis_genes_expanded.csv，请先运行 00c_download_ferrdb_genes.R")
}

ferro_genes <- read.csv(ferro_file)$Gene
valid_ferro <- intersect(ferro_genes, rownames(seu))

if (length(valid_ferro) < 3) {
  warning("在单细胞数据中找到的铁死亡基因太少，AddModuleScore 结果可能不稳定。")
}

seu <- AddModuleScore(seu, features = list(valid_ferro), name = "Ferroptosis_Score")

p_umap_ferro <- FeaturePlot(seu, features = "Ferroptosis_Score1") +
  ggtitle("Ferroptosis module score")

ggsave(file.path(plot_dir, "scCHB_umap_ferroptosis_score.pdf"), p_umap_ferro, width = 7, height = 6)

p_vln_ferro <- VlnPlot(seu, features = "Ferroptosis_Score1", group.by = "cell_type", pt.size = 0) +
  ggtitle("Ferroptosis score across cell types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "scCHB_violin_ferroptosis_score.pdf"), p_vln_ferro, width = 8, height = 6)

# 如 TFRC 等核心基因存在，额外画 FeaturePlot
core_genes_file <- file.path("../results", "diagnostic_core_genes_CHB.csv")
if (file.exists(core_genes_file)) {
  core_genes <- read.csv(core_genes_file, header = TRUE)[, 1]
  core_in_sc <- intersect(core_genes, rownames(seu))
  if (length(core_in_sc) > 0) {
    for (g in core_in_sc) {
      p_g <- FeaturePlot(seu, features = g) + ggtitle(paste0(g, " expression"))
      ggsave(file.path(plot_dir, paste0("scCHB_feature_", g, ".pdf")), p_g, width = 7, height = 6)
    }
  }
}

# ---------------------------
# 5. 拟时序分析（以某一细胞类型为例）
# ---------------------------
# 这里示范将当前全部细胞做 Monocle2 拟时序，实际可只选肝细胞或巨噬细胞子集。

expr_mat <- as(as.matrix(GetAssayData(seu, slot = "counts")), "sparseMatrix")
pheno    <- new("AnnotatedDataFrame", data = seu@meta.data)
feature  <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
feature  <- new("AnnotatedDataFrame", data = feature)

cds <- newCellDataSet(expr_mat,
                      phenoData       = pheno,
                      featureData     = feature,
                      lowerDetectionLimit = 0.5,
                      expressionFamily    = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 简化：用变异最大的基因作为排序基因
var_genes <- apply(expr_mat, 1, var)
ordering_genes <- names(sort(var_genes, decreasing = TRUE))[1:min(1000, length(var_genes))]

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)

pdf(file.path(plot_dir, "scCHB_monocle_trajectory.pdf"), width = 7, height = 6)
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

# 如果 TFRC 或其他核心基因在单细胞中存在，绘制随拟时序变化的趋势
ferro_interest <- intersect(c("TFRC", ferro_genes), rownames(cds))

if (length(ferro_interest) > 0) {
  pdf(file.path(plot_dir, "scCHB_monocle_ferroptosis_genes_trend.pdf"), width = 7, height = 6)
  plot_genes_in_pseudotime(cds[ferro_interest, ], color_by = "cell_type")
  dev.off()
}

saveRDS(seu, file.path(proc_dir, "seurat_liver_CHB_like.rds"))

message("单细胞（肝脏图谱）分析完成：UMAP / 铁死亡评分 / 拟时序图已生成。")
