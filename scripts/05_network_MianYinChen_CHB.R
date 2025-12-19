#!/usr/bin/env Rscript

# 05_network_MianYinChen_CHB.R
# 绵茵陈网络药理与 CHB 铁死亡/免疫核心基因的整合分析
# 使用方法: Rscript scripts/05_network_MianYinChen_CHB.R

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(VennDiagram)
  library(grid)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pheatmap)
})

# 禁用VennDiagram的日志输出
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

ref_dir  <- "data/references"
res_dir  <- "results"
plot_dir <- "plots"

dir.create(res_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1. 读取基因集：HCC DEGs / Ferroptosis / MianYinChen 靶点
# ---------------------------

# (1) HCC 差异基因（来自 02b_multi_cohort_DEG.R，GSE14520主队列）
# 按照调整后的故事线，以GSE14520(HCC)为主队列

deg_file <- file.path(res_dir, "deg_GSE14520_sig.csv")
if (!file.exists(deg_file)) {
  # 备选：使用旧的CHB DEG文件
  deg_file <- file.path(res_dir, "bulk_deg_sig_CHB.csv")
  if (!file.exists(deg_file)) {
    stop("找不到 deg_GSE14520_sig.csv 或 bulk_deg_sig_CHB.csv，请先运行 02b_multi_cohort_DEG.R")
  }
  message("注意: 使用备选文件 bulk_deg_sig_CHB.csv")
}

deg_df <- read.csv(deg_file, check.names = FALSE)
# 获取基因名
if ("Gene" %in% colnames(deg_df)) {
  deg_genes <- deg_df$Gene
} else if ("X" %in% colnames(deg_df)) {
  deg_genes <- deg_df$X
} else {
  deg_genes <- rownames(deg_df)
}
deg_genes <- unique(deg_genes[!is.na(deg_genes) & deg_genes != ""])
message("加载DEG基因: ", length(deg_genes), " 个 (来自 ", basename(deg_file), ")")

# (2) 铁死亡基因

# 使用扩展版铁死亡基因集 (108个基因，来自FerrDb+KEGG+MSigDB)
ferro_file <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")
if (!file.exists(ferro_file)) {
  stop("找不到 ferroptosis_genes_expanded.csv，请先运行 00c_download_ferrdb_genes.R")
}
ferro_genes <- read.csv(ferro_file)$Gene

# (3) 绵茵陈靶点（需要你从 TCMSP 下载并整理为 Ingredient,Target 两列）

# 绵茵陈靶点数据文件
# 来源: TCMSP数据库 (https://tcmsp-e.com/) 或 HERB数据库 (http://herb.ac.cn/)
# 用户需要从数据库下载绵茵陈(Artemisia capillaris)的活性成分及其靶点
# 文件格式要求: CSV文件，包含 Ingredient 和 Target 两列
# 
# 重要: 根据项目spec规则(需求11.2, 11.4)，禁止使用虚拟数据！
# 必须从权威数据库获取真实数据。
tcm_file <- file.path(ref_dir, "tcm_targets_CHB.csv")
if (!file.exists(tcm_file)) {
  stop(paste0(
    "\n",
    "================================================================================\n",
    "错误: 未找到绵茵陈靶点数据文件\n",
    "================================================================================\n",
    "\n",
    "根据项目规范(需求11.2, 11.4)，必须从权威数据库获取真实数据。\n",
    "禁止使用虚拟数据或硬编码数据！\n",
    "\n",
    "请按以下步骤准备数据:\n",
    "\n",
    "【方法1: TCMSP数据库 (推荐)】\n",
    "1. 访问 https://tcmsp-e.com/\n",
    "2. 搜索 '茵陈' 或 'Artemisia capillaris' 或 'Yinchen'\n",
    "3. 筛选活性成分: OB ≥ 30%, DL ≥ 0.18\n",
    "4. 下载成分-靶点数据\n",
    "5. 整理为CSV格式\n",
    "\n",
    "【方法2: HERB数据库】\n",
    "1. 访问 http://herb.ac.cn/\n",
    "2. 搜索 '茵陈' 或 'Herba Artemisiae Scopariae'\n",
    "3. 下载成分-靶点数据\n",
    "\n",
    "【文件格式要求】\n",
    "CSV文件，包含两列:\n",
    "  - Ingredient: 活性成分英文名 (如 Scoparone, Capillarisin)\n",
    "  - Target: 靶点基因符号 (HGNC标准，如 PTGS2, TNF)\n",
    "\n",
    "【保存位置】\n",
    tcm_file, "\n",
    "\n",
    "【引用要求】\n",
    "论文中需引用: Ru J, et al. TCMSP: a database of systems pharmacology\n",
    "for drug discovery from herbal medicines. J Cheminform. 2014;6:13.\n",
    "\n",
    "详细指南请参考: data/references/README_DATA_SOURCES.md\n",
    "================================================================================\n"
  ))
}

tcm_targets <- read.csv(tcm_file)
if (!all(c("Ingredient", "Target") %in% colnames(tcm_targets))) {
  stop("tcm_targets_CHB.csv 必须包含列: Ingredient, Target")
}

tcm_genes <- unique(tcm_targets$Target)

# ---------------------------
# 2. 三集合交集分析 & Venn 图
# ---------------------------

set_deg    <- unique(deg_genes)
set_tcm    <- unique(tcm_genes)
set_ferro  <- unique(ferro_genes)

inter_deg_tcm      <- intersect(set_deg, set_tcm)
inter_tcm_ferro    <- intersect(set_tcm, set_ferro)
inter_deg_ferro    <- intersect(set_deg, set_ferro)
inter_all_three    <- Reduce(intersect, list(set_deg, set_tcm, set_ferro))

message("Disease ∩ TCM targets: ", length(inter_deg_tcm))
message("TCM targets ∩ Ferroptosis: ", length(inter_tcm_ferro))
message("CHB DEGs ∩ Ferroptosis: ", length(inter_deg_ferro))
message("三者交集 (Drug–Disease–Ferroptosis hub genes): ", length(inter_all_three))

pdf(file.path(plot_dir, "venn_CHB_DEG_TCM_Ferroptosis.pdf"))
draw.triple.venn(
  area1 = length(set_deg),
  area2 = length(set_tcm),
  area3 = length(set_ferro),
  n12   = length(intersect(set_deg, set_tcm)),
  n23   = length(intersect(set_tcm, set_ferro)),
  n13   = length(intersect(set_deg, set_ferro)),
  n123  = length(inter_all_three),
  category = c("CHB DEGs", "MianYinChen targets", "Ferroptosis"),
  fill = c("skyblue", "pink", "lightgreen")
)
dev.off()

# ---------------------------
# 3. 富集分析：Drug–Disease 交集 & 三者交集
# ---------------------------

# 先对 (CHB DEGs ∩ TCM) 做 GO/KEGG，展示“中药作用靶点”的整体功能

if (length(inter_deg_tcm) > 0) {
  eg_deg_tcm <- bitr(inter_deg_tcm, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  if (!is.null(eg_deg_tcm) && nrow(eg_deg_tcm) > 0) {
    ego <- enrichGO(gene         = eg_deg_tcm$ENTREZID,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    readable      = TRUE)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      pdf(file.path(plot_dir, "GO_enrichment_CHB_DEG_TCM.pdf"), width = 8, height = 6)
      print(dotplot(ego, showCategory = 20) + ggtitle("GO (BP) of CHB DEGs ∩ MianYinChen targets"))
      dev.off()
    }

    ekegg <- tryCatch({
      enrichKEGG(gene         = eg_deg_tcm$ENTREZID,
                 organism     = "hsa",
                 pvalueCutoff = 0.05)
    }, error = function(e) {
      message("KEGG富集分析失败（可能是网络问题）: ", e$message)
      NULL
    })
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      pdf(file.path(plot_dir, "KEGG_enrichment_CHB_DEG_TCM.pdf"), width = 8, height = 6)
      print(dotplot(ekegg, showCategory = 20) + ggtitle("KEGG of CHB DEGs and MianYinChen targets"))
      dev.off()
    }
  }
}

# 对三者交集(Drug–Disease–Ferroptosis hub genes)再做一次富集，突出铁死亡/免疫通路

if (length(inter_all_three) > 0) {
  eg_all <- bitr(inter_all_three, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (!is.null(eg_all) && nrow(eg_all) > 0) {
    ego_all <- enrichGO(gene = eg_all$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont   = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        readable      = TRUE)
    if (!is.null(ego_all) && nrow(as.data.frame(ego_all)) > 0) {
      pdf(file.path(plot_dir, "GO_enrichment_hub_threeway.pdf"), width = 8, height = 6)
      print(dotplot(ego_all, showCategory = 20) + ggtitle("GO (BP) of Drug–Disease–Ferroptosis hub genes"))
      dev.off()
    }

    ekegg_all <- tryCatch({
      enrichKEGG(gene = eg_all$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff = 0.05)
    }, error = function(e) {
      message("KEGG富集分析失败（可能是网络问题）: ", e$message)
      NULL
    })
    if (!is.null(ekegg_all) && nrow(as.data.frame(ekegg_all)) > 0) {
      pdf(file.path(plot_dir, "KEGG_enrichment_hub_threeway.pdf"), width = 8, height = 6)
      print(dotplot(ekegg_all, showCategory = 20) + ggtitle("KEGG of Drug-Disease-Ferroptosis hub genes"))
      dev.off()
    }
  }
}

# ---------------------------
# 4. 导出 Cytoscape 网络文件：Herb–Ingredient–Target
# ---------------------------

# 只保留与 CHB DEGs 相关的靶点，尤其是三者交集的 hub genes

network_edges <- tcm_targets %>%
  filter(Target %in% inter_deg_tcm)

if (nrow(network_edges) == 0) {
  warning("绵茵陈靶点与 CHB DEGs 交集为空，network_edges 将为示例数据。")
}

write.csv(network_edges,
          file.path(res_dir, "cytoscape_edges_MianYinChen_CHB.csv"),
          row.names = FALSE)

nodes <- data.frame(
  Node = c(unique(network_edges$Ingredient), unique(network_edges$Target)),
  Type = c(rep("Ingredient", length(unique(network_edges$Ingredient))),
           rep("Target",     length(unique(network_edges$Target))))
)

# 标记三者交集 hub genes
nodes$Is_Hub <- ifelse(nodes$Node %in% inter_all_three, "Yes", "No")

write.csv(nodes,
          file.path(res_dir, "cytoscape_nodes_MianYinChen_CHB.csv"),
          row.names = FALSE)

message("绵茵陈网络药理整合完成：Venn、富集分析及 Cytoscape 网络文件已生成。")
