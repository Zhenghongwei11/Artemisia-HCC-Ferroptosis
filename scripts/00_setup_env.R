# 00_setup_env.R
# 统一安装与加载项目所需 R 包（在高配置电脑上执行一次即可）

# ============================================================================
# 禁用交互式图形设备 (防止XQuartz弹出)
# ============================================================================
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

cran_repo <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

base_pkgs <- c(
  # 基础数据处理与可视化
  "tidyverse", "ggplot2", "pheatmap", "ggpubr", "RColorBrewer", "VennDiagram", "cowplot",
  
  # GEO 数据与批量转录组
  "GEOquery", "limma", "edgeR", "sva",
  
  # 生存分析 / 临床模型 / 可视化
  "survival", "survminer", "rms", "timeROC", "gridExtra",

  # 注释数据库与映射
  "AnnotationDbi", "hgu133plus2.db",

  # WGCNA 与机器学习
  "WGCNA", "caret", "randomForest", "glmnet", "pROC",
  
  # 富集分析
  "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "DOSE", "enrichplot",
  
  # 单细胞 & 细胞通讯
  "Seurat", "patchwork", "monocle", "Matrix",

  # 外部验证（TCGA）
  "TCGAbiolinks", "SummarizedExperiment",

  # DCA / 日志等（部分脚本使用）
  "ggDCA", "futile.logger",
  
  # 热图等
  "ComplexHeatmap"
)

for (pkg in base_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    message("Installing: ", pkg)
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }, error = function(e) {
      message("[WARN] Failed to install ", pkg, ": ", e$message)
    })
  }
}

if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

message("如需使用 CellChat，请手动运行: devtools::install_github('sqjin/CellChat') 并 library(CellChat)")

message("环境准备完成（可能有个别包需手动处理依赖）")
