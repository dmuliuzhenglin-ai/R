## ======================================================================
## 耐药TLE多组学分析流程 · 脚本14/17：补充图S11：辛伐他汀相关靶点（TLR4/LY96/TNF/NFKB1/CCL2/HMGCR）DotPlot
## 对应论文：补充图S11
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

########################################################################
## 补充图S11 终版：修正 TLE/Control 配色方向 + 缺行诊断
## 用法：RStudio 打开本文件 → 点 Source 按钮
## 输出：
##   file.path(cfg$dir_sc, "结果图表面板/补充图S11_辛伐他汀靶点_DotPlot.pdf")
########################################################################
library(Seurat)
library(ggplot2)

setwd(cfg$dir_sc)
sc_data <- readRDS("Seurat_final_annotation_corrected.rds")
DefaultAssay(sc_data) <- "RNA"
Idents(sc_data) <- "celltype.2nd"

stopifnot("Meningeal_Fibroblast" %in% sc_data$celltype.2nd)

## —— 诊断：各细胞类型 × 分组的细胞数（含 NA）——
cat("\n=== 诊断：细胞数表（看 Oligodendrocyte 的 Control 列与 NA 列）===\n")
print(table(sc_data$celltype.2nd, sc_data$condition, useNA = "always"))
cat("\ncondition 原始取值：\n")
print(unique(as.character(sc_data$condition)))

sc_data$condition <- factor(sc_data$condition, levels = c("Control", "TLE"))

statin_genes <- c("TLR4", "LY96", "TNF", "NFKB1", "CCL2", "HMGCR")
statin_genes <- intersect(statin_genes, rownames(sc_data[["RNA"]]))
stopifnot(length(statin_genes) >= 5)

## 上一版实测配色反了（Seurat 把 cols 逆序映射），这里对调：目标 Control=蓝、TLE=红
p_dot <- DotPlot(sc_data, features = statin_genes,
                 group.by = "celltype.2nd", split.by = "condition",
                 cols = c("#E64B35", "#3C5488"), assay = "RNA") +
  RotatedAxis() +
  labs(x = NULL, y = NULL)
ggsave("结果图表面板/补充图S11_辛伐他汀靶点_DotPlot.pdf",
       p_dot, width = 12, height = 9)
cat("\n补充图S11 终版已保存。请检查：TLE 行应为红色点、Control 行应为蓝色点。\n")
cat("如颜色仍反，把 cols 两个颜色再对调一次即可。\n")
