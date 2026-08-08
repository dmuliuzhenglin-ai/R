## ======================================================================
## 耐药TLE多组学分析流程 · 脚本11/17：补充图S3：5个Hub+EDIL3 细胞类型DotPlot（Control蓝/TLE红）
## 对应论文：补充图S3
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

# ============================================================
# 补充图S3：枢纽基因细胞类型DotPlot（重画 = 分组配色修正版）
# 旧版问题：TLE=浅灰几乎不可见、Control=蓝，与全套色板相反
# 新版：Control=蓝(#3C5488)、TLE=红(#E64B35)，与补充图S11一致
# 用法：RStudio 打开本文件 → 点 "Source" 按钮（不要粘贴进Console）
# 输出：补充图S3_枢纽基因DotPlot_final.pdf（保存在rds所在文件夹）
# ============================================================
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})

## 定位16类型终版对象（路径见 config/config.R）
rds_file <- file.path(cfg$dir_sc, "Seurat_final_annotation_corrected.rds")
if (!file.exists(rds_file)) {
  stop("找不到 ", rds_file, "：请按 README 将 Seurat_final_annotation_corrected.rds 放入 data/singlecell/")
}
message("读取对象: ", rds_file)
sc <- readRDS(rds_file)
stopifnot("Meningeal_Fibroblast" %in% sc$celltype.2nd)  # 保险：必须是16类型终版对象
DefaultAssay(sc) <- "RNA"
Idents(sc) <- "celltype.2nd"
sc <- tryCatch(JoinLayers(sc), error = function(e) sc)
sc$condition <- factor(sc$condition, levels = c("Control", "TLE"))

genes <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2", "EDIL3")
genes <- genes[genes %in% rownames(sc)]
stopifnot(length(genes) == 6)

p <- DotPlot(sc, features = genes,
             group.by = "celltype.2nd",
             split.by = "condition",
             cols = c("#E64B35", "#3C5488"),   # Seurat实测按levels逆序映射 → Control蓝、TLE红
             assay = "RNA") +
  RotatedAxis() +
  labs(x = NULL, y = NULL, title = NULL)

out <- file.path(dirname(rds_file), "补充图S3_枢纽基因DotPlot_final.pdf")
ggsave(out, p, width = 8.5, height = 6.5)
message("已保存: ", out)

## 目检要点：
## 1) 所有 _TLE 行的点为红色、_Control 行的点为蓝色（若相反请告诉我）
## 2) Plasma_Cell_Control 与 Endothelial_Activated_Control 行缺如属正常（对照组0个细胞）
