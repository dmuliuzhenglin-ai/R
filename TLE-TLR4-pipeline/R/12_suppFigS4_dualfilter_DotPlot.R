## ======================================================================
## 耐药TLE多组学分析流程 · 脚本12/17：补充图S4：小胶质双重过滤后标志物DotPlot
## 对应论文：补充图S4
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

# ============================================================
# 补充图S4：小胶质细胞双重过滤后的标志物DotPlot
# 证明：过滤后4种小胶质状态表达P2RY12/TMEM119/CX3CR1，
#       基本不表达浆细胞标志JCHAIN/MZB1/IGKC；浆细胞行相反
# 用法：RStudio 打开本文件 → 点 "Source" 按钮
# 输出：补充图S4_双重过滤DotPlot.pdf（保存在rds所在文件夹）
# ============================================================
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})

rds_file <- file.path(cfg$dir_sc, "Seurat_final_annotation_corrected.rds")
if (!file.exists(rds_file)) {
  stop("找不到 ", rds_file, "：请按 README 将 Seurat_final_annotation_corrected.rds 放入 data/singlecell/")
}
message("读取对象: ", rds_file)
sc <- readRDS(rds_file)
stopifnot("Meningeal_Fibroblast" %in% sc$celltype.2nd)  # 必须是16类型终版对象
DefaultAssay(sc) <- "RNA"

micro_types <- c("Microglia_Homeostatic", "Microglia_Activated",
                 "Microglia_DAM", "Microglia_DAM_Prolif")
keep <- c(micro_types, "Plasma_Cell")
stopifnot(all(keep %in% sc$celltype.2nd))

sub <- subset(sc, subset = celltype.2nd %in% keep)
sub <- tryCatch(JoinLayers(sub), error = function(e) sub)
sub$celltype.2nd <- factor(as.character(sub$celltype.2nd), levels = keep)
Idents(sub) <- "celltype.2nd"

feats <- c("P2RY12", "TMEM119", "CX3CR1", "JCHAIN", "MZB1", "IGKC")
feats <- feats[feats %in% rownames(sub)]
stopifnot(length(feats) == 6)

p <- DotPlot(sub, features = feats, group.by = "celltype.2nd", assay = "RNA") +
  RotatedAxis() +
  labs(x = NULL, y = NULL, title = NULL)

out <- file.path(dirname(rds_file), "补充图S4_双重过滤DotPlot.pdf")
ggsave(out, p, width = 6.5, height = 3.5)
message("已保存: ", out)

## 质控打印：小胶质4行的 JCHAIN/MZB1/IGKC 的 pct.exp 应接近0；Plasma_Cell行相反
print(p$data[, c("id", "features.plot", "pct.exp", "avg.exp.scaled")])
