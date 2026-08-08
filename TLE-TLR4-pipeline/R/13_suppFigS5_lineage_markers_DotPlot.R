## ======================================================================
## 耐药TLE多组学分析流程 · 脚本13/17：补充图S5：16种细胞类型谱系标志物DotPlot（注释质量核查）
## 对应论文：补充图S5
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

# ============================================================
# 补充图S5：16种细胞类型谱系标志物DotPlot（注释质量核查）
# 替代上传的"注释检查_气泡图"（该图y轴为簇编号、缺小胶质/少突/免疫
# 标志物、顶部分组标签被裁、带占位标题，不宜直接使用）
# 用法：RStudio 打开本文件 → 点 "Source" 按钮
# 输出：补充图S5_注释标志物DotPlot.pdf（保存在rds所在文件夹）
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
sc <- tryCatch(JoinLayers(sc), error = function(e) sc)

type_order <- c("Excitatory_Neuron", "Inhibitory_Neuron",
                "Oligodendrocyte", "OPC",
                "Astrocyte", "Astrocyte_Reactive",
                "Microglia_Homeostatic", "Microglia_Activated",
                "Microglia_DAM", "Microglia_DAM_Prolif",
                "Endothelial", "Endothelial_Activated",
                "Pericyte", "Meningeal_Fibroblast",
                "T_Cell", "Plasma_Cell")
stopifnot(all(type_order %in% sc$celltype.2nd))
sc$celltype.2nd <- factor(as.character(sc$celltype.2nd), levels = type_order)
Idents(sc) <- "celltype.2nd"

markers <- c("SLC17A7", "SATB2", "RBFOX3",      # 兴奋性神经元
             "GAD1", "GAD2",                    # 抑制性神经元
             "MBP", "MOG",                      # 少突胶质
             "PDGFRA", "VCAN",                  # OPC
             "AQP4", "ALDH1L1",                 # 星形胶质（稳态）
             "GFAP", "CD44",                    # 星形胶质（反应性）
             "P2RY12", "TMEM119", "CX3CR1",     # 小胶质
             "CLDN5", "VWF",                    # 内皮
             "VCAM1", "SELE",                   # 激活内皮
             "PDGFRB", "RGS5",                  # 周细胞
             "COL1A1", "DCN",                   # 脑膜成纤维细胞
             "CD3D", "CD8A",                    # T细胞
             "JCHAIN", "MZB1")                  # 浆细胞
missing <- markers[!markers %in% rownames(sc)]
if (length(missing) > 0) message("以下标志物不在对象中，已自动跳过: ", paste(missing, collapse = ", "))
markers <- markers[markers %in% rownames(sc)]

p <- DotPlot(sc, features = markers, group.by = "celltype.2nd", assay = "RNA") +
  RotatedAxis() +
  labs(x = NULL, y = NULL, title = NULL)

out <- file.path(dirname(rds_file), "补充图S5_注释标志物DotPlot.pdf")
ggsave(out, p, width = 13, height = 6.5)
message("已保存: ", out)

## 目检要点：每类细胞应在其对应标志物列呈现大而深的点，对角线模式清晰
