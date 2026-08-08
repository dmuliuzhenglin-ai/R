## ======================================================================
## 耐药TLE多组学分析流程 · 脚本16/17：EDIL3配体-受体定位（整合素αvβ3/αvβ5、LFA-1）+ CellChat组间差异互作
## 对应论文：3.4节相关补充分析
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

########################################################################
# 图11 EDIL3 配体-受体分析脚本
# 目的：EDIL3 是646个共有DEG中仅6个反向变化基因之一（兴奋↑ 小胶↓ -1.92）
#       且为分泌型整合素配体，本脚本定位其受体在哪类细胞上表达
# 已知生物学：
#   EDIL3 (DEL-1) 结合整合素 αvβ3 (ITGAV+ITGB3)、αvβ5 (ITGAV+ITGB5)
#     → 促进巨噬细胞/小胶质 efferocytosis（吞噬清除凋亡细胞）
#   EDIL3 拮抗 LFA-1 (ITGAL+ITGB2)
#     → 抑制中性粒细胞/白细胞跨内皮浸润
# 结构：
#   模块0  环境准备 + 读取v2总对象
#   模块1  EDIL3 本身在各细胞类型×分组的表达（验证反向变化）
#   模块2  候选受体在各细胞类型的表达定位（核心结果，快速）
#   模块3  【可选】完整 CellChat 全局互作分析（耗时较长）
#   模块4  【可选】NicheNet 配体活性分析（需先装 nichenetr）
# 输出文件夹：图11_EDIL3配体受体分析
########################################################################


# ==================== 模块0 环境准备 ====================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 读取小胶质重分析时写回的 v2 总对象（含 celltype.2nd、condition 列）
work_dir <- cfg$dir_sc  # 工作目录（参数化，见 config/config.R）
setwd(work_dir)                                   # 应用工作目录
dir.create("图11_EDIL3配体受体分析", showWarnings = FALSE)  # ★必须在setwd之后
sc_data <- readRDS("Seurat_final_annotation_v2.rds")


# 确认分组列存在
table(sc_data$condition)          # 应显示 Control / TLE 两组细胞数
table(sc_data$celltype.2nd)       # 确认细胞类型注释列

Idents(sc_data) <- "celltype.2nd"


# ==================== 模块1 EDIL3本身：各细胞类型×分组表达 ====================
# 验证散点图中的反向变化，并看 EDIL3 在全部细胞类型里的分布

# --- 1.1 全细胞类型 DotPlot（按分组拆分） ---
p1 <- DotPlot(sc_data, features = "EDIL3",
              group.by = "celltype.2nd",
              split.by = "condition") +
  RotatedAxis() +
  ggtitle("EDIL3 expression: celltype x condition")
ggsave("图11_EDIL3配体受体分析/图11_1_EDIL3_各细胞类型分组DotPlot.pdf",
       p1, width = 10, height = 5)

# --- 1.2 定量表：各细胞类型中 EDIL3 的 Control/TLE 平均表达与检出率 ---
avg_exp <- AverageExpression(sc_data, features = "EDIL3",
                             group.by = c("celltype.2nd", "condition"),
                             assays = "RNA", slot = "data")$RNA
edil3_tbl <- as.data.frame(avg_exp)
edil3_tbl$gene <- "EDIL3"
write.csv(edil3_tbl,
          "图11_EDIL3配体受体分析/图11_1_EDIL3_各细胞类型分组平均表达.csv")

# --- 1.3 各细胞类型内 EDIL3 的 TLE vs Control 差异检验 ---
# 对每种主要细胞类型单独做 Wilcoxon，看反向变化在哪些细胞类型显著
celltypes <- unique(sc_data$celltype.2nd)
deg_list <- list()
for (ct in celltypes) {
  sub <- subset(sc_data, subset = celltype.2nd == ct)
  # 两组细胞数都够才检验
  if (min(table(sub$condition)) >= 30) {
    Idents(sub) <- "condition"
    m <- FindMarkers(sub, ident.1 = "TLE", ident.2 = "Control",
                     features = "EDIL3", logfc.threshold = 0,
                     min.pct = 0, test.use = "wilcox")
    m$celltype <- ct
    m$gene <- "EDIL3"
    deg_list[[ct]] <- m
  }
}
edil3_deg <- bind_rows(deg_list) %>%
  dplyr::select(celltype, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  arrange(avg_log2FC)
write.csv(edil3_deg,
          "图11_EDIL3配体受体分析/图11_1_EDIL3_各细胞类型组间差异检验.csv",
          row.names = FALSE)
print(edil3_deg)   # 重点看：兴奋性神经元是否↑、小胶质是否↓且显著
# --- 1.4 可视化：EDIL3 组间 log2FC 条形图（按细胞类型排序） ---
# 柱高 = log2FC(TLE vs Control)，颜色 = 方向，星号 = 显著性
edil3_deg <- edil3_deg %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up in TLE", "Down in TLE"),
         sig = case_when(p_val_adj < 0.001 ~ "***",
                         p_val_adj < 0.01  ~ "**",
                         p_val_adj < 0.05  ~ "*",
                         TRUE              ~ "ns"),
         celltype = factor(celltype, levels = celltype[order(avg_log2FC)]))
p_fc <- ggplot(edil3_deg, aes(x = celltype, y = avg_log2FC, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sig, hjust = ifelse(avg_log2FC > 0, -0.2, 1.2)), size = 5) +
  coord_flip() +
  scale_fill_manual(values = c("Down in TLE" = "#2166AC", "Up in TLE" = "#B2182B")) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(title = "EDIL3: TLE vs Control log2FC by cell type",
       x = NULL, y = "avg_log2FC (TLE vs Control)") +
  theme_bw(base_size = 13)
ggsave("图11_EDIL3配体受体分析/图11_1_EDIL3_各细胞类型log2FC条形图.pdf",
       p_fc, width = 8, height = 6)

# --- 1.5 可视化：EDIL3 检出率对比（pct.1=TLE vs pct.2=Control 哑铃图） ---
# 检出率变化比FC更直观（内皮36%→4%这类断崖一目了然）
pct_long <- edil3_deg %>%
  dplyr::select(celltype, pct.1, pct.2) %>%
  tidyr::pivot_longer(cols = c(pct.1, pct.2),
                      names_to = "group", values_to = "pct") %>%
  mutate(group = recode(group, "pct.1" = "TLE", "pct.2" = "Control"))
p_dumbbell <- ggplot(pct_long) +
  geom_line(aes(x = pct, y = celltype, group = celltype), color = "grey70") +
  geom_point(aes(x = pct, y = celltype, color = group), size = 3.5) +
  scale_color_manual(values = c("Control" = "#2166AC", "TLE" = "#B2182B")) +
  labs(title = "EDIL3 detection rate: Control vs TLE",
       x = "Percent expressed", y = NULL) +
  theme_bw(base_size = 13)
ggsave("图11_EDIL3配体受体分析/图11_1_EDIL3_检出率哑铃图.pdf",
       p_dumbbell, width = 8, height = 6)


# ==================== 模块2 候选受体表达定位（核心） ====================
# EDIL3 的已知/候选受体亚基
#   ITGAV + ITGB3 = αvβ3  → efferocytosis（巨噬细胞/小胶质吞噬清除）
#   ITGAV + ITGB5 = αvβ5  → efferocytosis
#   ITGAL + ITGB2 = αLβ2 (LFA-1) → 白细胞黏附/浸润（EDIL3是其拮抗剂）
#   ITGB1 作为泛整合素β链对照一并检查

receptors <- c("ITGAV", "ITGB3", "ITGB5", "ITGAL", "ITGB2", "ITGB1")

# --- 2.1 受体×细胞类型 DotPlot（总览，回答"受体在哪类细胞上"） ---
p2 <- DotPlot(sc_data, features = receptors,
              group.by = "celltype.2nd") +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("EDIL3 candidate receptors across cell types")
ggsave("图11_EDIL3配体受体分析/图11_2_EDIL3受体_各细胞类型DotPlot.pdf",
       p2, width = 9, height = 6)

# --- 2.2 受体×细胞类型×分组 DotPlot（看受体本身是否也有组间变化） ---
p3 <- DotPlot(sc_data, features = receptors,
              group.by = "celltype.2nd",
              split.by = "condition") +
  RotatedAxis() +
  ggtitle("EDIL3 receptors: celltype x condition")
ggsave("图11_EDIL3配体受体分析/图11_2_EDIL3受体_各细胞类型分组DotPlot.pdf",
       p3, width = 11, height = 6)


# --- 2.3 定量表：受体的各细胞类型检出率 ---
# DotPlot 里的 pct.exp 数字导出，方便写结果段引用
pct_mat <- DotPlot(sc_data, features = receptors,
                   group.by = "celltype.2nd")$data
write.csv(pct_mat,
          "图11_EDIL3配体受体分析/图11_2_EDIL3受体_检出率定量表.csv",
          row.names = FALSE)

# --- 2.4 受体组间差异检验（小胶质内重点） ---
# 看小胶质自身受体是否也下调（受体+配体双降 = 通路整体关闭，故事更强）
sort(table(sc_data$celltype.2nd), decreasing = TRUE)
all_labels <- unique(as.character(sc_data$celltype.2nd))
micro_labels <- grep("microglia", all_labels, ignore.case = TRUE, value = TRUE)
cat("检测到的细胞类型标签:\n"); print(sort(all_labels))
cat("识别为小胶质的标签:\n"); print(micro_labels)
if (length(micro_labels) == 0) {
  stop("未找到小胶质标签！请查看上面列出的全部标签，把 micro_labels 手动改成实际名称，如: micro_labels <- c('你的小胶质标签名')")
}
micro <- subset(sc_data, subset = celltype.2nd %in% micro_labels)
if (min(table(micro$condition)) >= 30) {
  Idents(micro) <- "condition"
  rec_deg <- FindMarkers(micro, ident.1 = "TLE", ident.2 = "Control",
                         features = receptors, logfc.threshold = 0,
                         min.pct = 0, test.use = "wilcox")
  rec_deg$gene <- rownames(rec_deg)
  write.csv(rec_deg,
            "图11_EDIL3配体受体分析/图11_2_EDIL3受体_小胶质内组间差异.csv",
            row.names = FALSE)
  print(rec_deg)
  # --- 2.5 可视化：小胶质内受体 log2FC 条形图 ---
  rec_plot_df <- rec_deg %>%
    mutate(direction = ifelse(avg_log2FC > 0, "Up in TLE", "Down in TLE"),
           sig = case_when(p_val_adj < 0.001 ~ "***",
                           p_val_adj < 0.01  ~ "**",
                           p_val_adj < 0.05  ~ "*",
                           TRUE              ~ "ns"),
           gene = factor(gene, levels = gene[order(avg_log2FC)]))
  p_rec <- ggplot(rec_plot_df, aes(x = gene, y = avg_log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig, hjust = ifelse(avg_log2FC > 0, -0.2, 1.2)), size = 5) +
    coord_flip() +
    scale_fill_manual(values = c("Down in TLE" = "#2166AC", "Up in TLE" = "#B2182B")) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(title = "EDIL3 receptors in Microglia: TLE vs Control",
         subtitle = "Ligand EDIL3 is DOWN; antagonized LFA-1 (ITGAL/ITGB2) is UP",
         x = NULL, y = "avg_log2FC (TLE vs Control)") +
    theme_bw(base_size = 13)
  ggsave("图11_EDIL3配体受体分析/图11_2_EDIL3受体_小胶质log2FC条形图.pdf",
         p_rec, width = 8, height = 5)
  
  # --- 2.6 可视化：配体-受体对仗汇总图（主图，一图讲完故事） ---
  ligand_row <- edil3_deg %>% filter(grepl("microglia", celltype, ignore.case = TRUE)) %>%
    summarise(avg_log2FC = mean(avg_log2FC), p_val_adj = min(p_val_adj)) %>%
    mutate(gene = "EDIL3 (ligand)")
  pair_df <- bind_rows(
    ligand_row,
    rec_deg %>% dplyr::select(gene, avg_log2FC, p_val_adj)) %>%
    mutate(direction = ifelse(avg_log2FC > 0, "Up in TLE", "Down in TLE"),
           sig = case_when(p_val_adj < 0.001 ~ "***",
                           p_val_adj < 0.01  ~ "**",
                           p_val_adj < 0.05  ~ "*",
                           TRUE              ~ "ns"),
           gene = factor(gene, levels = gene[order(avg_log2FC)]))
  p_pair <- ggplot(pair_df, aes(x = gene, y = avg_log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig, hjust = ifelse(avg_log2FC > 0, -0.2, 1.2)), size = 5) +
    coord_flip() +
    scale_fill_manual(values = c("Down in TLE" = "#2166AC", "Up in TLE" = "#B2182B")) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(title = "EDIL3-integrin axis in Microglia (TLE vs Control)",
         subtitle = "Ligand down + adhesion receptor up + efferocytosis receptor down",
         x = NULL, y = "avg_log2FC (TLE vs Control)") +
    theme_bw(base_size = 13)
  ggsave("图11_EDIL3配体受体分析/图11_2_EDIL3配体受体对仗图.pdf",
         p_pair, width = 8, height = 5.5)
  }
list.files(file.path(cfg$dir_sc, "图11_EDIL3配体受体分析"))   # 输出图位置
getwd()                                                          # 确认当前工作目录
# ★ 判读要点：
# ① ITGAV/ITGB3/ITGB5 若主要在小胶质/巨噬细胞检出
#    → "神经元/内皮释放EDIL3 → 小胶质通过αvβ3/αvβ5接收，执行吞噬清除"
#    → TLE中小胶质EDIL3自身也↓ = 自分泌清除信号减弱，与DAM表型互证
# ② ITGAL/ITGB2 若在浸润免疫细胞检出
#    → EDIL3↓ = LFA-1拮抗解除 = 外周免疫细胞浸润门槛降低
#    → 直接接通你的"肠-免疫-BBB-脑"级联主题


# ==================== 模块3 【可选】完整 CellChat 分析 ====================
# 注意：CellChat 安装较重（依赖多），分析需30分钟~2小时，按机器而定
# 安装（跑一次即可）：
# install.packages("BiocManager")
# BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"))
# install.packages("devtools")
# devtools::install_github("sqjin/CellChat")
# 西湖大学镜像：CRAN
options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/cran"))
# 西湖大学镜像：Bioconductor
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "BiocNeighbors", "ggalluvial"),
                     update = FALSE, ask = FALSE)
install.packages(c("devtools", "NMF", "circlize", "mnormt", "ggrepel",
                   "cowplot", "igraph", "reshape2", "svglite"))

options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/cran"))

# 一个一个装，哪个崩了一眼就能定位
install.packages("devtools")
install.packages("NMF")
install.packages("circlize")
install.packages("mnormt")
install.packages("ggrepel")
install.packages("cowplot")
install.packages("igraph")
install.packages("reshape2")
install.packages("svglite")
# codeload 是 GitHub 的下载通道，不吃 API 限流
download.file(
  "https://codeload.github.com/jinworks/CellChat/zip/refs/heads/main",
  destfile = "CellChat-main.zip", mode = "wb")
devtools::install_local("CellChat-main.zip")

  library(CellChat)
  packageVersion("CellChat")
  # --- 3.1 构建 CellChat 对象 ---
  # 用 log-normalized 表达矩阵 + 细胞类型注释
  # --- 诊断：找出缺失的亚基基因 ---
  db <- CellChatDB.human
  sub_cols <- grep("subunit", colnames(db$complex), value = TRUE)
  subunits <- unique(unlist(strsplit(unlist(db$complex[, sub_cols]), "_")))
  subunits <- subunits[subunits != "" & !is.na(subunits)]
  missing_genes <- setdiff(subunits, rownames(data.input))
  cat("缺失的亚基基因数:", length(missing_genes), "\n")
  print(missing_genes)     # 看看缺的是哪些（预期是免疫/T细胞相关基因）
  
  # --- 修复：补全零行 ---
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(0,
                               nrow = length(missing_genes),
                               ncol = ncol(data.input),
                               dimnames = list(missing_genes, colnames(data.input)),
                               sparse = TRUE)
    data.input <- rbind(data.input, zero_mat)
    cat("补零后矩阵维度:", dim(data.input), "\n")
  }
  
  # --- 重建 CellChat 对象并完整重跑 3.2 ---
  # 用补零后的 data.input 重建
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  download.file(
    "https://codeload.github.com/immunogenomics/presto/zip/refs/heads/master",
    destfile = "presto.zip", mode = "wb")
  devtools::install_local("presto.zip")
  library(presto)   # 验证
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)          # 现在缺基因已补上，不会再越界
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, "图11_EDIL3配体受体分析/cellchat_全部细胞.rds")
  
  # --- 3.2 预处理与推断 ---
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  dim(cellchat@data.signaling)          # 行数（信号基因数）是多少？为0或很小就是匹配失败
  sum(is.na(cellchat@idents))           # 分组里有没有NA
  table(cellchat@idents)                # 各细胞类型细胞数，有没有异常
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, "图11_EDIL3配体受体分析/cellchat_全部细胞.rds")
  # --- 3.2 预处理与推断（标准流程） ---
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)          # 最耗时的一步
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  saveRDS(cellchat, "图11_EDIL3配体受体分析/cellchat_全部细胞.rds")

  # --- 3.3 全局互作图 ---
  groupSize <- as.numeric(table(cellchat@idents))
  pdf("图11_EDIL3配体受体分析/图11_3_CellChat_全局互作网络.pdf", 9, 9)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Number of interactions")
  dev.off()

  # --- 3.4 查 EDIL3 相关互作 ---
  # CellChatDB 中检索 EDIL3（若收录，会列出其互作对）
  subsetDB <- cellchat@DB$interaction
  edil3_hits <- subsetDB[grepl("EDIL3", subsetDB$ligand) |
                         grepl("EDIL3", subsetDB$receptor), ]
  write.csv(edil3_hits,
            "图11_EDIL3配体受体分析/图11_3_CellChatDB中EDIL3条目.csv",
            row.names = FALSE)
  print(edil3_hits)   # 若为空 = CellChatDB未收录EDIL3，以模块2结果为准

  # 若数据库收录了，画小胶质作为接收方的相关信号
  df.net <- subsetCommunication(cellchat)
  edil3_net <- df.net[df.net$ligand == "EDIL3" | df.net$receptor == "EDIL3", ]
  write.csv(edil3_net,
            "图11_EDIL3配体受体分析/图11_3_EDIL3相关互作对.csv",
            row.names = FALSE)

  # ==================== 模块5 【核心】Control vs TLE 差异互作分析 ====================
  # 模块3是全细胞混合推断（只看结构）；模块5按分组各自推断再对比（才能进正文）
  # 前置条件：模块0已跑（sc_data存在）、模块3已跑通（data.input已补零、meta存在）
  # 耗时：两次computeCommunProb，约2~4小时；建议晚上挂跑
  
  # --- 5.1 按分组拆分细胞 ---
  cells_control <- rownames(meta)[meta$samples == "Control"]
  cells_tle     <- rownames(meta)[meta$samples == "TLE"]
  cat("Control:", length(cells_control), "细胞; TLE:", length(cells_tle), "细胞\n")
  
  # --- 5.2 封一个函数：建对象+完整推断（两组各跑一次） ---
  run_cellchat_pipeline <- function(cells_use, tag) {
    cc <- createCellChat(object = data.input[, cells_use],
                         meta = meta[cells_use, ],
                         group.by = "labels")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc)                    # 最耗时
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    saveRDS(cc, paste0("图11_EDIL3配体受体分析/cellchat_", tag, ".rds"))
    return(cc)
  }
  
  cellchat_control <- run_cellchat_pipeline(cells_control, "Control")
  cellchat_tle     <- run_cellchat_pipeline(cells_tle,     "TLE")
  
  # --- 5.3 合并两组对象 ---
  object.list <- list(Control = cellchat_control, TLE = cellchat_tle)
  cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
  saveRDS(cellchat_merged, "图11_EDIL3配体受体分析/cellchat_merged_Control_TLE.rds")
  
  # --- 5.4 互作总数与总强度对比（条形图，最概览的一张） ---
  gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE,
                             group = c(1, 2), measure = "count")
  gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE,
                             group = c(1, 2), measure = "weight")
  p_count <- gg1 + gg2
  ggsave("图11_EDIL3配体受体分析/图11_5_Control_vs_TLE_互作数量与强度.pdf",
         p_count, width = 8, height = 4)
  
  # --- 5.5 差异互作网络图（红=TLE增强, 蓝=TLE减弱） ---
  pdf("图11_EDIL3配体受体分析/图11_5_差异互作网络_数量.pdf", 9, 9)
  netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
  dev.off()
  pdf("图11_EDIL3配体受体分析/图11_5_差异互作网络_强度.pdf", 9, 9)
  netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE,
                            measure = "weight")
  dev.off()
  
  # --- 5.6 差异互作热图（细胞类型×细胞类型，适合正文） ---
  pdf("图11_EDIL3配体受体分析/图11_5_差异互作热图_数量.pdf", 7, 6)
  netVisual_heatmap(cellchat_merged)
  dev.off()
  pdf("图11_EDIL3配体受体分析/图11_5_差异互作热图_强度.pdf", 7, 6)
  netVisual_heatmap(cellchat_merged, measure = "weight")
  dev.off()
  
  # --- 5.7 信号通路层面差异排序（rankNet：哪些通路在TLE整体增强/减弱） ---
  p_rank <- rankNet(cellchat_merged, mode = "comparison",
                    stacked = TRUE, do.stat = TRUE)
  ggsave("图11_EDIL3配体受体分析/图11_5_差异信号通路_rankNet.pdf",
         p_rank, width = 9, height = 7)
  
  # --- 5.8 小胶质作为接收方的差异互作气泡图（与DAM故事对接，正文级） ---
  micro_targets <- grep("Microglia", levels(cellchat_merged@idents$joint), value = TRUE)
  p_bubble <- netVisual_bubble(cellchat_merged,
                               targets.use = micro_targets,
                               comparison = c(1, 2),
                               angle.x = 45,
                               title.name = "Incoming signaling to Microglia: Control vs TLE")
  ggsave("图11_EDIL3配体受体分析/图11_5_小胶质接收信号_组间气泡图.pdf",
         p_bubble, width = 12, height = 8)
  
  # --- 5.9 定量表：所有互作的两组对比（写结果段引用数字用） ---
  df_control <- subsetCommunication(cellchat_control)
  df_tle     <- subsetCommunication(cellchat_tle)
  write.csv(df_control, "图11_EDIL3配体受体分析/图11_5_互作明细_Control.csv", row.names = FALSE)
  write.csv(df_tle,     "图11_EDIL3配体受体分析/图11_5_互作明细_TLE.csv",     row.names = FALSE)
  
  # ★ 判读要点：
  # 5.4 → TLE总互作数/强度是否上升（炎症脑预期上升）
  # 5.5/5.6 → 哪些细胞类型对之间变化最大（重点看指向Microglia_*的）
  # 5.7 → 红色(TLE富集)通路里找 MIF/SPP1/TNF/CCL/BAFF/APOE 等炎症-DAM相关通路
  # 5.8 → 谁在给DAM发信号、TLE中增强还是减弱 = 级联的"上游"
  # EDIL3不在CellChatDB中，本模块不会有EDIL3结果，其证据以模块2为准
  
  # 5.8-改良版：只看TLE增强的免疫相关通路 + 小胶质接收
  immune_pathways <- c("COMPLEMENT", "MIF", "SPP1", "TNF", "BAFF",
                       "GALECTIN", "ANNEXIN", "CX3C", "CSF", "APP", "TGFb")
  # 先确认这些通路名在数据里的实际写法
  print(cellchat_merged@netP$pathways)   # 按输出微调上面的名单（大小写要一致）
  
  p_bubble2 <- netVisual_bubble(cellchat_merged,
                                targets.use = micro_targets,
                                signaling = immune_pathways,
                                comparison = c(1, 2),
                                angle.x = 45,
                                remove.isolate = TRUE,
                                title.name = "Immune-related incoming signaling to Microglia")
  ggsave("图11_EDIL3配体受体分析/图11_5_小胶质接收信号_免疫通路气泡图.pdf",
         p_bubble2, width = 10, height = 6)
 
   p_axis <- netVisual_bubble(cellchat_merged,
                             sources.use = "Plasma_Cell",
                             targets.use = grep("Microglia", levels(cellchat_merged@idents$joint), value = TRUE),
                             comparison = c(1, 2), angle.x = 45,
                             remove.isolate = TRUE,
                             title.name = "Plasma cell to Microglia: Control vs TLE")
  ggsave("图11_EDIL3配体受体分析/图11_5_浆细胞到小胶质_组间对比.pdf",
         p_axis, width = 9, height = 6)

