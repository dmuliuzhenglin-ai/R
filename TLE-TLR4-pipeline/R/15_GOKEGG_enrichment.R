## ======================================================================
## 耐药TLE多组学分析流程 · 脚本15/17：GO/KEGG富集分析（5个细胞类型基因集，上/下调分开富集）
## 对应论文：3.4节（细胞类型特异性病理程序）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

################################################################################
# 文件名 : 图8_GOKEGG富集分析.R
# 用  途 : TLE vs Control 组间差异基因的 GO/KEGG 富集分析(文件夹8)
#          五个细胞类型基因集: ① 兴奋性神经元 ② 小胶质 ③ 星形胶质
#          ④ 抑制性神经元 ⑤ 兴奋∩小胶质共有DEG(主线)
#          均分开上/下调分别富集, 最后出多集合对比图
# 前  提 : 工作目录下存在 Seurat_final_annotation_v2.rds(修正后的总对象)
# 输  出 : 图8_ 开头的pdf图和csv表
# 包安装 : 若缺包, 先运行(只需一次):
#          if (!require("BiocManager")) install.packages("BiocManager")
#          BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
# 备  注 : enrichKEGG 需要联网; 若想改用亚型间DEG(如Stressed vs Normal),
#          只需替换模块1的输入数据框, 模块2-5完全通用
################################################################################


# ============================================================================
# 模块0 | 前置准备
# ============================================================================
# 第1步: 确认BiocManager已安装
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 第2步: 安装全部所需包(逐个装, 方便定位哪个出错)
BiocManager::install("Biostrings")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("DOSE")
BiocManager::install(c("GO.db", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
# ---- 0.1 加载R包 ----
library(Seurat)            # 读取对象/差异分析
library(dplyr)             # 数据处理
library(ggplot2)           # 绘图
library(clusterProfiler)   # 富集分析主包
library(org.Hs.eg.db)      # 人基因注释(ID转换)
library(enrichplot)        # 富集结果可视化

# ---- 0.2 路径与对象 ----
work_dir <- cfg$dir_sc  # 工作目录（参数化，见 config/config.R）
setwd(work_dir)
sc_data <- readRDS("Seurat_final_annotation_v2.rds")  # 修正后的总对象
DefaultAssay(sc_data) <- "RNA"
sc_data$group <- factor(sc_data$condition, levels = c("Control", "TLE"))  # 组间比较分组


# ============================================================================
# 模块1 | 计算三个DEG列表(TLE vs Control, 组间)
# ============================================================================

# ---- 1.1 兴奋性神经元(合并两个亚型) TLE vs Control ----
exc <- subset(sc_data, subset = celltype.2nd %in% c("Excitatory_Neuron",
                                                    "Excitatory_Neuron_Stressed"))
Idents(exc) <- "group"                                  # 以疾病分组为比较身份
deg_exc <- FindMarkers(exc, ident.1 = "TLE", ident.2 = "Control",  # 正值=TLE高表达
                       test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
deg_exc$gene <- rownames(deg_exc)                       # 行名转gene列

# ---- 1.2 小胶质(合并三个亚群, 排除Other) TLE vs Control ----
micro <- subset(sc_data, subset = celltype.2nd %in% c("Microglia_Homeostatic",
                                                      "Microglia_Activated",
                                                      "Microglia_DAM"))
Idents(micro) <- "group"
deg_micro <- FindMarkers(micro, ident.1 = "TLE", ident.2 = "Control",
                         test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
deg_micro$gene <- rownames(deg_micro)
# 注: 合并亚群做组间比较, 结果会部分反映亚群比例变化; 若要更纯的通路信号,
#     可对 Microglia_DAM 等单个亚群单独跑(替换subset条件即可)

# ---- 1.2b 星形胶质 / 抑制性神经元 TLE vs Control(扩展分析) ----
astro <- subset(sc_data, subset = celltype.2nd == "Astrocyte_Reactive")  # 星形胶质
Idents(astro) <- "group"
deg_astro <- FindMarkers(astro, ident.1 = "TLE", ident.2 = "Control",
                         test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
deg_astro$gene <- rownames(deg_astro)

inh <- subset(sc_data, subset = celltype.2nd == "Inhibitory_Neuron")     # 抑制性神经元
Idents(inh) <- "group"
deg_inh <- FindMarkers(inh, ident.1 = "TLE", ident.2 = "Control",
                       test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
deg_inh$gene <- rownames(deg_inh)

# ---- 1.3 显著DEG筛选与上下调拆分 ----
sig_exc   <- subset(deg_exc,   p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)  # 兴奋显著DEG
sig_micro <- subset(deg_micro, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)  # 小胶质显著DEG

exc_up     <- sig_exc$gene[sig_exc$avg_log2FC > 0]      # 兴奋-TLE上调
exc_down   <- sig_exc$gene[sig_exc$avg_log2FC < 0]      # 兴奋-TLE下调(功能崩溃)
micro_up   <- sig_micro$gene[sig_micro$avg_log2FC > 0]  # 小胶质-TLE上调(异常活化)
micro_down <- sig_micro$gene[sig_micro$avg_log2FC < 0]  # 小胶质-TLE下调

sig_astro <- subset(deg_astro, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)  # 星形胶质显著DEG
sig_inh   <- subset(deg_inh,   p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)  # 抑制性显著DEG
astro_up   <- sig_astro$gene[sig_astro$avg_log2FC > 0]   # 星形胶质-TLE上调
astro_down <- sig_astro$gene[sig_astro$avg_log2FC < 0]   # 星形胶质-TLE下调
inh_up     <- sig_inh$gene[sig_inh$avg_log2FC > 0]       # 抑制性-TLE上调
inh_down   <- sig_inh$gene[sig_inh$avg_log2FC < 0]       # 抑制性-TLE下调

# ---- 1.4 二者共有DEG(级联候选) ----
common_genes <- intersect(sig_exc$gene, sig_micro$gene) # 两细胞类型共同显著
common_df <- data.frame(
  Gene         = common_genes,
  Exc_log2FC   = sig_exc$avg_log2FC[match(common_genes, sig_exc$gene)],     # 兴奋中FC
  Micro_log2FC = sig_micro$avg_log2FC[match(common_genes, sig_micro$gene)]  # 小胶质中FC
)
common_df$同向变化 <- sign(common_df$Exc_log2FC) == sign(common_df$Micro_log2FC)  # 方向是否一致
write.csv(common_df, "图8_共有DEG_组间.csv", row.names = FALSE)

cat("兴奋DEG: 上调", length(exc_up), " 下调", length(exc_down), "\n")
cat("小胶质DEG: 上调", length(micro_up), " 下调", length(micro_down), "\n")
cat("星形胶质DEG: 上调", length(astro_up), " 下调", length(astro_down), "\n")
cat("抑制性DEG: 上调", length(inh_up), " 下调", length(inh_down), "\n")
cat("共有DEG:", length(common_genes), "个\n")
if (length(common_genes) < 30) {
  cat("提示: 共有基因<30个, 富集统计功效低, 建议仅做展示性分析(hub/网络)\n")
}


# ============================================================================
# 模块2 | 基因ID转换(SYMBOL→ENTREZID)与背景基因集
# ============================================================================

# ---- 2.1 待富集基因列表(9个集合) ----
gene_sets <- list(
  Exc_Up = exc_up,     Exc_Down = exc_down,       # 兴奋-上调/下调
  Micro_Up = micro_up, Micro_Down = micro_down,   # 小胶质-上调/下调
  Astro_Up = astro_up, Astro_Down = astro_down,   # 星形胶质-上调/下调
  Inh_Up = inh_up,     Inh_Down = inh_down,       # 抑制性神经元-上调/下调
  Common = common_genes                           # 共有(仍取兴奋∩小胶质, 主线不变)
)
gene_sets <- gene_sets[sapply(gene_sets, length) >= 10]  # 去掉基因数<10的集合(无法富集)

# ---- 2.2 背景基因集: 四类细胞检测到的基因取并集 ----
# (FindMarkers输出的行名即为参与检验的基因, 比全基因组背景更合理)
universe_genes <- union(deg_exc$gene,
                        union(deg_micro$gene, union(deg_astro$gene, deg_inh$gene)))

# ---- 2.3 SYMBOL → ENTREZID 转换(clusterProfiler的KEGG需要ENTREZ) ----
univ_map <- bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)                    # 背景ID映射表
universe_entrez <- univ_map$ENTREZID                      # 背景ENTREZ向量

to_entrez <- function(symbols) {                          # 单集合转换函数
  df <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  df$ENTREZID                                             # 返回ENTREZ向量
}
gene_sets_entrez <- lapply(gene_sets, to_entrez)          # 全部集合转ENTREZ
print(sapply(gene_sets_entrez, length))                   # 转换后各集合基因数(会有少量丢失)


# ============================================================================
# 模块3 | GO富集分析(BP, 生物学过程)
# ============================================================================

# ---- 3.1 逐个集合富集并导出 ----
go_results <- list()                                       # 存储各集合结果
for (nm in names(gene_sets_entrez)) {
  ego <- enrichGO(gene          = gene_sets_entrez[[nm]],  # 目标基因
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",                    # 生物学过程(可改"MF"/"CC")
                  pAdjustMethod = "BH",                    # BH校正
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  universe      = universe_entrez,         # 背景基因集
                  readable      = TRUE)                    # ENTREZ转回SYMBOL便于阅读
  go_results[[nm]] <- ego
  write.csv(as.data.frame(ego), paste0("图8_GO_BP_", nm, ".csv"), row.names = FALSE)  # 导出表格
}

# ---- 3.2 单集合气泡图(每个集合一张) ----
for (nm in names(go_results)) {
  if (!is.null(go_results[[nm]]) && nrow(as.data.frame(go_results[[nm]])) > 0) {
    p <- dotplot(go_results[[nm]], showCategory = 12) +    # 展示前12条通路
         ggtitle(paste0("GO BP - ", nm))
    ggsave(paste0("图8_GO_BP_", nm, "_气泡图.pdf"), p, width = 9, height = 7)
  }
}

# ---- 3.3 五集合对比气泡图(一图看清三类基因集各富集什么) ----
ck_go <- compareCluster(geneCluster = gene_sets_entrez,    # 多集合比较
                        fun         = "enrichGO",
                        OrgDb       = org.Hs.eg.db,
                        keyType     = "ENTREZID",
                        ont         = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2,
                        universe    = universe_entrez,
                        readable    = TRUE)
p <- dotplot(ck_go, showCategory = 6) +                    # 每集合展示前6条
     theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
     ggtitle("GO BP comparison: Excitatory vs Microglia vs Common")
ggsave("图8_GO_BP_五集合对比气泡图.pdf", p, width = 12, height = 10)
write.csv(as.data.frame(ck_go), "图8_GO_BP_对比汇总.csv", row.names = FALSE)


# ============================================================================
# 模块4 | KEGG富集分析(需要联网)
# ============================================================================

# ---- 4.1 逐个集合富集并导出 ----
kegg_results <- list()
for (nm in names(gene_sets_entrez)) {
  kk <- enrichKEGG(gene          = gene_sets_entrez[[nm]],
                   organism      = "hsa",                  # 人
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   universe      = universe_entrez)
  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")  # 转回SYMBOL
  kegg_results[[nm]] <- kk
  write.csv(as.data.frame(kk), paste0("图8_KEGG_", nm, ".csv"), row.names = FALSE)
}

# ---- 4.2 单集合气泡图 ----
for (nm in names(kegg_results)) {
  if (!is.null(kegg_results[[nm]]) && nrow(as.data.frame(kegg_results[[nm]])) > 0) {
    p <- dotplot(kegg_results[[nm]], showCategory = 12) +
         ggtitle(paste0("KEGG - ", nm))
    ggsave(paste0("图8_KEGG_", nm, "_气泡图.pdf"), p, width = 9, height = 7)
  }
}

# ---- 4.3 五集合对比气泡图 ----
ck_kegg <- compareCluster(geneCluster = gene_sets_entrez,
                          fun         = "enrichKEGG",
                          organism    = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          universe    = universe_entrez)
ck_kegg <- setReadable(ck_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
p <- dotplot(ck_kegg, showCategory = 6) +
     theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
     ggtitle("KEGG comparison across cell types")
ggsave("图8_KEGG_多集合对比气泡图.pdf", p, width = 13, height = 10)
write.csv(as.data.frame(ck_kegg), "图8_KEGG_对比汇总.csv", row.names = FALSE)

# ---- 4.4 (可选)基因-通路网络图: 看hub基因连接哪些通路 ----
# 以共有基因集的KEGG结果为例(基因数太少时图会很稀疏, 属正常)
if ("Common" %in% names(kegg_results) &&
    !is.null(kegg_results[["Common"]]) &&
    nrow(as.data.frame(kegg_results[["Common"]])) > 0) {
  p <- cnetplot(kegg_results[["Common"]], showCategory = 8)  # 通路-基因网络(新版enrichplot无colorEdge参数)
  ggsave("图8_KEGG_Common_基因通路网络图.pdf", p, width = 12, height = 10)
}


# ============================================================================
# 运行结束汇总
# ============================================================================
cat("\n========== GO/KEGG 富集分析完成 ==========\n")
cat("输出文件(图8_开头):\n")
cat("  图8_GO_BP_*.csv / 气泡图.pdf      —— 各集合GO结果\n")
cat("  图8_KEGG_*.csv / 气泡图.pdf       —— 各集合KEGG结果\n")
cat("  图8_GO_BP_五集合对比气泡图.pdf    —— 三类基因集GO对比(论文主图候选)\n")
cat("  图8_KEGG_五集合对比气泡图.pdf     —— 三类基因集KEGG对比\n")
cat("  图8_共有DEG_组间.csv              —— 共有基因及两细胞类型FC\n")
