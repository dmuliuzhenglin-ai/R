## ======================================================================
## 耐药TLE多组学分析流程 · 脚本01/17：bulk发现队列：DESeq2 + gutMGene三重交集 + WGCNA + PPI/cytoHubba枢纽筛选
## 对应论文：3.1节（图1、补充图S1/S2、补充表S2/S5）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

################################################################################
## 课题：基于三重递进筛选与优先级排序解析肠道代谢-神经免疫-血脑屏障
##       交互调控耐药颞叶癫痫的纯生信分析
##
## 版本 v3（2026-07）：
##   - 核心数据集仅 GSE134697，使用用户本地下载的 CountMatrix.tsv
##     （featureCounts格式：前6列注释 + 36个样本计数列，Ensembl ID）
##   - Step1 改为 DESeq2（计数数据），不再下载GEO、不再使用芯片RMA流程
##   - 不使用小鼠KD数据集与GSE143272，无KD加权层
##   - gutMGene 使用本地两个CSV，筛选 human + causally
##
## 分析流程：
##   Step1: 读取count matrix -> DESeq2差异分析（耐药TLE vs Control）
##   Step2: gutMGene预处理 + 免疫基因集 -> 三重硬交集（主候选池）
##   Step3: WGCNA共表达网络 -> 表型相关模块 -> 与主候选池二次交集
##   Step4: STRING PPI -> 4种cytoHubba算法 -> Composite Score -> Top5 Hub
##   Step5: Hub GO/KEGG富集 + gutMGene反向查询（精确三联体）
################################################################################


#===============================================================================
# 0. 环境准备
#===============================================================================

## ---- 0.0 工作目录（参数化，见 config/config.R）--------------------------------
setwd(cfg$dir_bulk)

## ---- 0.1 加载包（计数矩阵流程：无需GEOquery/oligo/affy/hgu133plus2.db）--------
suppressPackageStartupMessages({
  library(DESeq2)          # RNA-seq差异分析
  library(WGCNA)           # 加权共表达网络
  library(STRINGdb)        # STRING蛋白互作
  library(igraph)          # 网络分析（cytoHubba四种算法的R实现）
  library(clusterProfiler) # GO/KEGG富集
  library(org.Hs.eg.db)    # 人类基因注释（Ensembl->Symbol）
  library(msigdbr)         # MSigDB免疫相关基因集
  library(data.table)      # 快速读取大表
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})
allowWGCNAThreads()

## ---- 0.2 目录结构 -------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("rds",     showWarnings = FALSE)

save_rds <- function(obj, name) saveRDS(obj, file.path("rds", paste0(name, ".rds")))


#===============================================================================
# Step 1. GSE134697 差异分析（RNA-seq count matrix，DESeq2）
#===============================================================================

## ---- 1.1 读取计数矩阵（featureCounts格式）--------------------------------------
## 文件结构：第1-6列 = Geneid/Chr/Start/End/Strand/Length（注释），
##           第7列起 = 36个样本的整数计数
counts_raw <- fread("data/GSE134697_CountMatrix.tsv/GSE134697_CountMatrix.tsv",
                    data.table = FALSE)
counts <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(counts) <- counts_raw$Geneid
storage.mode(counts) <- "integer"
rm(counts_raw)  # 释放内存
cat("计数矩阵维度：", nrow(counts), "基因 ×", ncol(counts), "样本\n")

## ---- 1.2 Ensembl ID -> Gene Symbol ----------------------------------------------
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = sub("\\..*$", "", rownames(counts)),  # 去版本号
                       column = "SYMBOL", keytype = "ENSEMBL",
                       multiVals = "first")
keep <- !is.na(gene_symbols) & gene_symbols != ""
counts <- counts[keep, ]
gene_symbols <- gene_symbols[keep]
## 同一Symbol对应多个Ensembl -> 计数合并
counts_gene <- rowsum(counts, group = gene_symbols, reorder = FALSE)
cat("注释合并后基因数：", nrow(counts_gene), "\n")

## ---- 1.3 样本分组（配对设计：17对患者海马H vs 配对皮层C）---------------------
sample_names <- colnames(counts_gene)

## 只保留17对配对样本（CHKJ-C/H），剔除2个无配对的CTRL
paired_idx <- grepl("CHKJ-[CH]\\d+", sample_names)
counts_paired <- counts_gene[, paired_idx]
sn <- sub("^.*/(CHKJ-[CH])(\\d+)-Brain.*$", "\\1-\\2", colnames(counts_paired))
colnames(counts_paired) <- sn

## 分组：H=患者海马(TLE)，C=同一患者新皮层(Control)
group <- factor(ifelse(grepl("CHKJ-H", sn), "TLE", "Control"),
                levels = c("Control", "TLE"))

## 患者编号：CHKJ-C-10 和 CHKJ-H-10 是同一个人 -> patient = "10"
patient <- factor(sub("^CHKJ-[CH]-?(\\d+)$", "\\1", sn))

cat("配对样本分组：\n"); print(table(group))          # 应输出 Control 17 / TLE 17
cat("患者数：", nlevels(patient), "\n")               # 应输出 17
print(data.frame(sample = sn, patient = patient, group = group))
stopifnot(all(table(patient, group) == 1))
cat("配对核对通过：17个患者，每人1个Control + 1个TLE\n")
## ---- 1.4 低表达过滤 + DESeq2配对差异分析 ----------------------------------------
keep_expr <- rowSums(counts_paired >= 10) >= 17   # 至少在一组中出现
counts_filt <- counts_paired[keep_expr, ]
cat("低表达过滤后基因数：", nrow(counts_filt), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = data.frame(group = group, patient = patient,
                       row.names = colnames(counts_filt)),
  design = ~ patient + group)        # ★ 配对设计：校正患者个体差异
dds <- DESeq(dds)
resultsNames(dds)                    # 确认最后一项是 group_TLE_vs_Control
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("apeglm", update = FALSE)
library(apeglm)
res_shrink <- lfcShrink(dds, coef = "group_TLE_vs_Control", type = "apeglm")
deg_all <- as.data.frame(res_shrink) %>%
  tibble::rownames_to_column("Gene") %>%
  filter(!is.na(padj)) %>%
  dplyr::select(Gene, log2FoldChange, pvalue, padj, baseMean) %>%
  dplyr::rename(logFC = log2FoldChange, adj.P.Val = padj) %>%
  arrange(adj.P.Val)

## 阈值：|log2FC| > 0.58 且 FDR < 0.05
deg_sig <- deg_all %>% filter(abs(logFC) > 0.58, adj.P.Val < 0.05)
cat("显著差异基因数：", nrow(deg_sig),
    "（上调", sum(deg_sig$logFC > 0), "，下调", sum(deg_sig$logFC < 0), "）\n")
fwrite(deg_all, "results/Step1_DEG_all_GSE134697.csv")
fwrite(deg_sig, "results/Step1_DEG_significant_GSE134697.csv")
save_rds(deg_all, "Step1_DEG_all"); save_rds(deg_sig, "Step1_DEG_sig")
## ---- 1.5 标准化表达矩阵（VST，供WGCNA/热图使用）-------------------------------------
vsd       <- assay(vst(dds, blind = FALSE))
expr_gene <- vsd
rownames(expr_gene) <- toupper(rownames(expr_gene))  # 与gutMGene/免疫集大小写对齐
save_rds(expr_gene, "Step1_expression_genelevel")

## ---- 1.6 质控：PCA ---------------------------------------------------------------
pdf("figures/Step1_QC_PCA.pdf", width = 6, height = 5)
pca <- prcomp(t(expr_gene), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = 19, col = as.numeric(group) + 1,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
     main = "PCA of samples (VST)")
legend("topright", legend = levels(group), col = 2:3, pch = 19)
dev.off()

## ---- 1.7 火山图 + 差异基因热图 --------------------------------------------------------
deg_all <- deg_all %>%
  mutate(sig = case_when(adj.P.Val < 0.05 & logFC >  0.58 ~ "Up",
                         adj.P.Val < 0.05 & logFC < -0.58 ~ "Down",
                         TRUE ~ "NS"))
p_volcano <- ggplot(deg_all, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NS" = "grey80")) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = deg_all %>% filter(sig != "NS") %>%
                    arrange(adj.P.Val) %>% head(15),
                  aes(label = Gene), size = 3, max.overlaps = 20) +
  theme_bw() + labs(title = "GSE134697: drug-resistant TLE vs Control (DESeq2)")
ggsave("figures/Step1_volcano.pdf", p_volcano, width = 7, height = 6)

top_genes <- deg_sig %>% arrange(adj.P.Val) %>% head(50) %>% pull(Gene) %>% toupper()
top_genes <- intersect(top_genes, rownames(expr_gene))
pdf("figures/Step1_heatmap_top50DEG.pdf", width = 8, height = 10)
pheatmap(expr_gene[top_genes, ], scale = "row", show_colnames = TRUE,
         annotation_col = data.frame(Group = group, row.names = colnames(expr_gene)),
         main = "Top 50 DEGs (row-scaled VST)")
dev.off()

#===============================================================================
# Step 2. 三重硬交集（策略A主候选池）
#   主候选池 = 人海马DEG ∩ gutMGene因果靶基因(human+causally) ∩ 神经免疫基因集
#   【本版本无小鼠KD加权层，主候选池直接进入后续步骤】
#===============================================================================

## ---- 2.1 gutMGene预处理（用户已下载的两个CSV） ---------------------------------
## 文件1：菌群-宿主基因（Microbe-Gene）；文件2：菌源代谢物-宿主基因（Metabolite-Gene）
gm_microbe    <- fread("data/Gut Microbe-Host Gene.csv")
gm_metabolite <- fread("data/Microbial metabolite-Host Gene.csv")

## 筛选条件：人源 + 因果关联
##   列 "human/mouse" 取值 human/mouse；列 "Associative mode" 取值 causally/correlatively
gm_microbe_f <- gm_microbe %>%
  filter(`human/mouse` == "human",
         grepl("causally", `Associative mode`, ignore.case = TRUE))
gm_metabolite_f <- gm_metabolite %>%
  filter(`human/mouse` == "human",
         grepl("causally", `Associative mode`, ignore.case = TRUE))

cat("gutMGene 菌群-基因因果条目：", nrow(gm_microbe_f),
    "；代谢物-基因因果条目：", nrow(gm_metabolite_f), "\n")

## 提取唯一Gene Symbol，取两类关联的并集作为gutMGene基因集
gutmgene_genes <- union(gm_microbe_f$Gene, gm_metabolite_f$Gene) %>%
  toupper() %>% unique() %>% na.omit()
cat("gutMGene 去重后靶基因数：", length(gutmgene_genes), "\n")
save_rds(gutmgene_genes, "Step2_gutMGene_geneset")
save_rds(list(microbe = gm_microbe_f, metabolite = gm_metabolite_f),
         "Step2_gutMGene_filtered_tables")

## ---- 2.2 神经免疫/炎症/血脑屏障基因集（msigdbr） --------------------------------
# 1. 装remotes（如果没有）
install.packages("remotes")

# 2. 卸载新版msigdbr
remove.packages("msigdbr")

# 3. 从CRAN镜像装旧版7.5.1（数据内置，无需联网下载MSigDB）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
remotes::install_version("msigdbr", version = "7.5.1", upgrade = "never")

# 4. 验证：重启R会话后（重要，避免已加载的旧dll干扰）
library(msigdbr)
library(dplyr)
packageVersion("msigdbr")   # 应显示 '7.5.1'
## 从MSigDB Hallmark + C5(GO BP) + C7(免疫) 中按关键词筛选
msig_all <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "H"),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
  msigdbr(species = "Homo sapiens", category = "C7")
)
kw <- paste0("IMMUN|INFLAM|CYTOKINE|TLR|NLRP|NFKB|MICROGLIA|ASTROCYTE|",
             "BLOOD_BRAIN_BARRIER|COMPLEMENT|INTERLEUKIN|TNF|INTERFERON|",
             "CHEMOKINE|NEUROINFLAM|GLIAL")
immune_sets  <- msig_all %>% filter(grepl(kw, gs_name, ignore.case = TRUE))
immune_genes <- unique(immune_sets$gene_symbol) %>% toupper()
cat("神经免疫基因集基因数：", length(immune_genes),
    "（来自", length(unique(immune_sets$gs_name)), "个基因集）\n")
save_rds(immune_genes, "Step2_immune_geneset")

## ---- 2.3 三重硬交集 ---------------------------------------------------------------
deg_genes   <- unique(toupper(deg_sig$Gene))
main_pool   <- Reduce(intersect, list(deg_genes, gutmgene_genes, immune_genes))
cat("====== 三重硬交集结果 ======\n")
cat("人海马DEG:", length(deg_genes),
    " ∩ gutMGene:", length(gutmgene_genes),
    " ∩ Immune:", length(immune_genes),
    " => 主候选池基因数:", length(main_pool), "\n")
print(main_pool)
writeLines(main_pool, "results/Step2_main_pool_triple_intersection.txt")
save_rds(main_pool, "Step2_main_pool")

## 交集可视化（Venn）
if (requireNamespace("ggvenn", quietly = TRUE)) {
  library(ggvenn)
  p_venn <- ggvenn(list(DEG_GSE134697 = deg_genes,
                        gutMGene_causal = gutmgene_genes,
                        Neuroimmune = immune_genes),
                   fill_color = c("#E64B35", "#4DBBD5", "#00A087"),
                   stroke_size = 0.5, set_name_size = 4)
  ggsave("figures/Step2_venn_triple_intersection.pdf", p_venn, width = 6, height = 6)
} else {
  message("提示：安装 ggvenn 包可绘制三元韦恩图")
}


#===============================================================================
# Step 3. WGCNA 共表达网络（基于GSE134697全基因组表达矩阵）
#   -> 找与耐药TLE表型最相关模块 -> 模块基因 ∩ 主候选池（二次交集）
#===============================================================================

## ---- 3.1 输入数据：方差最大的前5000个基因（降低计算量、提高稳定性） -------------
datExpr <- t(expr_gene)                     # WGCNA要求行=样本、列=基因
gsv     <- apply(datExpr, 2, var)
datExpr <- datExpr[, order(gsv, decreasing = TRUE)[1:min(5000, ncol(datExpr))]]

## 样本聚类检查离群样本
pdf("figures/Step3_sample_clustering.pdf", width = 10, height = 5)
plot(hclust(dist(datExpr), method = "average"),
     main = "Sample clustering (check outliers)", xlab = "", sub = "")
dev.off()  # 如有明显离群样本，请手动剔除后再继续

## 表型向量：TLE=1, Control=0
trait <- data.frame(TLE = as.numeric(group == "TLE"))
rownames(trait) <- rownames(datExpr)

## ---- 3.2 软阈值选择 --------------------------------------------------------------
library(WGCNA)
allowWGCNAThreads()
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                         corFnc = "bicor", verbose = 0)
pdf("figures/Step3_soft_threshold.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold power", ylab = "Scale-free fit (signed R^2)", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red"); abline(h = 0.8, col = "blue", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
par(mfrow = c(1, 1)); dev.off()
## 自动选取：第一个使signed R^2 >= 0.8的power；若无则取12（signed网络经验值）
soft_power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq >= 0.8)[1]]
if (is.na(soft_power)) soft_power <- 12
cat("选定软阈值 power =", soft_power, "\n")

## ---- 3.3 一步法网络构建 + 模块检测 -------------------------------------------------
net <- blockwiseModules(datExpr, power = soft_power, networkType = "signed",
                        corType = "bicor", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25,
                        numericLabels = TRUE, saveTOMs = FALSE, verbose = 0)
module_colors <- labels2colors(net$colors)
save_rds(net, "Step3_WGCNA_network")

## 模块-表型相关性
MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
mod_trait_cor  <- cor(MEs, trait, use = "p", method = "pearson")
mod_trait_pval <- corPvalueStudent(mod_trait_cor, nrow(datExpr))

pdf("figures/Step3_module_trait_heatmap.pdf", width = 6, height = 8)
labeledHeatmap(Matrix = mod_trait_cor, xLabels = "TLE", yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = paste0(signif(mod_trait_cor, 2), "\n(",
                                   signif(mod_trait_pval, 1), ")"),
               setStdMargins = FALSE, cex.text = 0.7,
               main = "Module-trait relationships")
dev.off()

## ---- 3.4 选择最相关模块，与主候选池二次交集 -----------------------------------------
## 规则：|cor|最大且P<0.05的模块（排除grey模块）
cor_df <- data.frame(module = rownames(mod_trait_cor),
                     cor = mod_trait_cor[, 1], p = mod_trait_pval[, 1]) %>%
  filter(module != "MEgrey") %>% arrange(desc(abs(cor)))
print(cor_df)
top_module <- cor_df$module[1]
cat("与耐药TLE最相关模块：", top_module,
    "（cor =", round(cor_df$cor[1], 3), ", p =", signif(cor_df$p[1], 3), "）\n")

top_color      <- sub("ME", "", top_module)
module_genes   <- colnames(datExpr)[module_colors == top_color]
module_genes   <- toupper(module_genes)
cat("模块基因数：", length(module_genes), "\n")

final_candidates <- intersect(main_pool, module_genes)
cat("====== WGCNA二次交集 ======\n")
cat("主候选池:", length(main_pool), " ∩ 模块基因:", length(module_genes),
    " => 进入PPI的基因数:", length(final_candidates), "\n")
print(final_candidates)
writeLines(final_candidates, "results/Step3_final_candidates_WGCNA_intersect.txt")
save_rds(final_candidates, "Step3_final_candidates")
## 兜底：若二次交集为空，则回退使用主候选池（保证流程可继续，并在报告中注明）
if (length(final_candidates) == 0) {
  warning("WGCNA模块与主候选池无交集！回退使用主候选池进入PPI。")
  final_candidates <- main_pool
  save_rds(final_candidates, "Step3_final_candidates")
}


#===============================================================================
# Step 4. PPI网络 + cytoHubba四种算法Hub锁定（无KD加权，直接Composite Score排序）
#===============================================================================
library(STRINGdb)
## ---- 4.1 STRING网络（score_threshold = 700，高置信度） ----------------------------
string_db <- STRINGdb$new(version = "12.0", species = 9606,
                          score_threshold = 700, input_directory = "data")
cand_df   <- data.frame(gene = final_candidates)
mapped    <- string_db$map(cand_df, "gene", removeUnmappedRows = TRUE)
interactions <- string_db$get_interactions(mapped$STRING_id)
cat("STRING互作边数：", nrow(interactions), "\n")
save_rds(list(mapped = mapped, interactions = interactions), "Step4_STRING")

## 构建igraph对象（边权 = STRING combined score）
id2sym <- setNames(mapped$gene, mapped$STRING_id)
edges  <- interactions %>%
  mutate(from = id2sym[from], to = id2sym[to]) %>%
  filter(!is.na(from), !is.na(to))
library(igraph)
g <- graph_from_data_frame(edges[, c("from", "to", "combined_score")],
                           directed = FALSE)
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
E(g)$weight <- E(g)$combined_score
nrow(interactions)   # STRING返回的互作边数
ecount(g)            # igraph里的边数
vcount(g)            # 节点数
nrow(mapped)         # 成功映射到STRING的基因数
## 重建网络：先赋权重，再简化
g <- graph_from_data_frame(edges[, c("from", "to", "combined_score")],
                           directed = FALSE)
E(g)$weight <- E(g)$combined_score
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE,
              edge.attr.comb = "first")

## 确认一下
ecount(g); vcount(g); edge_attr_names(g)   # 应能看到 weight 属性
## ---- 4.2 cytoHubba四种算法的R实现 -------------------------------------------------
## (1) Degree：节点度数
score_degree <- degree(g)

## (2) MCC (Maximal Clique Centrality)：节点在所有极大团中的得分之和
##     MCC(v) = sum over maximal cliques C containing v of (|C|-1)!
all_cliques <- max_cliques(g)
score_mcc <- setNames(rep(0, vcount(g)), V(g)$name)
for (cl in all_cliques) {
  w <- factorial(length(cl) - 1)
  for (v in names(cl)) score_mcc[v] <- score_mcc[v] + w
}

## (3) BottleNeck：介数中心性（betweenness）
score_bottleneck <- betweenness(g, directed = FALSE, normalized = TRUE)

## (4) EPC (Edge Percolated Component)：边渗透分量
##     简化实现：随机删除一定比例的边，重复多次，
##     计算每个节点与网络最大连通分量的连接保持率（近似EPC得分）
set.seed(123)
n_sim <- 1000; remove_frac <- 0.5
nodes <- V(g)$name
epc_mat <- matrix(0, nrow = length(nodes), ncol = n_sim,
                  dimnames = list(nodes, NULL))
for (i in seq_len(n_sim)) {
  g_sub <- delete_edges(g, sample(E(g), floor(ecount(g) * remove_frac)))
  comp  <- components(g_sub)
  giant <- names(comp$membership)[comp$membership ==
                                    which.max(comp$csize)]
  epc_mat[giant, i] <- 1
}
score_epc <- rowMeans(epc_mat)  # 节点留在最大连通分量中的频率

## ---- 4.3 标准化到0-1并计算Composite Score ------------------------------------------
norm01 <- function(x) if (max(x) == min(x)) rep(0.5, length(x)) else
  (x - min(x)) / (max(x) - min(x))

hub_scores <- data.frame(
  Gene       = nodes,
  MCC        = norm01(score_mcc[nodes]),
  Degree     = norm01(score_degree[nodes]),
  BottleNeck = norm01(score_bottleneck[nodes]),
  EPC        = norm01(score_epc[nodes])
) %>%
  mutate(Composite_Score = (MCC + Degree + BottleNeck + EPC) / 4) %>%
  arrange(desc(Composite_Score))
library(data.table)
fwrite(hub_scores, "results/Step4_hub_scores_all.csv")
save_rds(hub_scores, "Step4_hub_scores")

## 取Top 5为最终Hub（本版本无KD加权，直接按Composite Score降序）
top_n <- min(5, nrow(hub_scores))
hub_genes <- hub_scores$Gene[seq_len(top_n)]
cat("====== 最终Hub基因（Top", top_n, "）======\n")
print(hub_scores[seq_len(top_n), ])
fwrite(hub_scores[seq_len(top_n), ], "results/Step4_Top5_hub_genes.csv")
save_rds(hub_genes, "Step4_hub_genes")

## ---- 4.4 可视化：PPI网络图 + Hub得分条形图 ------------------------------------------
## PPI网络图（Top Hub高亮红色）
V(g)$color <- ifelse(V(g)$name %in% hub_genes, "#E64B35", "#4DBBD5")
V(g)$size  <- 4 + 8 * norm01(score_degree[V(g)$name])
pdf("figures/Step4_PPI_network.pdf", width = 9, height = 9)
set.seed(42)
plot(g, layout = layout_with_fr(g), vertex.label.cex = 0.7,
     vertex.label.color = "black", edge.width = E(g)$weight / 200,
     main = "STRING PPI network (red = Top5 hubs)")
dev.off()

## Hub得分条形图（四种算法+综合分，长格式）
hub_long <- hub_scores[seq_len(top_n), ] %>%
  tidyr::pivot_longer(cols = c(MCC, Degree, BottleNeck, EPC, Composite_Score),
                      names_to = "Method", values_to = "Score")
p_hub <- ggplot(hub_scores[seq_len(top_n), ],
                aes(x = reorder(Gene, Composite_Score), y = Composite_Score)) +
  geom_bar(stat = "identity", fill = "#00A087", width = 0.6) +
  coord_flip() + theme_bw() +
  labs(title = "Top hub genes by Composite Score",
       x = NULL, y = "Composite Score (mean of 4 normalized algorithms)")
ggsave("figures/Step4_hub_score_barplot.pdf", p_hub, width = 6, height = 4)


#===============================================================================
# Step 5. Hub功能注释（GO/KEGG）+ gutMGene反向查询（精确三联体）
#===============================================================================

## ---- 5.1 GO/KEGG富集 ---------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
hub_entrez <- bitr(hub_genes, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_bp <- enrichGO(gene = hub_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, readable = TRUE)
kegg <- enrichKEGG(gene = hub_entrez$ENTREZID, organism = "hsa",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05)
save_rds(list(GO_BP = go_bp, KEGG = kegg), "Step5_enrichment")
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  fwrite(go_bp@result, "results/Step5_GO_BP_hub.csv")
  pdf("figures/Step5_GO_BP_dotplot.pdf", width = 8, height = 6)
  print(dotplot(go_bp, showCategory = 15) + ggtitle("GO BP enrichment of hub genes"))
  dev.off()
}
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  fwrite(kegg@result, "results/Step5_KEGG_hub.csv")
  pdf("figures/Step5_KEGG_dotplot.pdf", width = 8, height = 6)
  print(dotplot(kegg, showCategory = 15) + ggtitle("KEGG enrichment of hub genes"))
  dev.off()
}

## ---- 5.2 gutMGene反向查询：精确三联体 -------------------------------------------------
## 对每个Hub，从筛选后的人源因果表中回查：
##   代谢物表：Substrate -> Metabolite -> Gene -> Alteration（调控方向）
##   菌群表：  Gut Microbiota(Strain) -> Gene -> Alteration
rev_met <- gm_metabolite_f %>%
  filter(toupper(Gene) %in% toupper(hub_genes)) %>%
  transmute(Type = "Metabolite-Gene",
            Microbe = NA_character_,
            Substrate, Metabolite, Gene,
            Regulation = Alteration,
            Sample, Condition, PMID)

rev_mic <- gm_microbe_f %>%
  filter(toupper(Gene) %in% toupper(hub_genes)) %>%
  transmute(Type = "Microbe-Gene",
            Microbe = ifelse(is.na(Strain) | Strain == "",
                             `Gut Microbiota`,
                             paste0(`Gut Microbiota`, " (", Strain, ")")),
            Substrate = NA_character_, Metabolite = NA_character_,
            Gene, Regulation = Alteration,
            Sample, Condition, PMID)

gutmgene_reverse <- bind_rows(rev_met, rev_mic) %>% arrange(Gene, Type)
cat("====== Hub的gutMGene反向查询条目数：", nrow(gutmgene_reverse), "======\n")
print(gutmgene_reverse)
fwrite(gutmgene_reverse, "results/Step5_gutMGene_reverse_query_hub.csv")
save_rds(gutmgene_reverse, "Step5_gutMGene_reverse")

## ---- 5.3 汇总报告 ---------------------------------------------------------------------
sink("results/FINAL_SUMMARY.txt")
cat("课题：肠道代谢-神经免疫-血脑屏障交互调控耐药颞叶癫痫（纯生信）\n")
cat("版本：v2（仅GSE134697；无小鼠KD加权；gutMGene限human+causally）\n\n")
cat("Step1 显著DEG数：", nrow(deg_sig), "\n")
cat("Step2 三重硬交集主候选池：", length(main_pool), "个：",
    paste(main_pool, collapse = ", "), "\n")
cat("Step3 WGCNA最相关模块：", top_module,
    "；二次交集后基因：", length(final_candidates), "个：",
    paste(final_candidates, collapse = ", "), "\n")
cat("Step4 最终Hub（Top", top_n, "，按Composite Score）：",
    paste(hub_genes, collapse = ", "), "\n")
cat("Step5 gutMGene反向查询条目：", nrow(gutmgene_reverse), "\n")
sink()
cat("\n分析完成！结果见 results/，图片见 figures/，中间对象见 rds/\n")
deg_sig <- readRDS("rds/Step1_DEG_sig.rds")
cat("DEG:", nrow(deg_sig), "上调:", sum(deg_sig$logFC > 0),
    "下调:", sum(deg_sig$logFC < 0), "\n")

net <- readRDS("rds/Step3_WGCNA_network.rds")
table(labels2colors(net$colors))   # 各模块基因数
print(sft$fitIndices[, c("Power", "SFT.R.sq")])
## ============================================================
## Bulk-scRNA整合图：5个Hub × 细胞类型 dotplot
## 点大小 = TLE组表达细胞比例(pct.1)，颜色 = 组间log2FC
## ============================================================
