## ======================================================================
## 耐药TLE多组学分析流程 · 脚本08/17：GSE143272 外周血PBMC验证：5Hub表达 + TLR4分组差异 + ROC
## 对应论文：3.7节（图5b）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

################################################################################
## 模块四：外周PBMC验证（GSE143272）
## 目标：检验5个Hub（TGFB1/ICAM1/TLR4/TLR2/CCL2）在癫痫患者外周血中的表达
## 输入：file.path(cfg$dir_pbmc, "GSE143272_non-normalized (1).txt.gz")
################################################################################
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(tidyr)
})

workdir <- cfg$dir_pbmc
setwd(workdir)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#===============================================================================
# 第0步：结构侦察（先跑这段，把输出发我）
#===============================================================================
pbmc_raw <- fread("GSE143272_non-normalized (1).txt.gz", data.table = FALSE, nrows = 5)
print(dim(pbmc_raw))
print(pbmc_raw[1:5, 1:min(8, ncol(pbmc_raw))])
## 重点看：
## (1) 第一列是什么ID（Symbol / Ensembl / Entrez / 探针ID）
## (2) 数值是整数(counts)还是小数(FPKM/TPM/信号值)
## (3) 列名长什么样（GSM号？样本名？含不含分组信息）

#===============================================================================
# 第1步：读取 + 提取信号列（剔除Detection Pval列）
#===============================================================================
pbmc <- fread("GSE143272_non-normalized (1).txt.gz", data.table = FALSE)
probe_ids <- pbmc$PROBE_ID
sig_cols <- !grepl("Detection Pval", colnames(pbmc)) & colnames(pbmc) != "PROBE_ID"
expr <- as.matrix(pbmc[, sig_cols])
rownames(expr) <- probe_ids
cat("探针数：", nrow(expr), "；样本数：", ncol(expr), "\n")
print(colnames(expr))   # 确认全部样本名（应该还有VA Responder组）

#===============================================================================
# 第2步：探针注释（ILMN -> Symbol）+ 分组
#===============================================================================
## 2.1 注释：需要illuminaHumanv4.db（HT-12 v4平台）
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("illuminaHumanv4.db", update = FALSE)
library(illuminaHumanv4.db)
sym_map <- mapIds(illuminaHumanv4.db, keys = probe_ids,
                  column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
keep <- !is.na(sym_map)
expr <- expr[keep, ]; sym_map <- sym_map[keep]
## 同Symbol多探针 -> 取均值表达最高的探针
avg <- rowMeans(expr)
keep_probe <- tapply(seq_along(sym_map), sym_map, function(idx) idx[which.max(avg[idx])])
expr_gene <- expr[unlist(keep_probe), ]
rownames(expr_gene) <- names(keep_probe)
cat("注释去重后基因数：", nrow(expr_gene), "\n")

## 2.2 标准化：log2 + 分位数标准化（芯片非标准化数据的标准流程）
library(limma)
expr_norm <- normalizeBetweenArrays(log2(expr_gene + 1), method = "quantile")

## 2.3 分组（完整版：7个组）
sample_names <- colnames(expr_norm)
group <- case_when(
  grepl("Healthy Control", sample_names) ~ "Control",
  substr(sample_names, 1, 4) == "Drug"   ~ "DrugNaive",     # 用前缀匹配，避开编码问题
  grepl("Non-responder", sample_names)   ~ "NonResponder",
  grepl("Responder", sample_names)       ~ "Responder",
  TRUE ~ NA_character_)
group <- factor(group, levels = c("Control", "DrugNaive", "Responder", "NonResponder"))
print(table(group))   # 应输出 Control 51 / DrugNaive 34 / Responder 33 / NonResponder 24

## 剔除技术重复样本（"Healthy Control-3 (Replicate)"）
dup_idx <- grepl("Replicate", sample_names)
if (any(dup_idx)) {
  expr_norm <- expr_norm[, !dup_idx]
  group <- group[!dup_idx]
}
print(table(group))   # Control应变50

#===============================================================================
# 第3步：全基因组差异分析（NonResponder vs Control，主对比）
#===============================================================================
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_norm, design)
contr <- makeContrasts(
  NR_vs_C  = NonResponder - Control,
  NR_vs_R  = NonResponder - Responder,
  DN_vs_C  = DrugNaive - Control,
  levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

for (cf in c("NR_vs_C", "NR_vs_R", "DN_vs_C")) {
  deg <- topTable(fit2, coef = cf, number = Inf, adjust.method = "BH") %>%
    tibble::rownames_to_column("Gene")
  fwrite(deg, paste0("results/PBMC_DEG_", cf, ".csv"))
  cat(cf, " |logFC|>0.58 & FDR<0.05:",
      deg %>% filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% nrow(), "\n")
}
#===============================================================================
# 第4步：5个Hub的组间比较（核心分析）
#===============================================================================
hubs <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2")
deg_nr <- fread("results/PBMC_DEG_NR_vs_C.csv")
deg_dn <- fread("results/PBMC_DEG_DN_vs_C.csv")

hub_check <- data.frame(Gene = hubs) %>%
  left_join(deg_nr %>% dplyr::select(Gene, logFC_NR = logFC, pval_NR = P.Value, fdr_NR = adj.P.Val), by = "Gene") %>%
  left_join(deg_dn %>% dplyr::select(Gene, logFC_DN = logFC, pval_DN = P.Value, fdr_DN = adj.P.Val), by = "Gene") %>%
  mutate(bulk_direction = "Up",   # bulk海马里5个Hub全部上调
         consistent_NR = ifelse(logFC_NR > 0, "一致", "不一致"),
         consistent_DN = ifelse(logFC_DN > 0, "一致", "不一致"))
print(hub_check)
fwrite(hub_check, "results/PBMC_hub_validation_summary.csv")
## ★ 与bulk海马方向核对：bulk里5个Hub全是上调，这里direction=Up即为方向一致


#===============================================================================
# 第5步：可视化
#===============================================================================
## 5.1 Hub表达箱线图
hubs <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2")
hubs_in <- intersect(hubs, rownames(expr_norm))

hub_df <- lapply(hubs_in, function(g) {
  data.frame(Gene = g, Expression = expr_norm[g, ], Group = group)
}) %>% bind_rows()

p_box <- ggplot(hub_df, aes(Group, Expression, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Control" = "#3C5488", "DrugNaive" = "#E8A33D",
                               "Responder" = "#0E8A7D", "NonResponder" = "#E64B35")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Hub gene expression in PBMC (GSE143272)", x = NULL)
ggsave("figures/PBMC_hub_boxplot.pdf", p_box, width = 10, height = 4.5)
## 5.2 ROC曲线（生物标志物潜力）
if (requireNamespace("pROC", quietly = TRUE)) {
  library(pROC)
  roc_list <- lapply(hubs_in, function(g) {
    roc(group, expr_norm[g, ], quiet = TRUE)
  })
  aucs <- sapply(roc_list, function(r) round(as.numeric(auc(r)), 3))
  names(aucs) <- hubs_in
  print(aucs)
  pdf("figures/PBMC_hub_ROC.pdf", width = 6, height = 6)
  plot(roc_list[[1]], col = 1, main = "ROC of hub genes (PBMC)")
  for (i in seq_along(roc_list)[-1]) plot(roc_list[[i]], col = i, add = TRUE)
  legend("bottomright", legend = paste0(hubs_in, " (AUC=", aucs, ")"),
         col = seq_along(roc_list), lwd = 2, cex = 0.8)
  dev.off()
  fwrite(data.frame(Gene = hubs_in, AUC = aucs), "results/PBMC_hub_AUC.csv")
} else install.packages("pROC")  # 装完重跑5.2

cat("\n完成！结果见", workdir, "/results 和 /figures\n")
deg_dn <- fread("results/PBMC_DEG_DN_vs_C.csv") %>%
  mutate(sig = case_when(adj.P.Val < 0.05 & logFC >  0.58 ~ "Up",
                         adj.P.Val < 0.05 & logFC < -0.58 ~ "Down", TRUE ~ "NS"))
p_vol <- ggplot(deg_dn, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NS" = "grey80")) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(data = deg_dn %>% filter(Gene %in% c("TLR4","TLR2","CCL2","ICAM1")),
                           aes(label = Gene), size = 4, color = "black") +
  theme_bw() + labs(title = "PBMC: drug-naive epilepsy vs control")
ggsave("figures/PBMC_volcano_DN_vs_C.pdf", p_vol, width = 6.5, height = 5.5)
## ============================================================
## TLR4（及其余Hub）ROC分析：DrugNaive vs Control
## 在 PBMC 脚本环境中运行（需要 expr_norm 和 group 对象）
## ============================================================
install.packages("pROC")   # 没装过的话先装
library(pROC)

## 只用 DrugNaive vs Control（未用药组，最干净的疾病信号）
idx <- group %in% c("Control", "DrugNaive")
group_roc <- droplevels(group[idx])
expr_roc  <- expr_norm[, idx]

roc_res <- lapply(hubs_in, function(g) {
  r <- roc(group_roc, expr_roc[g, ], levels = c("Control", "DrugNaive"), quiet = TRUE)
  ci <- ci.auc(r)
  data.frame(Gene = g,
             AUC = round(as.numeric(auc(r)), 3),
             CI_low = round(ci[1], 3), CI_high = round(ci[3], 3))
})
roc_df <- bind_rows(roc_res) %>% arrange(desc(AUC))
print(roc_df)
fwrite(roc_df, "results/PBMC_hub_ROC_AUC.csv")

## 画图：4个Hub的ROC曲线
cols <- c("TLR4" = "#E64B35", "CCL2" = "#E8A33D", "ICAM1" = "#3C5488", "TLR2" = "#0E8A7D")
pdf("figures/PBMC_hub_ROC.pdf", width = 6, height = 6)
first <- TRUE
for (g in hubs_in) {
  r <- roc(group_roc, expr_roc[g, ], levels = c("Control", "DrugNaive"), quiet = TRUE)
  a <- round(as.numeric(auc(r)), 3)
  if (first) { plot(r, col = cols[g], lwd = 2.5,
                    main = "ROC: hubs in PBMC (drug-naive vs control)")
    first <- FALSE
  } else plot(r, col = cols[g], lwd = 2.5, add = TRUE)
}
abline(0, 1, lty = 2, col = "grey60")
legend("bottomright",
       legend = paste0(roc_df$Gene, " (AUC = ", roc_df$AUC, ")"),
       col = cols[roc_df$Gene], lwd = 2.5, cex = 0.9, bty = "n")
dev.off()
cat("ROC完成\n")
