## ======================================================================
## 耐药TLE多组学分析流程 · 脚本03/17：SMR + HEIDI：5个Hub基因顺式eQTL（BrainMeta）与局灶性癫痫
## 对应论文：3.6节（补充表S3-C）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

## ============================================================================
## SMR + HEIDI 分析完整流程
## 目的：检验 5 个 Hub 基因（TGFB1/ICAM1/TLR4/TLR2/CCL2）在 TLE 中的上调
##       是否由生殖系顺式调控变异驱动
## 数据：结局  ILAE 2018 focal epilepsy GWAS（OpenGWAS: ieu-b-10, hg19）
##       暴露  BrainMeta v2 cis-eQTL（n=2,865 皮层，BESD 格式）
##             GTEx v8 Brain Hippocampus cis-eQTL（敏感性分析）
##       LD    1000 Genomes Phase 3（n=2,490，plink 格式）
## 目录：cfg$dir_smr（见 config/config.R）
## ============================================================================

library(data.table)
library(ggplot2)

setwd(cfg$dir_smr)

## ============================================================================
## 第 1 部分：ieu-b-10 VCF -> 5 个 Hub 基因区域的 SMR .ma 文件
## ============================================================================

N <- 39348   # ieu-b-10 总样本量（9,671 cases + 29,677 controls）

## Hub 基因区域（GRCh37/hg19 坐标，±1 Mb）
regions <- data.frame(
  gene  = c("TGFB1",  "ICAM1",   "TLR4",     "TLR2",     "CCL2"),
  chr   = c(19,        19,        9,          4,          17),
  start = c(41347791,  10271120,  120463118,  154612652,  32581971) - 1e6,
  end   = c(41356991,  10281807,  120493118,  154633652,  32591971) + 1e6
)

## 读取 OpenGWAS VCF（仅保留必需列；FORMAT = ES:SE:LP:AF）
vcf <- fread("ieu-b-10.vcf.gz", skip = "#CHROM",
             select = c("#CHROM", "POS", "ID", "REF", "ALT", "ieu-b-10"))
cat("VCF 总行数:", nrow(vcf), "\n")

parts <- tstrsplit(vcf[["ieu-b-10"]], ":", fixed = TRUE)
ma <- data.frame(
  CHR  = vcf[["#CHROM"]],
  BP   = vcf$POS,
  SNP  = vcf$ID,
  A1   = toupper(vcf$ALT),          # 效应等位基因 = ALT
  A2   = toupper(vcf$REF),
  freq = as.numeric(parts[[4]]),
  b    = as.numeric(parts[[1]]),
  se   = as.numeric(parts[[2]]),
  p    = 10^(-as.numeric(parts[[3]])),
  N    = N
)
rm(vcf, parts); gc()
ma <- as.data.table(ma)

## 裁剪到 Hub 区域
ma_hub <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
  r <- regions[i, ]
  ma[CHR == r$chr & BP >= r$start & BP <= r$end][, gene := r$gene]
}))
cat("Hub 区域 SNP 数:", nrow(ma_hub), "\n")
print(table(ma_hub$gene))

## 写出 SMR .ma 格式：SNP A1 A2 freq b se p N
fwrite(ma_hub[, .(SNP, A1, A2, freq, b, se, p, N)],
       "epilepsy_focal_hub.ma", sep = "\t")
cat("已生成 epilepsy_focal_hub.ma\n")

## ============================================================================
## 第 2 部分：SMR + HEIDI（Windows cmd 命令，非 R 代码，记录备查）
## 注意：SMR 的 --plot 数据文件中也含探针版本号（如 ENSG00000136869.16）
## ============================================================================
##
## 在 cfg$dir_smr 目录下执行 SMR 命令行（Windows 示例：cd /d <cfg$dir_smr 实际路径>）
##
## --- 主分析：BrainMeta v2（按染色体） ---
## smr-1.3.1-win.exe --bfile 1000G_phase3_common_norel/1000G_phase3_common_norel ^
##   --gwas-summary epilepsy_focal_hub.ma ^
##   --beqtl-summary BrainMeta_cis_eqtl_summary/.../BrainMeta_cis_eQTL_chr19 ^
##   --maf 0.01 --diff-freq-prop 0.9 --out smr_chr19 --thread-num 4
## （chr19 / chr9 / chr4 / chr17 各跑一条）
##
## --- 敏感性分析：GTEx v8 海马 ---
## smr-1.3.1-win.exe --bfile 1000G_phase3_common_norel/1000G_phase3_common_norel ^
##   --gwas-summary epilepsy_focal_hub.ma ^
##   --beqtl-summary Brain_Hippocampus/Brain_Hippocampus/Brain_Hippocampus ^
##   --maf 0.01 --diff-freq-prop 0.9 --out smr_hippo --thread-num 4
##
## --- 基因座图数据（TLR4） ---
## smr-1.3.1-win.exe --bfile 1000G_phase3_common_norel/1000G_phase3_common_norel ^
##   --gwas-summary epilepsy_focal_hub.ma ^
##   --beqtl-summary BrainMeta_cis_eqtl_summary/.../BrainMeta_cis_eQTL_chr9 ^
##   --out plot_TLR4 --plot --probe ENSG00000136869.16 --probe-wind 500 ^
##   --gene-list plot/plot/glist_hg19_refseq.txt --diff-freq-prop 0.9 --thread-num 4

## ============================================================================
## 第 3 部分：SMR 结果合并与判读
## ============================================================================

smr_all <- rbindlist(lapply(paste0("smr_chr", c(19, 9, 4, 17), ".smr"), fread))
cat("BrainMeta 总探针数:", nrow(smr_all), "\n")

hubs <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2")

hub_smr <- smr_all[Gene %in% hubs]
hub_smr[, verdict := ifelse(p_SMR < 0.05 & p_HEIDI > 0.05, "阳性(共享因果变异)",
                     ifelse(p_SMR < 0.05, "连锁,不可信", "阴性"))]
print(hub_smr[, .(probeID, Gene, topSNP, b_SMR, se_SMR, p_SMR, p_HEIDI, verdict)])

fwrite(smr_all, "SMR_BrainMeta_focalepi_all_probes.csv")
fwrite(hub_smr, "SMR_BrainMeta_focalepi_hubs.csv")

## 敏感性分析：GTEx 海马
hippo <- fread("smr_hippo.smr")
hub_hippo <- hippo[Gene %in% hubs]
hub_hippo[, verdict := ifelse(p_SMR < 0.05 & p_HEIDI > 0.05, "阳性(共享因果变异)",
                       ifelse(p_SMR < 0.05, "连锁,不可信", "阴性"))]
print(hub_hippo[, .(probeID, Gene, topSNP, b_SMR, p_SMR, p_HEIDI, verdict)])
fwrite(hub_hippo, "SMR_GTEx_hippocampus_hubs.csv")

## ============================================================================
## 第 4 部分：可视化
## ============================================================================

## ---- Fig 1: 5 Hub 的 SMR 状态森林图 ----
bm2 <- merge(data.frame(Gene = hubs),
             hub_smr[, .(Gene, b_SMR, se_SMR, p_SMR)],
             by = "Gene", all.x = TRUE)

p1 <- ggplot(bm2, aes(b_SMR, Gene)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = b_SMR - 1.96 * se_SMR,
                     xmax = b_SMR + 1.96 * se_SMR), height = 0.2, na.rm = TRUE) +
  geom_point(size = 3.5, na.rm = TRUE, color = "#4DBBD5") +
  geom_text(aes(label = ifelse(is.na(p_SMR), "no cis-eQTL",
                        sprintf("p_SMR = %.2f", p_SMR)), x = 0.09),
            hjust = 0, size = 3.5) +
  xlim(-0.08, 0.16) +
  theme_classic(base_size = 12) +
  labs(x = "b_SMR (brain expression -> focal epilepsy)", y = NULL,
       title = "SMR: hub gene expression vs focal epilepsy (BrainMeta v2)")
ggsave("Fig_SMR_hubs_forest.pdf", p1, width = 7, height = 3.5)

## ---- Fig 2/3: TLR4 基因座图 + 效应对比图（自绘版，不依赖官方 plot_SMR.r）----
## 解析 SMR --plot 生成的多段文本文件（$probe/$SNP/$GWAS/$eQTL/$Gene）
parse_smr_plot <- function(file, target_probe) {
  txt <- readLines(file)
  sec <- grep("^\\$", txt)
  ## $SNP 块：SNP Chr BP A1 A2
  i_snp <- grep("^\\$SNP", txt)
  n_snp <- as.integer(strsplit(txt[i_snp], " ")[[1]][2])
  snp <- fread(text = txt[(i_snp + 1):(i_snp + n_snp)], header = FALSE,
               col.names = c("SNP", "Chr", "BP", "A1", "A2"))
  ## $GWAS 块：SNP b se（p 缺失，由 z = b/se 重算）
  i_g <- grep("^\\$GWAS", txt)
  n_g <- as.integer(strsplit(txt[i_g], " ")[[1]][2])
  gwas <- fread(text = txt[(i_g + 1):(i_g + n_g)], header = FALSE,
                col.names = c("SNP", "b_GWAS", "se_GWAS"))
  gwas[, p_GWAS := 2 * pnorm(-abs(b_GWAS / se_GWAS))]
  ## $eQTL 块：按探针分块，取目标探针
  i_e <- grep("^\\$eQTL", txt)
  e_lines <- txt[(i_e + 1):grep("^\\$Gene", txt) - 1]
  probe_head <- grep("^ENSG", e_lines)
  eq <- NULL
  for (h in probe_head) {
    pid <- strsplit(e_lines[h], " ")[[1]][1]
    n   <- as.integer(strsplit(e_lines[h], " ")[[1]][2])
    if (pid == target_probe) {
      eq <- fread(text = e_lines[(h + 1):(h + n)], header = FALSE,
                  col.names = c("SNP", "b_eQTL", "se_eQTL", "p_na"))
      eq[, p_eQTL := 2 * pnorm(-abs(b_eQTL / se_eQTL))]
      break
    }
  }
  pos <- setNames(snp$BP, snp$SNP)
  gwas[, BP := pos[SNP]]; eq[, BP := pos[SNP]]
  list(gwas = gwas, eqtl = eq)
}

pd <- parse_smr_plot("plot/plot_TLR4.ENSG00000136869.16.txt",
                     "ENSG00000136869.16")
top <- "rs1329057"

## Fig 2: 双轨基因座图
library(gridExtra)
g1 <- ggplot(pd$gwas, aes(BP / 1e6, -log10(p_GWAS))) +
  geom_point(size = 0.8, color = "#4DBBD5", alpha = 0.7) +
  geom_point(data = pd$gwas[SNP == top], color = "#E64B35", size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60") +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "GWAS -log10(p)",
       title = "TLR4 locus: focal epilepsy GWAS vs BrainMeta eQTL (SMR)")
g2 <- ggplot(pd$eqtl, aes(BP / 1e6, -log10(p_eQTL))) +
  geom_point(size = 0.8, color = "#00A087", alpha = 0.7) +
  geom_point(data = pd$eqtl[SNP == top], color = "#E64B35", size = 2.5) +
  geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "grey60") +
  theme_classic(base_size = 11) +
  labs(x = "Chromosome 9 position (Mb)", y = "eQTL -log10(p)")
pdf("Fig_SMR_locus_TLR4.pdf", width = 8, height = 7)
grid.arrange(g1, g2, ncol = 1)
dev.off()

## Fig 3: eQTL-GWAS 效应对比图
m <- merge(pd$gwas, pd$eqtl[, .(SNP, b_eQTL, se_eQTL)], by = "SNP")
p3 <- ggplot(m, aes(b_eQTL, b_GWAS)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(size = 1.2, alpha = 0.35, color = "#4DBBD5") +
  geom_point(data = m[SNP == top], color = "#E64B35", size = 3) +
  theme_classic(base_size = 12) +
  labs(x = "eQTL effect (TLR4 expression)",
       y = "GWAS effect (focal epilepsy)",
       title = "No concordance between eQTL and GWAS effects at TLR4 locus")
ggsave("Fig_SMR_effect_TLR4.pdf", p3, width = 6, height = 6)

cat("SMR 全部分析与作图完成\n")
