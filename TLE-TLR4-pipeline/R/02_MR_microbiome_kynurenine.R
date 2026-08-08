## ======================================================================
## 耐药TLE多组学分析流程 · 脚本02/17：TwoSampleMR：菌群MR扫描 + 血浆犬尿氨酸->局灶性癫痫
## 对应论文：3.6节（图4a/b、补充图S6、补充表S3）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

################################################################################
## 完整版：TwoSampleMR 环境配置 + 菌群MR分析
## 使用方法：分段运行（第①段配置好后必须重启R再跑第②段）
################################################################################

## ===================== 第①段：装包 + 配token（跑一次即可） =====================

# 1. 装包（已装过会自动跳过）
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  install.packages("TwoSampleMR",
                   repos = c("https://mrcieu.r-universe.dev",
                             "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
}

# 2. token配置（安全版，勿将token写入代码或上传GitHub）：
#    在 api.opengwas.io 生成token后写入本地 ~/.Renviron：
#      OPENGWAS_JWT="你的token"
#    重启R生效；下面的检查用于确认配置
if (Sys.getenv("OPENGWAS_JWT") == "") {
  message("未检测到 OPENGWAS_JWT：请在 ~/.Renviron 配置后重启R（见 config/config.R 注释与 README）")
}

# 3. ★★★ 重启R会话（Session → Restart R）★★★


## ===================== 第②段：验证 + 捞目标菌（重启后运行） =====================

suppressPackageStartupMessages({
  library(TwoSampleMR); library(ieugwasr)
  library(data.table); library(dplyr); library(ggplot2)
})
setwd(cfg$dir_bulk)
dir.create("results/MR", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/MR", showWarnings = FALSE, recursive = TRUE)

# 4. 测试token是否生效（不报401即成功）
mb_ids <- paste0("ebi-a-GCST900", 16908:17118)
mb_info <- ieugwasr::gwasinfo(mb_ids)
cat("成功获取", nrow(mb_info), "个MiBioGen数据集信息\n")

# 5. 筛选目标菌（你的gutMGene三联体相关菌）
target <- mb_info[grepl("Akkermansia|Faecalibacterium|Bacteroides|Escherichia|Clostridiales",
                        mb_info$trait, ignore.case = TRUE), ]
print(target[, c("id", "trait", "sample_size")])
fwrite(target, "results/MR/target_taxa_info.csv")

# 6. 检索癫痫结局GWAS，人工选定
epi_search <- available_outcomes()
epi_list <- subset(epi_search, grepl("epilep", trait, ignore.case = TRUE))
print(epi_list[, c("id", "trait", "year", "sample_size", "ncase", "author")])

cat("实际返回：", nrow(mb_info), "\n")
print(mb_info$trait)   # 把完整名单发我，我帮你写准确的筛选正则
print(epi_list[, c("id", "trait", "year", "sample_size", "ncase", "author")], n = 27)
# 或导出来看
fwrite(epi_list[, c("id","trait","year","sample_size","ncase","author")],
       "results/MR/epilepsy_gwas_candidates.csv")
## ===================== 第③段：MR主分析（确认ID后运行） =====================
OUTCOME_ID <- "ieu-b-10"      # 主分析：focal epilepsy
## 跑完主分析后，把 OUTCOME_ID 分别换成 "ieu-b-8" 和 "ieu-b-14" 再各跑一遍，
## 作为敏感性验证和探索性分析（结果文件夹改名区分）
P_EXPOSURE <- 1e-5; CLUMP_R2 <- 0.001; CLUMP_KB <- 10000

save_mr_plots <- function(dat, res, tag) {
  ggsave(paste0("figures/MR/", tag, "_scatter.pdf"), mr_scatter_plot(res, dat)[[1]], width = 7, height = 6)
  ggsave(paste0("figures/MR/", tag, "_forest.pdf"),  mr_forest_plot(mr_singlesnp(dat))[[1]], width = 7, height = 6)
  ggsave(paste0("figures/MR/", tag, "_funnel.pdf"),  mr_funnel_plot(mr_singlesnp(dat))[[1]], width = 6, height = 6)
}

run_mr_pair <- function(exp_id, out_id, p_exp, tag) {
  message("\n========== ", tag, " ==========")
  tryCatch({
    exp_dat <- extract_instruments(outcomes = exp_id, p1 = p_exp,
                                   clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB)
    if (is.null(exp_dat) || nrow(exp_dat) < 3) { message("工具变量不足，跳过"); return(NULL) }
    out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id)
    if (is.null(out_dat) || nrow(out_dat) < 3) { message("结局SNP不足，跳过"); return(NULL) }
    dat <- harmonise_data(exp_dat, out_dat)
    res <- mr(dat)
    fwrite(res, paste0("results/MR/", tag, "_results.csv"))
    fwrite(mr_heterogeneity(dat),  paste0("results/MR/", tag, "_heterogeneity.csv"))
    fwrite(mr_pleiotropy_test(dat), paste0("results/MR/", tag, "_pleiotropy.csv"))
    save_mr_plots(dat, res, tag)
    print(res)
    list(dat = dat, res = res)
  }, error = function(e) { message("失败: ", e$message); NULL })
}

## 7. 正向MR：菌群 -> 癫痫
mr_forward <- list()
for (i in seq_len(nrow(target))) {
  tx <- target$trait[i]; tid <- target$id[i]
  tag <- gsub("[^A-Za-z0-9]+", "_", paste0("taxa2epi_", tx))
  mr_forward[[tx]] <- run_mr_pair(tid, OUTCOME_ID, P_EXPOSURE, tag)
}
ivw_sum <- bind_rows(lapply(names(mr_forward), function(tx) {
  r <- mr_forward[[tx]]; if (is.null(r)) return(NULL)
  r$res %>% filter(method == "Inverse variance weighted") %>% mutate(taxon = tx, .before = 1)
}))
print(ivw_sum); fwrite(ivw_sum, "results/MR/ALL_taxa2epi_IVW_summary.csv")
all_taxa <- mb_info  # 全部211个
mr_discovery <- list()
for (i in seq_len(nrow(all_taxa))) {
  tx <- all_taxa$trait[i]; tid <- all_taxa$id[i]
  tag <- gsub("[^A-Za-z0-9]+", "_", paste0("disc_", substr(tx, 1, 40)))
  mr_discovery[[tx]] <- run_mr_pair(tid, "ieu-b-10", 1e-5, tag)
}
disc_sum <- bind_rows(lapply(names(mr_discovery), function(tx) {
  r <- mr_discovery[[tx]]; if (is.null(r)) return(NULL)
  r$res %>% filter(method == "Inverse variance weighted") %>% mutate(taxon = tx, .before = 1)
}))
disc_sum$p_fdr <- p.adjust(disc_sum$pval, method = "BH")
fwrite(disc_sum, "results/MR/discovery_211taxa_IVW.csv")
print(disc_sum %>% arrange(pval) %>% head(20))
## 8. 反向MR：癫痫 -> 菌群
epi_exp <- extract_instruments(outcomes = OUTCOME_ID, p1 = 5e-8,
                               clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB)
mr_reverse <- list()
if (!is.null(epi_exp) && nrow(epi_exp) >= 3) {
  for (i in seq_len(nrow(target))) {
    tx <- target$trait[i]; tid <- target$id[i]
    tag <- gsub("[^A-Za-z0-9]+", "_", paste0("epi2taxa_", tx))
    mr_reverse[[tx]] <- tryCatch({
      out_dat <- extract_outcome_data(snps = epi_exp$SNP, outcomes = tid)
      if (is.null(out_dat) || nrow(out_dat) < 3) return(NULL)
      dat <- harmonise_data(epi_exp, out_dat)
      res <- mr(dat); fwrite(res, paste0("results/MR/", tag, "_results.csv")); res
    }, error = function(e) NULL)
  }
}

## 9. 代谢物 -> 癫痫（乙酸/BHBA/犬尿氨酸）
metabolites <- c("Acetate" = "ebi-a-GCST90092803",
                 "BHBA"    = "ebi-a-GCST90092811",
                 "Kynurenine" = "met-a-375")   # ID不符时按第6步检索修正
mr_metab <- list()
for (m in names(metabolites)) {
  mr_metab[[m]] <- run_mr_pair(metabolites[[m]], OUTCOME_ID, 5e-8, paste0("metab2epi_", m))
}

## 10. 两步中介MR：菌群 -> 代谢物 -> 癫痫
mediation <- list()
for (i in seq_len(nrow(target))) {
  tx <- target$trait[i]; tid <- target$id[i]
  for (m in names(metabolites)) {
    step1 <- run_mr_pair(tid, metabolites[[m]], P_EXPOSURE,
                         gsub("[^A-Za-z0-9]+", "_", paste0("taxa2metab_", tx, "_", m)))
    if (is.null(step1) || is.null(mr_metab[[m]]) || is.null(mr_forward[[tx]])) next
    get_b <- function(r) r$res %>% filter(method == "Inverse variance weighted") %>% pull(b)
    b1 <- get_b(step1); b2 <- get_b(mr_metab[[m]]); bt <- get_b(mr_forward[[tx]])
    if (length(c(b1, b2, bt)) < 3 || any(is.na(c(b1, b2, bt)))) next
    mediation[[paste(tx, m, sep = " -> ")]] <- data.frame(
      taxon = tx, metabolite = m, beta1 = b1, beta2 = b2, beta_total = bt,
      indirect = b1 * b2, prop_mediated = (b1 * b2) / bt)
  }
}
if (length(mediation) > 0) { med_df <- bind_rows(mediation); print(med_df)
fwrite(med_df, "results/MR/mediation_summary.csv") }

cat("\n全部完成！结果见 results/MR/，图片见 figures/MR/\n")
epi_exp <- extract_instruments(outcomes = "ieu-b-10", p1 = 5e-8,
                               clump = TRUE, r2 = 0.001, kb = 10000)
nrow(epi_exp)   # 大概率 <3
epi_exp8 <- extract_instruments(outcomes = "ieu-b-8", p1 = 5e-8,
                                clump = TRUE, r2 = 0.001, kb = 10000)
cat("ieu-b-8 工具变量数：", nrow(epi_exp8), "\n")
mr_reverse8 <- list()
if (!is.null(epi_exp8) && nrow(epi_exp8) >= 3) {
  for (i in seq_len(nrow(target))) {
    tx <- target$trait[i]; tid <- target$id[i]
    tag <- gsub("[^A-Za-z0-9]+", "_", paste0("epi8_2taxa_", tx))
    mr_reverse8[[tx]] <- tryCatch({
      out_dat <- extract_outcome_data(snps = epi_exp8$SNP, outcomes = tid)
      if (is.null(out_dat) || nrow(out_dat) < 3) return(NULL)
      dat <- harmonise_data(epi_exp8, out_dat)
      res <- mr(dat); fwrite(res, paste0("results/MR/", tag, "_results.csv")); res
    }, error = function(e) NULL)
  }
}
################################################################################
## MR结果可视化汇总图
## 输入：results/MR/ 下的汇总csv；输出：figures/MR_summary/
################################################################################
suppressPackageStartupMessages({library(data.table); library(dplyr); library(ggplot2)})
setwd(cfg$dir_bulk)
dir.create("figures/MR_summary", showWarnings = FALSE, recursive = TRUE)

shorten <- function(x) gsub("Gut microbiota abundance \\(|\\)", "", x)

## ============ 图1：发现型扫描火山图（203菌） ============
disc <- fread("results/MR/discovery_211taxa_IVW.csv")
disc <- disc %>%
  mutate(taxon_short = shorten(taxon),
         neglogp = -log10(pval),
         sig = case_when(pval < 0.05 & b > 0 ~ "Risk (P<0.05)",
                         pval < 0.05 & b < 0 ~ "Protective (P<0.05)",
                         TRUE ~ "NS"))
top8 <- disc %>% arrange(pval) %>% head(8)

p1 <- ggplot(disc, aes(b, neglogp, color = sig)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = c("Risk (P<0.05)" = "#E64B35",
                                "Protective (P<0.05)" = "#3C5488", "NS" = "grey80")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, color = "grey60") +
  ggrepel::geom_text_repel(data = top8, aes(label = taxon_short),
                           size = 2.8, max.overlaps = 15, show.legend = FALSE) +
  theme_bw() +
  labs(x = "Effect size (β, IVW) on focal epilepsy", y = "-log10(P)",
       color = NULL, title = "Discovery MR: gut microbiota → focal epilepsy")
ggsave("figures/MR_summary/Fig1_discovery_volcano.pdf", p1, width = 7.5, height = 6)
ggsave("figures/MR_summary/Fig1_discovery_volcano.png", p1, width = 7.5, height = 6, dpi = 300)

## ============ 图2：名义显著菌森林图 ============
top9 <- disc %>% filter(pval < 0.05) %>%
  mutate(or = exp(b), lo = exp(b - 1.96*se), hi = exp(b + 1.96*se),
         taxon_short = reorder(taxon_short, b))
p2 <- ggplot(top9, aes(or, taxon_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = b > 0), height = 0, linewidth = 1) +
  geom_point(aes(color = b > 0), size = 3) +
  geom_text(aes(x = 4.2, label = sprintf("P=%.3f", pval)), size = 3, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
  scale_x_log10(limits = c(0.5, 6), breaks = c(0.5, 1, 2, 4)) +
  theme_bw() +
  labs(x = "Odds ratio (95% CI, log scale)", y = NULL,
       title = "Nominally significant taxa (P<0.05)")
ggsave("figures/MR_summary/Fig2_nominal_forest.pdf", p2, width = 7, height = 5)
ggsave("figures/MR_summary/Fig2_nominal_forest.png", p2, width = 7, height = 5, dpi = 300)

## ============ 图3：目标菌 + 代谢物 综合效应图（先验层） ============
target_mr <- fread("results/MR/ALL_taxa2epi_IVW_summary.csv") %>%
  mutate(group = "Gut taxa", exposure_short = shorten(taxon))
metab_mr <- bind_rows(lapply(c("Acetate","BHBA","Kynurenine"), function(m) {
  f <- paste0("results/MR/metab2epi_", m, "_results.csv")
  if (!file.exists(f)) return(NULL)
  fread(f) %>% filter(method == "Inverse variance weighted") %>%
    mutate(group = "Metabolites", exposure_short = m)
}))
prior <- bind_rows(target_mr, metab_mr) %>%
  mutate(or = exp(b), lo = exp(b - 1.96*se), hi = exp(b + 1.96*se),
         exposure_short = reorder(exposure_short, b),
         sig_star = ifelse(pval < 0.05, "*", ""))
p3 <- ggplot(prior, aes(or, exposure_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = b > 0), height = 0, linewidth = 1) +
  geom_point(aes(color = b > 0), size = 3) +
  geom_text(aes(label = sig_star), vjust = -0.8, size = 6, color = "#E64B35") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"), guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
  theme_bw() +
  labs(x = "Odds ratio (95% CI, log scale)", y = NULL,
       title = "Prioritized exposures → focal epilepsy (IVW)",
       caption = "* P<0.05")
ggsave("figures/MR_summary/Fig3_prioritized_forest.pdf", p3, width = 7, height = 7)
ggsave("figures/MR_summary/Fig3_prioritized_forest.png", p3, width = 7, height = 7, dpi = 300)

## ============ 图4：犬尿氨酸重点验证（散点+leave-one-out） ============
## 注意：需要重新在线取数（token有效），或在TwoSampleMR对象还在时用现成函数
if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
  library(TwoSampleMR)
  exp_dat <- extract_instruments("met-a-375", p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
  out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = "ieu-b-10")
  dat_k  <- harmonise_data(exp_dat, out_dat)
  res_k  <- mr(dat_k)
  p4a <- mr_scatter_plot(res_k, dat_k)[[1]] +
    theme_bw() + ggtitle("Kynurenine → focal epilepsy (scatter)")
  ggsave("figures/MR_summary/Fig4a_kynurenine_scatter.pdf", p4a, width = 7, height = 6)
  ggsave("figures/MR_summary/Fig4a_kynurenine_scatter.png", p4a, width = 7, height = 6, dpi = 300)
  p4b <- mr_leaveoneout_plot(mr_leaveoneout(dat_k))[[1]] +
    theme_bw() + ggtitle("Leave-one-out sensitivity")
  ggsave("figures/MR_summary/Fig4b_kynurenine_loo.pdf", p4b, width = 7, height = 5)
  ggsave("figures/MR_summary/Fig4b_kynurenine_loo.png", p4b, width = 7, height = 5, dpi = 300)
}

cat("图片已输出到 figures/MR_summary/\n")
# 在 R 中运行
library(TwoSampleMR)

# 获取版本号——这就是文章中【待补】的内容
packageVersion("TwoSampleMR")
# 输出示例：[1] '0.5.6'
# 填入文章："TwoSampleMR（版本 0.5.6）"

packageVersion("ieugwasr")
# 输出示例：[1] '0.1.5'

library(TwoSampleMR)
library(data.table)


# 计算F统计量
calculate_F <- function(dat) {
  R2 <- sum(2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure) /
              (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure) + 
                 2 * (dat$se.exposure^2) * dat$samplesize.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure)))
  k <- nrow(dat)
  n <- mean(dat$samplesize.exposure)
  F <- (R2/(k-1)) / ((1-R2)/(n-k))
  return(round(F, 2))
}

# 取数
exp_k <- extract_instruments("met-a-375", p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
out_k <- extract_outcome_data(snps = exp_k$SNP, outcomes = "ieu-b-10")
dat_k <- harmonise_data(exp_k, out_k)
res_k <- mr(dat_k, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

# 先取IVW结果，避免长代码
ivw_res <- res_k[res_k$method == "Inverse variance weighted", ]
b_ivw <- ivw_res$b
se_ivw <- ivw_res$se
p_ivw <- ivw_res$pval
or_ivw <- exp(b_ivw)
or_lo <- exp(b_ivw - 1.96*se_ivw)
or_hi <- exp(b_ivw + 1.96*se_ivw)

# 统计量
F_k <- calculate_F(dat_k)
hetero_k <- mr_heterogeneity(dat_k)
p_hetero <- hetero_k$Q_pval[hetero_k$method == "Inverse variance weighted"]
pleio_k <- mr_pleiotropy_test(dat_k)
p_pleio <- pleio_k$pval
loo_k <- mr_leaveoneout(dat_k)
stable <- all(sign(loo_k$b) == sign(b_ivw))

# 汇总表
kynurenine_summary <- data.frame(
  exposure = "Kynurenine",
  outcome = "Focal epilepsy",
  n_SNP = nrow(dat_k),
  F_stat = F_k,
  OR = round(or_ivw, 2),
  OR_lo = round(or_lo, 2),
  OR_hi = round(or_hi, 2),
  P_IVW = signif(p_ivw, 3),
  P_heterogeneity = signif(p_hetero, 3),
  P_pleiotropy = signif(p_pleio, 3),
  leave_one_out_stable = ifelse(stable, "Yes", "No")
)
# 先创建文件夹
dir.create("results/MR", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/MR", showWarnings = FALSE, recursive = TRUE)

# 保存结果，用英文文件名避免路径错误
fwrite(kynurenine_summary, "results/MR/kynurenine_MR_summary.csv", row.names = FALSE)
p_loo <- mr_leaveoneout_plot(loo_k)[[1]] + theme_bw()
ggsave("figures/MR/kynurenine_leave_one_out.pdf", p_loo, width = 7, height = 6)
ggsave("figures/MR/kynurenine_leave_one_out.pdf", p_loo, width = 7, height = 6)
fwrite(kynurenine_summary, "results/MR/犬尿氨酸_MR完整结果.csv", row.names = FALSE)

# 留一法图
p_loo <- mr_leaveoneout_plot(loo_k)[[1]] + theme_bw()
ggsave("figures/MR/犬尿氨酸_留一法.pdf", p_loo, width = 7, height = 6)

# 打印结果
cat("=== 犬尿氨酸MR结果 ===\n")
cat("OR =", round(or_ivw,2), "95%CI:", round(or_lo,2), "-", round(or_hi,2), "P =", signif(p_ivw,3), "\n")
cat("F统计量:", F_k, "\n")
cat("异质性P:", signif(p_hetero,3), "\n")
cat("多效性P:", signif(p_pleio,3), "\n")
cat("留一法稳定:", stable, "\n")
library(ggplot2)

# ====================== 读取MR结果表（已参数化，见 config/config.R）======================
#（路径已由 cfg$dir_mr 统一管理）
mr_taxa <- read.csv(file.path(cfg$dir_mr, "犬尿氨酸_MR完整结果.csv"), 
                    stringsAsFactors = F, encoding = "UTF-8") 
# 如果你用中文文件名的那个结果，把下面这行文件名改成"犬尿氨酸_MR完整结果.csv"即可
mr_kyn <- read.csv(file.path(cfg$dir_mr, "犬尿氨酸_MR完整结果.csv"), 
                   stringsAsFactors = F, encoding = "UTF-8") 
# 直接用文件里现成的OR/CI/P值，不用任何计算
p_kyn_forest <- ggplot(mr_kyn, aes(x = OR, y = "Kynurenine")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.15, color = "black", linewidth = 0.8) +
  geom_point(shape = 18, size = 6, color = "#E64B35") +
  scale_x_log10(limits = c(0.8, 4), breaks = c(1, 2, 4), labels = c("1", "2", "4")) +
  geom_text(aes(x = OR_hi + 0.2, label = "*"), size = 7, color = "black", fontface = "bold") +
  geom_text(aes(x = 3.7, label = paste0("P=", round(P_IVW, 3))), size = 3.5, color = "black") +
  labs(x = "Odds ratio (95% CI, log scale)", y = "", title = "IVW: Kynurenine -> Focal epilepsy") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

# 保存
ggsave(file.path(cfg$dir_mr, "Fig4_Kynurenine_forest.pdf"), 
       plot = p_kyn_forest, width = 4, height = 1.8, dpi = 300, device = cairo_pdf)

cat("✅ 成功！OR =", mr_kyn$OR, "，95%CI：", mr_kyn$OR_lo, "-", mr_kyn$OR_hi, "，P =", mr_kyn$P_IVW, "\n")
cat("文件保存在：", file.path(cfg$dir_mr, "Fig4_Kynurenine_forest.pdf"), "\n")
