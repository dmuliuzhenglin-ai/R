## ======================================================================
## 耐药TLE多组学分析流程 · 脚本09/17：bulk免疫浸润反卷积：MCP-counter + quanTIseq 双算法互证
## 对应论文：3.5节末段
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

## ============================================================
## 免疫浸润反卷积：GSE134697（TLE vs Control 海马）
## ============================================================
options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN"))

options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("quantiseqr", update = FALSE, ask = FALSE)
library(quantiseqr)
## MCP-counter 安装（一次性，公开来源）：
## install.packages("remotes")
## remotes::install_github("ebecht/MCPcounter", subdir = "Source")
library(MCPcounter)
## ============================================================
## 免疫浸润反卷积：GSE134697（TLE vs Control 海马）
## 算法：MCP-counter + quanTIseq（双算法互证）
## 前提：mat_train（symbol×34样本 vst 矩阵）、y_train 在会话中
## ============================================================
## quanTIseq：第一列是样本名，其余列是细胞
d_q <- as.data.frame(res_qseq)
colnames(d_q)[1] <- "sample"
long_q <- pivot_longer(d_q, -sample, names_to = "cell", values_to = "score")
long_q$method <- "quanTIseq"

all_res <- rbind(long_m[, c("sample","cell","score","method")],
                 long_q[, c("sample","cell","score","method")])
all_res$group <- factor(
  sub(".*CHKJ[.-]([CH]).*", "\\1", all_res$sample),
  levels = c("C", "H"), labels = c("Control", "TLE"))

print(head(all_res[all_res$method == "quanTIseq", ]))
print(table(all_res$method, all_res$group))
## 重跑统计 + 画图（之前那段 stats 循环 + p1/p2）
print(head(all_res[all_res$method == "quanTIseq", ]))
print(head(all_res[all_res$method == "MCP-counter", ]))
## ---- 统计 ----
combos <- unique(all_res[, c("method", "cell")])
res_list <- list()
for (i in seq_len(nrow(combos))) {
  x <- all_res[all_res$method == combos$method[i] & all_res$cell == combos$cell[i], ]
  w <- wilcox.test(score ~ group, data = x)
  res_list[[i]] <- data.frame(
    method = combos$method[i], cell = combos$cell[i],
    p = w$p.value,
    med_TLE  = median(x$score[x$group == "TLE"]),
    med_Ctrl = median(x$score[x$group == "Control"]))
}
stats <- do.call(rbind, res_list)
stats$FDR <- p.adjust(stats$p, "BH")
stats$direction <- ifelse(stats$med_TLE > stats$med_Ctrl, "TLE升高", "TLE降低")
sig <- stats[stats$p < 0.05, ]
print(stats[order(stats$p), ])
fwrite(stats, "results/Deconv/deconv_group_stats.csv")

## ---- 画图 ----
p1 <- ggplot(all_res[all_res$method == "MCP-counter", ],
             aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  facet_wrap(~ cell, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(Control = "#4DBBD5", TLE = "#E64B35")) +
  theme_classic(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "Abundance score",
       title = "MCP-counter deconvolution of TLE hippocampus (GSE134697)")
ggsave("results/Deconv/Fig_deconv_mcp.pdf", p1, width = 11, height = 5.5)
print(p1)

if (nrow(sig) > 0) {
  d_sig <- all_res[all_res$cell %in% unique(sig$cell), ]
  p2 <- ggplot(d_sig, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.5, width = 0.55) +
    facet_grid(cell ~ method, scales = "free_y") +
    scale_fill_manual(values = c(Control = "#4DBBD5", TLE = "#E64B35")) +
    theme_classic(base_size = 10) + theme(legend.position = "none") +
    labs(x = NULL, y = NULL, title = "Significant cell types (cross-method)")
  ggsave("results/Deconv/Fig_deconv_sig.pdf", p2, width = 7,
         height = max(2.5, 1 + 1.3 * length(unique(sig$cell))))
  print(p2)
}
cat("完成\n")
## 图保存到指定文件夹
out_dir <- file.path(cfg$dir_supp, "bulk 免疫浸润反卷积/结果")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 图1：MCP-counter 全细胞类型（主图）
p1 <- ggplot(all_res[all_res$method == "MCP-counter", ],
             aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  facet_wrap(~ cell, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(Control = "#4DBBD5", TLE = "#E64B35")) +
  theme_classic(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "Abundance score",
       title = "MCP-counter deconvolution of TLE hippocampus (GSE134697)")
ggsave(file.path(out_dir, "Fig_deconv_MCPcounter.pdf"), p1, width = 11, height = 5.5)

## 图2：显著细胞类型 × 双算法交叉验证
sig_cells <- unique(stats$cell[stats$p < 0.05])
d_sig <- all_res[all_res$cell %in% sig_cells, ]
p2 <- ggplot(d_sig, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.55) +
  facet_grid(cell ~ method, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4DBBD5", TLE = "#E64B35")) +
  theme_classic(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = NULL,
       title = "Significant cell types (cross-method validation)")
ggsave(file.path(out_dir, "Fig_deconv_significant_2methods.pdf"), p2,
       width = 7, height = max(2.5, 1 + 1.3 * length(sig_cells)))

## 统计表也存一份过去
fwrite(stats, file.path(out_dir, "deconv_group_stats.csv"))
cat("已保存到:", out_dir, "\n")
