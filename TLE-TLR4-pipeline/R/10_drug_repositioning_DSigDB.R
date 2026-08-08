## ======================================================================
## 耐药TLE多组学分析流程 · 脚本10/17：DSigDB 双输入药物重定位（enrichR）
## 对应论文：3.8节（补充图S10）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

## ============================================================
## 药物重定位：Enrichr DSigDB
## ============================================================
options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN"))
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)
library(data.table)

setwd(cfg$dir_bulk)
dir.create("results/Drug", showWarnings = FALSE)

## ---- 输入1：5 个 Hub ----
hubs <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2")

## ---- 输入2：上调 DEG Top150 ----
deg <- fread("results/Step1_DEG_all_GSE134697.csv")
deg_up <- deg[logFC > 0 & adj.P.Val < 0.05][order(-logFC)]
genes_up150 <- head(deg_up$Gene, 150)

## ---- DSigDB 查询 ----
dbs <- "DSigDB"
res_hub  <- enrichr(hubs, dbs)[[1]]
res_up   <- enrichr(genes_up150, dbs)[[1]]

print(head(res_hub[, c("Term", "P.value", "Adjusted.P.value", "Overlap", "Genes")], 15))
print(head(res_up[, c("Term", "P.value", "Adjusted.P.value", "Overlap", "Genes")], 15))

fwrite(res_hub, "results/Drug/DSigDB_5hubs.csv")
fwrite(res_up,  "results/Drug/DSigDB_upDEG150.csv")

## ---- 画 Top10 条形图 ----
library(ggplot2)
plot_top <- function(res, tag) {
  d <- head(res[order(res$Adjusted.P.value), ], 11)
  d$Term <- factor(d$Term, levels = rev(d$Term))
  ggplot(d, aes(-log10(Adjusted.P.value), Term)) +
    geom_col(fill = "#4DBBD5") +
    theme_classic(base_size = 11) +
    labs(x = "-log10(FDR)", y = NULL,
         title = paste0("DSigDB drug repurposing (", tag, ")"))
}
p1 <- plot_top(res_hub, "5 hub genes")
p2 <- plot_top(res_up, "top150 upregulated DEGs")
ggsave("results/Drug/Fig_DSigDB_hubs.pdf", p1, width = 7, height = 4)
ggsave("results/Drug/Fig_DSigDB_up150.pdf", p2, width = 7, height = 4)
print(p1); print(p2)
