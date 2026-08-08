## ======================================================================
## 耐药TLE多组学分析流程 · 脚本06/17：TGFB1+TLR4 签名 LOOCV 内部验证（34例留出预测汇总）
## 对应论文：3.7节（补充图S8、补充表S4）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

## ============================================================
## LOOCV 内部验证（独立运行版，不依赖其他脚本）
## 数据：GSE134697 训练集 17 对配对样本（34个样本）
## 模型：TGFB1 + TLR4 双基因签名 + 逻辑回归
## 依赖：edgeR, org.Hs.eg.db, AnnotationDbi, pROC
## ============================================================
suppressPackageStartupMessages({
  library(edgeR)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(pROC)
})

out_dir <- file.path(cfg$dir_supp, "验证")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. 读训练集 count 矩阵 ----
train_counts <- read.table(
  file.path(cfg$dir_bulk, "data/GSE134697_CountMatrix.tsv/GSE134697_CountMatrix.tsv"),
  header = TRUE, row.names = 1, sep = "\t",
  check.names = FALSE, skip = 1)          # 跳过featureCounts命令行注释行
train_counts <- as.matrix(train_counts[, -c(1:5)])   # 去掉Chr/Start/End/Strand/Length

## 样本名提取（与原脚本一致）
colnames(train_counts) <- sapply(colnames(train_counts), function(x) {
  if (grepl("CHKJ-", x)) {
    gsub(".*CHKJ-([CH][0-9]+)-Brain.*", "\\1", x)
  } else if (grepl("CTRL-", x)) {
    gsub(".*CTRL-0*([0-9]+)-Brain.*", "CTRL\\1", x)
  } else x
})

## 分组 + 剔除额外CTRL样本，保留17对
train_group <- data.frame(
  sample = colnames(train_counts),
  group  = ifelse(grepl("^H", colnames(train_counts)), "TLE", "Control"),
  row.names = colnames(train_counts))
keep_sample <- !grepl("^CTRL", colnames(train_counts))
train_counts <- train_counts[, keep_sample]
train_group  <- train_group[keep_sample, ]
cat("训练集分组：\n"); print(table(train_group$group))   # 应为 Control 17 / TLE 17

## ---- 2. 基因注释（Ensembl -> Symbol）----
ens <- sub("\\..*", "", rownames(train_counts))
sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
keep <- !is.na(sym)
train_counts <- train_counts[keep, ]; sym <- sym[keep]
o <- order(rowMeans(train_counts), decreasing = TRUE)
train_counts <- train_counts[o, ]; sym <- sym[o]
train_counts <- train_counts[!duplicated(sym), ]
rownames(train_counts) <- sym[!duplicated(sym)]
cat("注释后基因数:", nrow(train_counts), "\n")

## ---- 3. log2 CPM 标准化 + 提取双基因 ----
expr_tr <- edgeR::cpm(train_counts, log = TRUE, prior.count = 1)
sig_genes <- c("TGFB1", "TLR4")
stopifnot(all(sig_genes %in% rownames(expr_tr)))
Xall <- t(expr_tr[sig_genes, ])                       # 34 x 2
yall <- factor(train_group$group, levels = c("Control", "TLE"))

## ---- 4. LOOCV 主循环 ----
n <- nrow(Xall); loo_prob <- numeric(n)
set.seed(2024)
for (i in seq_len(n)) {
  mu  <- colMeans(Xall[-i, , drop = FALSE])
  sd_ <- apply(Xall[-i, , drop = FALSE], 2, sd)
  Xtr <- scale(Xall[-i, , drop = FALSE], center = mu, scale = sd_)
  Xte <- scale(Xall[i,  , drop = FALSE], center = mu, scale = sd_)
  fit_i <- glm(yall[-i] ~ ., data = data.frame(Xtr), family = binomial)
  loo_prob[i] <- as.numeric(predict(fit_i, newdata = data.frame(Xte), type = "response"))
}

## ---- 5. 效能：AUC + 95%CI + 阈值0.48下的敏感度/特异度 ----
roc_loo <- roc(yall, loo_prob, levels = c("Control", "TLE"), direction = "<", quiet = TRUE)
ci_loo  <- ci.auc(roc_loo)

pred_cls <- ifelse(loo_prob > 0.48, "TLE", "Control")
cm <- table(Prediction = factor(pred_cls, levels = c("Control", "TLE")),
            Reference  = factor(ifelse(yall == "TLE", "TLE", "Control"),
                                levels = c("Control", "TLE")))
sens <- cm["TLE","TLE"] / sum(cm[,"TLE"])
spec <- cm["Control","Control"] / sum(cm[,"Control"])
acc  <- sum(diag(cm)) / sum(cm)

cat("\n========== LOOCV 结果 ==========\n")
cat(sprintf("LOOCV AUC = %.3f (95%%CI: %.3f - %.3f)\n",
            as.numeric(auc(roc_loo)), ci_loo[1], ci_loo[3]))
cat(sprintf("阈值0.48下：敏感度 %.1f%%，特异度 %.1f%%，准确率 %.1f%%\n",
            sens*100, spec*100, acc*100))

## ---- 6. 保存结果 ----
write.csv(data.frame(sample = rownames(Xall), group = yall, loo_prob = loo_prob),
          file.path(out_dir, "LOOCV_逐样本预测概率.csv"), row.names = FALSE)
write.csv(data.frame(metric = c("AUC","CI_low","CI_high","sensitivity","specificity","accuracy"),
                     value  = c(round(as.numeric(auc(roc_loo)),3), round(ci_loo[1],3),
                                round(ci_loo[3],3), round(sens,3), round(spec,3), round(acc,3))),
          file.path(out_dir, "LOOCV_效能汇总.csv"), row.names = FALSE)

## ---- 7. ROC图（pROC自带画图，无需ggplot2）----
pdf(file.path(out_dir, "Fig_LOOCV_ROC_GSE134697.pdf"), width = 6, height = 6)
plot(roc_loo, col = "#3C5488", lwd = 2.5,
     main = "LOOCV internal validation (GSE134697, n = 34)")
abline(0, 1, lty = 2, col = "grey60")
legend("bottomright", bty = "n", lwd = 2.5, col = "#3C5488",
       legend = sprintf("LOOCV AUC = %.3f (95%%CI: %.3f-%.3f)",
                        as.numeric(auc(roc_loo)), ci_loo[1], ci_loo[3]))
dev.off()
png(file.path(out_dir, "Fig_LOOCV_ROC_GSE134697.png"), width = 1800, height = 1800, res = 300)
plot(roc_loo, col = "#3C5488", lwd = 2.5,
     main = "LOOCV internal validation (GSE134697, n = 34)")
abline(0, 1, lty = 2, col = "grey60")
legend("bottomright", bty = "n", lwd = 2.5, col = "#3C5488",
       legend = sprintf("LOOCV AUC = %.3f (95%%CI: %.3f-%.3f)",
                        as.numeric(auc(roc_loo)), ci_loo[1], ci_loo[3]))
dev.off()

cat("\n✅ 全部完成！结果文件在:", out_dir, "\n")
cat("   - LOOCV_效能汇总.csv（投稿用数字）\n")
cat("   - Fig_LOOCV_ROC_GSE134697.pdf/png（可入补充图）\n")

