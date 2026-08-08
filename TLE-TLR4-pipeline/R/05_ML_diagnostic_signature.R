## ======================================================================
## 耐药TLE多组学分析流程 · 脚本05/17：诊断签名：LASSO + SVM-RFE 交集 -> TGFB1+TLR4 逻辑回归 + 三队列ROC
## 对应论文：3.7节（图5c、补充表S4）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

## ============================================================
## 第 0 段：训练集表达矩阵（GSE134697, vst 标准化）
## ============================================================
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN"))

library(data.table)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(glmnet)
library(e1071)
library(pROC)
library(ggplot2)

setwd(cfg$dir_bulk)

## ---- 读 count matrix ----
counts <- fread("data/GSE134697_CountMatrix.tsv/GSE134697_CountMatrix.tsv",
                data.table = FALSE)
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
counts <- as.matrix(counts)
cat("counts 维度:", dim(counts), "\n")
print(colnames(counts)[1:6])
print(colnames(counts))


## 看后面 7 列是什么
print(colnames(counts)[35:41])

## ---- 分组（从路径中提取样本名）----
sn <- basename(colnames(counts))            # 去掉路径
sn <- sub("-Brain-RNA.*$", "", sn)          # CHKJ-C10 / CHKJ-H10 / ...
print(sn)

group <- ifelse(grepl("^CHKJ-C", sn), "Control",
                ifelse(grepl("^CHKJ-H", sn), "TLE", NA))
keep <- !is.na(group)
counts <- counts[, keep]
sn     <- sn[keep]
group  <- factor(group[keep], levels = c("Control", "TLE"))
patient <- factor(sub("^CHKJ-[CH]-?(\\d+)$", "\\1", sn))

print(table(group))
print(data.frame(sample = sn, patient, group))


## ---- 分组（配对设计）----
sn <- colnames(counts)
group <- ifelse(grepl("CHKJ-C", sn), "Control",
                ifelse(grepl("CHKJ-H", sn), "TLE", NA))
keep <- !is.na(group)
counts <- counts[, keep]
group  <- factor(group[keep], levels = c("Control", "TLE"))
patient <- factor(sub("^CHKJ-[CH]-?(\\d+)$", "\\1", colnames(counts)))
print(table(group))

## 检查数值范围和被剔除的列都是什么
cat("全矩阵范围:", range(counts, na.rm = TRUE), "\n")

all_sn <- basename(colnames(counts))
print(all_sn)   # 看现在还剩哪些列（keep 之后的）

## 如果之前有 5 列被 NA 规则剔掉了，看看它们叫什么
## （在原始 41 列里）
counts_full <- fread("data/GSE134697_CountMatrix.tsv/GSE134697_CountMatrix.tsv",
                     data.table = FALSE)
print(colnames(counts_full)[35:42])
print(head(counts_full[, 37:41], 3))

print(dim(counts))
print(length(group)); print(length(patient))
print(table(patient, group))
print(levels(patient))
## 直接从 counts 列名重建 sn / group / patient
sn <- sub("-Brain-RNA.*$", "", basename(colnames(counts)))
print(sn)

group   <- factor(ifelse(grepl("^CHKJ-C", sn), "Control", "TLE"),
                  levels = c("Control", "TLE"))
patient <- factor(sub("^CHKJ-[CH](\\d+)$", "\\1", sn))

print(levels(patient))            # 必须是 "1".."17"
print(table(patient, group))      # 每个 patient 必须 Control 1 / TLE 1

dds <- DESeqDataSetFromMatrix(counts, data.frame(group, patient),
                              ~ patient + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
vsd <- vst(dds, blind = FALSE)
mat_train <- assay(vsd)
cat("训练集矩阵:", dim(mat_train), "\n")   # 列数必须是 34
## ---- DESeq2 + vst ----
dds <- DESeqDataSetFromMatrix(counts, data.frame(group, patient),
                              ~ patient + group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
vsd <- vst(dds, blind = FALSE)
mat_train <- assay(vsd)
cat("训练集矩阵:", dim(mat_train), "\n")

## ---- Ensembl -> symbol（行名已是 symbol 则自动跳过）----
if (grepl("^ENSG", rownames(mat_train)[1])) {
  sym <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(mat_train)),
                "SYMBOL", "ENSEMBL", multiVals = "first")
  k <- !is.na(sym); mat_train <- mat_train[k, ]; sym <- sym[k]
  o <- order(rowMeans(mat_train), decreasing = TRUE)
  mat_train <- mat_train[o, ]; sym <- sym[o]
  mat_train <- mat_train[!duplicated(sym), ]
  rownames(mat_train) <- sym[!duplicated(sym)]
}
cat("注释后基因数:", nrow(mat_train), "\n")

## ---- 转置成 样本×基因 ----
x_train <- t(mat_train)
y_train <- group[colnames(mat_train)]
stopifnot(all(rownames(x_train) == names(y_train)))
cat("x_train:", dim(x_train), "；分组:", table(y_train), "\n")

dir.create("results/ML", showWarnings = FALSE)
## ---- 1. 候选基因 ----
hubs <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2")
hubs <- intersect(hubs, colnames(x_train))
X <- x_train[, hubs]
cat("候选基因:", paste(hubs, collapse = ", "), "\n")
## 按位置对齐（mat_train 列顺序 = counts 列顺序 = group 顺序）
stopifnot(ncol(mat_train) == length(group))
y_train <- group
print(table(y_train))   # 必须 Control 17 / TLE 17

x_train <- t(mat_train)

## 重跑第 1-4 段
hubs <- c("TGFB1", "ICAM1", "TLR4", "TLR2", "CCL2")
hubs <- intersect(hubs, colnames(x_train))
X <- x_train[, hubs]

set.seed(123)
cvfit <- cv.glmnet(X, y_train, family = "binomial", alpha = 1, nfolds = 10)
pdf("results/ML/Fig_ML_LASSO.pdf", width = 10, height = 4)
par(mfrow = c(1, 2))
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE)
plot(cvfit)
dev.off()
lasso_genes <- rownames(coef(cvfit, s = "lambda.min"))[
  which(as.vector(coef(cvfit, s = "lambda.min")) != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
cat("LASSO 选中:", paste(lasso_genes, collapse = ", "), "\n")

svmRFE <- function(X, y) {
  feats <- colnames(X); ranking <- c()
  while (length(feats) > 0) {
    fit <- svm(X[, feats, drop = FALSE], y, kernel = "linear",
               type = "C-classification", scale = TRUE)
    w <- t(fit$coefs) %*% fit$SV
    rk <- sort(w[1, ]^2, decreasing = TRUE)
    worst <- names(rk)[length(rk)]
    ranking <- c(worst, ranking)
    feats <- setdiff(feats, worst)
  }
  rev(ranking)
}


print(dim(X))        # 必须是 34 × 5
print(colnames(X))   # 必须是 TGFB1 ICAM1 TLR4 TLR2 CCL2
print(table(y_train))



## 输入 Q 退出 Browse 后，跑这个：

debuggingState(FALSE)

set.seed(123)
fit_svm <- svm(X, y_train, kernel = "linear",
               type = "C-classification", scale = TRUE)
w <- t(fit_svm$coefs) %*% fit_svm$SV
rank_svm <- names(sort(w[1, ]^2, decreasing = TRUE))
cat("SVM 重要性排序:", paste(rank_svm, collapse = " > "), "\n")

svm_genes <- rank_svm[seq_len(max(1, length(lasso_genes)))]
sig_genes <- intersect(lasso_genes, svm_genes)
if (length(sig_genes) == 0) sig_genes <- lasso_genes
cat("★ 最终诊断签名:", paste(sig_genes, collapse = " + "), "\n")

df_train <- data.frame(y = y_train, X[, sig_genes, drop = FALSE])
fit_lr <- glm(y ~ ., data = df_train, family = binomial)
print(summary(fit_lr))
saveRDS(list(model = fit_lr, genes = sig_genes), "results/ML/signature_model.rds")

## ============================================================
## ML 验证全套：签名重拟合 -> 三队列 ROC -> 校准曲线
## 前提：X, y_train, sig_genes 已在会话中
## ============================================================
library(glmnet); library(pROC); library(ggplot2)

## ---- 5. 用最终签名重拟合 + 训练集 ROC ----
set.seed(123)
cvfit_sig <- cv.glmnet(X[, sig_genes], y_train,
                       family = "binomial", alpha = 1, nfolds = 10)

pred_train <- as.vector(predict(cvfit_sig, newx = X[, sig_genes],
                                s = "lambda.min", type = "response"))
roc_train <- roc(y_train, pred_train, quiet = TRUE)
cat("训练集 GSE134697 AUC:", round(auc(roc_train), 3), "\n")

## ---- 5a. GSE256068 外部验证 ----
## expr/meta 不在会话里时自动重建
if (!exists("expr") || !exists("meta")) {
  library(GEOquery); library(org.Hs.eg.db); library(AnnotationDbi)
  gse <- getGEO(filename = file.path(cfg$dir_supp, "data/GSE256068_series_matrix.txt.gz"))
  ph  <- pData(gse)
  expr <- as.matrix(read.table(gzfile(file.path(cfg$dir_supp, "data/GSE256068_raw_data.csv.gz")),
                               header = TRUE, row.names = 1, check.names = FALSE))
  ph$sample_id <- gsub(" ", "", ph$description)
  meta <- data.frame(sample = ph$sample_id,
                     diag = sub("disease state: ", "", ph$characteristics_ch1.3),
                     row.names = ph$sample_id)
  meta <- meta[meta$diag %in% c("TLE-HS", "ControlHippocampus"), ]
  meta$group <- factor(ifelse(meta$diag == "TLE-HS", "TLE_HS", "Control"),
                       levels = c("Control", "TLE_HS"))
  common <- intersect(colnames(expr), meta$sample)
  meta <- meta[common, ]; expr <- expr[, common]
  sym <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(expr)),
                "SYMBOL", "ENSEMBL", multiVals = "first")
  k <- !is.na(sym); expr <- expr[k, ]; sym <- sym[k]
  o <- order(rowMeans(expr), decreasing = TRUE)
  expr <- expr[o, ]; sym <- sym[o]
  expr <- expr[!duplicated(sym), ]; rownames(expr) <- sym[!duplicated(sym)]
  cat("GSE256068 重建完成:", dim(expr), "\n")
}

X_ext <- t(expr[sig_genes, , drop = FALSE])
pred_ext <- as.vector(predict(cvfit_sig, newx = X_ext, s = "lambda.min",
                              type = "response"))
roc_ext <- roc(meta$group, pred_ext, quiet = TRUE)
cat("GSE256068 AUC:", round(auc(roc_ext), 3), "\n")
print(colnames(pb_sig))
## ===== PBMC 验证：从零重建（一个块跑完）=====
library(data.table); library(limma)
library(illuminaHumanv4.db); library(AnnotationDbi)
library(pROC)

pb <- fread(file.path(cfg$dir_pbmc, "GSE143272_non-normalized (1).txt.gz"),
            data.table = FALSE)

## 信号列
mat0 <- pb[, !grepl("Detection", colnames(pb))]
rownames(mat0) <- mat0$PROBE_ID
mat0$PROBE_ID <- NULL
## 在分组之前加这行转码
sn0 <- iconv(colnames(mat0), from = "latin1", to = "UTF-8")
print(sn0[grepl("Drug", sn0)][1:3])   # 确认能打出 Drug-naïve

grp0 <- ifelse(grepl("Control", sn0, ignore.case = TRUE), "Control",
               ifelse(grepl("Drug-na", sn0, ignore.case = TRUE), "DrugNaive",
                      ifelse(grepl("Non-Responder|NonResponder", sn0, ignore.case = TRUE), "NonResponder",
                             ifelse(grepl("Responder", sn0, ignore.case = TRUE), "Responder", NA))))
print(table(grp0, useNA = "always"))

k <- grp0 %in% c("Control", "DrugNaive") & !grepl("Replicate", sn0, ignore.case = TRUE)
mat1 <- mat0[, k]
grp1 <- factor(grp0[k], levels = c("Control", "DrugNaive"))
print(table(grp1))     ## 必须 Control 50 / DrugNaive 34

## 注释 + 标准化
sym1 <- mapIds(illuminaHumanv4.db, rownames(mat1),
               "SYMBOL", "PROBEID", multiVals = "first")
kk <- !is.na(sym1); mat1 <- mat1[kk, ]; sym1 <- sym1[kk]
expr_pb <- normalizeBetweenArrays(log2(as.matrix(mat1) + 1), "quantile")
o <- order(rowMeans(expr_pb), decreasing = TRUE)
expr_pb <- expr_pb[o, ]; sym1 <- sym1[o]
expr_pb <- expr_pb[!duplicated(sym1), ]
rownames(expr_pb) <- sym1[!duplicated(sym1)]
cat("expr_pb:", dim(expr_pb), "；TLR4 在:", "TLR4" %in% rownames(expr_pb), "\n")

## TLR4 单基因 ROC
roc_pb <- roc(grp1, as.numeric(expr_pb["TLR4", ]), quiet = TRUE)
cat("PBMC TLR4 AUC:", round(auc(roc_pb), 3), "\n")


## ---- 6. 三队列 ROC 合图 ----
setwd(cfg$dir_bulk)
pdf("results/ML/Fig_ML_ROC_3cohorts.pdf", width = 5.4, height = 5)
plot(roc_train, col = "#E64B35", legacy.axes = TRUE,
     main = "TGFB1 + TLR4 diagnostic signature")
plot(roc_ext, col = "#4DBBD5", add = TRUE)
plot(roc_pb,  col = "#00A087", add = TRUE)
legend("bottomright",
       legend = c(paste0("GSE134697 hippocampus, 2-gene (train)  AUC = ", round(auc(roc_train), 3)),
                  paste0("GSE256068 hippocampus, 2-gene (val)    AUC = ", round(auc(roc_ext), 3)),
                  paste0("GSE143272 PBMC, TLR4 (val)             AUC = ", round(auc(roc_pb), 3))),
       col = c("#E64B35", "#4DBBD5", "#00A087"), lwd = 2, cex = 0.75)
dev.off()

## ---- 7. 校准曲线（训练集）----
library(ggplot2)
cal_df <- data.frame(pred = pred_train,
                     obs = ifelse(y_train == "TLE", 1, 0))
cal_df$bin <- cut(cal_df$pred, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
cal_plot <- aggregate(obs ~ bin, cal_df, mean)
cal_plot$pred <- aggregate(pred ~ bin, cal_df, mean)$pred

p_cal <- ggplot(cal_plot, aes(pred, obs)) +
  geom_abline(linetype = 2, color = "grey50") +
  geom_line(color = "#E64B35") +
  geom_point(size = 3, color = "#E64B35") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic(base_size = 12) +
  labs(x = "Predicted probability", y = "Observed proportion",
       title = "Calibration curve (GSE134697)")
ggsave("results/ML/Fig_ML_calibration.pdf", p_cal, width = 4.5, height = 4)
print(p_cal)

## ---- 8. 模型系数（论文写评分公式用）----
print(coef(cvfit_sig, s = "lambda.min"))
saveRDS(list(model = cvfit_sig, genes = sig_genes),
        "results/ML/signature_model_glmnet.rds")

## ---- 9. 混淆矩阵（训练集最佳阈值下的敏感度/特异度）----
best <- coords(roc_train, "best", ret = c("threshold", "sensitivity", "specificity"))
print(best)
