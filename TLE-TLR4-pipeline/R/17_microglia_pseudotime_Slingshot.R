## ======================================================================
## 耐药TLE多组学分析流程 · 脚本17/17：小胶质 Slingshot 拟时序（Homeostatic→DAM 激活轨迹，单图版）
## 对应论文：3.3节（图2d）
## 【参数化版本】所有路径来自 config/config.R（勿在脚本内写绝对路径）
## 运行：RStudio 打开 TLE-TLR4-pipeline.Rproj -> 打开本脚本 -> Source
## ======================================================================
source(here::here("config", "config.R"))

########################################################################
# 图12 小胶质细胞拟时序分析（Slingshot 版）
# 与 Monocle3 版目标一致：验证 Homeostatic → DAM 的激活轨迹
# 优势：直接用 Seurat 已有的 UMAP + microglia_state，无需重新降维聚类
# 安装（Bioconductor，西湖大学镜像可装）：
#   options(repos = BiocManager::repositories())
#   options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
#   BiocManager::install("slingshot", update = FALSE, ask = FALSE)
# 前置：True_Microglia_RPCA.rds（含 microglia_state、umap）
# 输出文件夹：图12_小胶质拟时序分析
########################################################################
# 第1步：装包（Bioconductor，这次是真的在仓库里）
options(repos = BiocManager::repositories())
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("slingshot", update = FALSE, ask = FALSE)

# 第2步：验证后直接跑脚本
library(slingshot)
library(Seurat)
library(slingshot)
library(dplyr)
library(ggplot2)

work_dir <- cfg$dir_sc
setwd(work_dir)
dir.create("图12_小胶质拟时序分析", showWarnings = FALSE)

true_microglia <- readRDS("True_Microglia_RPCA.rds")
# 排除 Other（污染/技术簇不参与激活轨迹）
micro <- subset(true_microglia, subset = microglia_state != "Other")
micro$microglia_state <- droplevels(micro$microglia_state)
table(micro$microglia_state, micro$condition)


# ==================== 模块1 跑 Slingshot ====================

rd <- Embeddings(micro, "umap")                    # 直接用重分群的UMAP
cl <- as.character(micro$microglia_state)

# start.clus 指定轨迹起点为 Homeostatic（生物学依据：稳态是激活起点）
set.seed(42)
sds <- slingshot(rd, clusterLabels = cl,
                 start.clus = "Homeostatic",
                 reducedDim = "UMAP")

saveRDS(sds, "图12_小胶质拟时序分析/sds_小胶质.rds")

# 查看有几条谱系（可能多条；我们主要用第1条主谱系）
print(slingLineages(sds))
pt_mat <- slingPseudotime(sds)                     # 每细胞×每谱系的拟时序
pt <- pt_mat[, 1]                                  # 主谱系拟时序（NA=不在该谱系上）


# ==================== 模块2 轨迹可视化 ====================

umap_df <- data.frame(UMAP_1 = rd[, 1], UMAP_2 = rd[, 2],
                      state = cl,
                      condition = micro$condition,
                      pseudotime = pt)

# 提取轨迹曲线坐标（供 ggplot 叠加）
curves <- slingCurves(sds)
curve_list <- lapply(seq_along(curves), function(i) {
  cr <- curves[[i]]
  data.frame(UMAP_1 = cr$s[cr$ord, 1],
             UMAP_2 = cr$s[cr$ord, 2],
             lineage = paste0("L", i))
})
curve_df <- bind_rows(curve_list)

state_cols <- c("Homeostatic" = "#4393C3", "Activated" = "#F4A582",
                "DAM" = "#B2182B", "DAM_Prolif" = "#7B3294")

# 2.1 按状态着色 + 轨迹线（核心图）
p1 <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = state), size = 0.6, alpha = 0.7) +
  geom_path(data = curve_df, aes(group = lineage),
            color = "black", linewidth = 1.2, arrow = arrow(length = unit(0.3, "cm"))) +
  scale_color_manual(values = state_cols) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Microglial trajectory (Slingshot): states") +
  theme_bw(base_size = 13)
ggsave("图12_小胶质拟时序分析/图12_1_轨迹_按状态着色.pdf", p1, width = 8, height = 6)

# 2.2 按分组着色 + 轨迹线
p2 <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = condition), size = 0.6, alpha = 0.7) +
  geom_path(data = curve_df, aes(group = lineage), color = "black", linewidth = 1.2) +
  scale_color_manual(values = c("Control" = "#2166AC", "TLE" = "#B2182B")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Microglial trajectory: condition") +
  theme_bw(base_size = 13)
ggsave("图12_小胶质拟时序分析/图12_2_轨迹_按分组着色.pdf", p2, width = 8, height = 6)

# 2.3 按拟时序着色 + 轨迹线
p3 <- ggplot(umap_df %>% filter(!is.na(pseudotime)), aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = pseudotime), size = 0.6, alpha = 0.7) +
  geom_path(data = curve_df, aes(group = lineage), color = "black", linewidth = 1.2) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Microglial trajectory: pseudotime") +
  theme_bw(base_size = 13)
ggsave("图12_小胶质拟时序分析/图12_3_轨迹_按拟时序着色.pdf", p3, width = 8, height = 6)


# ==================== 模块3 拟时序验证与统计 ====================

pt_df <- umap_df %>% filter(!is.na(pseudotime))

# 3.1 拟时序 × 状态箱线图（验证 Homeostatic < DAM 递增）
p4 <- ggplot(pt_df, aes(x = state, y = pseudotime, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.3) +
  scale_fill_manual(values = state_cols) +
  labs(title = "Pseudotime by microglial state", y = "Pseudotime", x = NULL) +
  theme_bw(base_size = 13)
ggsave("图12_小胶质拟时序分析/图12_4_拟时序_按状态箱线图.pdf", p4, width = 7, height = 5)

# 3.2 拟时序密度 × 分组（预期 TLE 峰偏右）
p5 <- ggplot(pt_df, aes(x = pseudotime, fill = condition)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Control" = "#2166AC", "TLE" = "#B2182B")) +
  labs(title = "Pseudotime density: Control vs TLE",
       x = "Pseudotime (Homeostatic -> DAM)", y = "Density") +
  theme_bw(base_size = 13)
ggsave("图12_小胶质拟时序分析/图12_5_拟时序密度_按分组.pdf", p5, width = 7, height = 5)

# 3.3 统计：两组拟时序中位数比较
pt_stat <- wilcox.test(pseudotime ~ condition, data = pt_df)
med_ctl <- median(pt_df$pseudotime[pt_df$condition == "Control"])
med_tle <- median(pt_df$pseudotime[pt_df$condition == "TLE"])
cat("拟时序中位数: Control", round(med_ctl, 3), "vs TLE", round(med_tle, 3),
    ", Wilcoxon p =", format.pval(pt_stat$p.value), "\n")
sink("图12_小胶质拟时序分析/拟时序_组间统计.txt")
cat("median Control:", med_ctl, "\nmedian TLE:", med_tle, "\n")
print(pt_stat)
sink()

# 3.4 各状态拟时序中位数（验证轨迹方向）
print(pt_df %>% group_by(state) %>%
        summarise(median_pseudotime = median(pseudotime), n = n()))


# ==================== 模块4 DAM marker 沿拟时序表达 ====================

dam_markers <- intersect(c("TREM2", "APOE", "SPP1", "TYROBP", "P2RY12", "SALL1"),
                         rownames(micro))
expr_mat <- GetAssayData(micro, assay = "RNA", layer = "data")[dam_markers, , drop = FALSE]

expr_df <- data.frame(pseudotime = pt,
                      t(as.matrix(expr_mat)),
                      check.names = FALSE) %>%
  filter(!is.na(pseudotime)) %>%
  tidyr::pivot_longer(-pseudotime, names_to = "gene", values_to = "expr")

p6 <- ggplot(expr_df, aes(x = pseudotime, y = expr)) +
  geom_smooth(se = TRUE, color = "#B2182B", linewidth = 1, method = "loess", span = 0.6) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(title = "Marker expression along pseudotime",
       x = "Pseudotime (Homeostatic -> DAM)", y = "Expression") +
  theme_bw(base_size = 12)
ggsave("图12_小胶质拟时序分析/图12_6_marker沿拟时序表达.pdf", p6, width = 10, height = 6)

########################################################################
# 结果解读备忘：
# 图12_4 → 若拟时序按 Homeostatic < (Activated) < DAM 递增 = 轨迹成立
# 图12_5 → 若 TLE 密度峰偏右 = TLE 细胞位于轨迹激活端
# 图12_6 → TREM2/APOE/SPP1 应上升，P2RY12/SALL1 应下降
# 写作句式（数字出来后填入）:
#   "Pseudotime analysis (Slingshot) reconstructed a Homeostatic-to-DAM
#   trajectory along which TLE nuclei were skewed toward the activated end
#   (median pseudotime x.xx vs x.xx, Wilcoxon p = x.x)"
# 若有多条谱系：报告主谱系(L1)结果即可；全部谱系可用 slingPseudotime 各行
########################################################################
