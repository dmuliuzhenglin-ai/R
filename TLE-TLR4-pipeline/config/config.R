## ======================================================================
## 全局配置（整个流程唯一需要按需修改的文件）
## 使用方式：RStudio 打开 TLE-TLR4-pipeline.Rproj，here::here() 自动定位项目根目录
## ======================================================================
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

cfg <- new.env()

## ---- 数据目录（相对项目根目录；按 README.md 第 3 节下载数据后放入）----
cfg$dir_bulk <- here::here("data", "bulk")        # GSE134697 CountMatrix、gutMGene CSV、Step1 输出
cfg$dir_sc   <- here::here("data", "singlecell")  # snRNA-seq 矩阵、Seurat_final_annotation_corrected.rds
cfg$dir_supp <- here::here("data", "supplement")  # GSE256068 series matrix / raw data、验证输出
cfg$dir_pbmc <- here::here("data", "pbmc")        # GSE143272 非标准化矩阵
cfg$dir_smr  <- here::here("data", "SMR")         # ieu-b-10 VCF、BrainMeta BESD、1000G LD(plink)
cfg$dir_mr   <- here::here("data", "MR")          # MR 结果表与森林图
cfg$dir_refdata <- here::here("data", "refdata")  # NicheNet/MultiNicheNet 先验矩阵（Zenodo，约262MB）

## ---- OpenGWAS token（安全要求：勿写入任何代码文件！）----
## 1) 在 https://api.opengwas.io 生成 token
## 2) 写入本地 ~/.Renviron 一行：  OPENGWAS_JWT="你的token"
## 3) 重启 R；TwoSampleMR/ieugwasr 会自动读取
## .Renviron 已在 .gitignore 屏蔽

cfg$seed <- 1234
