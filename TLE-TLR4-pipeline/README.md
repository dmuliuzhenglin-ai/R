# 耐药颞叶癫痫多组学分析流程（TLE-TLR4-pipeline）

**论文**：The TLR4-driven microglia–vascular endothelial inflammatory axis in drug-resistant temporal lobe epilepsy: multi-omics dissection and translational exploration
（TLR4驱动的小胶质细胞–血管内皮炎症轴：耐药颞叶癫痫的多组学机制解析与转化探索）

**作者**：刘政霖（Zhenglin Liu）、尹剑（Jian Yin，通讯作者）
**单位**：大连医科大学附属第二医院（The Second Affiliated Hospital of Dalian Medical University）
**投稿期刊**：IJMS
本仓库为论文全部分析代码的参数化公开版：所有路径集中于 `config/config.R`，**不含任何个人路径与密钥**。

---

## 1. 目录结构

```
TLE-TLR4-pipeline/
├── TLE-TLR4-pipeline.Rproj     # 双击打开（RStudio），here::here() 自动定位根目录
├── config/
│   └── config.R                # ★ 唯一需要按需修改的文件（数据目录定义）
├── R/
│   ├── 01_bulk_discovery_gutMGene.R
│   ├── 02_MR_microbiome_kynurenine.R
│   ├── 03_SMR_HEIDI_BrainMeta.R
│   ├── 04_singlecell_snRNAseq_full.R
│   ├── 05_ML_diagnostic_signature.R
│   ├── 06_LOOCV_internal_validation.R
│   ├── 07_external_validation_GSE256068.R
│   ├── 08_PBMC_validation_GSE143272.R
│   ├── 09_bulk_deconvolution_MCPcounter_quanTIseq.R
│   ├── 10_drug_repositioning_DSigDB.R
│   ├── 11_suppFigS3_hubgenes_DotPlot.R
│   ├── 12_suppFigS4_dualfilter_DotPlot.R
│   ├── 13_suppFigS5_lineage_markers_DotPlot.R
│   ├── 14_suppFigS11_statin_targets_DotPlot.R
│   ├── 15_GOKEGG_enrichment.R
│   ├── 16_EDIL3_ligand_receptor.R
│   └── 17_microglia_pseudotime_Slingshot.R
├── data/                       # 自行下载（见第 3 节；.gitignore 已屏蔽）
└── README.md
```

## 2. 运行环境

- Windows 10/11 + R ≥ 4.3（开发环境 R 4.4.3）+ RStudio
- 主要 R 包：DESeq2、WGCNA、STRINGdb、igraph、clusterProfiler、org.Hs.eg.db、msigdbr、
  Seurat、harmony、slingshot、miloR、CellChat、nichenetr、multinichenetr、
  TwoSampleMR、ieugwasr、glmnet、e1071、pROC、edgeR、limma、GEOquery、
  MCPcounter（`remotes::install_github("ebecht/MCPcounter", subdir = "Source")`）、
  quantiseqr、enrichR、illuminaHumanv4.db、data.table、tidyverse
- SMR 命令行工具（脚本 03 需外部 `smr` 可执行文件）

## 3. 数据下载清单（放入 `data/` 对应子目录）

| 数据 | 来源 | 放置目录 |
|---|---|---|
| GSE134697 CountMatrix（17对耐药TLE硬化海马/新皮层） | GEO 附件 | `data/bulk/data/GSE134697_CountMatrix.tsv/` |
| gutMGene 两个 CSV（human + causally 筛选） | gutMGene 数据库 | `data/bulk/` |
| GSE256068 series matrix + raw data（64 TLE-HS vs 13 对照） | GEO | `data/supplement/data/` |
| GSE143272 非标准化矩阵（141例外周血） | GEO | `data/pbmc/` |
| MiBioGen 菌群 GWAS | OpenGWAS API（脚本02在线获取） | — |
| ILAE 局灶性癫痫 GWAS（ieu-b-10 VCF） | OpenGWAS | `data/SMR/` |
| 血浆犬尿氨酸 GWAS（met-a-375，Shin 等） | OpenGWAS API | — |
| BrainMeta 脑组织顺式 eQTL（BESD 格式） | SMR 数据资源页 | `data/SMR/` |
| 1000 Genomes Phase 3 LD 参考（plink 格式） | SMR 数据资源页 | `data/SMR/` |
| snRNA-seq 矩阵与 Seurat_final_annotation_corrected.rds | 脚本04前半部分自行生成 | `data/singlecell/` |
| NicheNet ligand_target_matrix / lr_network | Zenodo（record/10229222 等） | `data/refdata/` |

## 4. OpenGWAS token 配置（脚本02需要）

1. 在 https://api.opengwas.io 注册并生成 token；
2. 写入本地 `~/.Renviron` 一行：`OPENGWAS_JWT="你的token"`；
3. 重启 R 生效。
**切勿**把 token 写入任何脚本或提交到 GitHub（`.Renviron` 已在 `.gitignore` 屏蔽）。

## 5. 运行顺序与论文对应

| 脚本 | 内容 | 论文章节/图表 |
|---|---|---|
| 01 | GSE134697 DESeq2 + gutMGene/MSigDB三重交集 + WGCNA + PPI/cytoHubba Top5枢纽 | 3.1节，图1，补充图S1/S2，补充表S2/S5 |
| 02 | 菌群MR扫描 + 犬尿氨酸→局灶性癫痫两样本MR | 3.6节，图4a/b，补充图S6，补充表S3 |
| 03 | SMR + HEIDI（BrainMeta 顺式eQTL，5个Hub） | 3.6节，补充表S3-C |
| 04 | snRNA-seq 全流程（整合/注释/双重过滤/拟时序/Milo/pseudobulk/CellChat/MultiNicheNet/LINCS扰动GSEA） | 3.2–3.5、3.8节，图2/3/6a，补充图S7 |
| 05 | LASSO + SVM-RFE → TGFB1+TLR4 签名 + 三队列ROC | 3.7节，图5c，补充表S4 |
| 06 | 签名 LOOCV 内部验证 | 3.7节，补充图S8 |
| 07 | GSE256068 外部验证 + 双队列DEG重叠 | 3.1/3.7节，图5a，补充图S9 |
| 08 | GSE143272 外周血验证 | 3.7节，图5b |
| 09 | MCP-counter + quanTIseq 双算法反卷积 | 3.5节末段 |
| 10 | DSigDB 双输入药物重定位 | 3.8节，补充图S10 |
| 11–14 | 补充图 S3/S4/S5/S11 重画脚本 | 对应补充图 |
| 15 | GO/KEGG 富集（5个细胞类型基因集，上/下调分开） | 3.4节 |
| 16 | EDIL3 配体-受体定位 + CellChat 组间差异互作 | 3.4节相关补充分析 |
| 17 | 小胶质 Slingshot 拟时序（单图版） | 3.3节，图2d |

**前置依赖**：脚本15、16 读取 `Seurat_final_annotation_v2.rds`（脚本04 小胶质重分群模块写回的对象）；
脚本17 读取 `True_Microglia_RPCA.rds`（脚本04 小胶质重分群产物）。均放入 `data/singlecell/`。

用法：RStudio 打开 `TLE-TLR4-pipeline.Rproj` → 打开脚本 → 点 **Source** 按钮（勿粘贴进 Console）。

## 6. 关键分析参数（与论文一致）

| 参数 | 取值 |
|---|---|
| bulk DEG 阈值 | \|log₂FC\| > 0.58 且 FDR < 0.05（DESeq2 配对模型） |
| WGCNA | power = 14，signed R² ≈ 0.80 |
| 单细胞 | Harmony 校正；前 40 个主成分；resolution = 0.6 |
| PPI | STRING combined score ≥ 700；cytoHubba MCC/Degree/BottleNeck/EPC 四算法综合 |
| MR | P < 5×10⁻⁸，r² < 0.001，kb = 10,000；IVW 主分析 |
| 诊断模型 | 训练阈值 0.48；L1 逻辑回归 |
| MultiNicheNet | logFC 阈值 0.50，P < 0.05 |
| 随机种子 | 1234（config.R） |

## 7. 重要声明

- **虚拟敲除（virtual KO）探索段已移除**：脚本04早期版本含一段 R 启发式"虚拟KO"模拟（非真实
  CellOracle 运行，部分响应率为手工设定常数），其定量结果无法复核，已从论文撤除（2026-08-08）。
  论文扰动分析以脚本04中 LINCS L1000 响应 panel 的 GSEA 为准。
- 本研究为公共数据集二次分析，所有数据来自已公开的去标识化数据集。
- 本研究未接受任何基金资助；作者声明无利益冲突。

## 8. 引用

如使用本代码，请引用上述论文（IJMS，投稿中）。
