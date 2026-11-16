# 单细胞分析报告 · demo-pbmc

- **数据集**：`ds1`（demo-pbmc）
- **来源**：内置 demo 生成器（test/make_sample.py）
- **规模**：2208 个细胞 × 1200 个基因
- **引擎**：builtin（Python 3.14.2；dsh-singlecell v0.1.0）
- **生成时间**：2026-11-06T07:08:06.514Z

## 数据摘要

| 项目 | 值 |
| --- | --- |
| 细胞数 | 2208 |
| 基因数 | 1200 |
| 已完成步骤 | qc → filter → normalize → hvg → pca → neighbors → cluster → embed → markers |
| 过滤保留 | 2208 / 2400 细胞，1200 / 1200 基因 |
| 聚类 | 10 个簇 |

## Methods

Single-cell RNA-seq data were analysed with **dsh-singlecell** v0.1.0 (builtin NumPy/SciPy/scikit-learn implementation; Python 3.14.2; NumPy 2.5.1; SciPy 1.17.1; scikit-learn 1.9.0; AnnData 0.12.19; matplotlib 3.11.0).

Cells were retained when they had at least 200 detected genes per cell, mitochondrial counts below 20% of total counts. Genes detected in fewer than 3 cells were discarded. This retained 2208 of 2400 cells and 1200 of 1200 genes.

Counts were scaled to a target library size of 10000 per cell and log1p-transformed (log(1 + x)).

Highly variable genes were selected using the Seurat-style mean-binned normalised dispersion method (20 equal-frequency mean bins, top 2000 genes; 1200 genes selected).

Principal component analysis was performed on the scaled highly variable genes with 30 components (randomized SVD, random seed 0); PC1 explained 1.2% of the variance.

A k-nearest-neighbour graph (k = 15, distance metric = cosine) was built in PCA space and symmetrised.

Cells were clustered with Louvain community detection (resolution = 1) on the kNN graph; communities were numbered by decreasing size (10 clusters).

Two-dimensional visualisation was produced with t-SNE (perplexity = 30, PCA initialisation, random seed 0).

Cluster marker genes were identified with one-vs-rest Welch t-tests on log-normalised expression (top 10 genes per cluster); p-values were adjusted with the Benjamini-Hochberg procedure.

All random seeds were fixed at 0, and all figures were generated with matplotlib.

## 图表

![聚类嵌入图](fig-ds1-embedding-2084529dc992.png)

![QC 指标分布](fig-ds1-qc_violin-a8dbedd4ee05.png)

## 各簇 top markers

| cluster | gene | logFC | adj.P | pct_in |
| --- | --- | --- | --- | --- |
| 0 | CD6 | 5.104 | 0 | 1 |
| 0 | TRBC1 | 5.067 | 0 | 0.9918 |
| 0 | TRBC2 | 4.69 | 0 | 0.9891 |
| 0 | LTB | 3.362 | 0 | 0.9946 |
| 0 | LCK | 3.311 | 0 | 0.9891 |
| 1 | KLRF1 | 4.802 | 0 | 1 |
| 1 | GZMB | 4.798 | 0 | 1 |
| 1 | GZMA | 4.226 | 0 | 0.9973 |
| 1 | SH2D1B | 4.193 | 0 | 0.9891 |
| 1 | NKG7 | 3.655 | 0 | 0.9783 |
| 2 | CLU | 3.65 | 0 | 1 |
| 2 | ITGB3 | 3.528 | 0 | 1 |
| 2 | TUBB1 | 2.956 | 0 | 0.9966 |
| 2 | TREML1 | 3.039 | 0 | 0.9966 |
| 2 | NRGN | 2.356 | 0 | 1 |
| 3 | CD14 | 3.621 | 0 | 1 |
| 3 | S100A9 | 3.31 | 0 | 0.9964 |
| 3 | CTSD | 3.273 | 0 | 0.9928 |
| 3 | S100A8 | 2.892 | 0 | 1 |
| 3 | FCN1 | 3.31 | 0 | 0.9892 |
| 4 | TNFRSF13C | 2.894 | 0 | 1 |
| 4 | VPREB3 | 2.91 | 0 | 0.9951 |
| 4 | CD79A | 3.028 | 0 | 0.9901 |
| 4 | FCRLA | 2.536 | 0 | 0.9951 |
| 4 | CD37 | 2.154 | 0 | 1 |
| 5 | CLEC10A | 2.965 | 0 | 1 |
| 5 | CD1C | 2.927 | 0 | 1 |
| 5 | BATF3 | 2.807 | 0 | 1 |
| 5 | ITGAX | 2.442 | 0 | 1 |
| 5 | HLA-DPA1 | 2.893 | 0 | 0.9948 |
| 6 | CLEC10A | 2.845 | 0 | 1 |
| 6 | BATF3 | 2.651 | 0 | 1 |
| 6 | CCR7 | 2.631 | 0 | 1 |
| 6 | FLT3 | 2.927 | 0 | 0.9943 |
| 6 | HLA-DPA1 | 2.71 | 0 | 1 |
| 7 | VPREB3 | 2.732 | 0 | 1 |
| 7 | CD79A | 2.832 | 0 | 1 |
| 7 | TNFRSF13C | 2.678 | 0 | 1 |
| 7 | FCRLA | 2.462 | 0 | 1 |
| 7 | CD37 | 2.068 | 0 | 1 |
| 8 | CD14 | 2.674 | 0 | 1 |
| 8 | S100A9 | 2.294 | 0 | 1 |
| 8 | CTSD | 2.44 | 0 | 1 |
| 8 | FCN1 | 2.376 | 0 | 1 |
| 8 | CSF1R | 1.929 | 0 | 1 |
| 9 | ITGB3 | 2.484 | 0 | 1 |
| 9 | TUBB1 | 2.185 | 0 | 1 |
| 9 | GP1BA | 1.873 | 0 | 1 |
| 9 | CLU | 2.458 | 0 | 0.9865 |
| 9 | TREML1 | 2.1 | 0 | 1 |

## 自动注释建议

参考面板：`builtin-marker-panel-v1`（仅建议，未写入 obs）

| cluster | 建议类型 | score | margin | evidence |
| --- | --- | --- | --- | --- |
| 0 | T cell | 0.2397 | 0.181 | ZAP70, TRAC, CD3E, CD3D, TRBC1 |
| 1 | NK cell | 0.2115 | 0.1212 | NCAM1, KLRF1, PRF1, CTSW, TYROBP |
| 2 | Platelet | 0.2219 | 0.1622 | SELP, NRGN, CLU, GP1BA, PPBP |
| 3 | Monocyte | 0.2295 | 0.1582 | CSF1R, CST3, FCGR3A, S100A9, CTSD |
| 4 | B cell | 0.2111 | 0.1515 | CD37, CD19, MS4A1, TNFRSF13C, CD74 |
| 5 | Dendritic cell | 0.2294 | 0.1652 | FCER1A, CST3, CD74, CLEC10A, CD1E |
| 6 | Dendritic cell | 0.2077 | 0.1437 | FCER1A, CST3, CLEC10A, CD74, CD1E |
| 7 | B cell | 0.1944 | 0.1338 | CD37, CD19, TNFRSF13C, CD74, MS4A1 |
| 8 | Monocyte | 0.2083 | 0.1323 | CSF1R, CST3, FCGR3A, S100A9, CD14 |
| 9 | Platelet | 0.2004 | 0.1381 | SELP, NRGN, CLU, GP1BA, GP9 |

## 自检 Audit

**score = 90 / 100（ok）**

- [OK] **先质控再聚类**：已完成 qc/filter，聚类在过滤后的数据上运行
- [OK] **线粒体比例**：过滤后 pct_counts_mt 中位数 1.40%
- [OK] **细胞数**：过滤后 2208 个细胞 × 1200 个基因
- [OK] **簇数**：10 个簇
- [OK] **簇大小分布**：最大簇占比 17%，分布均衡
- [OK] **标记基因**：10 个簇都有显著标记基因
- [OK] **二维嵌入**：已生成 TSNE 二维嵌入
- [WARN] **细胞类型注释**：尚未写回 cell_type（仍是簇号）；建议：先运行 autoannotate 看参考面板建议，再用 annotate 写回
- [OK] **双细胞风险**：n_genes 极高（> Q3+1.5·IQR）的细胞仅占 0.50%
- [OK] **分辨率稳健性**：已扫描分辨率：0.4, 0.8, 1.2

## 分析历史

| 时间 | 步骤 | 参数 | 结果摘要 | 耗时 |
| --- | --- | --- | --- | --- |
| 2026-11-06T07:08:02.988Z | qc | `{"min_genes": 200, "max_genes": null, "min_counts": null, "max_pct_mt": 20.0}` | 2400 个细胞 · 1200 个基因；每细胞基因数中位数 570，总计数中位数 997.5 | 2 ms |
| 2026-11-06T07:08:02.999Z | filter | `{"min_genes": 200, "max_genes": null, "min_counts": null, "max_pct_mt": 20.0,...` | 保留 2208/2400 个细胞（过滤掉 192 个）；保留 1200/1200 个基因（0 个基因在 < 3 个细胞中检出） | 8 ms |
| 2026-11-06T07:08:03.012Z | normalize | `{"target_sum": 10000.0}` | 每细胞缩放到目标总计数 10000，再做 log1p；标准化矩阵非零元素 1285371 个 | 11 ms |
| 2026-11-06T07:08:03.019Z | hvg | `{"n_top_genes": 2000, "flavor": "seurat"}` | 筛选出 1200 个高变基因（n_top_genes=2000，flavor=seurat）；可用基因只有 1200 个（少于 n_top_genes=2000），已取全部 | 5 ms |
| 2026-11-06T07:08:03.039Z | pca | `{"n_pcs": 30}` | 取 30 个主成分；PC1 解释方差 1.2% | 17 ms |
| 2026-11-06T07:08:03.085Z | neighbors | `{"n_neighbors": 15, "metric": "cosine"}` | k=15（metric=cosine）的 kNN 图；邻接图 22742 条边，4 个连通分量 | 44 ms |
| 2026-11-06T07:08:03.205Z | cluster | `{"resolution": 1.0, "method": "louvain"}` | resolution=1 → 10 个簇；簇大小（前 10）：0:368，1:368，2:294，3:279，4:203，5:194，6:174，7:165，8:89，9:74 | 117 ms |
| 2026-11-06T07:08:05.444Z | embed | `{"layout": "tsne", "perplexity": 30.0, "n_iter": 1000}` | TSNE 二维嵌入完成（2208 个细胞）；perplexity=30，max_iter=1000 | 2235 ms |
| 2026-11-06T07:08:05.473Z | markers | `{"groupby": "clusters", "top_n": 10, "method": "t-test"}` | 10 个分组 × top 10 = 100 个标记基因；显著正向标记（adj.P<0.05 且 logFC>0）：348 个 | 25 ms |
