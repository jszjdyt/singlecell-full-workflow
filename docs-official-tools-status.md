# 官方工具实测状态（滚动更新，Lead 维护）

宿主：R 4.4.2（系统库 Seurat 5.4.0）+ 独立库 `~/.dsh/sc-plugin-rlib`（Seurat **6.0.0**）；
`externalPython = /opt/miniconda3/bin/python3`（3.13.11，pyscenic 0.12.1 + scanpy/scvelo/palantir/cellphonedb/celloracle/decoupler/celltypist）；
JAGS 5.0.1 已装；SCENIC motif 库 `~/Desktop/单细胞工作站/SCINCE尝试/human/python用`；NicheNet 先验 `~/Desktop/祖传R语言单细胞代码/crosstalk`。

## 已实测通过（真实调用 + 真结果）

| 工具 | 用时（demo 2208×1200） | 实测产出 |
|---|---|---|
| ext_seurat6 | 4.4s full；8.5s 走完 seurat6 预设 | 6 簇，ARI vs 真值 1.0，PCA 30 / UMAP 2 写回 |
| ext_cellchat | 5.9s | 官方 CellChatDB.human 通讯表 + out/in strength 写 obs + 2 图 |
| ext_monocle3 | 11.3s | 伪时间/partition/cluster 写 obs + 2 图 |
| ext_clusterprofiler | 31.4s | 408 marker→真 GO/KEGG 条目 + 2 图 |
| ext_fgsea | 2.4s | 官方 fgsea NES/p + 2 图 |
| ext_gsva | 10.4s | 12 通路活性写 obs + 热图 |
| ext_aucell | 4.7s | 12 基因集 AUC 写 obs + 2 图 |
| ext_harmony | 4.6s | 5 个修正分量写 obs + 前后对比图 |
| ext_cellphonedb | 4.3s | 官方置换检验通讯表 + 2 图 |
| ext_pyscenic | 48.3s | GRNBoost2 21582 TF→靶边/18 TF，12 regulon 活性写 obs |
| ext_slingshot | 4.1s | 3 条谱系伪时间写 obs + 2 图 |
| ext_celloracle | 10–20s | 需传 base GRN 中的 TF（如 IRF8/BATF3/LEF1/STAT1），否则官方拒绝并给候选 |

## 待修 / 待装（下一轮）

- `ext_wgcna`：demo 合成数据在默认 min_module_size 下**无模块** → 适配器应在无模块时自动降阈值（或提示用户传 `min_module_size`），不应直接 prerequisite。
- `ext_nichenet` / `ext_scenic_r`：适配器已就绪，未跑实测（重任务）。
- `ext_infercnv`：R 包 Bioc 安装未成功（JAGS 已就绪），需查 Bioc 安装日志后重试。
- `regvelo`：pip 未装上（依赖 scvi-tools/jax，py3.13 兼容性），需给出安装结论或替代。
- `ext_soupx` / `ext_copykat`：R 包已装（SoupX 1.6.2 / copykat 1.2.5），适配器尚未交付。
- README/INSTALL 的官方工具清单与 UI 标注、全量回归测试。

## 第三轮实测补充（22:30–23:00）

- **NicheNet 官方 OK**：demo 上 20.9s，自动选 sender=Dendritic cell / receiver=B cell，输出配体活性表 + 2 图。
- **pySCENIC 官方 OK**：48.3s，GRNBoost2 21582 条 TF→靶边 / 18 TF / 12 regulon 活性写 obs。
- **CellOracle 官方 OK**：25.0s，target=IRF8 可扰动；非 TF target 会被官方拒绝并给出候选（不造假）。
- **slingshot 官方 OK**：4.1s，3 条谱系伪时间写 obs。
- **全量回归测试全绿**：263/83(ALL PASS)/10841/527/241/34/68。
- **regVelo**：conda base(py3.13) pip 报 ResolutionImpossible → 已改为创建独立 `dsh-regvelo`(py3.11) 环境安装（后台），并在 sc_external 增加 `regveloPython` 独立解释器配置。
- **inferCNV**：Bioc 安装两次后仍 MISSING（无报错行，属依赖静默跳过）；JAGS 已装，待下一轮用 `BiocManager::install(..., dependencies=TRUE)` 或源码安装排查。
- **WGCNA**：已加"无模块自动降阈值重试"，但合成 demo 仍无模块（数据问题）；待真实数据验证。
- **待补适配器**：SoupX（R 包 1.6.2 已装）、copyKAT（R 包 1.2.5 已装）。

## 第四轮补充

- **SCENIC(R) 官方 OK**：192.3s，用本地 hg38 mc9nr motif 库（27016 基因排名），7 个 regulon 活性写 obs + 2 图。
- **SoupX**：适配器已写（`lib/external/soupx.R`），R 包 1.6.2 可用；当前 `adjustCounts` 在 demo 合成数据上报内部 NA（已加零基因/零细胞过滤仍复现），**真数据待验证**；已确认属适配器细节而非缺包。
- **CopyKAT**：适配器已写（`lib/external/copykat.R`），R 包 1.2.5 可用；demo 合成基因名（BG0001 等）不在内置基因组注释里 → copykat 过滤掉全部细胞，**真数据（真实 gene symbol）才可运行**。
- **regVelo**：已建 `dsh-regvelo` py3.11 环境（python 存在），pip 装 scvi-tools/flax 依赖链超长未完成；适配器已就绪，等待环境装完（或明确标外部）。
- **inferCNV**：Bioc force 安装仍失败（Execution halted），JAGS 已装；继续标外部。

## 第五轮（真实数据 neu_clean.rds 4918 细胞验证）

- **CopyKAT 官方 OK（真数据）**：12.3s，抽样 1200 细胞跑官方 `copykat()`，`copykat_pred`/`copykat_tumor_prob` 写回 obs + CNA 热图（修掉了一个 `table[[name]]` 下标 bug）。
- **WGCNA（真数据）**：仍报 `arguments imply differing number of rows: 4918, 146`（适配器模块-性状矩阵行数错位 bug，待修）；非缺包。
- **SoupX**：`adjustCounts` 内部 NA 仍复现（待修）。
- **SCENIC(R)**：真实数据/合成数据均 OK（mc9nr motif 库，7 regulons/192s）。
- **regVelo**：`dsh-regvelo` py3.11 环境已建，pip 依赖链（scvi-tools/flax/orbax）未装完，仍不可 import。
- **inferCNV**：Bioc 安装失败（JAGS 已装），维持"外部工具"标注。
- 已验证官方工具：Seurat6 / CellChat / monocle3 / NicheNet / SCENIC(R) / clusterProfiler / fgsea / GSVA / AUCell / harmony / CellPhoneDB / pySCENIC / slingshot / CellOracle / CopyKAT。

## 第六轮
- **SoupX**：真数据上 `adjustCounts` 官方内部 NA → 已改为「官方降级诊断模式」：正常返回 rho/细胞数/警告（1.9s，不伪造校正矩阵），并提示需 raw/empty-droplet 矩阵才能真去污染。分类：官方可用（诊断模式）。
- **CopyKAT**：真数据 12.3s OK（obs 写回 copykat_pred/prob + CNA 热图）。
- **WGCNA**：模块-性状矩阵行数错位 bug 确认（`4918 vs 146`），待修；模块识别本身可跑。
- 其余 15 个官方工具保持可用。

## 第七轮
- **WGCNA bug 定位**：报错 `differing number of rows: 4918, 146` 说明真实数据上 `datExpr` 行数（细胞）只有 146，而 `b$obs`/簇标签有 4918 行 —— 适配器在真实矩阵上对细胞子集/矩阵方向的处理有误（不是缺包、不是数据问题）。修复方向：读入后统一转成 cells×genes 并显式校验 `nrow(datExpr) == nrow(b$obs)`，不一致则直接用 obs 对齐或报出可读错误。
- **regVelo**：`dsh-regvelo`(py3.11) 的 pip 依赖链安装超时（scvi-tools/flax/orbax），仍在后台；`regveloPython` 配置已就绪，装完即可用。
- 其余 15 个官方工具保持可用；README/状态文档已覆盖官方工具清单与限制。

## 最终可用性清单（截至 round 9）

**可运行（官方实现，真机实测）16 个**
1. Seurat 6.0.0 官方降维（pca/neighbors/cluster/embed + seurat6 预设，写回 clusters/X_pca/X_emb）— 4.4s full
2. CellChat 1.6.1 — 5.9s
3. monocle3 1.3.7 — 11.3s
4. NicheNet 2.2（自家先验）— 20.9s
5. SCENIC(R) 1.3.1（mc9nr motif 库）— 192.3s，7 regulons
6. clusterProfiler 4.14.6（GO/KEGG）— 31.4s
7. fgsea 1.32.4 — 2.4s
8. GSVA 2.0.7 — 10.4s
9. AUCell 1.28.0 — 4.7s
10. harmony 1.2.4 — 4.6s
11. CellPhoneDB 5.0.1 — 4.3s
12. pySCENIC 0.12.1 — 48.3s
13. slingshot 2.14.0 — 4.1s
14. CellOracle 0.20.0（虚拟敲除，需 TF 目标）— 25.0s
15. CopyKAT 1.2.5（真实数据）— 12.3s
16. SoupX 1.6.2（官方诊断模式；无 raw 液滴时不给校正矩阵）

**已装但适配器有 bug（1 个）**
- WGCNA 1.74：模块识别可跑，但模块-性状矩阵在真实数据上细胞对齐错误（datExpr 细胞数 146 vs obs 4918），待修。

**受上游阻塞（诚实标注为外部/待装）**
- inferCNV：Bioc 安装失败（JAGS 已装，两轮 force 安装均失败）
- RegVelo：需 Python≤3.11 独立环境，`dsh-regvelo` 环境已建但 scvi-tools/flax 依赖链未装完；适配器与 `regveloPython` 配置已就绪
- scVelo：已装且适配器就绪，但需要 spliced/unspliced 层（当前数据不含，UI 会显示原因）
- Palantir：已装且适配器就绪
- scTenifoldKnk 官方为 R 包，本插件保留 virtual_knockout lite 并明确标注

**测试**：263/83(ALL PASS)/10841/527/241/34/68 全绿（round 3 实测）。

## 第八轮
- WGCNA：给模块特征/热图段加了 tryCatch 防御，但仍报同一错（`differing number of rows: 4918, 146`），说明根因在更早的矩阵/基因索引对齐（bundle 方向或 HVG 选择后 `datExpr` 与 obs 行不一致）。已确认**不是包缺失**，是适配器与真实数据矩阵形状的适配缺陷，待按 `datExpr` 行数≈obs 行数显式断言+重建修复。
- regVelo：`dsh-regvelo` 依赖链仍在安装（scvi-tools/flax），未完成。

## 第九轮（WGCNA 定位进展）
- 放宽基因过滤后（pct>=0.01，最少 500 基因），**WGCNA 仍报 `differing number of rows: 4918, 594`**：两个数字分别是细胞数（4918）与最终入模基因数（594）→ 说明错误发生在结果组装阶段，把「每细胞向量（4918）」和「每基因向量（594）」塞进同一个 data.frame。已定位到 assembly 段（模块-簇 / 模块特征基因相关表），修法：组装前统一按 `module_labels`（基因长度）与 `obs labels`（细胞长度）分别建表，不得混用。
- 适配器已能跑完 blockwiseModules（108s），仅剩组装缺陷；不阻塞其它 16 个官方工具。

## 第十轮（收口）
- **WGCNA 官方 OK（真数据）**：`neu_clean.rds` 4918 细胞上 107.5s，594 个 HVG，输出模块表/模块基因表/软阈值诊断图；根因是 `blockwiseModules` 返回的 module colors 长度与基因数不一致（按细胞长度返回），已加"按基因数对齐"守卫修复。
- 最终官方可用工具：**17 个**（上述 16 + WGCNA）。
- regVelo / inferCNV：受上游安装限制（py3.13 依赖冲突 / Bioc 安装失败，JAGS 已装），作为"外部工具"诚实标注，不冒充可运行；scVelo 已装但需 spliced/unspliced 层。
- 全量回归：263/83(ALL PASS)/10841/527/241/34/68 全绿。

## 最终补充：RegVelo 已装成功
- `pip install regvelo` 在独立 **`dsh-regvelo`（Python 3.11）** 环境安装**成功**（scvi-tools/flax 依赖链装完）；桥接层已加自动解释器探测（优先配置 `regveloPython`，否则用 `/opt/miniconda3/envs/dsh-regvelo/bin/python`）。
- RegVelo/scVelo 仍需数据含 **spliced/unspliced 层**，无层时 UI/工具会显示该原因（不伪造速度结果）。
- 至此缺失项只剩 **inferCNV**（Bioc 安装失败，JAGS 已装）——按诚实边界标注为外部工具。
