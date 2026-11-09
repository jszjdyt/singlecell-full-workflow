# dsh-singlecell 参赛介绍

> **仓库**：https://github.com/jszjdyt/singlecell-full-workflow · **插件包名**：dsh-singlecell
>
> **一句话**：面向生物信息学小白和代码苦手的单细胞分析工作台 —— 在 DeepSeek Harness 里，GEO 数据一键读、流程点着跑、官方工具零配置调用，还能直接对 Agent 说人话；同时通过四个带契约的扩展点，让会写代码的人不断把它变大。

---

## 一、背景：生信入门的五道墙

单细胞分析真正劝退新手的，从来不是算法，而是这五件事：

1. **格式墙**：一个 GEO 数据集能同时给出 10x 三件套、`.h5`、`.loom`、`series_matrix.txt.gz`、CSV、Seurat `.rds`、`.RData`、`.tar.gz`……每个都要单独查教程。
2. **环境墙**：R 要装 Seurat/SeuratObject/Bioconductor，Python 要装 scanpy/anndata，版本冲突一天就没了。
3. **分组墙**：GEO 常常「一个 GSM 一个文件」，样本名与大组信息只存在于文件名里，用界面根本没法指定分组，组间比较直接卡死。
4. **工具墙**：CellChat、monocle3、SCENIC、NicheNet、inferCNV……每个一套 API、一篇教程、一堆依赖。
5. **Agent 墙**：大模型能写代码，但它不知道你本机装了什么、工作台里有什么数据，结果就是一遍遍「重新写脚本」。

**dsh-singlecell 把五道墙都拆进 DeepSeek Harness Web GUI**：给不写代码的人一个可视化工作台，给想偷懒的人一套 Agent 工具，两者共用同一套实现。

## 二、目标用户

| 用户 | 他们的痛 | 插件给什么 |
|---|---|---|
| **生物信息学小白** | 不懂 R/Python，只会点网页；GEO 数据下载完就懵 | 可视化导入、点按钮跑流程、自动出图出表 |
| **代码苦手（湿实验同学）** | 能看懂图，不会调包；被环境配置反复折磨 | 一键安装依赖检测、报错说人话、导出 methods/report 直接用 |
| **干湿结合的研究者** | 要自己选样本分组、做组间比较、写 Methods | 样本&分组编辑器、可审计 protocol、英文 Methods 一键生成 |
| **生信老手** | 不想重复写导入/画图胶水代码 | 把重复劳动交给工作台，自己专注写新的分析模块（可拓展性） |

## 三、我们做了什么

### 1）GEO 上能下到的格式，全都能读

h5ad、10x `.h5`、10x MatrixMarket、loom、CSV/TSV/TXT（含 `.gz`）、GEO `series_matrix.txt.gz`、**Seurat `.rds`（v3/v4/v5）**、**`.RData/.rda`（多对象自动列候选）**、zip/tar 压缩包、URL 直链、GEO `GSExxxxx` 号。
**多文件合并导入**：一次选多个 `GSM*.csv.gz`（或目录），每个文件=一个样本，自动生成 `sample` 列，基因取交集。

### 2）样本 & 大组编辑器（本次参赛核心场景，不需要写代码）

- 表格逐行改**样本名**、指定**大组**（young/old、GC/LM…可新建）；
- 按文件名正则**批量自动分组**（如 `(GC|LM)`）；
- 导入「样本→大组」对照 CSV；
- 导入细胞级 metadata CSV（按 barcode join，显示匹配率）；
- 设置分组列后，**所有组间分析自动使用它**。

### 3）可视化 10 步流程 + 交互结果

QC → 过滤 → 标准化 → HVG → PCA → 邻接图 → 聚类 → 二维嵌入 → 标记基因 → 细胞类型注释。
交互式嵌入图（缩放/拖拽/悬停/按任意基因或元数据着色）、QC 图、marker dotplot/heatmap、簇大小、桑基图、火山图。
另有：确定性**自检**（10 项，0–100 分）、**自动注释建议**（内置参考 marker 面板）、**分辨率扫描**（ARI 对比）、**预设流程**（standard/quick/qc/**seurat6**）。

### 4）高级分析：27 项原生 + 17 项官方工具

原生（本插件实现，装好 Python 就能跑）：细胞周期、亚型细分、批次校正、组间比例差异、组间差异基因、火山图、桑基图、GO/KEGG、GSEA、GSVA、AUCell、cNMF、WGCNA、伪时间、PAGA、熵干性、CytoTRACE、CNV/ITH、配体-受体通讯、SCENIC-lite、虚拟敲除/过表达/筛选。

**官方工具桥**：本机装了的官方 R/Python 工具会被**自动探测并调用官方实现**，界面显示真实包版本：
Seurat 6.0（独立 R 库，不污染系统版本）、CellChat、monocle3、NicheNet、SCENIC(R)、clusterProfiler、fgsea、GSVA、WGCNA、AUCell、harmony、CellPhoneDB、pySCENIC、slingshot、CellOracle（虚拟敲除）、CopyKAT、SoupX。
实测耗时 2–192 秒（demo 数据）。**官方版与 lite 版并存且明确标注，绝不拿 lite 冒充官方**；装不上/数据不满足的会直接显示原因。

### 5）可审计产物

每步 run 记录参数、耗时、引擎版本、随机种子；可导出 `protocol.json`（完整溯源）、`methods.md`（英文 Methods 可直接进论文）、`report.md`（摘要+图+marker 表+自检）、`script.py`（scanpy 风格复现脚本）、`.h5ad`。

### 6）Agent 双入口

页面按钮与 Agent 工具**共用同一套实现**。对话里可以直接说：

> 用单细胞工作台把这几个 GSM 文件合并导入，样本名改成 GC1/GC2/LM1，分组 young/old，跑 Seurat 6.0 流程，然后按 group 做组间差异基因和通讯分析

Agent 会调 `singlecell` 工具完成，页面同步显示状态与图。**这就是「AI 可操作」而不是「AI 帮你写个脚本」**。

## 四、可拓展性（重点）

插件不是黑盒，四个扩展点都有契约 + 测试，第三方可以按需长大：

1. **新增分析模块**（零 UI 改动）
   新建 `lib/sc_analyses_xxx.py`，实现：
   ```python
   ANALYSES = [{"id": "my_score", "title": "我的评分", "group": "基因集合分析",
                "method": "本插件实现", "params": [...], "needs": []}]
   def run(analysis_id, ctx, params, outdir):
       return {"summary": {"bullets": [...]}, "tables": [...],
               "figures": [...], "obs_values": {"my_score": [...]}}
   ```
   引擎自动发现 →「高级分析」目录和 Agent 工具里同时出现。

2. **新增官方工具适配器**
   新建 `lib/external/<tool>.R` 或 `.py`，统一 CLI `--bundle/--out/--params`，写 `result.json` 即被桥接层接管；探测、数据包导出、obs 回写、图注册全部自动。本仓库 17 个官方工具就是这么接的。

3. **新增导入格式**
   在 `lib/sc_io.py` 加一个分支即可（现有 10+ 格式与统一输出契约）。

4. **新增预设流程（recipe）**
   在 `lib/sc_engine.py` 的 `RECIPES` 里加一条 JSON，UI 下拉自动出现。

> 对插件事务所式的生态来说：**别人做算法模块，我们负责导入、可视化、Agent 接线、测试与产物**，互不重复造轮子。

## 五、技术亮点（工程实现）

- **纯 JS 客户端插件**：只 `require('react')`，不依赖任何 Harness 内部包，主题/深浅色自适应；
- **前缀路由**：`/singlecell/*` 一条进 Host，新增接口无需重载插件；
- **JSON-lines Python worker**：长驻进程、数据集常驻内存、请求排队/超时/崩溃自动重启；导入只读一次；
- **官方工具桥**：探测本机 R/Python 包、按需导出小数据包、统一 adapter 契约、obs 回写、失败原因可读；
- **不污染系统环境**：Seurat 6.0 装独立 R 库、regVelo 用独立 Python 3.11 环境、跨平台解释器探测（macOS/Linux/Windows）；
- **工具输出 lossless JSON 校验**：避免 `format` 这类 JSON Schema 保留字坑（已踩坑并写进契约）；
- **可审计**：protocol/methods/report/script 一键导出。

## 六、验证数据（都是本机实跑，不是演示视频）

| 项 | 结果 |
|---|---|
| 引擎主流程测试 | 263/263 PASS |
| 导入层测试 | ALL PASS（23 种格式/歧义/安全/URL/真实 rds） |
| 原生分析测试 | 10841 断言 + 527 checks ALL PASS |
| Seurat 6.0 测试 | 241 PASS / 0 FAIL |
| Host↔引擎联调 | 34/34 |
| 前端测试 | 68/68 |
| 真实数据 | 5 个 Seurat rds 全部读通（4918×26609 / 6103×26609） |

## 七、与现有方案的差异

| 方案 | 定位 | 我们的差异 |
|---|---|---|
| CellxGene | 数据浏览 | 我们做**分析全流程**，不只是看 |
| scanpy / Seurat | 代码库 | 我们用**界面与 Agent** 包装它们，并且直接调用你本机官方实现 |
| Galaxy / Cellenics | 通用流程平台 | 我们深度集成 **Agent**，且导入格式/分组编辑/中文场景更贴 GEO 现实 |
| 通用 AI 编程 | 写脚本 | 我们让 AI **操作一个真实工作台**，状态可审计、结果可复现 |

## 八、3 分钟体验路径

1. 侧边栏「单细胞」→ 顶栏「导入数据」→ 粘贴多个 `GSM*.csv.gz` 路径或选目录合并；
2. 「记录」→「样本 & 分组」→ 改样本名、选 young/old → 保存分组；
3. 「高级分析」→ 选「同细胞类型组间差异基因」→ 分组列已默认选中 → 运行。

## 九、诚实边界

- 27 项原生分析是本插件实现（`method` 标 lite/近似/代理），不等同官方工具结果；
- 官方工具能否用取决于本机装了什么，插件只探测与调用；装不上的（如 inferCNV）显示原因，不伪造；
- RNA velocity / RegVelo 需要 spliced/unspliced 层，没有就拒绝并提示；
- 自动注释是参考面板建议，自检是规则检查，不替代专家判断。

## 十、路线图

- 样本分组：支持按细胞名/文件交叉校验、分组方案保存为模板；
- 官方工具：inferCNV 换安装方式，补充 CellChat 官方多组对比、SCENIC+；
- 可拓展性：发布「插件生态 SDK」（模块模板 + adapter 模板 + 一键校验）；
- 产物：报告生成 Word/PDF、图一键进论文排版。

---

## 附：短版介绍（表单可直接粘贴，约 300 字）

> **dsh-singlecell** 是 DeepSeek Harness 上的单细胞分析工作台，面向**生物信息学小白和代码苦手**：GEO 数据不用转格式（10x/h5ad/loom/rds/RData/压缩包/URL 全支持，多文件自动合并成多样本），样本名与大组在界面表格里点选、组间分析自动可用；10 步流程可视化点击运行，交互嵌入图与出版级图即时生成；27 项分析内置，本机装了的 CellChat/monocle3/SCENIC/WGCNA/CellOracle 等 17 个官方工具会被自动探测并调用官方实现，官方与 lite 明确标注、不冒充。页面按钮与 Agent 工具共用同一套实现，对话里说一句「合并导入这些 GSM、设 young/old 分组、跑 Seurat 6.0 流程、做组间差异基因」即可完成，并导出可审计的 protocol/methods/report/script。插件提供四个带契约的扩展点（分析模块、官方工具适配器、导入格式、预设流程），便于社区持续扩展。

## 附：宣传语候选

1. 生信小白的第一个单细胞工作台，代码苦手的最后一篇教程。
2. 不再和格式、环境、分组三座大山搏斗 —— 点几下，或者说一句话。
3. 你本机的官方工具，终于有人替你接线。
4. AI 不只是帮你写脚本，它直接在真实工作台上替你操作。
