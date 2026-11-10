# dsh-singlecell 内部契约（Single-cell workbench）

> 本文件是三个并行开发者的**唯一事实来源**。任何实现与本契约冲突，以本契约为准；确需改协议，先在 `/singlecell` 群聊里说明并同步更新本文件。

## 0. 目标

在 DeepSeek Harness Web GUI 里提供一个**可视化、可操作、Agent 可驱动**的单细胞分析工作台：

- 左栏：数据导入 + 可点击的流程卡片（QC → 过滤 → 标准化 → HVG → PCA → 邻接图 → 聚类 → 二维嵌入 → 标记基因 → 注释）。
- 中栏：交互式嵌入图（canvas，可缩放/悬停/点聚类高亮）、QC 图、标记基因 dot plot、聚类柱状图、活动日志。
- 右栏：Agent 面板（当前状态、建议下一步 + 复制指令、Agent/用户操作时间线）。
- Host 同时注册 `singlecell` agent 工具：Agent 用**同一套** Host service 方法驱动全部操作（`lib/ops.mjs` 是唯一实现）。
- 计算引擎：Python worker 子进程（JSON-lines over stdio）。**自动探测 scanpy**；没有 scanpy 时使用内置 numpy/scipy/sklearn 实现（本机就是这种情况）。UI 必须显示当前引擎模式。

**硬性要求**：无 Python 依赖时 UI 必须给出可读的原因，不能白屏；单个步骤失败不能影响其它步骤与整页。

## 1. 目录与归属（写权限严格隔离）

```
plugins/dsh-singlecell/
  CONTRACT.md                    # 本文件（Lead 维护）
  README.md                      # Lead（用户文档）
  package.json / cordis.patch.yml / icon.svg / locale/*.json   # Lead
  index.js                       # Lead   Host 入口（单条前缀路由 + 工具 + prompt section）
  lib/impl.mjs                   # Lead   ★唯一操作实现（HTTP 与工具共用）+ 路由 + 工具 + prompt
  lib/engine.mjs                 # Lead   Python worker 管理（探测/spawn/队列/超时/重启）
  lib/sc_engine.py               # sc-engine  ★Python worker（协议 §2/§6 + 导入/分析分发）
  lib/sc_io.py                   # Lead   ★全格式导入层（§8；R 桥调用方）
  lib/r_bridge.R                 # sc-io  Rscript 桥（rds/RData → mtx/features/barcodes/obs/meta）
  lib/analysis_registry.py       # Lead   高级分析目录/ctx/分发 + 外部工具项登记（§9）
  lib/sc_analyses_sets.py        # sc-analysis-a  基因集/通路/矩阵类 14 项
  lib/gene_sets_builtin.py       # sc-analysis-a  内置 S/G2M 基因集 + 36 条通路 GMT
  lib/sc_analyses_dynamics.py    # sc-analysis-b  整合/轨迹/异质性/通讯/调控/虚拟敲除 13 项
  lib/biopriors.py               # sc-analysis-b  204 LR 对 + 644 TF + 家族映射
  client.js                      # sc-ui      ★浏览器 UI（§4/§9.5：导入弹层 + 7 tab）
  test/make_sample.py            # sc-engine  合成 h5ad + demo（含细胞周期背景基因）
  test/py_pipeline_test.py       # sc-engine  引擎主流程测试（263 项）
  test/engine.test.mjs           # Lead      node↔python 联调测试（34 项）
  test/make_formats.py           # sc-io      生成 23 种格式 fixture + truth.json
  test/py_import_test.py         # sc-io      导入层回归（83 项）
  test/py_analyses_sets_test.py  # sc-analysis-a  14 项分析测试（10841 断言）
  test/py_analyses_dynamics_test.py # sc-analysis-b 13 项分析测试（527 checks）
  test/client-logic.test.mjs     # sc-ui      纯函数测试（56 用例）
  examples/                      # Lead   示例产物（嵌入图 / dot plot / QC / 分析图 / report）
```

不要改别人作用域里的文件；需要对方改的东西，写进任务描述/消息里。**所有人写文件前先 `read` 目标文件**。

## 2. Python worker 协议（`lib/sc_engine.py`）

### 2.1 启动与传输

- Host 以 `python3 -u lib/sc_engine.py --workdir <abs_dir>` 启动；`--workdir` 已由 Host 创建，可写。
- **stdin**：一行一个 JSON 请求；**stdout**：一行一个 JSON 响应（UTF-8，`ensure_ascii=False`，`flush=True`）。
- stdout 只能有协议行；日志/警告一律写 stderr（Host 会收集到 `<workdir>/engine.log`）。
- 请求：`{"id":"<str>","cmd":"<name>","args":{...}}`
- 成功：`{"id":"<str>","ok":true,"result":{...}}`
- 失败：`{"id":"<str>","ok":false,"error":{"code":"bad_request|not_found|prerequisite|unsupported|internal","message":"中文可读原因","trace":"可选"}}`
- 可选进度：`{"id":"<str>","event":"progress","message":"...","pct":0-100}`（Host 只在同 id 请求未完成时缓存展示；可少发、可不发）。
- 引擎必须**逐条尽力容错**：单条命令异常只影响该条响应，进程不退出；只有 `shutdown` 或 stdin EOF 退出。
- 启动时不得向 stdout 打印任何东西（包括欢迎语）。

### 2.2 命令

> 命名/字段必须逐字一致。未列出的字段可以加（向前兼容），已列出字段不能少。

**`ping`** args `{}` →
```json
{"python":"3.14.2","platform":"darwin","mode":"builtin|scanpy",
 "packages":{"numpy":"2.5.1","pandas":"2.3.3","scipy":"1.17.1","sklearn":"1.9.0","anndata":"0.12.19","scanpy":null},
 "layouts":["tsne","pca"],"warnings":["..."]}
```
`mode` 决定 UI 徽章：检测到 scanpy 用 `"scanpy"`，否则 `"builtin"`；`layouts` 至少含 `"tsne","pca"`，有 umap-learn/scanpy 时加 `"umap"`。

**`open`** args `{"path"?:str,"demo"?:bool,"name"?:str}` → 见 §2.3 `dataset` 对象。
- `path` 为 `.h5ad`（必须支持；也读 `.h5` + `adata` key 余力支持）。路径不存在/解析失败 → `bad_request`。
- `demo:true` → 用 `test/make_sample.py` 的同一生成器在内存造数据（**不落盘 html**），名 `"demo-pbmc"`。
- 多个数据集可同时存在；`id` 由引擎生成（`ds1`,`ds2`,…）。同一路径重复 open 返回已存在的数据集（不重复加载）。

**`list`** args `{}` → `{"datasets":[<dataset 摘要>]}`

**`close`** args `{"dataset":id}` → `{"closed":true}`

**`preview`** args `{"dataset":id}` →
```json
{"genes":["CD3D", "...首 50 个"],"obs_columns":[{"key":"n_genes","kind":"numeric"},
 {"key":"leiden","kind":"categorical"}],"n_obs":2700,"n_vars":1200,"steps":["qc"],"current":"qc"}
```

**`run`** args `{"dataset":id,"step":str,"params":{}}` →
```json
{"step":"pca",
 "summary":{"title":"PCA","bullets":["取 30 个主成分","PC1 解释方差 12.3%"],
   "metrics":{"n_pcs":30,"pc1_variance_ratio":0.123}},
 "steps":["qc","filter","normalize","hvg","pca"],"current":"pca",
 "dataset":{"n_obs":2700,"n_vars":1200}}
```
步骤定义（顺序即依赖；缺前置 → `prerequisite`；重跑某步自动清空其下游状态）：

| step | 依赖 | params（默认值） | 语义 |
|---|---|---|---|
| `qc` | — | `min_genes=200,max_genes=null,min_counts=null,max_pct_mt=20` | 只算指标（`n_genes,n_counts,pct_counts_mt`），不过滤 |
| `filter` | qc | 同上 + `min_cells=3` | 过滤细胞/基因，`metrics` 返回前后数量 |
| `normalize` | filter | `target_sum=10000` | `log1p(X/row_sum*target_sum)` |
| `hvg` | normalize | `n_top_genes=2000,flavor="seurat"` | 记 `var["highly_variable"]` |
| `pca` | hvg | `n_pcs=30` | `obsm["X_pca"]` |
| `neighbors` | pca | `n_neighbors=15,metric="cosine"` | kNN 图 + 连通矩阵 |
| `cluster` | neighbors | `resolution=0.2,method="louvain"` | `obs["clusters"]`（字符串），`obs["cell_type"]` 同步为簇号 |
| `embed` | pca | `layout="umap",perplexity=30`（无 UMAP 时回退 t-SNE） | `obsm["X_emb"]`、`obs["_emb_x/_emb_y"]`；layout 不可用要回退并在 summary 说明 |
| `markers` | cluster | `groupby="clusters",top_n=10,method="t-test"` | 结果缓存，`metrics.n_markers` |
| `annotate` | cluster | `mapping={"0":"T cell",...}` | 写 `obs["cell_type"]`，`metrics.n_types` |

`qc`/`filter` 的阈值参数 UI 可改；默认值必须与上表一致。

**`embedding`** args `{"dataset":id,"color_by":{"kind":"cluster|obs|gene","key":"clusters"},"max_points":12000}` →
```json
{"layout":"tsne","n_points":2700,"downsampled":false,
 "x":[-1.24,...],"y":[...],
 "continuous":false,
 "legend_title":"cluster",
 "categories":[{"id":"0","label":"0","color":"#4e79a7","n":312}],
 "group_index":[0,0,1,...]}
```
- `kind="gene"`：`continuous:true`，`values`（log-normalized 表达，缺列=全 0 + warning），`categories` 换成 min/max 描述。
- `kind="obs"`：数值列→连续，类别列→分类。
- 坐标保留 3 位小数。`x/y/categories/group_index` 长度必须与 `n_points` 一致。

**`plot`** args `{"dataset":id,"kind":str,"params":{}}` →
```json
{"file":"/abs/path/fig-xxx.png","width":720,"height":480,"bytes":41233}
```
`kind` ∈ `qc_violin`（参数 `metrics:["n_genes","n_counts","pct_counts_mt"]`、可选 `thresholds`）、`embedding`（`color_by` 同 §embedding，画带聚类质心标签的散点）、`dotplot`（`top_n=5`，cluster × marker gene）、`heatmap`（`top_n=5`）、`cluster_bar`（各簇细胞数）、`pca_variance`（碎石图）、`qc_scatter`（`x_metric,y_metric,color_by_pct_mt`）。
- PNG，`dpi≈140`，白底，**图上文字一律英文**（避免缺 CJK 字体）。
- 文件写到 `<workdir>/figures/`，文件名含数据集 id + kind + 参数哈希（同参数复用）。
- 缺前置数据 → `prerequisite`。

**`markers`** args `{"dataset":id,"groupby":"clusters","top_n":10,"method":"t-test"}` →
```json
{"groups":[{"id":"0","label":"0","n_cells":312,
  "genes":[{"gene":"CD3D","score":8.4,"logfc":2.31,"pval_adj":1.2e-30,"pct_in":0.92,"pct_out":0.11}]}]}
```

**`annotate`** 由 `run` 的 `annotate` 步承担，无独立命令。

**`suggest`** args `{"dataset":id}` →
```json
{"suggestions":[{"id":"run:cluster","title":"运行聚类","detail":"邻接图已就绪，建议 resolution=0.2 聚类","step":"cluster","params":{"resolution":1.0}}]}
```
规则：无数据 → 建议 open/demo；有数据未 qc → qc；…；全流程完成 → 建议 markers/embedding/导出。**确定性**，不调用模型。

**`export`** args `{"dataset":id,"what":"h5ad|markers_csv|obs_csv|summary_json","path"?}` →
```json
{"file":"/abs/path/xxx","bytes":12345,"what":"h5ad"}
```
`path` 缺省 → `<workdir>/exports/<dataset>-<what>-<时间戳>`；`h5ad` 里要带 `obs`（含簇/细胞类型）、`obsm`（PCA/嵌入）、`X`（当前标准化矩阵）、`layers["counts"]`（原始计数）。

**`shutdown`** args `{}` → `{"bye":true}` 后退出进程。

### 2.3 `dataset` 对象

```json
{"id":"ds1","name":"demo-pbmc","path":"/abs/x.h5ad","n_obs":2700,"n_vars":1200,
 "obs_columns":[{"key":"clusters","kind":"categorical","n_unique":8}],
 "steps":["qc","filter"],"current":"filter",
 "warnings":["列 X 缺失，已按稠密矩阵处理"]}
```

### 2.4 内置引擎的计算要求（正确性，会被测试检查）

- 内部用 `scipy.sparse`（X 稀有）或 `numpy`（稠密）；`float32`。
- QC：`n_genes`=每细胞非零基因数；`n_counts`=总计数；`pct_counts_mt`=基因名以 `MT-` 开头（大小写不敏感）的计数占比。
- 标准化：`log1p(x / row_sum * target_sum)`，row_sum=0 的细胞置 0。
- HVG（builtin）：按均值分箱（20 箱）的标准化离散度，取 top N，结果布尔 mask。
- PCA：对 HVG、标准化、零均值（`std` 归一可选）后的矩阵做 `sklearn.decomposition.PCA(n_components=n_pcs, svd_solver="randomized")`。
- 邻接图：`sklearn.neighbors.NearestNeighbors`，k=`n_neighbors`。
- 聚类（builtin）：**纯 Python Louvain**（模块度优化 + 分辨率参数）在 kNN 图上跑；`leidenalg`/`scanpy` 存在时可用真 Leiden。簇 id 一律字符串 `"0","1",…`，按大小降序重排。
- 嵌入（builtin）：`sklearn.manifold.TSNE(perplexity=min(30,(n_obs-1)/3), init="pca", random_state=0)`；`n_obs>15000` 时提示并允许（可能慢）。无 umap 时 `layout="umap"` 回退到 tsne 并在 `summary.bullets` 写明。
- 标记基因（builtin）：one-vs-rest Welch t 检验（向量化）+ BH FDR + `pct_in/pct_out`；按 `score` 降序。
- 所有随机种子固定（`random_state=0` / `np.random.default_rng(0)`），同参数同结果（测试要求）。
- demo 数据：≥2000 细胞，≥1000 基因，6 个真实类型（每个类型有 12–20 个 marker 基因高表达），带线粒体基因与一个"低质量"亚群，保证 QC 过滤有意义；生成 < 3 秒。**用同一份生成器给 `test/make_sample.py` 复用**。

### 2.5 sc-engine 的测试要求

`test/py_pipeline_test.py`（不依赖 node，`python3 test/py_pipeline_test.py` 直接跑）：
1. demo open → run qc/filter/normalize/hvg/pca/neighbors/cluster/embed/markers 全链路；
2. 注入的 marker 基因必须出现在对应簇的 top markers（命中率 ≥ 60%）；
3. 同参数两次 `cluster` 完全一致（确定性）；
4. `embedding` 返回长度一致、坐标有限；
5. `plot` 生成 PNG 且首 8 字节是 PNG magic、宽高 > 200；
6. 缺前置时返回 `prerequisite` 错误码；
7. `open` 一个不存在的路径返回 `bad_request` 且进程不退出。

## 3. HTTP API（Host，前缀 `/singlecell`）

- 只接受回环地址（照抄 `dsh-paper-reader/lib/impl.mjs` 的 `hostIsLoopback` 逻辑）。
- 一律 `content-type: application/json; charset=utf-8` + `cache-control: no-store`（图片路由除外）。
- 请求体 JSON（upload 除外）。响应成功 `{"ok":true, ...}`，失败 `{"ok":false,"error":"中文可读原因","code":"..."}`（HTTP 400/404/500 视情况，前端按 body 判断）。
- 请求可带 `"actor":"user"`（UI 默认）；没有则 `"system"`。Agent 工具路径固定 `"agent"`。所有操作都要写进活动日志。

| method path | body | 说明 |
|---|---|---|
| `GET /singlecell/health` | — | `{"engine":{"ready":bool,"mode":"builtin\|scanpy\|none","python":str,"packages":{},"error":str?}}` |
| `GET /singlecell/state` | — | `{"active":id?, "datasets":[dataset摘要], "engine":{...}, "suggestions":[...], "activity":[last 50]}` |
| `POST /singlecell/dataset/open` | `{path}` | 打开 h5ad |
| `POST /singlecell/dataset/demo` | `{}` | 载入内置示例 |
| `POST /singlecell/dataset/upload` | raw bytes + header `x-dsh-filename` | 落盘到 workdir 后打开 |
| `POST /singlecell/dataset/select` | `{dataset}` | 切当前数据集 |
| `POST /singlecell/dataset/close` | `{dataset}` | 关闭 |
| `POST /singlecell/preview` | `{dataset}` | 基因/obs 列清单 |
| `POST /singlecell/run` | `{dataset,step,params}` | 跑一步 |
| `POST /singlecell/embedding` | `{dataset,color_by}` | 嵌入坐标 |
| `POST /singlecell/plot` | `{dataset,kind,params}` | → `{"figure":{"id","url","width","height","kind"}}` |
| `POST /singlecell/markers` | `{dataset,top_n,method}` | 标记基因 |
| `POST /singlecell/export` | `{dataset,what,path}` | 导出 |
| `GET /singlecell/figure?id=...` | — | 回 PNG，`id` 只允许 `[A-Za-z0-9._-]+\.png` |
| `GET /singlecell/activity` | — | `{"entries":[...]}` |
| `POST /singlecell/reset` | `{}` | 清空所有数据集与日志 |

活动日志条目：`{"at":"ISO","actor":"agent|user|system","action":"run:cluster","dataset":"ds1","detail":"resolution=0.8 → 8 簇","ok":true,"ms":812}`。

## 4. 浏览器端（`client.js`，sc-ui）

### 4.1 形态（照抄 `plugins/dsh-paper-reader/client.js` 的结构）

- 整文件 `;(function(){'use strict' ...})()`；顶部是**纯函数层**（不碰 DOM、不 require react），在 `typeof window==='undefined'` 时挂到 `globalThis.__SC_TEST__`，供 node 测试 import。
- 浏览器层用 `window.__ModuleLoader__.load({id:'dsh-singlecell', factory(require){...}})`，**只 `require('react')`**，绝不 require 任何 `@deepseek-ai/*` 包。
- `inject: ['slots','locale']`；`apply(ctx)` 里：
  - `locale.register('dsh-singlecell', {zh:DICT_ZH, en:DICT_EN})`（用 `ctx.effect` 包起来，照抄 paper-reader）
  - `slots.inject('main', () => slots.register({name:'main', key:'singlecell', locale:NS}, WorkbenchPage))`
  - `slots.inject('sidebar.panellist', () => slots.register({name:'sidebar.panellist', id:'singlecell', order:3, label:()=>T('nav'), locale:NS}, SidebarIcon))`
- 所有数据走**同源相对路径** `/singlecell/*`（`fetch`），不 import Host 模块。
- 必须包一层错误边界（照抄 paper-reader 的 `PluginErrorBoundary` 思路），单点异常不能把整个槽位打黑。

### 4.2 布局（整页，用 `height:100%`，min-height 72vh）

```
┌ 顶栏: 标题「单细胞分析工作台」·引擎徽章·数据集 chips·[载入示例] [打开 h5ad ▾]·[刷新]
├──────────────┬───────────────────────────────────┬──────────────┐
│ 左栏 268px   │ 中栏 flex                        │ 右栏 300px   │
│ 流程卡片      │ Tab: 嵌入图/QC/标记基因/聚类/日志  │ Agent 面板   │
└──────────────┴───────────────────────────────────┴──────────────┘
```

- **左栏流程卡片**（顺序）：qc→filter→normalize→hvg→pca→neighbors→cluster→embed→markers→annotate。每张卡：序号圆点 + 标题 + 状态（待运行/已完成/运行中/需前置）+ 可折叠参数区 + 运行按钮 + 完成后的一行结果摘要。已完成可"重跑"。参数控件自己写（number input / select / text），不要引第三方。
- **中栏 Tab**：
  - `嵌入图`：<canvas> 交互散点。支持：滚轮缩放（以鼠标为锚）、拖拽平移、双击复位、悬停 tooltip（细胞序号 + 分组名）、点击图例 chip 高亮该组（其余降到 15% 透明度）、color-by 选择器（cluster / cell_type / n_genes / pct_counts_mt / 基因搜索框）。数据来自 `POST /singlecell/embedding`；画布要处理 devicePixelRatio，缩放=对 x/y 做线性变换，不要每次缩放重新请求。
  - `QC`：`<img>` 显示 `qc_violin` 图 + 指标多选（勾选生成新图）；无图时给"生成 QC 图"按钮（用户点才请求，避免乱跑）。
  - `标记基因`：表格（每簇 top 基因，列：簇/基因/logFC/adj.P/pct_in）+ 顶部按需生成 `dotplot` 图。
  - `聚类`：`cluster_bar` 图 + 每个簇一行输入框可改细胞类型名（写回 `annotate`）。
  - `日志`：活动列表（时间/actor/动作/结果/耗时）。
- **右栏 Agent 面板**：引擎与数据状态摘要；"建议下一步"卡片（来自 `state.suggestions`，每张带 `复制指令` —— 复制一句可直接发给 Agent 的中文指令，如「用单细胞工作台把当前数据集跑到聚类，resolution=1.0」）；"Agent 活动"时间线（actor=agent 的条目用 🤖，user 用 👤）；底部一句提示"也可以在对话里直接让 Agent 操作这个工作台"。
- 顶部 `打开 h5ad`：可输入绝对路径（回车打开），也可 `<input type=file accept=".h5ad,.h5">` 上传（`fetch('/singlecell/dataset/upload', {method:'POST', body:file, headers:{'x-dsh-filename':file.name}})`)。
- 忙碌态：任何请求进行中按钮 disabled + 文案"运行中…"；错误用行内红色提示（`--dsw-alias-*` 里挑语义色，找不到就用带透明度的 `label-primary` + 文字），不要 alert。

### 4.3 样式约束

- 类名一律 `sc-` 前缀；**颜色只用 `--dsw-alias-*` 主题变量**（派生色用 `color-mix(in srgb, var(--dsw-alias-brand-primary) 12%, transparent)`），不要写死 #hex（聚类调色板例外：颜色由服务端下发或放在 JS 常量里，属于数据色）。
- 字号 12/13px 为主，和宿主页面一致；圆角 8–12px；边框 `--dsw-alias-border-l1`；卡片底 `--dsw-alias-bg-layer-1`，页面底 `--dsw-alias-bg-base`。
- 可见文案全部走 `T('key')`（`locale.register` 之后从 `ctx.locale` 拿翻译；照抄 paper-reader 的 `createTranslator`）。
- 不写 DOM 到组件外，不 append `document.body`，不读别的插件 DOM。

### 4.4 sc-ui 的测试要求

`test/client-logic.test.mjs`（`node --test` 或纯 node 脚本）：
1. 能 `import('../client.js')`（window 未定义时不炸）；
2. `STEPS` 顺序与依赖正确、`nextStep()` 逻辑正确；
3. 数值/百分比格式化、坐标变换（缩放平移后的屏幕坐标）、tooltip 文本正确；
4. API 客户端在注入的假 `fetch` 下生成正确 method/url/body，并能解析 `{ok:false,error}`；
5. 建议→中文指令文本生成正确。

## 5. 环境事实（所有人必读）

- macOS，工作目录 `/Users/dingyuntao/Desktop/deepseek harness`，插件目录 `plugins/dsh-singlecell`。
- Python：`python3` → `/Library/Frameworks/Python.framework/Versions/3.14/bin/python3`（3.14.2），有 numpy 2.5 / pandas 2.3 / scipy 1.17 / sklearn 1.9 / anndata 0.12 / matplotlib 3.11 / h5py；**无 scanpy、无 igraph/leidenalg、无 umap-learn**。
- `~/.matplotlib` 不可写：Python 侧必须 `os.environ["MPLBACKEND"]="Agg"` 且把 `MPLCONFIGDIR` 指到 `--workdir` 下（Host 已设，引擎里再兜底一次）。
- 宿主 Electron 进程的 PATH 可能很窄：Host 探测 Python 时按候选列表逐个 `-c "import numpy,anndata,sklearn"` 试。
- **禁止** `pip install`、禁止改用户 profile、禁止 `pnpm`、禁止写 `$DSH_HOME`。
- 插件安装由 Lead 用 `plugin_manager install_bundle` 完成；开发者只写包内文件并跑自己的测试。
- 提交前自检：`node --check client.js` / `node -e "import('./lib/impl.mjs')"` 之类能跑就跑；Python 用 `python3 -m py_compile lib/sc_engine.py`。

## 6. v2 增补：对齐 Claude Science 等成熟方案

参考要点（Anthropic Claude Science：<https://www.anthropic.com/news/claude-science-ai-workbench>、文档 <https://claude.com/docs/claude-science/overview>）：
每个产物带**可审计的历史与代码环境**、有**复核**环节、**skills/预设流程**可复用、**可视化产物可交互迭代**、参数可**对比**。
单细胞领域的成熟对照物：CellxGene 的交互嵌入、Cellenics/Galaxy 的流程卡片与历史、CellMarker/PanglaoDB 的参考 marker 面板、scanpy 最佳实践。

**总原则：§2/§3/§4 的既有字段一律不改，§6 只做追加。** Host 已按本段实现，缺字段会直接测试失败。

### 6.1 溯源历史（`dataset.history`）

每次 `run` 成功后追加一条（失败也记 `ok:false`，最多保留 50 条）：
```json
{"at":"2026-11-06T06:20:11.240Z","step":"cluster","params":{"resolution":1.0},
 "n_obs":2700,"n_vars":1200,"ok":true,"ms":812,
 "engine":{"mode":"builtin","python":"3.14.2","seed":0,
   "versions":{"numpy":"2.5.1","sklearn":"1.9.0"}},
 "summary":"resolution=1.0 → 8 簇"}
```
`dataset` 对象在所有返回里都带 `history`（可截断为最近 20 条给 UI）。`open` 时 `history:[]`。UI 与工具据此展示「分析记录」。

### 6.2 新命令

**`audit`** args `{"dataset":id}` → 确定性自检（reviewer 的轻量版，不调用模型）：
```json
{"score":78,"level":"warn","checks":[
  {"id":"qc_before_cluster","level":"ok|warn|fail","title":"先质控再聚类",
   "detail":"已完成 qc/filter","hint":"..."}]}
```
必须实现的检查项（id 固定）：`qc_before_cluster`、`mito_high`（过滤后 pct_counts_mt 中位数 >20%）、`cells_low`（<200 细胞）、`clusters_few`（簇数 <2）、`cluster_skew`（最大簇占比 >70%）、`markers_missing`（存在簇无显著 marker）、`embedding_missing`、`annotate_missing`、`doublet_risk`（n_genes 极高细胞占比 >5% 时 warn）、`resolution_untested`（未做分辨率扫描时 warn）。`score = clamp(100 - 40*Fail - 10*Warn, 0, 100)`。

**`autoannotate`** args `{"dataset":id,"groupby":"clusters"}` →
```json
{"panel":"builtin-marker-panel-v1",
 "suggestions":[{"cluster":"0","label":"T cell","score":0.82,"margin":0.19,
   "evidence":["CD3D","CD3E","IL7R"],"alternative":"NK cell"}]}
```
- 内置 ≥14 类参考面板（T、CD4 T、CD8 T、NK、B、Plasma、Monocyte、Macrophage、DC、Neutrophil、Mast、Platelet、Erythrocyte、Endothelial、Fibroblast、Epithelial）。
- 每类 6–12 个 marker 基因；簇得分 = 面板基因在该簇的平均 z-score；返回 top1/top2、`margin=top1-top2`、`evidence`=该簇表达最高的 5 个面板基因。**不写 obs**（由调用方确认后走 `annotate`）。

**`recipe`** args `{}` →
```json
{"recipes":[{"id":"standard","title":"标准流程","detail":"QC→聚类→标记基因→嵌入",
  "steps":[{"step":"qc","params":{}},{"step":"filter","params":{}}]}]}
```
至少 3 个：`standard`（推荐全流程）、`quick`（qc→filter→normalize→hvg→pca→neighbors→cluster→embed）、`qc`（qc→filter + 出 QC 图）。`recipe` 带 `{"id":"standard"}` 时返回该 recipe 的 steps。执行由 Host 顺序调 `run`（便于逐步展示与失败恢复）。

**`resolution_scan`** args `{"dataset":id,"resolutions":[0.4,0.8,1.2]}` →
```json
{"current":{"resolution":1.0,"n_clusters":8},
 "results":[{"resolution":0.4,"n_clusters":5,"ari_vs_current":0.61}]}
```
在既有 kNN 图上重跑聚类并用 `sklearn.metrics.adjusted_rand_score` 算 ARI；**不改变当前聚类**。输入为空时默认 `[0.4,0.6,0.8,1.0,1.2,1.6]`；最多 8 个值。

**`export` 扩展 `what`**：
- `protocol`：JSON（dataset 元信息、source path、n_obs/n_vars、完整 history、engine 版本、seed、图清单、插件版本）。
- `methods`：Markdown，英文 Methods 段落（步骤 + 参数 + 软件版本 + 阈值），可直接用于论文。
- `report`：Markdown 报告（标题、数据摘要、Methods、嵌入图/QC 图相对路径引用、各簇 top markers 表、自动注释建议、audit 结果）。同名 PNG 复制到报告同目录。
- `script`：`.py` 可复现脚本（scanpy 风格；注释说明 builtin 引擎的等价参数，scanpy 可用时为精确复现）。

### 6.3 demo 数据必须配合 autoannotate

demo 的 6 个类型必须取自内置参考面板（例如 T cell / B cell / NK / Monocyte / DC / Platelet），且每类至少 6 个面板基因强表达。测试要求：`autoannotate` 的 top1 建议与真实类型一致 ≥ 4/6。

### 6.4 sc-engine 测试增补

在 §2.5 之外追加：
1. `audit` 在只跑 `qc` 的数据上必须出现 warn/fail（如 `clusters_few`、`embedding_missing`）；全流程跑完后 score 明显更高；
2. `autoannotate` 在 demo 上 top1 命中 ≥ 4/6；
3. `resolution_scan` 返回条目数 = 输入数、ARI ∈ [-1,1]、当前簇未改变；
4. `export protocol/methods/report/script` 文件存在且非空，report 内含 PNG 引用与 markers 表；
5. `history` 随 run 增长且字段齐全（step/params/at/ok/engine.versions）。

### 6.5 sc-ui 的 v2 界面增补

- 顶栏加「预设流程」下拉（来自 `GET /singlecell/state` 的 `recipes`）+「一键运行」，运行中显示「第 n/m 步：xxx」。
- 中栏 tab 变成 6 个：`嵌入图 / QC / 标记基因 / 聚类 / 记录 / 日志`。新增的「记录」tab：
  - 数据集历史（`state.datasets[active].history`）：时间、步骤、参数摘要、结果摘要、耗时；
  - 「自检」→ `POST /singlecell/audit`，列表按 ok/warn/fail 着色（warn 黄、fail 红，用主题色派生）；
  - 「自动注释」→ `POST /singlecell/autoannotate`，每簇一行：建议 label + score + evidence 基因 + 「采用」按钮（把建议汇成 mapping 走 `run:annotate`）+ 「全部采用」；
  - 「分辨率扫描」：逗号分隔输入 → `POST /singlecell/resolution-scan` → 表格 resolution/n_clusters/ARI；
  - 导出按钮组：`protocol / methods / report / script / h5ad / markers_csv` → `POST /singlecell/export`，成功后显示文件绝对路径（可复制）。
- `GET /singlecell/state` 额外返回 `recipes` 与 `datasets[].history`。
- 新增路由：`POST /singlecell/recipe` `{dataset,id}`（Host 负责顺序执行并返回每步结果）、`POST /singlecell/audit` `{dataset}`、`POST /singlecell/autoannotate` `{dataset}`、`POST /singlecell/resolution-scan` `{dataset,resolutions}`。
- 右侧 Agent 面板不变；记录 tab 的「自检/自动注释/扫描」结果都要能通过 Agent 工具得到同样内容（Host 已保证同实现）。

## 8. 全格式导入层（v3，对标 GEO 实际可下格式 + R 生态）

**用户反馈**：只支持 h5ad 太窄。GEO 单细胞补充数据常见格式必须全部可读，R 的 `.rds`/`.RData` 也要读。

### 8.1 文件：`lib/sc_io.py`（sc-io 所有）+ `lib/r_bridge.R`（sc-io 所有）

```python
# 依赖只能是标准库 + numpy/scipy/pandas/h5py/anndata（不得新增第三方依赖）
SUPPORTED_FORMATS = ["h5ad","10x_h5","10x_mtx","loom","matrix_text","geo_series","rds","rdata","archive","auto"]
def inspect_source(source, *, workdir, options=None) -> dict
def load_source(source, *, workdir, options=None, progress=None) -> dict
```

`source` 是本地路径（文件或目录）或 `http(s)://`/`ftp://` URL；`options` 可含 `format`（强制格式）、`member`（压缩包内文件）、`object`（RData 内对象名）、`transpose`（矩阵方向）、`delimiter`、`name`。

`load_source` 返回（**键名固定**）：
```python
{"format":"10x_mtx","source":"/abs/x","name":"GSE123_10x",
 "X": scipy.sparse.csr_matrix(float32, cells×genes),   # 主矩阵（原始计数或已给矩阵）
 "counts": csr | None,                                 # 若格式里同时有 counts 层则给，否则 None
 "obs": pandas.DataFrame(index=cell names),             # 行数 = X.shape[0]
 "var_names": [str],                                    # 长度 = X.shape[1]
 "var": pandas.DataFrame | None,                        # 可选基因注释（行对齐 var_names）
 "warnings": [str], "notes": [str],
 "meta": {"n_obs":int,"n_vars":int,"species":str|None,"assays":[str],"has_counts":bool}}
```
- 失败抛 `EngineError(code, 中文消息)`；`code` ∈ `bad_request|unsupported|not_found|ambiguous`。
- `ambiguous` 时异常对象带 `candidates`（list）→ 引擎透出给前端选择。
- `inspect_source` 不加载完整矩阵（大文件只读头/元数据），返回 `{"format","candidates"?,"n_obs"?,"n_vars"?,"notes":[],"size":int}`。

### 8.2 必须支持的格式（逐条测试）

| format | 识别方式 | 说明 |
|---|---|---|
| `h5ad` | `.h5ad`/`.h5ad.gz` | 现有 anndata 路径（保留现有逻辑分支） |
| `10x_h5` | `.h5` 且 HDF5 里有 `/matrix` 组 | CellRanger v2/v3：`/matrix/{data,indices,indptr,shape}`，`barcodes`，`features`/`genes`（v2 用 `gene_names`），shape 为 genes×cells → 转置；`/matrix/features/feature_type` 存在时只保留 `Gene Expression` 行 |
| `loom` | `.loom` 且 HDF5 里有 `/matrix` dataset + `/row_attrs` | `/matrix` 是 genes×cells → 转置；基因名取 `/row_attrs/Gene`；`/col_attrs/*` 全部转成 obs 列 |
| `10x_mtx` | 目录含 `matrix.mtx(.gz)` + `barcodes.tsv(.gz)` + `features/genes.tsv(.gz)`；或单个 `.mtx`/`.mtx.gz`（旁边有 sidecar 就用） | 无 sidecar 时自动生成 `cell_1..`、`gene_1..` 并在 warnings 说明；`features.tsv` 第 2 列是基因名时用第 2 列 |
| `matrix_text` | `.csv/.tsv/.txt`（可 `.gz`） | 分隔符嗅探（`,`/`\t`/`;`/空白）；自动判方向：行名像基因（大写符号/无 barcode 特征）→ genes×cells 转置；行名像 barcode（`^[ACGT]{12,}` 或 `.*-\d+$`）→ cells×genes；不确定时用列数>行数等启发式并在 warnings 写明；**支持 `transpose` 参数强制** |
| `geo_series` | 文本含 `!series_matrix_table_begin`（通常在 `*_series_matrix.txt.gz`） | 取 table 段落，首列 `ID_REF` 为基因，其余列样本 → 转置成 cells×genes；附 `!Sample_*` 行作为 obs 列（同名去重）；warnings 说明"这是 GEO series matrix，本质是样本×基因汇总矩阵" |
| `rds` | `.rds`/`.RDS` | 走 R 桥（§8.3）；支持 Seurat v3/v4/v5、SingleCellExperiment、Matrix/dgCMatrix、data.frame |
| `rdata` | `.RData`/`.rda`（大小写不敏感） | 走 R 桥；多个对象时优先 Seurat/SCE，若仍不唯一 → `ambiguous` + candidates（含对象名/class/dim），可用 `object` 参数指定 |
| `archive` | `.zip/.tar/.tar.gz/.tgz/.tar.bz2` | 解压到 `<workdir>/unpack/<uuid>`（防路径穿越、跳过符号链接、解压后总大小上限 20GB）。解压后：① 找到 10x 三件套目录 → 按 `10x_mtx`；② 只有 1 个可识别文件 → 递归识别；③ 多个候选 → `ambiguous` + candidates（相对路径/大小/格式猜测）；`member` 参数可指定相对路径 |
| `gz` 单文件 | `*.h5ad.gz`/`*.csv.gz`/`*.mtx.gz`/`*.loom.gz` 等 | 解压到临时文件后按目标扩展名识别（`.mtx.gz` 交给 10x 分支） |

**URL 源**（`http(s)://`、`ftp://`）：下载到 `<workdir>/downloads/<hash>-<filename>`（文件名从 URL 或 `Content-Disposition` 取），带大小上限（默认 8GB，可 `options.max_bytes`），失败给可读错误；下载后可再走 archive/格式识别。支持 GEO 直链（`https://ftp.ncbi.nlm.nih.gov/geo/series/.../suppl/...`）。可选 `options.gse`：给 `GSE123456` 时用 `https://ftp.ncbi.nlm.nih.gov/geo/series/<GSE前缀>nnn/<GSE>/suppl/` 列表页解析出补充文件（`inspect_source` 返回 candidates，`load_source` 需 `member`）。

### 8.3 R 桥（`lib/r_bridge.R`）

- 本机 R：`/usr/local/bin/Rscript` 或 `/Library/Frameworks/R.framework/Resources/bin/Rscript`（4.4.2；已有 Seurat 5.4.0、SeuratObject 5.3.0、SingleCellExperiment 1.28.1、Matrix、SummarizedExperiment）。
- Python 侧：`subprocess.run(["Rscript","--vanilla",r_bridge,"--in",src,"--out",tmpdir,"--object",obj], timeout=600)`，捕获 stdout/stderr，非 0 退出 → `bad_request` + stderr 尾部。
- R 侧输出到 tmpdir：`matrix.mtx`（Matrix::writeMM，genes×cells 或 cells×genes 均可，但 meta.json 里写清 `"orientation"`）、`features.tsv`、`barcodes.tsv`、`obs.csv`（row.names=细胞名）、`meta.json`（`{class, assay, layers, dim, orientation, warnings, candidates?}`）。Python 侧据此构造返回 dict。
- 支持对象：Seurat（`SeuratObject::LayerData(obj,"counts")` → 回退 `GetAssayData(slot="counts")` → 回退 data 并在 warnings 说明）、SingleCellExperiment（`assay(sce,"counts")`/`assay(sce)` + `colData`）、matrix/Matrix（dimnames 必需，缺失则用生成名）、data.frame（全数值列，行=细胞）。
- 多 assay / 多 layer：默认取 `RNA`/第一个 assay 的 counts，其余写进 `meta.assays` 与 warnings。
- 中文/特殊字符：TSV 与 CSV 一律 `quote=FALSE`、UTF-8。

### 8.4 协议与 HTTP 扩展

**引擎 `open`** args 扩展：`{"path"|"url":str, "name"?, "format"?, "member"?, "object"?, "transpose"?, "max_bytes"?}`；`demo:true` 不变。返回的 `dataset` 对象新增：
```json
{"source_format":"10x_mtx","source":"/abs/x","notes":["..."],"warnings":["..."]}
```
**引擎新命令 `inspect`** args `{"path"|"url", ...同上}` → `inspect_source` 的结果 + `"ok":true`。

**HTTP**：
- `GET /singlecell/formats` → `{"formats":[{"id","title","extensions":[...],"available":bool,"note"}]}`
- `POST /singlecell/import/inspect` `{path|url,format?,member?,object?}` → inspect 结果
- `POST /singlecell/dataset/open` body 扩展 `{path?,url?,format?,member?,object?,transpose?}`
- 错误 `code:"ambiguous"` 时 body 带 `candidates`

**Agent 工具** `singlecell`：`open` 参数扩展 `url/data_format/member/object/transpose`；新增 action `inspect`。
注意：工具参数名不能叫 `format`。Host 把 `parameters` 原样投给模型侧，`format` 是 JSON Schema 保留关键字（值必须是 string），写了会让 DeepSeek API 直接返回 `Invalid schema for function 'singlecell'`。HTTP body 与 Python 引擎内部仍用 `format`，工具层做一次 `data_format → format` 映射。

### 8.5 sc-io 测试要求 `test/py_import_test.py` + `test/make_formats.py`

`test/make_formats.py` 在临时目录生成全套 fixture（用 numpy/anndata/h5py/zipfile/tarfile + Rscript）：
- h5ad、10x_h5、10x_mtx 目录、只给 `.mtx.gz`、loom、dense csv、tsv.gz、geo series matrix txt.gz、zip(10x 三件套)、tar.gz(csv)、Seurat rds、SCE rds、RData（含 2 个对象）、RDS 非 Seurat（matrix）、URL 用本地 `http.server` 临时起一个服务模拟下载。
- 每个 fixture 都写一份 `truth.json`（n_obs/n_vars/前几个基因名/首细胞某基因值）。
`test/py_import_test.py` 覆盖：每种格式 load 后 `truth.json` 一致；`inspect` 不加载矩阵且格式判断正确；`ambiguous` 分支（多对象 RData、多文件 zip）能拿到 candidates 且指定 member/object 后成功；坏文件/不存在路径给可读错误不崩；`.gz`/archive 路径穿越被拒；URL 下载（本地 http server）成功。全部通过退出码 0。
**真实数据验收**（只读，不改用户文件）：`/Users/dingyuntao/Desktop/deepseek harness/neu_analysis/objects/neu_clean.rds` 等 5 个 Seurat rds 必须能读，打印 n_obs/n_vars/obs 列数。

## 9. 高级分析层（v3，按用户给的能力图全量登记）

### 9.1 文件与归属

```
lib/sc_io.py                     sc-io
lib/r_bridge.R                   sc-io
lib/sc_analyses_sets.py          sc-analysis-a  # 基因集/通路/表达矩阵类
lib/sc_analyses_dynamics.py      sc-analysis-b  # 整合/轨迹/异质性/通讯/调控/敲除类
lib/analysis_registry.py         Lead           # 汇总两类模块 + 外部项登记 + ctx 构造 + 分发
test/make_formats.py             sc-io
test/py_import_test.py           sc-io
test/py_analyses_sets_test.py    sc-analysis-a
test/py_analyses_dynamics_test.py sc-analysis-b
```

### 9.2 模块 API（两个模块文件完全同构，互不依赖）

```python
ANALYSES = [{
  "id":"cell_cycle", "title":"细胞周期评分", "group":"基本分析",
  "method":"Seurat 式 S/G2M 打分（内置 Tirosh 基因集）",   # 必须写清实现方式（诚实标注 lite/代理）
  "params":[{"key":"...","label":"...","type":"number|string|select|boolean|column|file","default":...,"options":[...],"hint":"..."}],
  "needs":["groups:group"],   # 空=随时可跑；groups:<obs列>=需要某个分类列；spliced=需要剪接/未剪接层；gene_order=需要基因位置文件；gene_sets=需要基因集
}]
def run(analysis_id, ctx, params, outdir) -> dict   # 见 9.3
```

`ctx`（引擎构造，**只读**）：
```python
{"name":str, "n_obs":int, "n_vars":int,
 "X": csr_matrix float32, "X_kind":"counts"|"lognorm",   # 当前矩阵；已 normalize 则为 log 标准化
 "counts": csr_matrix float32,
 "obs": pandas.DataFrame,          # 索引=细胞名；含 clusters/cell_type/原有列
 "var_names":[str],
 "obsm":{"X_pca": ndarray|None, "X_emb": ndarray|None},
 "clusters":{"key":"clusters","labels":[str]}|None,
 "group_keys":[str],               # obs 里的分类列
 "engine_mode":"builtin"|"scanpy", "rng": numpy.random.Generator,
 "dataset_id":str}
```

`run` 返回（键固定）：
```python
{"id":..., "title":..., "method":...,
 "summary":{"bullets":[str],"metrics":{...}},
 "tables":[{"id","title","columns":[{"key","label","type":"string|number"}],"rows":[{...}],"note":str?}],
 "figures":[{"id","title","file":"/abs/*.png","width":int,"height":int}],
 "obs_columns":[{"key","kind":"categorical|numeric"}],   # 引擎按 obs.index 对齐合并
 "obs_values":{"key":[与 ctx.obs.index 等长的值列表]},   # 合并进 dataset.obs 的实际数值（JSON 可序列化，NaN→null）
 "warnings":[str]}
```
- 图：`outdir` 下写 PNG，白底 dpi≈140，**图上文字英文**，文件名 `<dataset>_<analysis>_<slug>.png`；返回绝对路径。
- 大数据：分析内部自行设 `max_cells`/`max_genes` 上限（默认 cells ≤ 20000、基因 ≤ 4000），超限随机下采样并在 warnings 写明。
- 确定性：同名参数结果一致（固定 rng/seed）。
- **缺失前置数据**（如 obs 里没有 group 列）→ 抛 `EngineError("prerequisite","...")`，不要编造假结果。

### 9.3 内置登记项（必须全部出现在 registry；`native` 必须可跑通，`external` 只登记+给出可行替代/安装指引）

**sc-analysis-a（`sc_analyses_sets.py`）**
| id | 标题 | 组 | 实现 | needs |
|---|---|---|---|---|
| `cell_cycle` | 细胞周期评分 | 基本分析 | 内置 S/G2M 基因集打分 + 相位 | — |
| `subcluster` | 亚型细分 | 基本分析 | 对指定簇子集重算 HVG/PCA/Louvain，写 obs 列 | clusters |
| `proportion_test` | 组间细胞类型比例差异 | 组间对比 | 各组占比 + 卡方/Dirichlet-multinomial 近似 + 反事实区间 | groups:group, clusters |
| `de_group` | 同细胞类型组间差异基因 | 组间对比 | 逐细胞类型 Welch t + BH，表+火山图 | groups:group, clusters |
| `volcano` | 差异基因火山图 | 数据可视化 | 读 de_group 结果或重新算 | groups:group, clusters |
| `proportion_sankey` | 细胞比例桑基图 | 数据可视化 | 两级 alluvial（组→细胞类型）matplotlib | groups:group, clusters |
| `featureplot` | 差异基因 Featureplot/Vlnplot/dotplot/heatmap | 数据可视化 | 复用现有嵌入/矩阵渲染逻辑 | clusters |
| `pathway_ora` | GO/KEGG 富集（ORA） | 通路分析 | 内置小 GMT + 用户 `.gmt` + 可选 Enrichr 在线（`source=enrichr`，需联网） | gene_sets |
| `gsea_prerank` | GSEA（预排序） | 通路分析 | 加权 KS running-sum + 置换 p | gene_sets |
| `gsva_like` | GSVA 通路活性 | 通路分析 | 每细胞基因集打分（ssGSEA 简化）写 obs | gene_sets |
| `auc_score` | 细胞基因集合分析（AUC） | 基因集合分析 | AUCell 式 rank-AUC 写 obs | gene_sets |
| `cnmf_programs` | cNMF 表达程序 | 基因集合分析 | sklearn NMF，k 可调，返回程序基因+细胞用法 | — |
| `wgcna_lite` | WGCNA 共表达模块 | 基因集合分析 | 相关矩阵层次聚类模块 + hub 基因 + 模块-簇相关 | — |
| `geneset_clinical` | 基因集分数联合临床因素 | 基因集合分析 | 分数 × obs 元数据关联（ANOVA/相关） | groups:group |

**sc-analysis-b（`sc_analyses_dynamics.py`）**
| id | 标题 | 组 | 实现 | needs |
|---|---|---|---|---|
| `integration` | 多样本整合 / 去除批次效应 | 基本分析 | ComBat 式位置-尺度校正 + 批次感知 kNN（BBKNN 式），返回校正后 PCA/嵌入、批次混合指标前后对比 | groups:batch |
| `trajectory_dpt` | 扩散伪时间（monocle-like） | 发育轨迹分析 | 扩散图 + DPT，自动/指定根细胞，分支+排序图 | clusters |
| `trajectory_paga` | 簇级轨迹图（PAGA-like） | 发育轨迹分析 | 簇间连通度图 + MST/ELG + 伪时间排序 | clusters |
| `entropy_stemness` | 基于熵的干性分析 | 发育轨迹分析 | 信号熵/表达熵打分写 obs | — |
| `cytotrace_lite` | CytoTRACE-like 分化潜能 | 发育轨迹分析 | 基因数-表达相关性平滑 + 近邻池化 | — |
| `cnv_score` | CNV 信号/ITH（inferCNV-like，需基因顺序） | 肿瘤异质性分析 | 参照细胞平滑 z 分数；有 `gene_order` 文件则按染色体排列并出区段图；无则不排序并标注"代理分数" | gene_order（可选） |
| `cnv_tree` | CNVscore/进化树 | 肿瘤异质性分析 | 基于 `cnv_score` 的簇级距离树（scipy hierarchy） | clusters |
| `ith_score` | ITH_score 肿瘤异质性 | 肿瘤异质性分析 | 簇内平均成对 CNV 距离 | clusters |
| `lr_communication` | 细胞通讯（CellChat/NATMI-lite） | 细胞通讯分析 | 内置人 LR 对表 + 表达乘积打分 + 置换 p + 通讯网络图 | clusters |
| `regulon_activity` | 转录调控（SCENIC-lite） | 转录调控分析 | 内置 TF 列表 + 共表达靶基因 + AUCell 式 regulon 活性写 obs | — |
| `virtual_knockout` | 虚拟敲除（scTenifoldKnk-lite） | 虚拟敲除 | HVG 相关/岭回归 GRN + KO 传播 + 受影响基因排序 | clusters |
| `virtual_overexpress` | 虚拟过表达 | 虚拟敲除 | 同上，方向取反 | clusters |
| `virtual_screen` | 虚拟筛选 | 虚拟敲除 | 对目标签名（如某簇 marker）逐一扰动排序 | clusters |

**外部项（`analysis_registry.py` 中登记，`native=False`，附 reason + 替代方案）**：`rna_velocity`（需 spliced/unspliced 层）、`monocle_exact`、`palantir_exact`、`nichenet`、`cellphonedb`、`cellchat_exact`、`pyscenic_exact`、`scenic_exact`、`aracne_ap`、`celloracle`、`regvelo`、`infercnv_exact`、`copykat`、`soupx_decontx`（污染校正）。每项：`{"id","title","group","native":false,"reason":"...","alternative":"本插件的 xxx 或安装参考 ..."}`。UI 显示为灰色"外部工具"，不得伪装可运行。

### 9.4 引擎/HTTP/工具扩展

- 引擎命令：`analyses {}` → `{"analyses":[{...,"available":bool,"reason":str?}]}`（读 registry + 当前数据集判断 needs）；
  `analysis {"dataset","id","params"}` → §9.2 结果 + `"obs_merged":[列名]`（引擎把 `obs_columns` 合并进 `dataset.obs`，供嵌入图着色）。
- HTTP：`GET /singlecell/analyses`、`POST /singlecell/analysis`。
- 工具：`singlecell` 新增 action `analyses`（列表+可用性）与 `analysis`（`id` + `params`）；`analysis` 的 message 给 bullets + 表行数/图路径，tables 只回前 20 行。
- 活动日志 action = `analysis:<id>`。
- `dataset.history` 里 `analysis` 不混入步骤历史（单独 `dataset.analyses_run` 记录 id/params/时间）。

### 9.5 UI（sc-ui）

- 顶栏「打开 h5ad」→ 「导入数据」：弹出导入面板（Dialog/抽屉）：
  - 本地路径输入（支持目录、`.h5ad/.h5/.loom/.mtx/.csv/.tsv/.txt` 及 `.gz`、`.zip/.tar.gz`、`.rds/.RData`）
  - URL 输入（GEO FTP/HTTP 直链）
  - 文件选择（accept 全格式）
  - 格式下拉：自动（默认）/ 各格式（来自 `GET /singlecell/formats`）
  - `ambiguous` 响应 → 展示 candidates 列表（文件名/大小/对象名/维度），点选后带 `member`/`object` 重试
  - 可选高级项：`transpose`（矩阵方向）、`object`（RData 对象名）
  - 忙碌计时 + 可读错误
- 中栏新增 tab「高级分析」：左列按 7 大类分组列出全部条目（native 可点，external 灰色 + reason/替代方案 tooltip）；点 native 条目 → 参数表单（按 `params` 的 type 生成：number/select/boolean/column/file/text）→ 运行 `POST /singlecell/analysis`；结果区：bullets、tables（表格，行数多时滚动）、figures（`<img>`）；新 obs 列出现后嵌入图 color-by 自动可选（下次刷新 state 即出现）。
- 分析面板不做自动跑；运行中禁用按钮 + 「运行中…」；失败行内红字。
- 现有流程卡/嵌入图/记录/Agent 面板保持不变；纯函数层补充 `ANALYSIS_GROUPS`、`groupAnalyses()`、`availableParams()`、`formatCandidates()` 等并加测试。


## 10. 官方工具桥 v4（Seurat 6.0 降维 + 官方 R/Python 工具 + 虚拟敲除）

用户机器实测已具备：R 4.4.2（CellChat 1.6.1、monocle3 1.3.7、monocle 2.34、nichenetr 2.2.1、SCENIC 1.3.1+GENIE3+RcisTarget+AUCell、GSVA 2.0.7、clusterProfiler 4.14.6+org.Hs.eg.db、fgsea、WGCNA、harmony、batchelor、slingshot）；conda base Python 3.13.11（pyscenic 0.12.1/arboreto/ctxcore/loompy/umap-learn）；SCENIC motif 库与 NicheNet 先验模型在用户目录。
要求：**官方实现全部接上；降维改用 Seurat 6.0；虚拟敲除的 Python 包由本插件负责安装**。本插件的 lite 实现保留并与官方版并存（目录里区分标注）。

### 10.1 文件与归属

```
lib/sc_external.py               Lead   工具注册/探测/bundle 导出/子进程分发
lib/external/seurat_pipeline.R   sc-seurat6  Seurat 6.0 降维
lib/external/cellchat.R          sc-ext-r   官方 CellChat
lib/external/monocle3.R          sc-ext-r   官方 monocle3
lib/external/nichenet.R          sc-ext-r   官方 NicheNet
lib/external/scenic_r.R          sc-ext-r   官方 SCENIC(R)
lib/external/clusterprofiler_ora.R sc-ext-r 官方 GO/KEGG
lib/external/fgsea.R             sc-ext-r   官方 fgsea
lib/external/gsva.R              sc-ext-r   官方 GSVA
lib/external/wgcna.R             sc-ext-r   官方 WGCNA
lib/external/aucell.R            sc-ext-r   官方 AUCell
lib/external/harmony.R           sc-ext-r   官方 harmony
lib/external/slingshot.R         sc-ext-r   官方 slingshot
lib/external/{soupx,copykat,infercnv}.R  sc-ext-r（装得上才交付）
lib/external/pyscenic.py         sc-ext-py  官方 pySCENIC
lib/external/cellphonedb.py      sc-ext-py  官方 CellPhoneDB
lib/external/scvelo.py           sc-ext-py  官方 scVelo（需 spliced/unspliced）
lib/external/palantir.py         sc-ext-py  官方 Palantir
lib/external/celloracle.py       sc-ext-py  官方 CellOracle（虚拟敲除）
lib/external/regvelo.py          sc-ext-py  官方 RegVelo（装不上则登记原因）
test/py_seurat6_test.py / py_external_r_test.py / py_external_py_test.py
```

### 10.2 bundle 格式（Lead 导出，所有适配器只读）

目录 `<workdir>/ext/<dataset>-<tool>/`：
- `counts.mtx` MatrixMarket，**genes × cells**（原始 counts）
- `data.mtx` 同形状，当前矩阵（已 normalize 则为 log-normalized，否则同 counts）
- `features.tsv` 行=基因名；`barcodes.tsv` 行=细胞名
- `obs.csv`：第一列细胞名 + 元数据列（值内逗号/换行替换为空格，UTF-8，无引号）
- `embed.csv`（可选）：`cell,dim1,dim2`（X_emb）
- `pca.csv`（可选）：`cell,PC1..PCn`
- `meta.json`：`{dataset, name, n_obs, n_vars, steps, X_kind, engine}`

R 读法：`m <- Matrix::readMM("counts.mtx"); m <- t(m)`（得到 cells×genes）。Python 读法：`scipy.io.mmread(...).T`。

### 10.3 适配器契约

- R：`Rscript --vanilla lib/external/<tool>.R --bundle <dir> --out <dir> --params '<json>'`
- Python：`<externalPython> lib/external/<tool>.py --bundle <dir> --out <dir> --params '<json>'`
- 退出码：`0` 成功；`3` 缺前置/数据不满足（stderr 给中文原因 → 引擎转 `prerequisite`）；其它非 0 = 失败（stderr 前 2KB 作为错误消息）。
- stdout 必须干净（日志走 stderr）。
- 必须写 `<out>/result.json`：
```json
{"title":"CellChat 细胞通讯","method":"CellChat 1.6.1（官方 R 包）",
 "summary":{"bullets":["..."],"metrics":{"n_significant":12}},
 "tables":[{"id":"interactions","title":"...","columns":[{"key":"source","label":"Source","type":"string"}],"rows":[{"source":"..."}]}],
 "figures":[{"id":"network","title":"...","file":"<out>/network.png"}],
 "obs_file":"<out>/obs_extra.csv",
 "warnings":[],"notes":[]}
```
- `obs_extra.csv`（可选）：第一列细胞名（与 barcodes 对齐）+ 新增列（数值或字符串）。
- 图：PNG，白底，dpi≥140，图上文字英文。
- 规模上限：cells ≤ 15000、genes ≤ 5000（超限下采样并写 warnings；CellChat/NicheNet/SCENIC 等更重的工具自行收紧到 ≤3000–5000 并在 notes 说明）。

### 10.4 注册、探测与配置

`lib/sc_external.py` 的 `TOOLS` 每项：`{"id","title","group","lang":"r|py","script","needs":[…],"params":[…],"probe":{"r":[…]} 或 {"py":[…]}, "method"}`。
- 探测：R `Rscript -e 'for (p in c(...)) cat(p, requireNamespace(p, quietly=TRUE), "\n")'`；Python `<externalPython> -c "import ..."`。结果缓存 5 分钟。
- 配置键（`cordis.patch.yml` 的 `config`）：
  - `externalPython`：默认探测 `/opt/miniconda3/bin/python3`（conda base）→ `python3`
  - `Rscript`：默认 `/usr/local/bin/Rscript` → 框架版 → PATH
  - `rLibs`：额外 R 库路径（Seurat 6 独立库 `~/.dsh/sc-plugin-rlib` 会由 sc-seurat6 创建）
  - `nichenetPriors`：默认 `/Users/dingyuntao/Desktop/祖传R语言单细胞代码/crosstalk`
  - `scenicMotifDb`：默认 `/Users/dingyuntao/Desktop/单细胞工作站/SCINCE尝试/human/python用`
  - 以上**全部"探测 + 配置覆盖"**：探测不到只让对应工具 `available:false` 并给 reason，不得崩。
- 引擎命令：`external_tools {}` → `{"tools":[…,"available":bool,"reason":str?]}`；`external_run {dataset,id,params}` → `{…result.json, obs_merged:[…], ms}`。
- `analyses` 目录合并 external：可用 external 项 `available:true, external:true`，Host 层把可用 external 的 `native` 置 true（UI 据此可运行），同时保留 `external:true` 供显示"官方实现"。不可用项 `available:false` + reason。

### 10.5 Seurat 6.0 降维（sc-seurat6）

- 安装到**独立 R 库** `~/.dsh/sc-plugin-rlib`（绝不覆盖系统 Seurat 5.4）；适配器脚本运行时 `.libPaths(c(<独立库>, .libPaths()))`。
- `seurat_pipeline.R` CLI 同上；`params.step` ∈ `pca|neighbors|cluster|embed|full`；官方 API：NormalizeData(LogNormalize) → FindVariableFeatures(vst, n=2000) → ScaleData → RunPCA(npcs) → FindNeighbors(dims) → FindClusters(resolution) → RunUMAP(dims)。
- 输出 `obs_extra.csv`（`seurat_clusters` 等）、`pca.csv`、`embed.csv` + `result.json`（含 Seurat 版本 metrics）。
- 引擎接入：`run` 步骤 `pca/neighbors/cluster/embed` 支持 `params.engine = "seurat"|"seurat6"`：走 `sc_external.seurat_step(dataset, step, params)`，把 clusters 写回 `obs["clusters"]`、PCA/UMAP 写回 `obsm["X_pca"]/["X_emb"]`（后续 markers/annotate/UI 全部可用）。
- 预设流程新增 `seurat6`：qc→filter→normalize→hvg→pca(engine=seurat)→neighbors(engine=seurat)→cluster(engine=seurat)→embed(engine=seurat)。
- UI：流程卡参数里 `engine`（自动/内置/Seurat 6）；预设下拉自动出现 seurat6。Seurat 6 不可用时该选项 disabled + reason。


---

## §12 Seurat 6 + Harmony 批次整合（v0.5.0，对齐「单细胞分析器」）

流程步骤可在 `params.engine="seurat"` 时走 Seurat 6.0 官方实现（`lib/external/seurat_pipeline.R`）：

- **默认参数对齐用户软件**：`resolution=0.2`、`n_top_genes=3000`、`n_pcs=30`、`n_neighbors(UMAP/kNN)=30`、`use_harmony=true`
- **Harmony 批次整合**：PCA 后按批次列（自动识别 `sample/orig.ident/group/batch`，可用 `params.batch_key` 指定）运行官方 harmony；Seurat v5/v6 上 `RunHarmony` 报 `non-conformable arguments` 时自动回退 `harmony::HarmonyMatrix` 矩阵版并挂回 DimReduc；批次列只有 1 个取值或无批次列时跳过并给 warning
- **元数据透传**：bundle 的 obs（sample/group/…）会附加到 Seurat 对象 meta.data，否则找不到批次列
- **聚类数提示**：簇数 >30 提示调低 resolution、<10 提示调高（参考 20–30）
- **不静默降采样**：`max_cells` 默认 200000（原 25000），导出 bundle 的细胞上限由参数控制
- **预设流程 `seurat6`** 已改为：QC→filter→normalize→hvg(3000)→Seurat PCA→Harmony→FindClusters(resolution=0.2)→UMAP(n.neighbors=30)
- 工具 schema 参数 `group_regex`（原 `regex`/`pattern` 均为 JSON Schema 保留字雷区，`TOOL_SCHEMA_VERSION=4`）


---

## §13 UI v6 与配套后端（v0.6.0）

### 13.1 按选中基因出图

`POST /plot`，`kind=dotplot|heatmap`，`params` 支持：
- `genes: ["CD3D","NKG7",...]`（也接受逗号/空格分隔字符串）→ 按选中基因 × 簇出图；不在表达矩阵中的基因忽略并写 warning；最多 120 个；
- 未传 `genes` 时保持原行为：每簇 top-N marker（`params.top_n`，默认 5）。

引擎实现：`_selected_gene_list` / `_expression_by_cluster`（均值 + 表达比例矩阵）。

### 13.2 按样本移除细胞

`POST /samples/remove {dataset, samples:[样本名或id]}` →
`{ok, dataset: 摘要, removed_cells, removed_samples}`。

- 真正从 `X/counts/obs/obsm/layers` 中剔除这些样本的细胞；
- 重置下游结果：清空 clusters/cell_type/_emb_*/obsm、markers/adj/knn/layout/邻接图等；
- `steps` 保留数据级前缀 `qc/filter/normalize/hvg`；`samples/groups` 同步更新；
- 全部样本都移除时返回 `bad_request`；
- Agent 工具同样可用：`action=samples, op=remove, samples=[...]`。

### 13.3 UI v6 入口

| 功能 | 位置 |
|---|---|
| 降维聚类控制卡（Seurat6/Harmony/resolution 0.2/UMAP 30） | 聚类 tab 顶部 |
| 一键降维聚类 / 只重聚类 | 同上 |
| 手动注释工作台（簇列表 / marker 表 / DotPlot / 热图 / FeaturePlot / 保存注释 / 粘贴 CSV） | 标记基因 tab 顶部 |
| 删除分组 / 清空分组 / 移除样本 | 记录 tab「样本 & 分组」 |
| Agent 输入框 + 建议卡 + 直接运行 | 右栏 Agent 顶部 |


---

## §14 对齐 pro.R 的 IntegrateLayers 整合（v0.7.0）

- `params.integration`：`"harmony" | "rpca" | "none"`（默认 harmony；`use_harmony=false` 等价 none）
- 使用 Seurat 官方 `IntegrateLayers(method = HarmonyIntegration / RPCAIntegration, orig.reduction="pca", new.reduction=...)`；单层对象在整合前按批次列 `split()` 并补齐 data/scale 层
- RPCA：`k.weight` 自适应并在失败时回退 Harmony（结果 metrics 标 `integration_used` 与真实 `reduction`）
- 默认参数：`dims=1:30`、`n_top_genes=3000`、`metric=euclidean`、resolution 默认 none→0.5 / 整合→0.8
- 预设：`seurat6`（Harmony 0.8）、`pro-rpca`（RPCA 0.8）、`pro-unintegrated`（PCA 0.5），均带 pro.R QC 参数
- 簇列：当前结果写 `clusters`；整合结果额外保留 `clusters_harmony` / `clusters_rpca`
- 修复：按样本 merge 构建多层对象会打乱细胞顺序 → 曾导致簇标签错位（ARI≈0）；已改回单对象 + split 路径


---

## §17 UI v8：Step3 三段式注释 + HVG 降权（v0.9.0）

Step3 按「先方案 → 再出图 → 最后分簇」重排：
1. 标记方案卡：细胞类型名 + marker 基因可编辑表格；内置人/鼠各 18 类参考方案；粘贴 CSV；校验空/重复/无基因/超限；基因大小写不敏感去重。
2. 看表达卡：[出 DotPlot]/[出热图] → `/plot {params:{genes, gene_groups, groupby:"clusters"}}`；[画 FeaturePlot] → `/embedding {color_by:{kind:"gene",key}}`。
3. cluster 分配卡：每簇下拉/手输 → [全部保存注释] → `/run annotate mapping`；注释表/覆盖率/按注释着色；保留小亚群提取。

Step2：`n_top_genes` 移入「进阶/参考参数」折叠区（默认 3000，对齐 pro.R，仅参考）；`defaultDimredParams()/buildDimredPayload()` 保证折叠不丢值。

### 17.1 gene_groups 出图（后端）

`POST /singlecell/plot`（dotplot/heatmap）`params.gene_groups = [{"label":"T cell","genes":[...]}, ...]`：按类型展平、组间分隔线 + 类型名；缺失基因忽略并 warning；同时给 `genes` 时以 gene_groups 为准。
