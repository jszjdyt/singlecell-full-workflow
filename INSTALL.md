# dsh-singlecell-workbench 安装说明

单细胞分析工作台：DeepSeek Harness 插件。导入 GEO 常见单细胞数据（h5ad / 10x h5 / 10x mtx / loom / 表达矩阵 / GEO series matrix / `.rds` / `.RData` / 压缩包 / URL），
可视化跑 QC→聚类→嵌入→标记基因，并内置一整套高级分析（通路/基因集/轨迹/异质性/通讯/调控/虚拟敲除）；Agent 与页面共用同一套操作。

## 0. 你会收到什么

一个压缩包 `dsh-singlecell-workbench-0.3.0.zip`。解压后是一个 `dsh-singlecell-workbench/` 目录，**这个目录就是插件本体**。

## 1. 环境要求

| 项 | 要求 | 说明 |
|---|---|---|
| DeepSeek Harness | 已安装并能正常使用 | 插件走 DSH 的 plugin_manager 安装 |
| Python | 3.10+ | 必装，见第 3 步 |
| R | 4.x + `Seurat`/`SeuratObject`/`SingleCellExperiment` | **可选**，只有读 `.rds`/`.RData` 才需要 |
| 网络 | 安装依赖时需要 | 运行分析不需要联网（Enrichr 在线富集是可选功能） |

## 2. 放好插件目录

把解压出来的 `dsh-singlecell-workbench/` 放到一个**不会被自动清理**的固定位置，例如：

- macOS / Linux：`~/dsh-plugins/dsh-singlecell-workbench`
- Windows：`C:\Users\你\dsh-plugins\dsh-singlecell-workbench`

> 注意：安装是「链接目录」方式，**这个目录以后不能删/不能挪**，否则插件失效。

## 3. 安装 Python 依赖

```bash
cd ~/dsh-plugins/dsh-singlecell-workbench
python3 -m pip install -r requirements.txt
```

Windows（如果 `python3` 不在 PATH）：

```bat
cd C:\Users\你\dsh-plugins\dsh-singlecell-workbench
python -m pip install -r requirements.txt
```

自检（应输出一行版本号）：

```bash
python3 -c "import numpy,pandas,scipy,sklearn,anndata,h5py,matplotlib;print('deps ok')"
```

插件会自动按顺序找解释器：DSH 配置里的 `python` → 环境变量 `DSH_SINGLECELL_PYTHON` → PATH 上的 `python3`/`python`/`py` → macOS 框架版 / Homebrew / 系统 Python。
如果你有多个 Python，建议在 DSH 配置里显式指定（见第 7 节）。

## 4. 在 DSH 里安装插件

打开 DeepSeek Harness，在任意会话里对 Agent 说：

> 帮我安装本地插件：plugin_manager install_bundle target=/你的绝对路径/dsh-singlecell-workbench

（Windows 路径示例：`C:\Users\你\dsh-plugins\dsh-singlecell-workbench`）

Agent 会调用 `plugin_manager` 并请求你批准；看到 `application: applied` 即安装成功。
命令行也可以（在 DSH 已启动的环境里由 Agent 代跑，不建议手改 profile 文件）。

## 5. 验证安装

1. **侧边栏**出现「单细胞」图标（点开是工作台整页）。
2. 顶栏点 **「载入示例数据」** → 左栏点几张流程卡（或顶栏「预设流程 → 标准流程」一键跑）→ 中栏「嵌入图」能看到彩色散点。
3. **Agent 工具**：新开/刷新会话后，让 Agent「列出可用工具」应看到 `singlecell`；或者直接说
   「用单细胞工作台载入示例数据并跑到聚类」。
4. 插件管理页能看到「单细胞分析工作台」及其图标、描述。

> 如果侧边栏没出现：刷新页面（Cmd+R / Ctrl+R）；仍没有就重启 DSH。

## 6. 可选：读 `.rds` / `.RData`

装 R 与 Seurat 生态：

```r
install.packages(c("Seurat", "SeuratObject", "SingleCellExperiment", "SummarizedExperiment", "Matrix"))
```

插件自动探测 Rscript：macOS `/usr/local/bin/Rscript`、R 框架版；Windows `C:\Program Files\R\R-*\bin\Rscript.exe`；其它情况找 PATH 上的 `Rscript`。
不装 R 也能用——只是打开 `.rds/.RData` 时会提示「没有找到 Rscript」，其它格式不受影响。

## 7. 可选配置

在 profile 的 `cordis.patch.yml` 里（或请 Agent 帮你改）：

```yaml
- id: dsh-singlecell-workbench
  name: dsh-singlecell-workbench
  config:
    python: /绝对路径/python3      # 指定解释器
    workdir: /tmp/sc-workdir       # 导入解包/图/导出目录（默认系统临时目录）
    timeoutMs: 300000              # 单请求超时
```

## 8. 数据与隐私

- 所有计算在本机进行：Python/R 都是本机子进程，数据不出本机。
- 页面只通过**本机回环地址**的 HTTP 路由与插件通信（非回环请求会被拒绝）。
- 只有主动使用「Enrichr 在线富集」或「URL 导入」时才会联网。

## 9. 常见问题

| 现象 | 处理 |
|---|---|
| 顶栏徽章「Python 环境未就绪」 | 按第 3 步装依赖；或看提示里缺哪个包；或在配置里指定 python |
| 打开 `.rds` 提示没有 Rscript | 按第 6 步装 R 与 Seurat，或先用 R 转成 `.h5ad` |
| 侧边栏没有「单细胞」 | 确认 `application: applied`；刷新页面/重启 DSH；看插件管理页是否启用 |
| Agent 工具里没有 `singlecell` | 刷新会话；插件未启用时工具不注册 |
| 导入大文件很慢/超时 | 改用本地路径导入（比网页上传快）；或调大 `timeoutMs` |
| 图显示不出来 | 图保存在临时目录，被系统清理后需重新生成 |
| 高级分析显示「不可用」 | 看原因：通常需要先聚类、或缺分组列、或数据缺基因位置文件 |

## 10. 更新与卸载

- **更新**：用新压缩包覆盖 `dsh-singlecell-workbench/` 目录 → 重启 DSH（或让 Agent 先 `remove_bundle dsh-singlecell-workbench` 再 `install_bundle`）。
- **卸载**：让 Agent 执行 `plugin_manager remove_bundle target=dsh-singlecell-workbench`，然后删除目录即可。卸载不会删除你导出的分析结果文件。

## 11. 限制（先说清楚）

内置分析均为本插件自行实现（numpy/scipy/sklearn/R 桥），`method` 字段会标注 "lite/近似/代理"，
**不等同于 pySCENIC/CellChat/NicheNet/inferCNV/CellOracle 等官方工具的结果**；需要这些官方工具的条目在分析目录里只登记并给出替代方案。
自动注释是基于内置参考 marker 面板的建议，自检是规则检查，请结合专业知识判断后再用于结论/发表。
