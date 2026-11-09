#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell — 官方工具桥（CONTRACT §10）

把用户机器上已安装 / 本插件代装的官方工具接进工作台：
  - R 侧：Seurat 6.0 降维、CellChat、monocle3、NicheNet、SCENIC、clusterProfiler、
          fgsea、GSVA、WGCNA、AUCell、harmony、slingshot（+可选 SoupX/copykat/inferCNV）
  - Python 侧：pySCENIC、CellPhoneDB、scVelo、Palantir、CellOracle、RegVelo

职责：
  1. 探测外部解释器与包（结果缓存 5 分钟），给出每个工具的 available/reason；
  2. 把当前数据集导出成适配器约定的 bundle（§10.2）；
  3. 以子进程方式调用 lib/external/<tool>.{R,py}，解析 result.json / obs_extra.csv；
  4. 提供 Seurat 6.0 分步（pca/neighbors/cluster/embed/full）给引擎的流程步骤用。

不做的事：不硬编码用户机器上的绝对路径到"没有就崩"——全部是"探测 + 配置覆盖"。
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import scipy.sparse as sp
    from scipy.io import mmwrite
except Exception:  # pragma: no cover
    sp = None
    mmwrite = None

LIB_DIR = Path(__file__).resolve().parent
EXTERNAL_DIR = LIB_DIR / "external"

PROBE_TTL = 300
_SEURAT_MIN = (6, 0)


class ExternalError(Exception):
    def __init__(self, code, message):
        self.code = code
        self.message = str(message)
        super().__init__(self.message)


# ---------------------------------------------------------------------------
# 配置（可通过 cordis config → DS_SC_CONFIG 环境变量覆盖）
# ---------------------------------------------------------------------------


def _config():
    raw = os.environ.get("DS_SC_CONFIG") or "{}"
    try:
        value = json.loads(raw)
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def _first_existing(paths):
    for item in paths:
        if not item:
            continue
        path = Path(os.path.expanduser(str(item)))
        if path.exists():
            return str(path)
    return None


def rscript_path():
    cfg = _config()
    explicit = cfg.get("Rscript")
    if explicit:
        return str(explicit)
    return _first_existing([
        "/usr/local/bin/Rscript",
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
        shutil.which("Rscript"),
        r"C:\\Program Files\\R\\R-*/bin/Rscript.exe",
    ]) or shutil.which("Rscript")


def external_python():
    cfg = _config()
    for key in ("externalPython", "external_python"):
        if cfg.get(key):
            return str(cfg[key])
    env = os.environ.get("DS_SC_EXTERNAL_PYTHON")
    if env:
        return env
    return _first_existing([
        "/opt/miniconda3/bin/python3",
        os.path.expanduser("~/miniconda3/bin/python3"),
        os.path.expanduser("~/anaconda3/bin/python3"),
        shutil.which("python3"),
    ]) or shutil.which("python3")


def _seurat_lib():
    cfg = _config()
    explicit = cfg.get("seuratLib") or cfg.get("rLibs")
    candidates = [explicit, os.path.expanduser("~/.dsh/sc-plugin-rlib")]
    for item in candidates:
        if item and Path(os.path.expanduser(str(item))).exists():
            return str(Path(os.path.expanduser(str(item))))
    return os.path.expanduser("~/.dsh/sc-plugin-rlib")


def _nichenet_priors():
    cfg = _config()
    if cfg.get("nichenetPriors"):
        return os.path.expanduser(str(cfg["nichenetPriors"]))
    return _first_existing([
        os.path.expanduser("~/nichenet_priors"),
    ])


def _scenic_db():
    cfg = _config()
    if cfg.get("scenicMotifDb"):
        return os.path.expanduser(str(cfg["scenicMotifDb"]))
    return _first_existing([
        os.path.expanduser("~/scenic_databases"),
        os.path.expanduser("~/cisTarget_databases"),
    ])


# ---------------------------------------------------------------------------
# 工具登记表
# ---------------------------------------------------------------------------


def _p(key, label, ptype="string", default=None, options=None, hint=None):
    item = {"key": key, "label": label, "type": ptype}
    if default is not None:
        item["default"] = default
    if options:
        item["options"] = options
    if hint:
        item["hint"] = hint
    return item


TOOLS = [
    # ---------------- R：官方实现 ----------------
    {
        "id": "ext_seurat6", "title": "Seurat 6.0 官方降维（全流程）", "group": "基本分析", "lang": "r",
        "script": "seurat_pipeline.R", "needs": [], "timeout": 900,
        "params": [_p("step", "步骤", "select", "full", ["full", "pca", "cluster", "embed"]),
                   _p("n_pcs", "主成分数", "number", 30), _p("resolution", "聚类分辨率", "number", 0.8),
                   _p("n_top_genes", "高变基因数", "number", 2000), _p("dims", "UMAP/HVG 维度", "number", 20)],
        "probe": {"r": ["Seurat"], "seurat_version": True},
        "method": "Seurat 6.0 官方流程（NormalizeData/FindVariableFeatures/ScaleData/RunPCA/FindNeighbors/FindClusters/RunUMAP）",
    },
    {
        "id": "ext_cellchat", "title": "CellChat 官方细胞通讯", "group": "细胞通讯分析", "lang": "r",
        "script": "cellchat.R", "needs": ["clusters"], "timeout": 1800,
        "export": {"max_cells": 4000, "max_genes": 15000},
        "params": [_p("min_cells", "每类最小细胞数", "number", 10), _p("thresh", "通讯 p 阈值", "number", 0.05)],
        "probe": {"r": ["CellChat"]},
        "method": "CellChat 1.6.1 官方实现（CellChatDB.human + computeCommunProb + aggregateNet）",
    },
    {
        "id": "ext_monocle3", "title": "monocle3 官方轨迹", "group": "发育轨迹分析", "lang": "r",
        "script": "monocle3.R", "needs": ["clusters"], "timeout": 1800,
        "params": [_p("resolution", "聚类分辨率", "number", 0.8), _p("n_top_genes", "高变基因数", "number", 2000),
                   _p("root_cell", "根细胞名（可选）")],
        "probe": {"r": ["monocle3"]},
        "method": "monocle3 1.3.7 官方实现（preprocess_cds/reduce_dimension/cluster_cells/learn_graph/order_cells）",
    },
    {
        "id": "ext_nichenet", "title": "NicheNet 官方配体-靶基因", "group": "细胞通讯分析", "lang": "r",
        "script": "nichenet.R", "needs": ["clusters"], "timeout": 3600,
        "export": {"max_cells": 3000, "max_genes": 15000},
        "params": [_p("cell_type_col", "细胞类型列", "column", "cell_type"),
                   _p("sender", "发送方细胞类型（留空=自动选一个）"),
                   _p("receiver", "接收方细胞类型（留空=自动选一个）",),
                   _p("top_n", "top 配体数", "number", 20)],
        "probe": {"r": ["nichenetr"], "nichenet_priors": True},
        "method": "NicheNet 2.2 官方实现（nichenetr + 用户本地先验模型）",
    },
    {
        "id": "ext_scenic_r", "title": "SCENIC(R) 官方转录调控", "group": "转录调控分析", "lang": "r",
        "script": "scenic_r.R", "needs": [], "timeout": 3600,
        "export": {"max_cells": 3000, "max_genes": 4000},
        "params": [_p("species", "物种", "select", "human", ["human", "mouse"]),
                   _p("n_top_genes", "高变基因数", "number", 2000),
                   _p("db_dir", "motif 数据库目录（留空=自动探测）")],
        "probe": {"r": ["SCENIC", "GENIE3", "RcisTarget", "AUCell"], "scenic_db": True},
        "method": "SCENIC 1.3.1 官方 R 实现（GENIE3 + RcisTarget motif 富集 + AUCell）",
    },
    {
        "id": "ext_clusterprofiler", "title": "clusterProfiler 官方 GO/KEGG 富集", "group": "通路分析", "lang": "r",
        "script": "clusterprofiler_ora.R", "needs": [], "timeout": 900,
        "params": [_p("source", "数据库", "select", "both", ["go", "kegg", "both"]),
                   _p("species", "物种", "select", "human", ["human", "mouse"]),
                   _p("gene_list", "基因列表（留空=用各簇 top markers）", "text", ""),
                   _p("top_n", "每簇取 top 基因", "number", 50)],
        "probe": {"r": ["clusterProfiler", "org.Hs.eg.db"]},
        "method": "clusterProfiler 4.14.6 官方实现（enrichGO/enrichKEGG + org.Hs.eg.db）",
    },
    {
        "id": "ext_fgsea", "title": "fgsea 官方 GSEA", "group": "通路分析", "lang": "r",
        "script": "fgsea.R", "needs": [], "timeout": 900,
        "params": [_p("ranking", "排序依据", "select", "pca1", ["pca1", "dispersion", "custom"]),
                   _p("gmt_file", "GMT 文件（留空=内置基因集）", "file", "")],
        "probe": {"r": ["fgsea"]},
        "method": "fgsea 1.32.4 官方实现（fgseaMultilevel）",
    },
    {
        "id": "ext_gsva", "title": "GSVA 官方通路活性", "group": "通路分析", "lang": "r",
        "script": "gsva.R", "needs": [], "timeout": 1200,
        "params": [_p("method", "算法", "select", "gsva", ["gsva", "ssgsea"]),
                   _p("gmt_file", "GMT 文件（留空=内置基因集）", "file", ""),
                   _p("n_sets", "通路数", "number", 12)],
        "probe": {"r": ["GSVA"]},
        "method": "GSVA 2.0.7 官方实现（gsva/ssgsea）",
    },
    {
        "id": "ext_wgcna", "title": "WGCNA 官方共表达模块", "group": "基因集合分析", "lang": "r",
        "script": "wgcna.R", "needs": [], "timeout": 1800,
        "params": [_p("n_top_genes", "高变基因数", "number", 800), _p("min_module_size", "最小模块大小", "number", 20)],
        "probe": {"r": ["WGCNA"]},
        "method": "WGCNA 1.74 官方实现（软阈值 + 模块识别 + 枢纽基因）",
    },
    {
        "id": "ext_aucell", "title": "AUCell 官方基因集打分", "group": "基因集合分析", "lang": "r",
        "script": "aucell.R", "needs": [], "timeout": 900,
        "params": [_p("gmt_file", "GMT 文件（留空=内置基因集）", "file", ""), _p("n_sets", "基因集数", "number", 12)],
        "probe": {"r": ["AUCell"]},
        "method": "AUCell 1.28.0 官方实现",
    },
    {
        "id": "ext_harmony", "title": "Harmony 官方批次校正", "group": "基本分析", "lang": "r",
        "script": "harmony.R", "needs": ["groups:batch"], "timeout": 900,
        "params": [_p("batch", "批次列", "column", "batch"), _p("dims", "使用维度", "number", 20)],
        "probe": {"r": ["harmony"]},
        "method": "harmony 1.2.4 官方实现",
    },
    {
        "id": "ext_slingshot", "title": "slingshot 官方谱系轨迹", "group": "发育轨迹分析", "lang": "r",
        "script": "slingshot.R", "needs": ["clusters"], "timeout": 1200,
        "params": [_p("cluster_col", "聚类列", "column", "clusters")],
        "probe": {"r": ["slingshot"]},
        "method": "slingshot 2.14.0 官方实现（谱系推断 + 平滑伪时间）",
    },
    {
        "id": "ext_soupx", "title": "SoupX 官方污染校正", "group": "基本分析", "lang": "r",
        "script": "soupx.R", "needs": [], "timeout": 1200,
        "params": [_p("contamination", "污染比例（留空=自动估计）", "number", 0.1)],
        "probe": {"r": ["SoupX"]},
        "method": "SoupX 官方实现（需要 raw/empty droplet 矩阵时请提供）",
    },
    {
        "id": "ext_copykat", "title": "CopyKAT 官方 CNV 推断", "group": "肿瘤异质性分析", "lang": "r",
        "script": "copykat.R", "needs": ["clusters"], "timeout": 2400,
        "export": {"max_cells": 5000, "max_genes": 15000},
        "params": [_p("n_cells", "抽样细胞数", "number", 3000)],
        "probe": {"r": ["copykat"]},
        "method": "CopyKAT 官方实现（基础 CNV + 亚克隆）",
    },
    {
        "id": "ext_infercnv", "title": "inferCNV 官方实现", "group": "肿瘤异质性分析", "lang": "r",
        "script": "infercnv.R", "needs": ["clusters"], "timeout": 2400,
        "export": {"max_cells": 5000, "max_genes": 15000},
        "params": [_p("reference_group", "参照细胞类型（可选）")],
        "probe": {"r": ["infercnv"], "jags": True},
        "method": "inferCNV 官方实现（需 JAGS）",
    },
    # ---------------- Python：官方实现 ----------------
    {
        "id": "ext_pyscenic", "title": "pySCENIC 官方转录调控", "group": "转录调控分析", "lang": "py",
        "script": "pyscenic.py", "needs": [], "timeout": 3600,
        "export": {"max_cells": 3000, "max_genes": 4000},
        "params": [_p("n_workers", "并行数", "number", 4), _p("n_top_genes", "高变基因数", "number", 2000),
                   _p("db_dir", "motif 数据库目录（留空=自动探测）")],
        "probe": {"py": ["pyscenic", "ctxcore", "arboreto"], "scenic_db": True},
        "method": "pySCENIC 0.12.1 官方实现（GRNBoost2 + ctx + AUCell）",
    },
    {
        "id": "ext_cellphonedb", "title": "CellPhoneDB 官方通讯", "group": "细胞通讯分析", "lang": "py",
        "script": "cellphonedb.py", "needs": ["clusters"], "timeout": 2400,
        "export": {"max_cells": 4000, "max_genes": 15000},
        "params": [_p("cell_type_col", "细胞类型列", "column", "clusters"), _p("n_perms", "置换次数", "number", 100)],
        "probe": {"py": ["cellphonedb"]},
        "method": "CellPhoneDB 官方统计（permutation test + 官方数据库）",
    },
    {
        "id": "ext_scvelo", "title": "scVelo 官方 RNA velocity", "group": "发育轨迹分析", "lang": "py",
        "script": "scvelo.py", "needs": ["spliced"], "timeout": 2400,
        "export": {"max_cells": 5000, "max_genes": 15000},
        "params": [_p("mode", "模式", "select", "stochastic", ["stochastic", "deterministic", "dynamical"])],
        "probe": {"py": ["scvelo"]},
        "method": "scVelo 官方实现（需要 spliced/unspliced 层）",
    },
    {
        "id": "ext_palantir", "title": "Palantir 官方轨迹", "group": "发育轨迹分析", "lang": "py",
        "script": "palantir.py", "needs": ["clusters"], "timeout": 2400,
        "params": [_p("root_cell", "根细胞名（可选）"), _p("n_components", "分量数", "number", 10)],
        "probe": {"py": ["palantir"]},
        "method": "Palantir 官方实现（diffusion maps + 终端状态 + 伪时间）",
    },
    {
        "id": "ext_celloracle", "title": "CellOracle 官方虚拟敲除", "group": "虚拟敲除", "lang": "py",
        "script": "celloracle.py", "needs": ["clusters"], "timeout": 3600,
        "params": [_p("target", "目标基因（如 CD3D）"), _p("mode", "扰动", "select", "knockout", ["knockout", "overexpress"]),
                   _p("n_top_genes", "高变基因数", "number", 2000)],
        "probe": {"py": ["celloracle"]},
        "method": "CellOracle 官方实现（base GRN + 表达 GRN + in silico perturbation）",
    },
    {
        "id": "ext_regvelo", "title": "RegVelo 官方调控速度", "group": "虚拟敲除", "lang": "py",
        "script": "regvelo.py", "needs": ["spliced"], "timeout": 3600, "python_key": "regveloPython",
        "params": [_p("target", "目标基因（可选）")],
        "probe": {"py": ["regvelo"]},
        "method": "RegVelo 官方实现（需要剪接/未剪接层 + 调控先验）",
    },
]


def _tool(tool_id):
    for item in TOOLS:
        if item["id"] == tool_id:
            return item
    return None


# ---------------------------------------------------------------------------
# 探测
# ---------------------------------------------------------------------------

_PROBE_CACHE = {"at": 0.0, "data": None}


def _run_probe(cmd, timeout=120):
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return completed.returncode, completed.stdout, completed.stderr
    except Exception as error:
        return 127, "", str(error)


def probe(force=False):
    now = time.time()
    if not force and _PROBE_CACHE["data"] and now - _PROBE_CACHE["at"] < PROBE_TTL:
        return _PROBE_CACHE["data"]
    data = {"r": {}, "py": {}, "rscript": None, "python": None, "r_versions": {}, "py_versions": {}, "errors": []}

    rscript = rscript_path()
    data["rscript"] = rscript
    packages = set()
    for item in TOOLS:
        for package in (item.get("probe") or {}).get("r") or []:
            packages.add(package)
    if rscript and packages:
        names = ", ".join('"%s"' % name for name in sorted(packages))
        code = (
            'for (p in c(%s)) { ok <- requireNamespace(p, quietly=TRUE); '
            'cat(p, if (ok) as.character(packageVersion(p)) else "MISSING", "\\n") }' % names
        )
        rc, out, err = _run_probe([rscript, "--vanilla", "-e", code], timeout=180)
        if rc == 0:
            for line in out.splitlines():
                parts = line.strip().split()
                if len(parts) >= 2:
                    data["r"][parts[0]] = parts[1] != "MISSING"
                    if parts[1] != "MISSING":
                        data["r_versions"][parts[0]] = parts[1]
        else:
            data["errors"].append("R 探测失败：%s" % (err or out)[-300:])

    # Seurat 版本（独立库优先）
    seurat_lib = _seurat_lib()
    if rscript:
        libs = [seurat_lib] if seurat_lib else []
        libs_code = ", ".join('"%s"' % item.replace("\\", "/") for item in libs)
        code = (
            'suppressWarnings(.libPaths(unique(c(c(%s), .libPaths())))); '
            'if (requireNamespace("Seurat", quietly=TRUE)) cat(as.character(packageVersion("Seurat"))) else cat("MISSING")'
            % libs_code
        )
        rc, out, err = _run_probe([rscript, "--vanilla", "-e", code], timeout=120)
        version = (out or "").strip().splitlines()[-1] if rc == 0 and out.strip() else "MISSING"
        data["seurat_lib"] = seurat_lib
        data["seurat_version"] = version
        try:
            major, minor = [int(part) for part in version.split(".")[:2]]
            data["seurat6"] = (major, minor) >= _SEURAT_MIN
        except Exception:
            data["seurat6"] = False

    python = external_python()
    data["python"] = python
    modules = set()
    for item in TOOLS:
        for module in (item.get("probe") or {}).get("py") or []:
            modules.add(module)
    if python and modules:
        code = (
            "import importlib, json\n"
            "out = {}\n"
            "for name in %r:\n"
            "    try:\n"
            "        mod = importlib.import_module(name)\n"
            "        out[name] = getattr(mod, '__version__', 'OK')\n"
            "    except Exception:\n"
            "        out[name] = None\n"
            "print(json.dumps(out))\n" % sorted(modules)
        )
        rc, out, err = _run_probe([python, "-c", code], timeout=180)
        if rc == 0:
            try:
                versions = json.loads(out.strip().splitlines()[-1])
                for name, version in versions.items():
                    data["py"][name] = bool(version)
                    if version and version != "OK":
                        data["py_versions"][name] = str(version)
            except Exception as error:
                data["errors"].append("Python 探测解析失败：%s" % error)
        else:
            data["errors"].append("Python 探测失败：%s" % (err or out)[-300:])

    # 文件类依赖
    data["scenic_db"] = _scenic_db()
    data["nichenet_priors"] = _nichenet_priors()
    data["jags"] = bool(shutil.which("jags") or Path("/usr/local/bin/jags").exists())

    _PROBE_CACHE["at"] = now
    _PROBE_CACHE["data"] = data
    return data


def _availability(tool, dataset_info=None):
    """返回 (available, reason)。dataset_info 是 None 或 {clusters:bool, group_keys:[...], has_spliced:bool}。"""
    script = EXTERNAL_DIR / str(tool.get("script") or "")
    if not script.exists():
        return False, "适配器脚本尚未就绪：%s" % script.name
    data = probe()
    spec = tool.get("probe") or {}
    reasons = []
    for package in spec.get("r") or []:
        if not data["r"].get(package):
            reasons.append("R 包 %s 未安装" % package)
    for module in spec.get("py") or []:
        if not data["py"].get(module):
            reasons.append("Python 包 %s 未安装（externalPython=%s）" % (module, data.get("python")))
    if spec.get("seurat_version"):
        version = data.get("seurat_version")
        if not data.get("seurat6"):
            reasons.append("需要 Seurat >= 6.0（独立库 %s 当前：%s）" % (data.get("seurat_lib"), version))
    if spec.get("scenic_db") and not data.get("scenic_db"):
        reasons.append("未找到 SCENIC motif 数据库目录（scenicMotifDb）")
    if spec.get("nichenet_priors") and not data.get("nichenet_priors"):
        reasons.append("未找到 NicheNet 先验模型目录（nichenetPriors）")
    if spec.get("jags") and not data.get("jags"):
        reasons.append("inferCNV 需要系统 JAGS（未找到 jags 可执行文件）")
    if dataset_info is not None:
        for need in tool.get("needs") or []:
            if need == "clusters" and not dataset_info.get("clusters"):
                reasons.append("需要先完成聚类（cluster 步骤）")
            elif need.startswith("groups:"):
                column = need.split(":", 1)[1]
                keys = dataset_info.get("group_keys") or []
                if column not in keys:
                    reasons.append("需要 obs 里的分类列 %s（现有：%s）" % (column, "、".join(keys[:6]) or "无"))
            elif need == "spliced" and not dataset_info.get("has_spliced"):
                reasons.append("需要 spliced/unspliced 层（scVelo/RegVelo 专用）")
    return (not reasons), ("；".join(reasons) if reasons else None)


def dataset_info(dataset):
    obs = dataset.obs
    group_keys = []
    for column in obs.columns:
        text = str(column)
        if text.startswith("_"):
            continue
        series = obs[column]
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            continue
        try:
            unique = int(series.nunique(dropna=True))
        except Exception:
            unique = 0
        if 1 < unique <= 200:
            group_keys.append(text)
    has_spliced = False
    try:
        layers = getattr(dataset, "layers", None) or {}
        has_spliced = any(str(key).lower() in ("spliced", "unspliced") for key in layers.keys())
    except Exception:
        has_spliced = False
    return {
        "clusters": "clusters" in obs.columns,
        "group_keys": group_keys,
        "has_spliced": has_spliced,
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
    }


def tools_payload(dataset=None):
    info = dataset_info(dataset) if dataset is not None else None
    items = []
    for tool in TOOLS:
        available, reason = _availability(tool, info)
        item = {
            "id": tool["id"], "title": tool["title"], "group": tool["group"],
            "lang": tool["lang"], "method": tool["method"],
            "params": tool.get("params") or [], "needs": tool.get("needs") or [],
            "external": True, "native": True,
            "available": bool(available) if dataset is not None else bool(
                _availability(tool, None)[0]
            ),
            "reason": None if dataset is None else (None if available else reason),
        }
        if dataset is None:
            base_available, base_reason = _availability(tool, None)
            item["available"] = base_available
            item["reason"] = base_reason
        items.append(item)
    data = probe()
    return {
        "tools": items,
        "rscript": data.get("rscript"),
        "python": data.get("python"),
        "seurat_version": data.get("seurat_version"),
        "seurat_lib": data.get("seurat_lib"),
        "versions": {"r": data.get("r_versions"), "py": data.get("py_versions")},
        "errors": data.get("errors") or [],
    }


# ---------------------------------------------------------------------------
# bundle 导出（§10.2）
# ---------------------------------------------------------------------------


def _clean_value(value):
    text = str(value)
    return text.replace(",", " ").replace("\n", " ").replace("\r", " ")


def _subset_matrix(dataset, options):
    options = options or {}
    max_cells = int(options.get("max_cells") or 15000)
    max_genes = int(options.get("max_genes") or 5000)
    X = dataset.X
    n_obs, n_vars = X.shape
    cell_idx = None
    gene_idx = None
    warnings = []
    if n_obs > max_cells:
        rng = np.random.default_rng(0)
        cell_idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))
        warnings.append("为适配外部工具对细胞下采样 %d → %d" % (n_obs, max_cells))
    if n_vars > max_genes:
        hvg = getattr(dataset, "var", None)
        picked = None
        if isinstance(hvg, pd.DataFrame) and "highly_variable" in hvg.columns:
            mask = hvg["highly_variable"].to_numpy(dtype=bool)
            if int(mask.sum()) >= 200:
                picked = np.where(mask)[0][:max_genes]
        if picked is None:
            sample = X if cell_idx is None else X[cell_idx]
            if sp is not None and sp.issparse(sample):
                mean = np.asarray(sample.mean(axis=0)).ravel()
                mean_sq = np.asarray(sample.multiply(sample).mean(axis=0)).ravel()
                var = mean_sq - mean * mean
            else:
                arr = np.asarray(sample, dtype=np.float32)
                var = arr.var(axis=0)
            picked = np.argsort(var)[::-1][:max_genes]
        gene_idx = np.sort(picked)
        warnings.append("为适配外部工具取高变/高方差基因 %d → %d" % (n_vars, len(gene_idx)))
    return cell_idx, gene_idx, warnings


def export_bundle(dataset, workdir, tag, options=None):
    if sp is None or mmwrite is None:
        raise ExternalError("internal", "缺少 scipy，无法导出外部工具数据包")
    root = Path(workdir) / "ext" / ("%s-%s" % (dataset.id, tag))
    if root.exists():
        shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    cell_idx, gene_idx, warnings = _subset_matrix(dataset, options)

    def matrix_of(source):
        matrix = source
        if cell_idx is not None:
            matrix = matrix[cell_idx]
        if gene_idx is not None:
            matrix = matrix[:, gene_idx]
        return matrix

    counts = matrix_of(dataset.counts)
    data = matrix_of(dataset.X)
    try:
        mmwrite(str(root / "counts.mtx"), sp.coo_matrix(counts.T))
        mmwrite(str(root / "data.mtx"), sp.coo_matrix(data.T))
    except Exception as error:
        raise ExternalError("internal", "导出 MatrixMarket 失败：%s" % error)

    var_names = list(dataset.var_names)
    if gene_idx is not None:
        var_names = [var_names[index] for index in gene_idx]
    obs = dataset.obs
    if cell_idx is not None:
        obs = obs.iloc[cell_idx]
    barcodes = [str(item) for item in obs.index]
    (root / "features.tsv").write_text("\n".join(str(item) for item in var_names) + "\n", encoding="utf-8")
    (root / "barcodes.tsv").write_text("\n".join(barcodes) + "\n", encoding="utf-8")

    obs_out = obs.copy()
    obs_out.insert(0, "cell", barcodes)
    for column in obs_out.columns:
        if column == "cell":
            continue
        obs_out[column] = [_clean_value(item) for item in obs_out[column].tolist()]
    obs_out.to_csv(root / "obs.csv", index=False, encoding="utf-8")

    obsm = getattr(dataset, "obsm", None) or {}
    if isinstance(obsm, dict):
        embed = obsm.get("X_emb")
        if embed is not None:
            embed = np.asarray(embed)
            if cell_idx is not None:
                embed = embed[cell_idx]
            frame = pd.DataFrame(embed[:, :2], columns=["dim1", "dim2"])
            frame.insert(0, "cell", barcodes)
            frame.to_csv(root / "embed.csv", index=False)
        pca = obsm.get("X_pca")
        if pca is not None:
            pca = np.asarray(pca)
            if cell_idx is not None:
                pca = pca[cell_idx]
            frame = pd.DataFrame(pca, columns=["PC%d" % (i + 1) for i in range(pca.shape[1])])
            frame.insert(0, "cell", barcodes)
            frame.to_csv(root / "pca.csv", index=False)

    meta = {
        "dataset": dataset.id, "name": dataset.name, "n_obs": len(barcodes), "n_vars": len(var_names),
        "steps": list(getattr(dataset, "steps", [])), "X_kind": "lognorm" if "normalize" in getattr(dataset, "steps", []) else "counts",
        "warnings": warnings,
    }
    (root / "meta.json").write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    return {"bundle": str(root), "warnings": warnings, "n_obs": len(barcodes), "n_vars": len(var_names)}


# ---------------------------------------------------------------------------
# 运行适配器
# ---------------------------------------------------------------------------


def _png_size(path):
    try:
        with open(path, "rb") as handle:
            head = handle.read(24)
        if len(head) >= 24 and head[:8] == b"\x89PNG\r\n\x1a\n":
            import struct
            width, height = struct.unpack(">II", head[16:24])
            return int(width), int(height)
    except Exception:
        pass
    return 0, 0


def _load_result(outdir, tool):
    result_path = Path(outdir) / "result.json"
    if not result_path.exists():
        raise ExternalError("internal", "适配器没有写 result.json（%s）" % tool["id"])
    try:
        result = json.loads(result_path.read_text(encoding="utf-8", errors="replace"))
    except Exception as error:
        raise ExternalError("internal", "result.json 解析失败：%s" % error)
    payload = dict(result)
    payload.setdefault("id", tool["id"])
    payload.setdefault("title", tool["title"])
    payload.setdefault("method", tool["method"])
    payload.setdefault("group", tool["group"])
    payload["external"] = True
    summary = payload.get("summary")
    if not isinstance(summary, dict):
        summary = {}
    summary.setdefault("bullets", [])
    summary.setdefault("metrics", {})
    payload["summary"] = summary
    figures = []
    for figure in payload.get("figures") or []:
        if not isinstance(figure, dict) or not figure.get("file"):
            continue
        file = str(figure["file"])
        if not os.path.isabs(file):
            file = str((Path(outdir) / file).resolve())
        if not os.path.exists(file):
            continue
        entry = dict(figure)
        entry["file"] = file
        if not entry.get("width") or not entry.get("height"):
            width, height = _png_size(file)
            entry["width"] = width
            entry["height"] = height
        figures.append(entry)
    payload["figures"] = figures
    payload["tables"] = [table for table in (payload.get("tables") or []) if isinstance(table, dict)]
    return payload


def run_external(dataset, tool_id, params, workdir, progress=None):
    tool = _tool(str(tool_id))
    if tool is None:
        raise ExternalError("bad_request", "未知外部工具：%s" % tool_id)
    info = dataset_info(dataset)
    available, reason = _availability(tool, info)
    if not available:
        raise ExternalError("unsupported", reason or "外部工具不可用")
    lang = tool["lang"]
    data = probe()
    if lang == "r":
        interpreter = data.get("rscript")
        script = EXTERNAL_DIR / tool["script"]
        if not interpreter:
            raise ExternalError("unsupported", "找不到 Rscript")
        libs = []
        if tool.get("probe", {}).get("seurat_version"):
            libs.append(_seurat_lib())
        cmd = [interpreter, "--vanilla", str(script)]
        env = dict(os.environ)
        if libs:
            env["R_LIBS_USER"] = libs[0]
    else:
        interpreter = data.get("python")
        python_key = tool.get("python_key")
        if python_key:
            configured = (_config() or {}).get(python_key)
            if configured:
                interpreter = str(configured)
            elif Path("/opt/miniconda3/envs/dsh-regvelo/bin/python").exists():
                interpreter = "/opt/miniconda3/envs/dsh-regvelo/bin/python"
        script = EXTERNAL_DIR / tool["script"]
        if not interpreter:
            raise ExternalError("unsupported", "找不到 externalPython")
        cmd = [interpreter, str(script)]
        env = dict(os.environ)
    if not script.exists():
        raise ExternalError("unsupported", "缺少适配器脚本：%s" % script.name)

    if progress:
        progress("导出数据包（%s）…" % tool["id"])
    bundle = export_bundle(dataset, workdir, tool["id"], tool.get("export"))
    outdir = Path(bundle["bundle"]) / "out"
    outdir.mkdir(parents=True, exist_ok=True)
    payload_json = json.dumps(params or {}, ensure_ascii=False)
    cmd += ["--bundle", bundle["bundle"], "--out", str(outdir), "--params", payload_json]
    env["DS_SC_SCENIC_DB"] = data.get("scenic_db") or ""
    env["DS_SC_NICHENET_PRIORS"] = data.get("nichenet_priors") or ""
    env["DS_SC_SEURAT_LIB"] = data.get("seurat_lib") or ""
    env["DS_SC_CONFIG"] = os.environ.get("DS_SC_CONFIG") or "{}"
    env.setdefault("OPENBLAS_NUM_THREADS", "2")
    env.setdefault("OMP_NUM_THREADS", "2")
    timeout = int(tool.get("timeout") or 1800)
    started = time.time()
    if progress:
        progress("运行 %s（%s）…" % (tool["title"], tool["id"]))
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        raise ExternalError("timeout", "%s 超时（>%d 秒）" % (tool["title"], timeout))
    elapsed = time.time() - started
    stderr = (completed.stderr or "")[-4000:]
    if completed.returncode == 3:
        raise ExternalError("prerequisite", stderr.strip()[-1500:] or "数据不满足该工具前置条件")
    if completed.returncode != 0:
        raise ExternalError("internal", "%s 运行失败（退出码 %d）：%s" % (tool["title"], completed.returncode, stderr.strip()[-1500:]))
    payload = _load_result(outdir, tool)
    payload["ms"] = int(round(elapsed * 1000))
    payload["warnings"] = list(payload.get("warnings") or []) + list(bundle.get("warnings") or [])
    obs_file = payload.get("obs_file") or str(outdir / "obs_extra.csv")
    obs_values = {}
    if obs_file and os.path.exists(str(obs_file)):
        try:
            frame = pd.read_csv(str(obs_file), index_col=0)
            index = [str(item) for item in dataset.obs.index]
            aligned = frame.reindex(index)
            for column in aligned.columns:
                values = aligned[column].tolist()
                obs_values[str(column)] = [None if (isinstance(value, float) and not np.isfinite(value)) else value for value in values]
            payload["obs_columns"] = [
                {"key": str(column), "kind": "numeric" if pd.api.types.is_numeric_dtype(aligned[column]) else "categorical"}
                for column in aligned.columns
            ]
        except Exception as error:
            payload.setdefault("warnings", []).append("obs_extra.csv 合并失败：%s" % error)
    payload["obs_values"] = obs_values
    # 清理大矩阵数据包，保留 out/ 里的结果图表
    try:
        for name in ("counts.mtx", "data.mtx", "pca.csv", "embed.csv"):
            path = Path(bundle["bundle"]) / name
            if path.exists():
                path.unlink()
    except Exception:
        pass
    return payload


# ---------------------------------------------------------------------------
# Seurat 6.0 分步（给引擎流程步骤用）
# ---------------------------------------------------------------------------


def seurat_available():
    data = probe()
    return bool(data.get("seurat6")), data.get("seurat_version"), data.get("seurat_lib")


def seurat_step(dataset, step, params, workdir, progress=None):
    available, version, seurat_lib = seurat_available()
    if not available:
        raise ExternalError("unsupported", "需要 Seurat >= 6.0（独立库 %s 当前：%s）" % (seurat_lib, version))
    tool = _tool("ext_seurat6")
    script = EXTERNAL_DIR / tool["script"]
    interpreter = probe().get("rscript")
    if not script.exists() or not interpreter:
        raise ExternalError("unsupported", "Seurat 6.0 适配器不可用（%s）" % script)
    bundle = export_bundle(dataset, workdir, "seurat6-%s" % step, {"max_cells": 25000, "max_genes": 8000})
    outdir = Path(bundle["bundle"]) / "out"
    outdir.mkdir(parents=True, exist_ok=True)
    step_params = dict(params or {})
    step_params["step"] = step
    env = dict(os.environ)
    env["R_LIBS_USER"] = seurat_lib
    env["DS_SC_SEURAT_LIB"] = seurat_lib
    env["DS_SC_CONFIG"] = os.environ.get("DS_SC_CONFIG") or "{}"
    cmd = [interpreter, "--vanilla", str(script), "--bundle", bundle["bundle"], "--out", str(outdir),
           "--params", json.dumps(step_params, ensure_ascii=False)]
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, env=env)
    except subprocess.TimeoutExpired:
        raise ExternalError("timeout", "Seurat 6.0 %s 超时" % step)
    stderr = (completed.stderr or "")[-3000:]
    if completed.returncode == 3:
        raise ExternalError("prerequisite", stderr.strip()[-1200:] or "Seurat 6.0 前置不满足")
    if completed.returncode != 0:
        raise ExternalError("internal", "Seurat 6.0 %s 失败（%d）：%s" % (step, completed.returncode, stderr.strip()[-1200:]))
    payload = _load_result(outdir, tool)
    result = {"payload": payload, "clusters": None, "pca": None, "embed": None, "warnings": bundle.get("warnings", [])}
    obs_file = outdir / "obs_extra.csv"
    if obs_file.exists():
        frame = pd.read_csv(obs_file, index_col=0)
        frame.index = [str(item) for item in frame.index]
        index = [str(item) for item in dataset.obs.index]
        frame = frame.reindex(index)
        for column in frame.columns:
            if column in ("seurat_clusters", "clusters"):
                result["clusters"] = [None if pd.isna(item) else str(item) for item in frame[column].tolist()]
            else:
                dataset.obs[str(column)] = frame[column]
    pca_file = outdir / "pca.csv"
    if pca_file.exists():
        frame = pd.read_csv(pca_file, index_col=0)
        frame.index = [str(item) for item in frame.index]
        frame = frame.reindex([str(item) for item in dataset.obs.index])
        result["pca"] = frame.to_numpy(dtype=np.float32)
    embed_file = outdir / "embed.csv"
    if embed_file.exists():
        frame = pd.read_csv(embed_file, index_col=0)
        frame.index = [str(item) for item in frame.index]
        frame = frame.reindex([str(item) for item in dataset.obs.index])
        result["embed"] = frame.to_numpy(dtype=np.float32)
    result["ms"] = int(round(0))
    return result


if __name__ == "__main__":  # 调试：python3 lib/sc_external.py --probe
    if "--probe" in sys.argv:
        print(json.dumps(probe(force=True), ensure_ascii=False, indent=2, default=str))
    elif "--tools" in sys.argv:
        print(json.dumps(tools_payload(), ensure_ascii=False, indent=2, default=str))
