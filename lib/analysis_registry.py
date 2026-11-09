"""
dsh-singlecell — 高级分析注册表（Host 侧，Lead 维护）

职责：
  1. 把 sc_analyses_sets.py / sc_analyses_dynamics.py 两个模块的 ANALYSES 合并成一份目录；
  2. 登记无法在本插件内原生实现、但用户能力图里要求的外部工具项（诚实标注 + 替代方案）；
  3. 从 Dataset（鸭子类型，不 import 引擎，避免循环依赖）构造分析 ctx；
  4. 按 analysis id 分发到所属模块，并规范化返回结构（补 figure 宽高、补齐固定键）。

模块加载是惰性且容错的：某个模块文件缺失或语法错误时，其它分析仍然可用，
目录里会带 diagnostics，UI/Agent 能看到原因。

契约见 CONTRACT.md §9（模块 API/ctx/结果结构/登记项）。
"""

from __future__ import annotations

import importlib.util
import shutil
import struct
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

LIB_DIR = Path(__file__).resolve().parent

# 能力图上的分组顺序（与用户给的图一致）
GROUPS = [
    "基本分析",
    "组间对比",
    "数据可视化",
    "通路分析",
    "基因集合分析",
    "发育轨迹分析",
    "肿瘤异质性分析",
    "细胞通讯分析",
    "转录调控分析",
    "虚拟敲除",
]

MODULE_FILES = [
    ("sets", "sc_analyses_sets.py"),
    ("dynamics", "sc_analyses_dynamics.py"),
]

# ---------------------------------------------------------------------------
# 导入格式清单（GET /singlecell/formats）
# ---------------------------------------------------------------------------


def _r_available():
    for candidate in ("/usr/local/bin/Rscript", "/Library/Frameworks/R.framework/Resources/bin/Rscript"):
        if Path(candidate).exists():
            return True
    return shutil.which("Rscript") is not None


def formats():
    r_ok = _r_available()
    return [
        {"id": "auto", "title": "自动识别", "extensions": [], "available": True,
         "note": "按内容与扩展名判断格式"},
        {"id": "h5ad", "title": "AnnData (.h5ad)", "extensions": [".h5ad", ".h5ad.gz"], "available": True,
         "note": "scverse/scanpy 标准格式"},
        {"id": "10x_h5", "title": "10x Genomics HDF5 (.h5)", "extensions": [".h5"], "available": True,
         "note": "CellRanger filtered_feature_bc_matrix.h5 / raw_feature_bc_matrix.h5"},
        {"id": "10x_mtx", "title": "10x Matrix Market (.mtx + barcodes/features)", "extensions": [".mtx", ".mtx.gz", ".tsv", ".tsv.gz"],
         "available": True, "note": "三件套目录，或单个 matrix.mtx（无 sidecar 时自动生成名）"},
        {"id": "loom", "title": "Loom (.loom)", "extensions": [".loom"], "available": True,
         "note": "GEO 上常见，用 h5py 原生解析"},
        {"id": "matrix_text", "title": "表达矩阵文本 (.csv/.tsv/.txt, 可 .gz)", "extensions": [".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz"],
         "available": True, "note": "自动判断行列方向；可用 transpose 强制"},
        {"id": "geo_series", "title": "GEO series matrix (…_series_matrix.txt.gz)", "extensions": [".txt.gz"],
         "available": True, "note": "解析 !series_matrix_table 并带上 !Sample_* 元数据"},
        {"id": "rds", "title": "R 对象 (.rds)", "extensions": [".rds", ".RDS"], "available": r_ok,
         "note": "Seurat v3/v4/v5、SingleCellExperiment、Matrix、data.frame" if r_ok else "未找到 Rscript"},
        {"id": "rdata", "title": "R 工作空间 (.RData/.rda)", "extensions": [".RData", ".rda", ".rdata"], "available": r_ok,
         "note": "自动挑选 Seurat/SCE 对象，多个候选时返回候选列表" if r_ok else "未找到 Rscript"},
        {"id": "archive", "title": "压缩包 (.zip/.tar/.tar.gz/.tgz)", "extensions": [".zip", ".tar", ".tar.gz", ".tgz"],
         "available": True, "note": "自动找 10x 三件套或唯一数据文件；多个候选时返回候选列表"},
        {"id": "url", "title": "URL 下载 (http/https/ftp)", "extensions": [], "available": True,
         "note": "GEO FTP/HTTP 直链，下载后自动识别"},
    ]


# ---------------------------------------------------------------------------
# 外部工具项（能力图里要求、但本插件不原生实现的分析）
# ---------------------------------------------------------------------------

EXTERNAL_ANALYSES = [
    {
        "id": "rna_velocity", "title": "RNA velocity", "group": "发育轨迹分析",
        "reason": "需要 spliced/unspliced 计数层（velocyto / scVelo 量化产物），普通表达矩阵无法推断方向。",
        "alternative": "本插件的 trajectory_dpt / trajectory_paga 可给伪时间与分支结构；提供 spliced/unspliced 层后可扩展接入 scVelo。",
    },
    {
        "id": "monocle_exact", "title": "Monocle（官方实现）", "group": "发育轨迹分析",
        "reason": "Monocle 3 是 R 包（含 UMAP/DDRTree/分支检验），本插件未内置。",
        "alternative": "用本插件 trajectory_dpt（扩散伪时间）或 trajectory_paga（簇级轨迹图）替代，结果口径不同。",
    },
    {
        "id": "palantir_exact", "title": "Palantir（官方实现）", "group": "发育轨迹分析",
        "reason": "Palantir 需要其专用 diffusion maps / 终端状态吸收概率实现。",
        "alternative": "本插件 trajectory_dpt 提供扩散图伪时间近似；终端状态概率未实现。",
    },
    {
        "id": "nichenet", "title": "NicheNet", "group": "细胞通讯分析",
        "reason": "需要官方配体-靶基因先验模型与排序框架，模型文件不随插件分发。",
        "alternative": "用本插件 lr_communication（内置 LR 对 + 表达乘积 + 置换检验）看通讯强度；靶基因层面需 NicheNet 先验。",
    },
    {
        "id": "cellphonedb", "title": "CellPhoneDB", "group": "细胞通讯分析",
        "reason": "依赖 CellPhoneDB 官方复合物数据库与置换统计框架。",
        "alternative": "用本插件 lr_communication（CellChat/NATMI-lite）；数据库更小但零配置可跑。",
    },
    {
        "id": "cellchat_exact", "title": "CellChat（官方实现）", "group": "细胞通讯分析",
        "reason": "CellChat 是 R 包，依赖其 curated 数据库与 trimean/置换检验细节。",
        "alternative": "本插件 lr_communication 为 CellChat-lite，打分口径更简单，需谨慎解释。",
    },
    {
        "id": "pyscenic_exact", "title": "pySCENIC（官方实现）", "group": "转录调控分析",
        "reason": "需要 cisTarget motif 数据库（数百 MB～GB）与 GRNBoost2/arboreto。",
        "alternative": "本插件 regulon_activity 为 SCENIC-lite（内置 TF 列表 + 共表达靶基因 + AUCell 式活性），不是 motif 驱动。",
    },
    {
        "id": "scenic_exact", "title": "SCENIC（R 版）", "group": "转录调控分析",
        "reason": "同 pySCENIC，需要 motif 排名数据库与 R 包。",
        "alternative": "regulon_activity（SCENIC-lite）。",
    },
    {
        "id": "aracne_ap", "title": "ARACNe-AP", "group": "转录调控分析",
        "reason": "Java 工具，需要单独安装与运行环境。",
        "alternative": "本插件 virtual_knockout 的 GRN 推断可作为轻量替代（相关/岭回归，非 ARACNe 互信息去卷积）。",
    },
    {
        "id": "celloracle", "title": "CellOracle", "group": "虚拟敲除",
        "reason": "需要 base GRN（染色质可及性先验）与 scATAC/velocity 数据。",
        "alternative": "本插件 virtual_knockout / virtual_screen 基于表达 GRN 做扰动传播，属近似。",
    },
    {
        "id": "regvelo", "title": "RegVelo", "group": "虚拟敲除",
        "reason": "同时需要 RNA velocity 与调控先验。",
        "alternative": "trajectory_dpt + regulon_activity + virtual_knockout 组合近似，但不是同一模型。",
    },
    {
        "id": "infercnv_exact", "title": "inferCNV（官方实现）", "group": "肿瘤异质性分析",
        "reason": "需要基因染色体位置与 R 包；官方实现含 HMM 与多种参照模式。",
        "alternative": "本插件 cnv_score 在提供 gene_order（gene, chr, start）时按染色体排序出热图；未提供时为代理分数，另有 cnv_tree / ith_score。",
    },
    {
        "id": "copykat", "title": "CopyKAT", "group": "肿瘤异质性分析",
        "reason": "R 包，含专门的全基因组 CNV 与亚克隆推断流程。",
        "alternative": "cnv_score + cnv_tree（近似）。",
    },
    {
        "id": "soupx_decontx", "title": "scRNA-seq 污染校正（SoupX/DecontX）", "group": "基本分析",
        "reason": "需要空液滴/原始（未过滤）10x 矩阵估计环境 RNA 分布，普通矩阵无法可靠去污染。",
        "alternative": "本插件提供 QC 过滤与线粒体比例检查；若有 raw_feature_bc_matrix 可先经宿主流程导入，但去污染仍需 SoupX/DecontX。",
    },
]


# ---------------------------------------------------------------------------
# 模块加载（惰性 + 容错）
# ---------------------------------------------------------------------------

_MODULE_CACHE = {}
_MODULE_ERRORS = {}


def _load_module(key, filename):
    cached = _MODULE_CACHE.get(key)
    if cached is not None:
        return cached
    name = "_dsh_sc_" + key
    if name in sys.modules:
        _MODULE_CACHE[key] = sys.modules[name]
        return sys.modules[name]
    path = LIB_DIR / filename
    if not path.exists():
        raise ImportError("缺少模块文件：%s" % path)
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError("无法加载模块：%s" % path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(name, None)
        raise
    _MODULE_CACHE[key] = module
    return module


def modules():
    """返回 [(key, module)]；单个模块失败只记录错误，不影响其它模块。"""
    loaded = []
    for key, filename in MODULE_FILES:
        try:
            loaded.append((key, _load_module(key, filename)))
            _MODULE_ERRORS.pop(key, None)
        except Exception as error:  # noqa: BLE001
            _MODULE_ERRORS[key] = "%s: %s" % (filename, error)
            _MODULE_CACHE.pop(key, None)
    return loaded


def diagnostics():
    return dict(_MODULE_ERRORS)


# ---------------------------------------------------------------------------
# ctx 构造（从 engine 的 Dataset 鸭子类型对象）
# ---------------------------------------------------------------------------


def build_ctx(dataset, engine_mode="builtin"):
    obs = dataset.obs.copy()
    clusters = None
    if "clusters" in obs.columns:
        labels = [str(value) for value in obs["clusters"].tolist()]
        clusters = {"key": "clusters", "labels": labels}

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

    obsm = {}
    raw_obsm = getattr(dataset, "obsm", None) or {}
    for key in ("X_pca", "X_emb"):
        value = raw_obsm.get(key) if isinstance(raw_obsm, dict) else None
        if value is not None:
            obsm[key] = np.asarray(value)

    x_kind = "lognorm" if "normalize" in getattr(dataset, "steps", []) else "counts"
    return {
        "dataset_id": str(getattr(dataset, "id", "ds")),
        "name": str(getattr(dataset, "name", "dataset")),
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
        "X": dataset.X,
        "X_kind": x_kind,
        "counts": dataset.counts,
        "obs": obs,
        "var_names": [str(gene) for gene in dataset.var_names],
        "obsm": obsm,
        "clusters": clusters,
        "group_keys": group_keys,
        "engine_mode": engine_mode,
        "rng": np.random.default_rng(0),
    }


# ---------------------------------------------------------------------------
# 目录与可用性
# ---------------------------------------------------------------------------


def _norm_params(definition):
    params = definition.get("params")
    if not isinstance(params, list):
        return []
    normalized = []
    for entry in params:
        if not isinstance(entry, dict):
            continue
        param = dict(entry)
        if "type" not in param:
            param["type"] = "string"
        normalized.append(param)
    return normalized


def _availability(needs, ctx, params):
    if ctx is None:
        return {"available": False, "reason": "先载入数据集（demo 或导入 h5ad/10x/loom/rds 等）"}
    missing = []
    obs = ctx["obs"]
    group_keys = ctx.get("group_keys") or []
    for need in needs or []:
        text = str(need)
        if text == "clusters":
            if "clusters" not in obs.columns:
                missing.append("需要先完成聚类（cluster 步骤）")
        elif text.startswith("groups:"):
            column = text.split(":", 1)[1]
            # 参数里的 column 类型可以覆盖默认列名（UI 会按当前 obs 列给出选项）
            override = None
            for key, value in (params or {}).items():
                if isinstance(value, str) and value in obs.columns and value in group_keys:
                    if key in (column, "group", "batch", "cell_type", "groupby", "column", "columns"):
                        override = value
                        break
            if column in obs.columns:
                try:
                    unique = int(obs[column].nunique(dropna=True))
                except Exception:
                    unique = 0
                if unique < 2:
                    missing.append("分组列 %s 只有一个取值" % column)
            elif override:
                pass
            elif group_keys:
                pass  # 有可用分类列，具体列由参数选择，模块会自行校验
            else:
                missing.append("数据没有可用的分组列（obs 里没有分类元数据）")
        elif text == "spliced":
            missing.append("需要 spliced/unspliced 层")
        elif text in ("gene_sets", "gene_order"):
            continue  # 内置基因集可用；基因顺序是可选项
        else:
            continue
    return {"available": not missing, "reason": "；".join(missing) if missing else None}


def catalog(ctx=None):
    items = []
    for key, module in modules():
        analyses = getattr(module, "ANALYSES", None) or []
        for definition in analyses:
            if not isinstance(definition, dict) or not definition.get("id"):
                continue
            item = {
                "id": str(definition["id"]),
                "title": str(definition.get("title") or definition["id"]),
                "group": str(definition.get("group") or "其它"),
                "method": str(definition.get("method") or ""),
                "native": True,
                "module": key,
                "params": _norm_params(definition),
                "needs": list(definition.get("needs") or []),
            }
            item.update(_availability(item["needs"], ctx, {}))
            items.append(item)
    for entry in EXTERNAL_ANALYSES:
        item = dict(entry)
        item["native"] = False
        item["available"] = False
        item["module"] = "external"
        item["params"] = []
        item["needs"] = []
        items.append(item)
    # 按 GROUPS 顺序稳定排序
    order = {name: index for index, name in enumerate(GROUPS)}
    items.sort(key=lambda item: (order.get(item.get("group"), 99), 0 if item.get("native") else 1, item.get("id")))
    columns = []
    if ctx is not None:
        columns = [{"key": str(col)} for col in ctx["obs"].columns if not str(col).startswith("_")]
    return {
        "groups": list(GROUPS),
        "analyses": items,
        "columns": columns,
        "group_keys": list(ctx.get("group_keys") or []) if ctx else [],
        "diagnostics": diagnostics(),
    }


# ---------------------------------------------------------------------------
# 分发
# ---------------------------------------------------------------------------


class AnalysisError(Exception):
    def __init__(self, code, message, details=None):
        self.code = str(code)
        self.message = str(message)
        self.details = details or None
        super().__init__(self.message)


def _png_size(path):
    try:
        with open(path, "rb") as handle:
            head = handle.read(24)
        if len(head) >= 24 and head[:8] == b"\x89PNG\r\n\x1a\n":
            width, height = struct.unpack(">II", head[16:24])
            return int(width), int(height)
    except Exception:
        pass
    return 0, 0


def _normalize_figures(figures):
    normalized = []
    for figure in figures or []:
        if not isinstance(figure, dict) or not figure.get("file"):
            continue
        entry = dict(figure)
        if not entry.get("width") or not entry.get("height"):
            width, height = _png_size(entry["file"])
            entry["width"] = width or entry.get("width") or 0
            entry["height"] = height or entry.get("height") or 0
        normalized.append(entry)
    return normalized


def _definition_for(analysis_id):
    for key, module in modules():
        for definition in getattr(module, "ANALYSES", None) or []:
            if str(definition.get("id")) == analysis_id:
                return key, module, definition
    return None, None, None


def run_analysis(analysis_id, ctx, params, outdir):
    key, module, definition = _definition_for(str(analysis_id))
    if module is None:
        known = [item.get("id") for item in catalog(ctx)["analyses"] if item.get("native")]
        raise AnalysisError(
            "bad_request",
            "未知分析：%s（可用：%s）" % (analysis_id, ", ".join(str(item) for item in known)),
        )
    available = _availability(definition.get("needs") or [], ctx, params)
    if not available["available"]:
        raise AnalysisError("prerequisite", available["reason"] or "缺少前置数据")
    try:
        result = module.run(str(analysis_id), ctx, params or {}, str(outdir))
    except Exception as error:  # noqa: BLE001
        code = getattr(error, "code", "internal")
        details = getattr(error, "details", None)
        message = getattr(error, "message", None) or str(error)
        if code == "internal" and not getattr(error, "code", None):
            traceback.print_exc()
            message = "分析执行失败：%s" % message
        raise AnalysisError(code, message, details)
    if not isinstance(result, dict):
        raise AnalysisError("internal", "分析模块返回了非对象结果（%s）" % analysis_id)

    payload = dict(result)
    payload.setdefault("id", str(analysis_id))
    payload.setdefault("title", definition.get("title") or str(analysis_id))
    payload.setdefault("group", definition.get("group"))
    payload.setdefault("method", definition.get("method"))
    payload["native"] = True
    summary = payload.get("summary")
    if not isinstance(summary, dict):
        summary = {"bullets": [], "metrics": {}}
    summary.setdefault("bullets", [])
    summary.setdefault("metrics", {})
    payload["summary"] = summary
    payload["tables"] = [table for table in (payload.get("tables") or []) if isinstance(table, dict)]
    payload["figures"] = _normalize_figures(payload.get("figures"))
    payload["obs_columns"] = [col for col in (payload.get("obs_columns") or []) if isinstance(col, dict)]
    payload.setdefault("warnings", [])
    return payload
