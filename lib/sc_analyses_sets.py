#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench — 基因集 / 通路 / 表达矩阵类高级分析模块（sc-analysis-a）。

契约：plugins/dsh-singlecell-workbench/CONTRACT.md §9.2（模块 API / ctx / 结果结构）、§9.3 表格。
本文件只依赖 numpy / pandas / scipy / scikit-learn / matplotlib（本机无 scanpy、
gseapy、decoupler），所有算法均在本文件内实现，不 import 引擎（EngineError 本地定义，
引擎按 ``.code`` / ``.message`` 属性识别）。

公共约定
--------
* 入口：``run(analysis_id, ctx, params, outdir)``；``ANALYSES`` 为元数据列表。
* ``ctx`` 为 dict（也兼容属性访问）：X(csr float32)、X_kind、counts、obs(DataFrame)、
  var_names、obsm{X_pca,X_emb}、clusters{key,labels}、group_keys、rng、n_obs。
* 返回固定键：id/title/group/method/summary/tables/figures/obs_columns/warnings；
  需要写回数据集 obs 的分析额外返回顶层 ``obs_values``（键→与 ctx.obs.index 等长的值列表）。
* 图：outdir 下写 ``<dataset>_<analysis>_<slug>.png``，Agg 后端、白底、dpi=140、
  图上文字一律英文；同一分析同参数结果确定（固定随机种子 0，不使用共享 ctx["rng"]，
  否则多次调用不可复现）。
* 大矩阵上限：cells ≤ 20000、genes ≤ 4000（相关性/NMF 等平方复杂度分析更低），
  超限做确定性随机下采样并在 warnings 写明。
* 缺前置数据抛 ``EngineError("prerequisite", 中文消息)``，绝不编造结果。
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import re
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# matplotlib：必须在 import pyplot 前设置 Agg + 可写缓存目录（~/.matplotlib 本机不可写）
# ---------------------------------------------------------------------------
os.environ.setdefault("MPLBACKEND", "Agg")
if not os.environ.get("MPLCONFIGDIR"):
    _MPL_FALLBACK = Path(tempfile.gettempdir()) / "dsh-singlecell-workbench-mpl"
    try:
        _MPL_FALLBACK.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(_MPL_FALLBACK)
    except Exception:  # pragma: no cover - 只影响缓存目录
        pass

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.stats as _stats  # noqa: E402
from scipy.cluster import hierarchy as _hierarchy  # noqa: E402
from scipy.spatial.distance import squareform as _squareform  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402
from matplotlib.path import Path as _MplPath  # noqa: E402
from matplotlib.patches import PathPatch as _PathPatch  # noqa: E402



# ===========================================================================
# 0. 错误与通用小工具
# ===========================================================================


class EngineError(Exception):
    """可预期的业务错误（引擎按 .code/.message 转协议错误）。"""

    CODES = {"bad_request", "not_found", "prerequisite", "unsupported", "internal"}

    def __init__(self, code, message):
        self.code = code if code in self.CODES else "internal"
        self.message = str(message)
        super().__init__(self.message)


_MODULE_CACHE = {}
_STYLE_READY = False
_TAB20 = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948",
    "#b07aa1", "#ff9da7", "#9c755f", "#bab0ab", "#1f77b4", "#ff7f0e",
    "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
]


def _builtin():
    """按文件路径加载 gene_sets_builtin.py（模块可能以任意模块名被宿主加载）。"""
    cached = _MODULE_CACHE.get("builtin")
    if cached is not None:
        return cached
    path = Path(__file__).resolve().parent / "gene_sets_builtin.py"
    try:
        spec = importlib.util.spec_from_file_location("_dsh_sc_gene_sets_builtin", str(path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception as error:  # pragma: no cover
        raise EngineError("internal", "无法加载内置基因集：%s" % error)
    _MODULE_CACHE["builtin"] = module
    return module


def _builtin_gmt():
    return _builtin().BUILTIN_GMT


def _cell_cycle_genes():
    return list(_builtin().CELL_CYCLE_S), list(_builtin().CELL_CYCLE_G2M)


def _rng():
    """固定种子随机数（契约要求同参数同结果，故不使用共享的 ctx["rng"]）。"""
    return np.random.default_rng(0)


def _ctx_get(ctx, key, default=None):
    if isinstance(ctx, dict):
        return ctx.get(key, default)
    return getattr(ctx, key, default)


def _ctx_obs(ctx):
    obs = _ctx_get(ctx, "obs")
    if not isinstance(obs, pd.DataFrame):
        raise EngineError("internal", "分析上下文缺少 obs 表")
    return obs


def _n_obs(ctx):
    obs = _ctx_obs(ctx)
    n_obs = _ctx_get(ctx, "n_obs", obs.shape[0])
    try:
        n_obs = int(n_obs)
    except Exception:
        n_obs = int(obs.shape[0])
    if n_obs != int(obs.shape[0]):
        raise EngineError("internal", "ctx.n_obs 与 obs 行数不一致")
    return n_obs


def _var_names(ctx):
    names = _ctx_get(ctx, "var_names")
    if not names:
        raise EngineError("internal", "分析上下文缺少 var_names（基因名）")
    return [str(name) for name in names]


def _num(value):
    """JSON 友好的数值：非有限值→None，极小值保留科学计数精度。"""
    if value is None:
        return None
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    try:
        number = float(value)
    except Exception:
        return None
    if not math.isfinite(number):
        return None
    if number == 0:
        return 0.0
    if abs(number) < 1e-4:
        return float("%.3e" % number)
    return round(number, 6)


def _num_list(values):
    return [_num(value) for value in np.asarray(values, dtype=np.float64).ravel().tolist()]


def _label_key(label):
    """簇/分组标签排序：数字优先按数值，其余按字符串。"""
    try:
        return (0, float(str(label)))
    except Exception:
        return (1, str(label))


def _ordered_labels(labels):
    return sorted({str(label) for label in labels}, key=_label_key)


def _slug(text, fallback="x"):
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", str(text or "")).strip("-._")
    return cleaned or fallback


def _shorten(text, limit=42):
    text = str(text)
    return text if len(text) <= limit else text[: max(1, limit - 1)] + "…"


def _bh_fdr(pvalues):
    """Benjamini-Hochberg FDR。"""
    values = np.asarray(pvalues, dtype=np.float64).ravel()
    total = values.size
    if total == 0:
        return values
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    ranks = np.arange(1, total + 1, dtype=np.float64)
    adjusted = np.minimum.accumulate((ranked * total / ranks)[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    output = np.empty_like(adjusted)
    output[order] = adjusted
    return output


def _zscore(values):
    values = np.asarray(values, dtype=np.float64).ravel()
    if values.size == 0:
        return values
    sd = float(np.std(values))
    if not math.isfinite(sd) or sd <= 0:
        return np.zeros_like(values)
    return (values - float(np.mean(values))) / sd


def _safe_p(pvalue):
    try:
        number = float(pvalue)
    except Exception:
        return 1.0
    if not math.isfinite(number):
        return 1.0
    return min(max(number, 0.0), 1.0)


def _minus_log10(pvalue):
    return float(-math.log10(max(_safe_p(pvalue), 1e-300)))


# ===========================================================================
# 1. 矩阵 / 基因 / obs 访问
# ===========================================================================


def _csr(matrix):
    if matrix is None:
        raise EngineError("prerequisite", "缺少表达矩阵（请先完成 normalize 步骤）")
    if sp.issparse(matrix):
        return sp.csr_matrix(matrix, dtype=np.float32)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float32))


def _log_normalize(matrix):
    """log1p(x / row_sum * 10000)，与 CONTRACT §2.4 的 normalize 一致。"""
    matrix = matrix.tocsr().astype(np.float64)
    row_sum = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.zeros_like(row_sum)
    np.divide(10000.0, row_sum, out=scale, where=row_sum > 0)
    result = sp.diags(scale) @ matrix
    result = result.tocsr()
    result.data = np.log1p(result.data)
    result.eliminate_zeros()
    return result.astype(np.float32)


def _expression(ctx, warnings=None):
    """取 log 标准化表达矩阵（若当前是原始计数则在分析内部补标准化）。"""
    matrix = _csr(_ctx_get(ctx, "X"))
    kind = str(_ctx_get(ctx, "X_kind", "lognorm") or "lognorm")
    if kind == "lognorm":
        return matrix
    counts = _ctx_get(ctx, "counts")
    source = _csr(counts) if counts is not None else matrix
    if warnings is not None:
        warnings.append("当前矩阵为原始计数，分析内部临时按 log1p(CP10K) 标准化")
    return _log_normalize(source)


def _gene_index(ctx):
    names = _var_names(ctx)
    exact = {name: index for index, name in enumerate(names)}
    upper = {}
    for index, name in enumerate(names):
        upper.setdefault(name.upper(), index)
    return exact, upper


def _gene_positions(index, names):
    """返回 (positions, found, missing)；大小写不敏感兜底匹配。"""
    exact, upper = index
    positions = []
    found = []
    missing = []
    seen = set()
    for name in names or []:
        text = str(name).strip()
        if not text:
            continue
        position = exact.get(text)
        if position is None:
            position = upper.get(text.upper())
        if position is None:
            missing.append(text)
        elif position not in seen:
            seen.add(position)
            positions.append(position)
            found.append(text)
    return positions, found, missing


def _limit_cells(ctx, matrix, max_cells, warnings, note=""):
    """对细胞做确定性随机下采样；返回 (matrix_sub, row_index, downsampled)。"""
    n_obs = int(matrix.shape[0])
    max_cells = int(max_cells)
    if n_obs <= max_cells:
        return matrix, np.arange(n_obs), False
    rows = np.sort(_rng().choice(n_obs, size=max_cells, replace=False))
    warnings.append(
        "细胞数 %d 超过%s上限 %d，已随机下采样到 %d（固定随机种子 0）"
        % (n_obs, note or "分析", max_cells, max_cells)
    )
    return matrix[rows], rows, True


def _gene_variance(matrix, rows, rng):
    """在抽样细胞上估计每个基因的方差（用于超限时的确定性基因筛选）。"""
    n_obs = int(matrix.shape[0])
    rows = np.asarray(rows, dtype=np.int64)
    if rows.size > 1500:
        rows = np.sort(rng.choice(rows, size=1500, replace=False))
    sub = matrix[rows]
    mean = np.asarray(sub.mean(axis=0)).ravel()
    mean_sq = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
    variance = mean_sq - mean * mean
    variance[~np.isfinite(variance)] = -np.inf
    return variance


def _limit_genes(ctx, matrix, rows, max_genes, warnings, note=""):
    """对基因做确定性筛选（超限时保留方差最大的基因）；返回 (matrix_sub, gene_index)。"""
    n_vars = int(matrix.shape[1])
    max_genes = int(max_genes)
    if n_vars <= max_genes:
        return matrix, np.arange(n_vars)
    variance = _gene_variance(matrix, rows, _rng())
    order = np.argsort(-variance, kind="stable")
    keep = np.sort(order[:max_genes])
    warnings.append(
        "基因数 %d 超过%s上限 %d，按方差保留前 %d 个基因" % (n_vars, note or "分析", max_genes, max_genes)
    )
    return matrix[:, keep], keep


def _gene_vector(matrix, position):
    column = matrix[:, int(position)]
    if sp.issparse(column):
        return np.asarray(column.todense()).ravel().astype(np.float64)
    return np.asarray(column, dtype=np.float64).ravel()


def _obs_column(ctx, name, role="obs 列"):
    obs = _ctx_obs(ctx)
    column = str(name or "").strip()
    if not column:
        raise EngineError("bad_request", "需要指定%s名称" % role)
    if column not in obs.columns:
        available = [str(col) for col in obs.columns]
        raise EngineError(
            "prerequisite",
            "缺少 obs %s %s；可用列：%s" % (role, column, "、".join(available[:20])),
        )
    return column, obs[column]


def _group_labels(ctx, column, warnings=None, role="分组列"):
    """返回 (values ndarray[str], levels list[str] 按细胞数降序, masks dict)。"""
    column, series = _obs_column(ctx, column, role=role)
    values = []
    for value in series.tolist():
        if value is None or (isinstance(value, float) and math.isnan(value)):
            values.append("NA")
        else:
            values.append(str(value))
    values = np.asarray(values, dtype=object)
    levels = _ordered_labels(values.tolist())
    if column == "clusters":
        levels = sorted(levels, key=_label_key)
    else:
        counts = pd.Series(values).value_counts()
        levels = [str(label) for label in counts.index.tolist()]
    if len(levels) < 2:
        raise EngineError("prerequisite", "%s %s 只有一个取值，无法做组间比较" % (role, column))
    if len(levels) > 60:
        raise EngineError("bad_request", "%s %s 取值过多（%d 个），不像分类变量" % (role, column, len(levels)))
    masks = {label: values == label for label in levels}
    return values, levels, masks


def _cluster_labels(ctx):
    """簇标签；缺失时抛 prerequisite。"""
    clusters = _ctx_get(ctx, "clusters")
    if isinstance(clusters, dict):
        labels = clusters.get("labels")
        if labels:
            return np.asarray([str(label) for label in labels], dtype=object)
    obs = _ctx_obs(ctx)
    if "clusters" in obs.columns:
        return np.asarray([str(label) for label in obs["clusters"].tolist()], dtype=object)
    raise EngineError("prerequisite", "需要先完成聚类（cluster 步骤）")


def _cell_type_labels(ctx, params, warnings):
    """细胞类型标签：params.cell_type → 默认 cell_type 列 → 回退 clusters。"""
    candidate = str(params.get("cell_type") or "").strip()
    obs = _ctx_obs(ctx)
    if candidate and candidate in obs.columns:
        return _group_labels(ctx, candidate, role="细胞类型列")
    if "cell_type" in obs.columns:
        return _group_labels(ctx, "cell_type", role="细胞类型列")
    if "clusters" in obs.columns:
        if warnings is not None:
            warnings.append("obs 没有 cell_type 列，改用 clusters 作为细胞类型")
        return _group_labels(ctx, "clusters", role="细胞类型列")
    raise EngineError("prerequisite", "需要 obs 里的细胞类型列（cell_type 或 clusters）")


# ===========================================================================
# 2. 结果构造（table / figure / obs）
# ===========================================================================


def _col(key, label, kind="string"):
    return {"key": str(key), "label": str(label), "type": kind}


def _table(table_id, title, columns, rows, note=None):
    table = {"id": str(table_id), "title": str(title), "columns": list(columns), "rows": list(rows)}
    if note:
        table["note"] = str(note)
    return table


def _cap_rows(rows, limit, warnings, what="表格"):
    limit = int(limit)
    if len(rows) > limit:
        warnings.append("%s行数 %d 超过上限，只保留前 %d 行" % (what, len(rows), limit))
        return rows[:limit]
    return rows


def _ensure_style():
    global _STYLE_READY
    if not _STYLE_READY:
        matplotlib.rcParams.update({
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "figure.dpi": 140,
        })
        _STYLE_READY = True


def _new_figure(width=7.0, height=4.6):
    _ensure_style()
    figure, axes = plt.subplots(figsize=(float(width), float(height)), dpi=140)
    return figure, axes


def _new_figure_grid(rows, columns, width=7.0, height=4.6):
    _ensure_style()
    figure, axes = plt.subplots(int(rows), int(columns), figsize=(float(width), float(height)),
                                dpi=140, squeeze=False)
    return figure, axes


def _png_size(path):
    try:
        with open(path, "rb") as handle:
            header = handle.read(24)
        if len(header) >= 24 and header[:8] == b"\x89PNG\r\n\x1a\n":
            return int.from_bytes(header[16:20], "big"), int.from_bytes(header[20:24], "big")
    except Exception:
        pass
    return 0, 0


def _save_figure(figure, ctx, outdir, analysis_id, slug, title):
    directory = Path(str(outdir)).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    dataset = _slug(_ctx_get(ctx, "name", "dataset"), fallback="dataset")
    path = directory / ("%s_%s_%s.png" % (dataset, _slug(analysis_id), _slug(slug)))
    try:
        figure.savefig(str(path), dpi=140, facecolor="white", bbox_inches="tight")
    finally:
        plt.close(figure)
    width, height = _png_size(path)
    return {
        "id": str(slug),
        "title": str(title),
        "file": str(path),
        "width": int(width),
        "height": int(height),
    }


def _obs_pair(ctx, key, kind, values):
    """构造 obs_columns 条目 + obs_values 值列表（长度必须等于细胞数）。"""
    values = list(values)
    n_obs = _n_obs(ctx)
    if len(values) != n_obs:
        raise EngineError("internal", "obs 列 %s 长度 %d 与细胞数 %d 不一致" % (key, len(values), n_obs))
    if kind == "categorical":
        cleaned = ["" if value is None else str(value) for value in values]
        return {"key": str(key), "kind": "categorical"}, {str(key): cleaned}
    cleaned = [_num(value) for value in values]
    return {"key": str(key), "kind": "numeric"}, {str(key): cleaned}


def _palette(count):
    """确定性调色板（前 20 个 Tableau 色，之后黄金角 HSV）。"""
    count = int(max(0, count))
    colors = list(_TAB20[: min(count, len(_TAB20))])
    index = len(colors)
    while len(colors) < count:
        hue = (index * 0.618033988749895) % 1.0
        import colorsys

        red, green, blue = colorsys.hsv_to_rgb(hue, 0.58, 0.86)
        colors.append("#%02x%02x%02x" % (int(red * 255), int(green * 255), int(blue * 255)))
        index += 1
    return colors


def _style_axis(axis, grid_axis="y"):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if grid_axis:
        axis.grid(axis=grid_axis, color="#e6e6e6", linewidth=0.6)
        axis.set_axisbelow(True)


def _spec(analysis_id):
    for item in ANALYSES:
        if item["id"] == analysis_id:
            return item
    raise EngineError(
        "not_found",
        "未知分析：%s（sc-analysis-a 支持：%s）" % (analysis_id, ", ".join(item["id"] for item in ANALYSES)),
    )


def _result(analysis_id, bullets, metrics, tables, figures, obs_columns=None,
            obs_values=None, warnings=None):
    spec = _spec(analysis_id)
    return {
        "id": spec["id"],
        "title": spec["title"],
        "group": spec["group"],
        "method": spec["method"],
        "summary": {"bullets": [str(item) for item in bullets or []],
                    "metrics": dict(metrics or {})},
        "tables": list(tables or []),
        "figures": list(figures or []),
        "obs_columns": list(obs_columns or []),
        "obs_values": dict(obs_values or {}),
        "warnings": [str(item) for item in warnings or []],
    }


# ===========================================================================
# 3. 基因集来源（内置 / 用户 GMT / Enrichr）
# ===========================================================================


def _read_gmt(path):
    """解析标准 GMT：term <TAB> description <TAB> gene1 <TAB> gene2 ..."""
    text = str(path or "").strip()
    if not text:
        raise EngineError("bad_request", "source=gmt 时必须提供 gmt_file 参数")
    file_path = Path(text).expanduser()
    if not file_path.exists() or not file_path.is_file():
        raise EngineError("bad_request", "GMT 文件不存在或不是文件：%s" % text)
    sets = {}
    try:
        with open(file_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\r\n")
                if not line.strip() or line.startswith("#"):
                    continue
                parts = [part.strip() for part in line.split("\t")]
                parts = [part for part in parts if part]
                if not parts:
                    continue
                term = parts[0]
                genes = parts[2:] if len(parts) >= 3 else parts[1:]
                genes = [gene for gene in genes if gene]
                if not term or not genes:
                    continue
                seen = set()
                ordered = []
                for gene in genes:
                    if gene not in seen:
                        seen.add(gene)
                        ordered.append(gene)
                sets[term] = ordered
    except EngineError:
        raise
    except Exception as error:
        raise EngineError("bad_request", "读取 GMT 文件失败：%s" % error)
    if not sets:
        raise EngineError("bad_request", "GMT 文件里没有解析到任何基因集：%s" % text)
    return sets


def _shorten_gmt_limit(sets, limit=200, warnings=None):
    """基因集过多时限制数量（保持文件/内置顺序，确定性）。"""
    if len(sets) <= limit:
        return sets
    if warnings is not None:
        warnings.append("基因集数量 %d 超过上限 %d，只使用前 %d 个" % (len(sets), limit, limit))
    return dict(list(sets.items())[:limit])


def _resolve_gene_sets(params, warnings):
    """按 params.source 返回 (sets: dict, source_label)。"""
    source = str(params.get("source") or "builtin").strip().lower()
    if source in ("builtin", "internal", "内置"):
        return _shorten_gmt_limit(_builtin_gmt(), 200, warnings), "builtin curated gene sets"
    if source in ("gmt", "file", "local"):
        return _shorten_gmt_limit(_read_gmt(params.get("gmt_file")), 200, warnings), "user GMT file"
    raise EngineError("bad_request", "未知基因集来源：%s（可用 builtin / gmt）" % source)


def _enrichr_ora(genes, library):
    """调用 Enrichr 在线 ORA；失败抛可读错误（超时 20s）。"""
    genes = [str(gene) for gene in genes if str(gene).strip()][:2000]
    if not genes:
        raise EngineError("bad_request", "Enrichr 在线富集需要非空基因列表")
    library = str(library or "KEGG_2021_Human").strip()
    try:
        payload = urllib.parse.urlencode({
            "list": "\n".join(genes),
            "description": "dsh-singlecell-workbench ORA",
        }).encode("utf-8")
        request = urllib.request.Request(
            "https://maayanlab.cloud/Enrichr/addList",
            data=payload,
            headers={"User-Agent": "dsh-singlecell-workbench/0.1", "Content-Type": "application/x-www-form-urlencoded"},
        )
        with urllib.request.urlopen(request, timeout=20) as response:
            added = json.loads(response.read().decode("utf-8"))
        user_list_id = added.get("userListId")
        if not user_list_id:
            raise EngineError("unsupported", "Enrichr addList 未返回 userListId")
        query = urllib.parse.urlencode({"userListId": user_list_id, "backgroundType": library})
        with urllib.request.urlopen(
            "https://maayanlab.cloud/Enrichr/enrich?" + query, timeout=20
        ) as response:
            enriched = json.loads(response.read().decode("utf-8"))
    except EngineError:
        raise
    except urllib.error.HTTPError as error:
        raise EngineError("unsupported", "Enrichr 在线服务返回 HTTP %s：可改用 source=builtin 或本地 .gmt" % error.code)
    except Exception as error:  # 网络不可用 / 超时 / 响应格式异常
        raise EngineError(
            "unsupported",
            "Enrichr 在线富集失败（%s）：请检查网络，或改用 source=builtin / 本地 .gmt" % error,
        )
    rows = enriched.get(library) or []
    output = []
    for entry in rows:
        try:
            rank, term, pvalue, zscore, combined, hits = (
                int(entry[0]), str(entry[1]), float(entry[2]), float(entry[3]),
                float(entry[4]), list(entry[5]),
            )
        except Exception:
            continue
        output.append({
            "rank": rank, "term": term, "pvalue": pvalue, "zscore": zscore,
            "combined_score": combined, "hits": [str(hit) for hit in hits],
        })
    if not output:
        raise EngineError("unsupported", "Enrichr 返回空结果（库 %s）：确认基因是否为人类 symbol" % library)
    pvalues = np.asarray([row["pvalue"] for row in output], dtype=np.float64)
    for row, adjusted in zip(output, _bh_fdr(pvalues)):
        row["pvalue_adj"] = _safe_p(adjusted)
    return output


# ===========================================================================
# 4. 统计分析小工具
# ===========================================================================


def _fig_label(label, fallback):
    """图上标签：非 ASCII 时退化为编号（本机无 CJK 字体，契约要求图上英文）。"""
    text = str(label)
    if text.isascii():
        return text
    return fallback


def _obs_key(prefix, term, used):
    """把通路名转成合法的 obs 列名（ASCII、唯一、不过长）。"""
    base = re.sub(r"[^A-Za-z0-9]+", "_", str(term)).strip("_")[:26] or "set"
    key = "%s_%s" % (prefix, base)
    candidate = key
    index = 2
    while candidate in used:
        candidate = "%s_%d" % (key, index)
        index += 1
    used.add(candidate)
    return candidate


def _group_stats(matrix, rows):
    """一组细胞在全部基因上的 (均值, 样本方差, 细胞数)，稀疏向量化。"""
    rows = np.asarray(rows, dtype=np.int64)
    n = int(rows.size)
    if n == 0:
        n_vars = int(matrix.shape[1])
        return np.zeros(n_vars), np.zeros(n_vars), 0
    sub = matrix[rows]
    total = np.asarray(sub.sum(axis=0)).ravel().astype(np.float64)
    total_sq = np.asarray(sub.multiply(sub).sum(axis=0)).ravel().astype(np.float64)
    mean = total / n
    if n > 1:
        variance = (total_sq - total * total / n) / (n - 1)
    else:
        variance = np.zeros_like(mean)
    return mean, np.maximum(variance, 0.0), n


def _pct_expressed(matrix, rows):
    rows = np.asarray(rows, dtype=np.int64)
    n = int(rows.size)
    if n == 0:
        return np.zeros(int(matrix.shape[1]))
    sub = matrix[rows]
    non_zero = np.asarray(sub.getnnz(axis=0)).ravel().astype(np.float64)
    return non_zero / n


def _welch_t(mean1, var1, n1, mean2, var2, n2):
    """Welch t 检验（向量化），返回 (t, df, p)。零方差基因给 t=0/p=1。"""
    se2 = var1 / max(n1, 1) + var2 / max(n2, 1)
    valid = se2 > 0
    t_stat = np.zeros_like(se2)
    t_stat[valid] = (mean1[valid] - mean2[valid]) / np.sqrt(se2[valid])
    numerator = np.square(se2)
    denominator = np.zeros_like(se2)
    if n1 > 1:
        denominator = denominator + np.square(var1 / n1) / (n1 - 1)
    if n2 > 1:
        denominator = denominator + np.square(var2 / n2) / (n2 - 1)
    df = np.full_like(se2, float(max(n1 + n2 - 2, 1)))
    good = denominator > 0
    df[good] = numerator[good] / denominator[good]
    df = np.maximum(df, 1.0)
    pvalue = 2.0 * _stats.t.sf(np.abs(t_stat), df)
    pvalue = np.where(valid, pvalue, 1.0)
    return t_stat, df, np.clip(pvalue, 0.0, 1.0)


def _beta_ci(count, total, level=0.95):
    """Beta(1,1) 先验下比例的后验可信区间（反事实区间的近似）。"""
    count = float(max(0.0, min(float(count), float(total))))
    total = float(total)
    if total <= 0:
        return 0.0, 1.0
    alpha = 1.0 + count
    beta = 1.0 + total - count
    low = float(_stats.beta.ppf((1.0 - level) / 2.0, alpha, beta))
    high = float(_stats.beta.ppf(1.0 - (1.0 - level) / 2.0, alpha, beta))
    return max(0.0, low), min(1.0, high)


def _rank_ascending(values):
    """0-based 升序排名（ordinal，确定性；平局按基因顺序）。"""
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=np.int64)
    ranks[order] = np.arange(len(values), dtype=np.int64)
    return ranks


def _safe_div(numerator, denominator, default=0.0):
    if denominator == 0 or not math.isfinite(float(denominator)):
        return default
    value = float(numerator) / float(denominator)
    return value if math.isfinite(value) else default


# ===========================================================================
# 5. 分析元数据（id / title / group / method / params / needs 与 §9.3 表格逐字一致）
# ===========================================================================

ANALYSES = [
    {
        "id": "cell_cycle",
        "title": "细胞周期评分",
        "group": "基本分析",
        "method": "Tirosh/Regev 内置 S/G2M 基因集均值打分 + 中心化相对阈值判相位（本插件实现，非 Seurat 官方实现）",
        "params": [
            {"key": "scoring", "label": "打分方式", "type": "select",
             "default": "mean", "options": ["mean", "zscore"],
             "hint": "mean=基因集平均表达；zscore=跨细胞 z-score"},
            {"key": "threshold", "label": "相对阈值", "type": "number", "default": 0.0,
             "hint": "中心化分数超过该值视为阳性（默认 0）"},
        ],
        "needs": [],
    },
    {
        "id": "subcluster",
        "title": "亚型细分",
        "group": "基本分析",
        "method": "簇内重算 HVG→PCA→kNN→纯 Python Louvain（本插件实现，非 scanpy/Leiden 官方实现）",
        "params": [
            {"key": "cluster", "label": "目标簇", "type": "string", "default": "0",
             "hint": "要细分的簇 id，如 0 或 3"},
            {"key": "resolution", "label": "分辨率", "type": "number", "default": 0.8,
             "hint": "Louvain 分辨率，越大亚簇越多"},
            {"key": "n_neighbors", "label": "近邻数", "type": "number", "default": 15},
            {"key": "n_top_genes", "label": "HVG 数", "type": "number", "default": 1000},
            {"key": "n_pcs", "label": "主成分数", "type": "number", "default": 20},
            {"key": "min_cells", "label": "最小组细胞数", "type": "number", "default": 5,
             "hint": "小于该数的亚簇并入邻居最多的亚簇"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "proportion_test",
        "title": "组间细胞类型比例差异",
        "group": "组间对比",
        "method": "各组列联表卡方检验 + 标准化残差富集方向 + Beta 后验区间（Dirichlet-multinomial 近似，本插件实现）",
        "params": [
            {"key": "group", "label": "分组列", "type": "column", "default": "group"},
            {"key": "cell_type", "label": "细胞类型列", "type": "column", "default": "cell_type"},
            {"key": "min_cells", "label": "最小组细胞数", "type": "number", "default": 10},
        ],
        "needs": ["groups:group", "clusters"],
    },
    {
        "id": "de_group",
        "title": "同细胞类型组间差异基因",
        "group": "组间对比",
        "method": "逐细胞类型 Welch t 检验 + BH FDR（logFC=log 标准化均值差；本插件实现）",
        "params": [
            {"key": "group", "label": "分组列", "type": "column", "default": "group"},
            {"key": "cell_type", "label": "细胞类型列", "type": "column", "default": "cell_type"},
            {"key": "group1", "label": "组 1", "type": "string", "default": "",
             "hint": "留空默认取细胞数最多的两个组"},
            {"key": "group2", "label": "组 2", "type": "string", "default": ""},
            {"key": "min_cells", "label": "每组最少细胞数", "type": "number", "default": 10},
            {"key": "top_n", "label": "每类返回基因数", "type": "number", "default": 50},
        ],
        "needs": ["groups:group", "clusters"],
    },
    {
        "id": "volcano",
        "title": "差异基因火山图",
        "group": "数据可视化",
        "method": "复用 de_group 的 Welch t + BH 结果绘制火山图（本插件实现）",
        "params": [
            {"key": "group", "label": "分组列", "type": "column", "default": "group"},
            {"key": "cell_type", "label": "细胞类型列", "type": "column", "default": "cell_type"},
            {"key": "group1", "label": "组 1", "type": "string", "default": ""},
            {"key": "group2", "label": "组 2", "type": "string", "default": ""},
            {"key": "top_labels", "label": "标注基因数", "type": "number", "default": 10},
            {"key": "padj_cutoff", "label": "FDR 阈值", "type": "number", "default": 0.05},
            {"key": "logfc_cutoff", "label": "logFC 阈值", "type": "number", "default": 0.25},
        ],
        "needs": ["groups:group", "clusters"],
    },
    {
        "id": "proportion_sankey",
        "title": "细胞比例桑基图",
        "group": "数据可视化",
        "method": "组→细胞类型两级 alluvial 桑基图（matplotlib 贝塞尔多边形实现，本插件实现）",
        "params": [
            {"key": "group", "label": "分组列", "type": "column", "default": "group"},
            {"key": "cell_type", "label": "细胞类型列", "type": "column", "default": "cell_type"},
            {"key": "min_cells", "label": "最小流向细胞数", "type": "number", "default": 5},
        ],
        "needs": ["groups:group", "clusters"],
    },
    {
        "id": "featureplot",
        "title": "差异基因特征图",
        "group": "数据可视化",
        "method": "嵌入坐标着色 / 小提琴图 / dotplot / 热图（本插件 matplotlib 实现，最多 6 个基因）",
        "params": [
            {"key": "genes", "label": "基因列表", "type": "text",
             "default": "CD3D,MS4A1,NKG7,LYZ,PPBP", "hint": "逗号分隔，最多 6 个"},
            {"key": "kind", "label": "图类型", "type": "select",
             "default": "embedding", "options": ["embedding", "vln", "dot", "heatmap"]},
            {"key": "groupby", "label": "分组列", "type": "column", "default": "clusters"},
            {"key": "max_genes", "label": "最多基因数", "type": "number", "default": 6},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "pathway_ora",
        "title": "GO/KEGG 富集（ORA）",
        "group": "通路分析",
        "method": "超几何检验 + BH FDR；内置精选 GMT / 用户 .gmt 文件 / 可选 Enrichr 在线（source=enrichr）",
        "params": [
            {"key": "source", "label": "基因集来源", "type": "select",
             "default": "builtin", "options": ["builtin", "gmt", "enrichr"]},
            {"key": "gene_list", "label": "前景基因列表", "type": "text", "default": "",
             "hint": "逗号/空格分隔；留空则用各簇 top 标记基因"},
            {"key": "cluster", "label": "指定簇（可选）", "type": "string", "default": "",
             "hint": "留空=所有簇标记基因并集"},
            {"key": "gmt_file", "label": "GMT 文件", "type": "file", "default": ""},
            {"key": "enrichr_library", "label": "Enrichr 库", "type": "select",
             "default": "KEGG_2021_Human",
             "options": ["KEGG_2021_Human", "GO_Biological_Process_2023", "Reactome_2022"]},
            {"key": "top_n", "label": "返回条数", "type": "number", "default": 15},
            {"key": "min_set_size", "label": "最小基因集大小", "type": "number", "default": 3},
        ],
        "needs": ["gene_sets"],
    },
    {
        "id": "gsea_prerank",
        "title": "GSEA（预排序）",
        "group": "通路分析",
        "method": "加权 KS running-sum + 置换检验算 NES/p（本插件实现，非 GSEA 官方实现）",
        "params": [
            {"key": "ranking", "label": "排序依据", "type": "select",
             "default": "pca1", "options": ["pca1", "dispersion", "custom"],
             "hint": "pca1=基因与 PC1 相关；dispersion=归一化离散度；custom=自定义基因分数"},
            {"key": "gene_scores", "label": "自定义排序", "type": "text", "default": "",
             "hint": "ranking=custom 时使用，格式 GENE:score,GENE:score"},
            {"key": "weight", "label": "权重指数", "type": "number", "default": 1.0,
             "hint": "0=不加权，1=经典加权 KS"},
            {"key": "n_perm", "label": "置换次数", "type": "number", "default": 1000},
            {"key": "top_n", "label": "返回条数", "type": "number", "default": 20},
            {"key": "min_set_size", "label": "最小基因集大小", "type": "number", "default": 5},
            {"key": "max_set_size", "label": "最大基因集大小", "type": "number", "default": 500},
            {"key": "gmt_file", "label": "GMT 文件", "type": "file", "default": ""},
        ],
        "needs": ["gene_sets"],
    },
    {
        "id": "gsva_like",
        "title": "GSVA 通路活性",
        "group": "通路分析",
        "method": "每细胞秩加权 KS 富集分数（ssGSEA 式简化，本插件实现，非 GSVA 官方实现）",
        "params": [
            {"key": "n_sets", "label": "通路数", "type": "number", "default": 12},
            {"key": "gmt_file", "label": "GMT 文件", "type": "file", "default": ""},
            {"key": "obs_prefix", "label": "obs 列前缀", "type": "string", "default": "gsva"},
            {"key": "max_cells", "label": "最大细胞数", "type": "number", "default": 5000},
        ],
        "needs": ["gene_sets"],
    },
    {
        "id": "auc_score",
        "title": "细胞基因集合分析（AUC）",
        "group": "基因集合分析",
        "method": "AUCell 式 rank-AUC（秩和公式，本插件实现；非 AUCell 官方实现）",
        "params": [
            {"key": "n_sets", "label": "基因集数", "type": "number", "default": 12},
            {"key": "gmt_file", "label": "GMT 文件", "type": "file", "default": ""},
            {"key": "obs_prefix", "label": "obs 列前缀", "type": "string", "default": "auc"},
            {"key": "max_cells", "label": "最大细胞数", "type": "number", "default": 5000},
            {"key": "threshold", "label": "活性阈值", "type": "number", "default": 0.0,
             "hint": "0=自动（均值+1SD）"},
        ],
        "needs": ["gene_sets"],
    },
    {
        "id": "cnmf_programs",
        "title": "cNMF 表达程序",
        "group": "基因集合分析",
        "method": "sklearn NMF 在 HVG 表达矩阵上分解 k 个表达程序（本插件实现，非 cNMF 官方实现）",
        "params": [
            {"key": "k", "label": "程序数 k", "type": "number", "default": 6},
            {"key": "n_top_genes", "label": "HVG 数", "type": "number", "default": 1000},
            {"key": "max_iter", "label": "最大迭代", "type": "number", "default": 300},
            {"key": "top_genes", "label": "每个程序基因数", "type": "number", "default": 15},
            {"key": "obs_prefix", "label": "obs 列前缀", "type": "string", "default": "nmf"},
        ],
        "needs": [],
    },
    {
        "id": "wgcna_lite",
        "title": "WGCNA 共表达模块",
        "group": "基因集合分析",
        "method": "WGCNA-lite（HVG 相关矩阵 + 层次聚类模块 + 模块特征基因/枢纽基因；本插件实现）",
        "params": [
            {"key": "n_top_genes", "label": "HVG 数", "type": "number", "default": 800},
            {"key": "min_corr", "label": "模块最小相关", "type": "number", "default": 0.5,
             "hint": "相关低于该值的基因分到不同模块（树高 = 1 - min_corr）"},
            {"key": "min_module_size", "label": "最小模块大小", "type": "number", "default": 20},
            {"key": "n_hub", "label": "每模块枢纽基因数", "type": "number", "default": 5},
        ],
        "needs": [],
    },
    {
        "id": "geneset_clinical",
        "title": "基因集分数联合临床因素",
        "group": "基因集合分析",
        "method": "基因集平均表达 z-score × obs 元数据关联（分类 ANOVA/Kruskal，数值 Spearman；本插件实现）",
        "params": [
            {"key": "gmt_file", "label": "GMT 文件", "type": "file", "default": ""},
            {"key": "n_sets", "label": "基因集数", "type": "number", "default": 8},
            {"key": "columns", "label": "元数据列", "type": "text", "default": "",
             "hint": "逗号分隔；留空自动检测 obs 里的分类/数值列"},
            {"key": "max_levels", "label": "分类列最多水平数", "type": "number", "default": 12},
        ],
        "needs": ["groups:group"],
    },
]


# ===========================================================================
# 6. 基本分析：cell_cycle / subcluster
# ===========================================================================


def _run_cell_cycle(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    index = _gene_index(ctx)
    s_genes, g2m_genes = _cell_cycle_genes()
    s_positions, s_found, _ = _gene_positions(index, s_genes)
    g_positions, g_found, _ = _gene_positions(index, g2m_genes)
    # 合成/非人类数据可能只命中很少（甚至单侧 0 个）细胞周期基因：此时照常打分并给出低置信度提示，
    # 只有两侧都完全未命中才无法计算。
    if len(s_positions) + len(g_positions) == 0:
        raise EngineError(
            "prerequisite",
            "内置细胞周期基因集在数据中完全未命中（S 0 个、G2M 0 个），无法打分；"
            "可能是基因符号不匹配或数据为合成/非人类数据",
        )
    if len(s_positions) == 0 or len(g_positions) == 0:
        warnings.append(
            "细胞周期基因集命中极少：S %d 个、G2M %d 个；空缺一侧的分数记为 0，相位判定置信度低"
            % (len(s_positions), len(g_positions))
        )
    elif len(s_positions) < 5 or len(g_positions) < 5:
        warnings.append(
            "细胞周期基因集在数据中命中较少（S %d 个、G2M %d 个），打分置信度低"
            % (len(s_positions), len(g_positions))
        )
    n_obs_cells = int(matrix.shape[0])
    if s_positions:
        s_score = np.asarray(matrix[:, s_positions].mean(axis=1)).ravel().astype(np.float64)
    else:
        s_score = np.zeros(n_obs_cells, dtype=np.float64)
    if g_positions:
        g2m_score = np.asarray(matrix[:, g_positions].mean(axis=1)).ravel().astype(np.float64)
    else:
        g2m_score = np.zeros(n_obs_cells, dtype=np.float64)
    scoring = str(params.get("scoring") or "mean").lower()
    if scoring == "zscore":
        s_score = _zscore(s_score)
        g2m_score = _zscore(g2m_score)
    else:
        scoring = "mean"
    try:
        threshold = float(params.get("threshold", 0.0))
    except Exception:
        threshold = 0.0
    n_obs = int(matrix.shape[0])
    s_relative = s_score - float(np.mean(s_score))
    g_relative = g2m_score - float(np.mean(g2m_score))
    s_high = s_relative > threshold
    g_high = g_relative > threshold
    phase = np.full(n_obs, "G1", dtype=object)
    phase[s_high & ~g_high] = "S"
    phase[g_high & ~s_high] = "G2M"
    both = s_high & g_high
    phase[both] = np.where(s_relative[both] >= g_relative[both], "S", "G2M")

    counts = {name: int(np.sum(phase == name)) for name in ("G1", "S", "G2M")}
    obs_columns = []
    obs_values = {}
    for key, kind, values in (
        ("S_score", "numeric", s_score),
        ("G2M_score", "numeric", g2m_score),
        ("phase", "categorical", phase),
    ):
        entry, values_map = _obs_pair(ctx, key, kind, values)
        obs_columns.append(entry)
        obs_values.update(values_map)

    # 图 1：相位分布 + S/G2M 相对分数散点
    palette = {"G1": "#b0b0b0", "S": "#4e79a7", "G2M": "#e15759"}
    figure, axes_grid = _new_figure_grid(1, 2, 9.6, 4.2)
    axes = axes_grid[0]
    bars = axes[0].bar(
        list(counts.keys()), [counts[key] for key in counts],
        color=[palette[key] for key in counts], width=0.62, edgecolor="#333333", linewidth=0.6,
    )
    for bar, key in zip(bars, counts):
        axes[0].text(
            bar.get_x() + bar.get_width() / 2.0, bar.get_height() + max(counts.values()) * 0.01,
            "%d (%.1f%%)" % (counts[key], 100.0 * counts[key] / max(n_obs, 1)),
            ha="center", va="bottom", fontsize=8,
        )
    axes[0].set_title("Cell cycle phase distribution")
    axes[0].set_ylabel("Cells")
    axes[0].set_xlabel("Phase")
    _style_axis(axes[0])

    sample = np.arange(n_obs)
    if n_obs > 3000:
        sample = np.sort(_rng().choice(n_obs, size=3000, replace=False))
    for name in ("G1", "S", "G2M"):
        mask = phase[sample] == name
        axes[1].scatter(
            s_relative[sample][mask], g_relative[sample][mask],
            s=9, alpha=0.55, linewidths=0, color=palette[name], label=name,
        )
    axes[1].axhline(threshold, color="#888888", linestyle="--", linewidth=0.8)
    axes[1].axvline(threshold, color="#888888", linestyle="--", linewidth=0.8)
    axes[1].set_xlabel("S score (centered)")
    axes[1].set_ylabel("G2M score (centered)")
    axes[1].set_title("S vs G2M relative scores")
    axes[1].legend(frameon=False, loc="best", markerscale=1.6)
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "cell_cycle", "phase_distribution", "Cell cycle phase distribution")]

    phase_rows = [
        {"phase": name, "n_cells": counts[name], "pct": _num(100.0 * counts[name] / max(n_obs, 1))}
        for name in ("G1", "S", "G2M")
    ]
    coverage_rows = [
        {"gene_set": "S", "n_genes": len(s_genes), "n_matched": len(s_positions),
         "pct_matched": _num(100.0 * len(s_positions) / max(len(s_genes), 1))},
        {"gene_set": "G2M", "n_genes": len(g2m_genes), "n_matched": len(g_positions),
         "pct_matched": _num(100.0 * len(g_positions) / max(len(g2m_genes), 1))},
    ]
    tables = [
        _table(
            "phase_summary", "Cell cycle phase summary",
            [_col("phase", "Phase"), _col("n_cells", "Cells", "number"), _col("pct", "Percent", "number")],
            phase_rows,
        ),
        _table(
            "gene_set_coverage", "Gene set coverage",
            [_col("gene_set", "Gene set"), _col("n_genes", "Genes in set", "number"),
             _col("n_matched", "Matched in data", "number"), _col("pct_matched", "Percent", "number")],
            coverage_rows,
        ),
    ]
    bullets = [
        "S 基因集命中 %d/%d，G2M 命中 %d/%d（内置 Tirosh/Regev 版）" % (
            len(s_positions), len(s_genes), len(g_positions), len(g2m_genes)),
        "相位划分：G1 %d、S %d、G2M %d（相对阈值 %.2f，打分=%s）" % (
            counts["G1"], counts["S"], counts["G2M"], threshold, scoring),
        "S_score/G2M_score/phase 已写回 obs（phase 可用于嵌入图着色）",
    ]
    if len(s_positions) < 5 or len(g_positions) < 5:
        bullets.append(
            "注意：细胞周期基因集命中较少（S %d、G2M %d），打分置信度低；"
            "建议使用人类基因符号的数据" % (len(s_positions), len(g_positions))
        )
    metrics = {
        "n_cells": n_obs,
        "n_s_genes": len(s_positions),
        "n_g2m_genes": len(g_positions),
        "pct_g1": _num(100.0 * counts["G1"] / max(n_obs, 1)),
        "pct_s": _num(100.0 * counts["S"] / max(n_obs, 1)),
        "pct_g2m": _num(100.0 * counts["G2M"] / max(n_obs, 1)),
    }
    return _result("cell_cycle", bullets, metrics, tables, figures,
                   obs_columns=obs_columns, obs_values=obs_values, warnings=warnings)


def _louvain_local_moving(adjacency, resolution=1.0, max_passes=50):
    """Louvain 局部移动阶段（模块度增益 = w_iC - resolution * tot_C * k_i / 2m）。"""
    n_obs = int(adjacency.shape[0])
    if n_obs == 0:
        return np.zeros(0, dtype=np.int64)
    indptr = adjacency.indptr
    indices = adjacency.indices
    data = adjacency.data
    degree = np.asarray(adjacency.sum(axis=1)).ravel().astype(np.float64)
    total_degree = float(degree.sum())
    labels = np.arange(n_obs, dtype=np.int64)
    if total_degree <= 0:
        return labels
    community_total = degree.copy()
    for _ in range(int(max_passes)):
        moved = 0
        for node in range(n_obs):
            ki = degree[node]
            if ki <= 0:
                continue
            start, end = int(indptr[node]), int(indptr[node + 1])
            neighbours = indices[start:end]
            if neighbours.size == 0:
                continue
            weights = data[start:end].astype(np.float64)
            neighbour_labels = labels[neighbours]
            unique_labels, inverse = np.unique(neighbour_labels, return_inverse=True)
            weight_sums = np.bincount(inverse, weights=weights, minlength=unique_labels.size)
            gains = weight_sums - resolution * community_total[unique_labels] * ki / total_degree
            best_index = int(np.argmax(gains))
            best_label = int(unique_labels[best_index])
            current_label = int(labels[node])
            stay = -resolution * community_total[current_label] * ki / total_degree
            current_position = np.where(unique_labels == current_label)[0]
            if current_position.size:
                stay += float(weight_sums[current_position[0]])
            if gains[best_index] > stay + 1e-12:
                labels[node] = best_label
                community_total[current_label] -= ki
                community_total[best_label] += ki
                moved += 1
        if moved == 0:
            break
    _, compact = np.unique(labels, return_inverse=True)
    return compact.astype(np.int64)


def _louvain(adjacency, resolution=1.0, max_levels=10):
    """紧凑版 Louvain（多层级聚合 + 按簇大小降序重排），返回 0..k-1 标签。"""
    adjacency = sp.csr_matrix(adjacency, dtype=np.float64)
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    adjacency = adjacency.maximum(adjacency.T).tocsr()
    n_obs = int(adjacency.shape[0])
    if n_obs == 0:
        return np.zeros(0, dtype=np.int64)
    mapping = np.arange(n_obs, dtype=np.int64)
    current = adjacency
    for _ in range(int(max_levels)):
        labels = _louvain_local_moving(current, resolution)
        count = int(labels.max()) + 1 if labels.size else 0
        if count <= 1:
            mapping = np.zeros(n_obs, dtype=np.int64)
            break
        if count == current.shape[0]:
            mapping = labels[mapping]
            break
        mapping = labels[mapping]
        indicator = sp.csr_matrix(
            (np.ones(current.shape[0], dtype=np.float64),
             (labels, np.arange(current.shape[0], dtype=np.int64))),
            shape=(count, current.shape[0]),
        )
        current = (indicator @ current @ indicator.T).tocsr()
        if current.shape[0] <= 2:
            break
    return _relabel_by_size(mapping)


def _relabel_by_size(labels):
    labels = np.asarray(labels, dtype=np.int64)
    if labels.size == 0:
        return labels
    unique, counts = np.unique(labels, return_counts=True)
    order = np.argsort(-counts, kind="stable")
    remap = np.zeros(int(unique.max()) + 1, dtype=np.int64)
    for new_id, old_index in enumerate(order):
        remap[unique[old_index]] = new_id
    return remap[labels]


def _merge_small_communities(adjacency, labels, min_size):
    """把小于 min_size 的簇并入连边权重最大的邻居簇。"""
    labels = np.asarray(labels, dtype=np.int64).copy()
    adjacency = sp.csr_matrix(adjacency)
    min_size = int(max(2, min_size))
    for _ in range(6):
        unique, counts = np.unique(labels, return_counts=True)
        small = [int(label) for label, count in zip(unique, counts) if count < min_size]
        if not small:
            break
        largest = int(unique[int(np.argmax(counts))])
        changed = False
        for label in small:
            members = np.where(labels == label)[0]
            weights = np.asarray(adjacency[members].sum(axis=0)).ravel()
            other = labels != label
            best_label = largest
            best_weight = -1.0
            for candidate in np.unique(labels[other]):
                candidate = int(candidate)
                mask = (labels == candidate) & other
                if not mask.any():
                    continue
                weight = float(weights[mask].sum())
                if weight > best_weight:
                    best_weight = weight
                    best_label = candidate
            labels[members] = best_label
            changed = True
        if not changed:
            break
    return _relabel_by_size(labels)


def _pca_representation(dense_matrix, n_components, seed=0):
    from sklearn.decomposition import PCA

    dense_matrix = np.asarray(dense_matrix, dtype=np.float64)
    n_obs, n_vars = dense_matrix.shape
    n_components = int(max(2, min(int(n_components), n_obs - 1, n_vars)))
    model = PCA(n_components=n_components, svd_solver="randomized", random_state=seed)
    representation = model.fit_transform(dense_matrix)
    return representation, np.asarray(model.explained_variance_ratio_, dtype=np.float64)


def _embed_2d(representation, seed=0):
    """2D t-SNE（样本过少时退化为 PCA 前两维）。"""
    from sklearn.manifold import TSNE

    representation = np.asarray(representation, dtype=np.float64)
    n_obs = representation.shape[0]
    if n_obs < 6 or representation.shape[1] < 2:
        padded = np.zeros((n_obs, 2))
        padded[:, : min(2, representation.shape[1])] = representation[:, : min(2, representation.shape[1])]
        return padded
    perplexity = float(min(max(2.0, (n_obs - 1) / 3.0), 30.0))
    try:
        model = TSNE(n_components=2, perplexity=perplexity, init="pca",
                     random_state=seed, max_iter=750)
        return np.asarray(model.fit_transform(representation), dtype=np.float64)
    except TypeError:  # pragma: no cover - 旧版 sklearn
        model = TSNE(n_components=2, perplexity=perplexity, init="pca",
                     random_state=seed, n_iter=750)
        return np.asarray(model.fit_transform(representation), dtype=np.float64)


def _subcluster_marker_matrix(matrix, labels, positions, top_per_cluster=5):
    """每个亚簇 top 基因（均值差），返回 (genes, 亚簇×基因 z-score 矩阵)。"""
    levels = sorted({int(label) for label in np.asarray(labels).tolist()})
    selected = []
    seen = set()
    for level in levels:
        rows = np.where(labels == level)[0]
        mask = labels == level
        mean_in, _, _ = _group_stats(matrix, rows)
        mean_out, _, _ = _group_stats(matrix, np.where(~mask)[0])
        diff = mean_in - mean_out
        order = np.argsort(-np.asarray(diff), kind="stable")
        for index in order[: int(top_per_cluster)]:
            index = int(index)
            if index not in seen:
                seen.add(index)
                selected.append(index)
    selected = selected[:60]
    if not selected:
        return [], np.zeros((len(levels), 0))
    sub = matrix[:, selected].toarray()
    rows = []
    for level in levels:
        mask = labels == level
        rows.append(sub[mask].mean(axis=0))
    values = np.asarray(rows, dtype=np.float64)
    sd = values.std(axis=0)
    sd[sd <= 0] = 1.0
    values = (values - values.mean(axis=0)) / sd
    genes = []
    names = _var_names  # placeholder 防止误用
    return selected, np.clip(values, -2.5, 2.5)


def _run_subcluster(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    parent_labels = _cluster_labels(ctx)
    parent_levels = _ordered_labels(parent_labels.tolist())
    target = str(params.get("cluster") or "").strip()
    if not target:
        raise EngineError("bad_request", "请指定要细分的簇（params.cluster，例如 0）")
    if target not in parent_levels:
        raise EngineError(
            "prerequisite",
            "簇 %s 不存在（现有簇：%s）" % (target, "、".join(parent_levels[:20])),
        )
    parent_mask = parent_labels == target
    rows = np.where(parent_mask)[0]
    n_parent = int(rows.size)
    if n_parent < 30:
        raise EngineError("prerequisite", "簇 %s 只有 %d 个细胞，太少无法细分（至少 30）" % (target, n_parent))
    if n_parent > 8000:
        keep = np.sort(_rng().choice(n_parent, size=8000, replace=False))
        rows = rows[keep]
        warnings.append("目标簇细胞数 %d 超过上限 8000，已随机下采样" % n_parent)
    resolution = float(params.get("resolution", 0.8) or 0.8)
    resolution = min(max(resolution, 0.05), 5.0)
    n_neighbors = int(params.get("n_neighbors", 15) or 15)
    n_neighbors = min(max(n_neighbors, 3), max(3, rows.size - 1))
    n_top_genes = int(params.get("n_top_genes", 1000) or 1000)
    n_top_genes = min(max(n_top_genes, 50), 2000, int(matrix.shape[1]))
    n_pcs = int(params.get("n_pcs", 20) or 20)
    min_cells = int(params.get("min_cells", 5) or 5)

    subset = matrix[rows]
    variance = _gene_variance(subset, np.arange(subset.shape[0]), _rng())
    hvg_order = np.argsort(-variance, kind="stable")[:n_top_genes]
    hvg = np.sort(hvg_order)
    dense = subset[:, hvg].toarray().astype(np.float64)
    mean = dense.mean(axis=0)
    sd = dense.std(axis=0)
    sd[sd <= 0] = 1.0
    dense = (dense - mean) / sd
    representation, variance_ratio = _pca_representation(dense, n_pcs, seed=0)
    from sklearn.neighbors import NearestNeighbors

    neighbours = NearestNeighbors(n_neighbors=int(n_neighbors) + 1, metric="euclidean")
    neighbours.fit(representation)
    graph = neighbours.kneighbors_graph(representation, mode="connectivity")
    graph = graph.maximum(graph.T).tocsr()
    raw_labels = _louvain(graph, resolution=resolution)
    labels = _merge_small_communities(graph, raw_labels, min_cells)
    n_sub = int(labels.max()) + 1 if labels.size else 0
    if n_sub < 2:
        raise EngineError(
            "prerequisite",
            "该簇在当前分辨率（%.2f）下没有分出新亚型；请提高 resolution 或检查簇内异质性" % resolution,
        )

    subcluster_names = np.array(["%s_%d" % (target, int(label)) for label in labels], dtype=object)
    full_values = [""] * _n_obs(ctx)
    for local_index, row_index in enumerate(rows):
        full_values[int(row_index)] = subcluster_names[local_index]
    obs_columns = []
    obs_values = {}
    entry, values_map = _obs_pair(ctx, "subcluster", "categorical", full_values)
    obs_columns.append(entry)
    obs_values.update(values_map)

    # 子集嵌入图
    embedding = _embed_2d(representation, seed=0)
    embedding_kind = "t-SNE"
    if rows.size > 3000:
        embedding = representation[:, :2]
        embedding_kind = "PCA (2D)"
        warnings.append("目标簇细胞数较多，子集嵌入图改用 PCA 前两维")
    palette = _palette(n_sub)
    figure, axis = _new_figure(6.6, 5.2)
    for label in range(n_sub):
        mask = labels == label
        axis.scatter(
            embedding[mask, 0], embedding[mask, 1], s=14, alpha=0.75, linewidths=0,
            color=palette[label], label="%s_%d (n=%d)" % (target, label, int(mask.sum())),
        )
        axis.text(
            float(np.median(embedding[mask, 0])), float(np.median(embedding[mask, 1])),
            "%s_%d" % (target, label), fontsize=8, weight="bold", ha="center", va="center",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.85},
        )
    axis.set_title("Subclusters of cluster %s (%s)" % (target, _fig_label(target, "target")))
    axis.set_xlabel("%s 1" % embedding_kind)
    axis.set_ylabel("%s 2" % embedding_kind)
    axis.legend(frameon=False, loc="best", fontsize=7, markerscale=1.4)
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "subcluster", "subcluster_embedding",
                            "Subcluster embedding")]

    # 亚簇 marker 热图
    selected_columns, zmatrix = _subcluster_marker_matrix(
        subset, labels, np.arange(hvg.size), top_per_cluster=5)
    marker_rows = []
    if selected_columns and zmatrix.size:
        genes = [_var_names(ctx)[int(column)] for column in selected_columns]
        levels = sorted({int(label) for label in np.asarray(labels).tolist()})
        figure, axis = _new_figure(8.6, 0.42 * len(levels) + 2.6)
        image = axis.imshow(zmatrix, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        axis.set_yticks(np.arange(len(levels)))
        axis.set_yticklabels([_fig_label(level, "S%d" % index) for index, level in enumerate(levels)])
        axis.set_xticks(np.arange(len(genes)))
        axis.set_xticklabels(genes, rotation=90, fontsize=7)
        axis.set_title("Top genes per subcluster (mean z-score)")
        figure.colorbar(image, ax=axis, fraction=0.024, pad=0.01, label="z-score")
        figure.tight_layout()
        figures.append(_save_figure(figure, ctx, outdir, "subcluster", "subcluster_markers",
                                    "Subcluster markers"))
        for level in levels:
            level_index = levels.index(level)
            mask = labels == level
            mean_in, _, _ = _group_stats(subset, np.where(mask)[0])
            mean_out, _, _ = _group_stats(subset, np.where(~mask)[0])
            diff = mean_in - mean_out
            order = np.argsort(-np.asarray(diff), kind="stable")[:10]
            for rank, column in enumerate(order, start=1):
                marker_rows.append({
                    "subcluster": "%s_%d" % (target, int(level)),
                    "rank": rank,
                    "gene": _var_names(ctx)[int(column)],
                    "mean_in": _num(mean_in[int(column)]),
                    "mean_out": _num(mean_out[int(column)]),
                    "diff": _num(diff[int(column)]),
                })
    counts_rows = [
        {"subcluster": "%s_%d" % (target, int(label)),
         "n_cells": int(np.sum(labels == label)),
         "pct_of_parent": _num(100.0 * int(np.sum(labels == label)) / max(rows.size, 1))}
        for label in range(n_sub)
    ]
    tables = [
        _table(
            "subcluster_sizes", "Subcluster sizes",
            [_col("subcluster", "Subcluster"), _col("n_cells", "Cells", "number"),
             _col("pct_of_parent", "% of parent", "number")],
            counts_rows,
        ),
    ]
    if marker_rows:
        tables.append(_table(
            "subcluster_markers", "Top genes per subcluster",
            [_col("subcluster", "Subcluster"), _col("rank", "Rank", "number"),
             _col("gene", "Gene"), _col("mean_in", "Mean in", "number"),
             _col("mean_out", "Mean out", "number"), _col("diff", "Difference", "number")],
            _cap_rows(marker_rows, 400, warnings, "亚簇 marker 表"),
        ))
    bullets = [
        "簇 %s 共 %d 个细胞，分辨率 %.2f 下细分出 %d 个亚簇" % (target, rows.size, resolution, n_sub),
        "亚簇命名 `%s_<n>`（按亚簇大小降序），已写回 obs 列 subcluster" % target,
        "子集嵌入使用 %s；HVG %d 个，PCA %d 维，kNN k=%d" % (
            embedding_kind, hvg.size, representation.shape[1], n_neighbors),
    ]
    metrics = {
        "parent_cluster": target,
        "n_parent_cells": int(rows.size),
        "n_subclusters": n_sub,
        "resolution": _num(resolution),
        "n_hvg": int(hvg.size),
        "n_pcs": int(representation.shape[1]),
        "pc1_variance_ratio": _num(variance_ratio[0] if variance_ratio.size else 0.0),
    }
    return _result("subcluster", bullets, metrics, tables, figures,
                   obs_columns=obs_columns, obs_values=obs_values, warnings=warnings)


# ===========================================================================
# 7. 组间对比：proportion_test / de_group / volcano
# ===========================================================================

_DE_CACHE = {}
_DE_CACHE_LIMIT = 6


def _obs_signature(ctx):
    obs = _ctx_obs(ctx)
    return (int(obs.shape[0]), tuple(str(column) for column in obs.columns))


def _de_cache_key(ctx, params):
    return (
        str(_ctx_get(ctx, "dataset_id", "ds")),
        _obs_signature(ctx),
        str(params.get("group") or "group"),
        str(params.get("cell_type") or ""),
        str(params.get("group1") or ""),
        str(params.get("group2") or ""),
        int(params.get("min_cells", 10) or 10),
    )


def _de_cache_put(key, payload):
    if key in _DE_CACHE:
        _DE_CACHE.pop(key, None)
    _DE_CACHE[key] = payload
    while len(_DE_CACHE) > _DE_CACHE_LIMIT:
        _DE_CACHE.pop(next(iter(_DE_CACHE)))


def _compute_de(ctx, params, warnings):
    """逐细胞类型 Welch t + BH；结果按 (group,cell_type,groups,min_cells) 缓存供 volcano 复用。"""
    cache_key = _de_cache_key(ctx, params)
    cached = _DE_CACHE.get(cache_key)
    if isinstance(cached, dict):
        warnings.append("复用 de_group 缓存结果（同分组/细胞类型/参数）")
        return cached
    matrix = _expression(ctx, warnings)
    group_column = str(params.get("group") or "group")
    values, levels, masks = _group_labels(ctx, group_column, warnings, role="分组列")
    if len(levels) < 2:
        raise EngineError("prerequisite", "分组列 %s 至少需要 2 个水平" % group_column)
    group1 = str(params.get("group1") or "").strip()
    group2 = str(params.get("group2") or "").strip()
    if not group1 and not group2:
        group1, group2 = levels[0], levels[1]
    elif not group1:
        group1 = next((label for label in levels if label != group2), "")
    elif not group2:
        group2 = next((label for label in levels if label != group1), "")
    for label in (group1, group2):
        if label not in masks:
            raise EngineError(
                "bad_request", "组 %s 不存在（可用组：%s）" % (label, "、".join(levels[:20])))
    cell_values, cell_levels, cell_masks = _cell_type_labels(ctx, params, warnings)
    min_cells = int(params.get("min_cells", 10) or 10)
    genes = _var_names(ctx)
    per_type = {}
    skipped = []
    for cell_type in cell_levels:
        in_group = cell_masks[cell_type] & masks[group1]
        out_group = cell_masks[cell_type] & masks[group2]
        n_in, n_out = int(in_group.sum()), int(out_group.sum())
        if n_in < min_cells or n_out < min_cells:
            skipped.append((cell_type, n_in, n_out))
            continue
        rows_in = np.where(in_group)[0]
        rows_out = np.where(out_group)[0]
        mean_in, var_in, _ = _group_stats(matrix, rows_in)
        mean_out, var_out, _ = _group_stats(matrix, rows_out)
        t_stat, _, pvalue = _welch_t(mean_in, var_in, n_in, mean_out, var_out, n_out)
        pvalue = np.asarray(pvalue, dtype=np.float64)
        padj = _bh_fdr(pvalue)
        pct_in = _pct_expressed(matrix, rows_in)
        pct_out = _pct_expressed(matrix, rows_out)
        per_type[cell_type] = {
            "n_in": n_in, "n_out": n_out,
            "mean_in": mean_in, "mean_out": mean_out,
            "logfc": mean_in - mean_out,
            "t_stat": np.asarray(t_stat, dtype=np.float64),
            "pvalue": pvalue, "padj": np.asarray(padj, dtype=np.float64),
            "pct_in": pct_in, "pct_out": pct_out,
        }
    if not per_type:
        raise EngineError(
            "prerequisite",
            "所有细胞类型在两个组内的细胞数都小于 min_cells=%d，无法做组间差异分析" % min_cells,
        )
    if skipped:
        warnings.append(
            "以下细胞类型因组内细胞数不足未检验：" + "、".join("%s(n=%d/%d)" % item for item in skipped[:8])
        )
    payload = {
        "matrix": matrix,
        "group_column": group_column,
        "group1": group1,
        "group2": group2,
        "levels": levels,
        "cell_type_column": "cell_type" if "cell_type" in _ctx_obs(ctx).columns else "clusters",
        "cell_types": list(per_type.keys()),
        "per_type": per_type,
        "genes": genes,
        "skipped": skipped,
        "min_cells": min_cells,
    }
    _de_cache_put(cache_key, payload)
    return payload


def _de_rows(data, cell_types, top_n=50):
    """把缓存结果转成表格行（每类按 p 值取前 top_n 个基因）。"""
    rows = []
    genes = data["genes"]
    for cell_type in cell_types:
        record = data["per_type"][cell_type]
        padj = record["padj"]
        order = np.argsort(padj, kind="stable")[: int(top_n)]
        for index in order:
            index = int(index)
            rows.append({
                "cell_type": cell_type,
                "gene": genes[index],
                "group1": data["group1"],
                "group2": data["group2"],
                "mean_in": _num(record["mean_in"][index]),
                "mean_out": _num(record["mean_out"][index]),
                "logfc": _num(record["logfc"][index]),
                "t_stat": _num(record["t_stat"][index]),
                "pval": _num(record["pvalue"][index]),
                "pval_adj": _num(padj[index]),
                "pct_in": _num(record["pct_in"][index]),
                "pct_out": _num(record["pct_out"][index]),
                "n_in": record["n_in"],
                "n_out": record["n_out"],
                "direction": "up" if record["logfc"][index] > 0 else "down",
            })
    rows.sort(key=lambda row: (row["pval_adj"] if row["pval_adj"] is not None else 1.0,
                               row["cell_type"], row["gene"]))
    return rows


def _de_summary_rows(data, padj_cutoff=0.05):
    rows = []
    for cell_type in data["cell_types"]:
        record = data["per_type"][cell_type]
        padj = record["padj"]
        significant = int(np.sum(padj < padj_cutoff))
        best = int(np.argmin(padj))
        up = int(np.sum((padj < padj_cutoff) & (record["logfc"] > 0)))
        down = int(np.sum((padj < padj_cutoff) & (record["logfc"] < 0)))
        rows.append({
            "cell_type": cell_type,
            "n_in": record["n_in"],
            "n_out": record["n_out"],
            "n_tested": int(len(padj)),
            "n_sig": significant,
            "n_up": up,
            "n_down": down,
            "top_gene": data["genes"][best],
            "top_pval_adj": _num(padj[best]),
        })
    return rows


def _draw_volcano(axis, record, title, padj_cutoff, logfc_cutoff, top_labels, label_prefix="",
                  fontsize=8):
    """单细胞类型火山图面板。"""
    logfc = np.asarray(record["logfc"], dtype=np.float64)
    padj = np.asarray(record["padj"], dtype=np.float64)
    yvalues = np.array([_minus_log10(value) for value in padj], dtype=np.float64)
    up = (padj < padj_cutoff) & (logfc > logfc_cutoff)
    down = (padj < padj_cutoff) & (logfc < -logfc_cutoff)
    neutral = ~(up | down)
    axis.scatter(logfc[neutral], yvalues[neutral], s=5, color="#c7c7c7", alpha=0.65, linewidths=0)
    axis.scatter(logfc[up], yvalues[up], s=8, color="#e15759", alpha=0.85, linewidths=0, label="up")
    axis.scatter(logfc[down], yvalues[down], s=8, color="#4e79a7", alpha=0.85, linewidths=0, label="down")
    axis.axhline(_minus_log10(padj_cutoff), color="#888888", linestyle="--", linewidth=0.7)
    axis.axvline(logfc_cutoff, color="#888888", linestyle="--", linewidth=0.7)
    axis.axvline(-logfc_cutoff, color="#888888", linestyle="--", linewidth=0.7)
    # 标注 top 基因，并做简单的贪心间距避让（避免文字叠在一起）
    order = np.argsort(padj, kind="stable")
    x_range = max(float(np.max(logfc) - np.min(logfc)) if logfc.size else 1.0, 1e-9)
    y_range = max(float(np.max(yvalues) - np.min(yvalues)) if yvalues.size else 1.0, 1e-9)
    placed = []
    labeled = 0
    gene_names = record.get("genes") or []
    for index in order:
        index = int(index)
        if labeled >= int(top_labels):
            break
        if padj[index] > 0.5:
            continue
        x_value = float(logfc[index])
        y_value = float(yvalues[index])
        if any(abs(x_value - previous_x) < 0.10 * x_range
               and abs(y_value - previous_y) < 0.09 * y_range
               for previous_x, previous_y in placed):
            continue
        label = str(gene_names[index]) if index < len(gene_names) else ""
        if not label:
            continue
        if x_value < 0:
            axis.annotate(label, (x_value, y_value), fontsize=fontsize, ha="right",
                          xytext=(-4, 3), textcoords="offset points", color="#222222")
        else:
            axis.annotate(label, (x_value, y_value), fontsize=fontsize, ha="left",
                          xytext=(4, 3), textcoords="offset points", color="#222222")
        placed.append((x_value, y_value))
        labeled += 1
    axis.set_xlabel("logFC (mean difference, log-normalized)")
    axis.set_ylabel("-log10(FDR)")
    axis.set_title(title, fontsize=9)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def _run_de_group(ctx, params, outdir):
    warnings = []
    data = _compute_de(ctx, params, warnings)
    top_n = int(params.get("top_n", 50) or 50)
    top_n = min(max(top_n, 1), 500)
    rows = _de_rows(data, data["cell_types"], top_n=top_n)
    rows = _cap_rows(rows, 800, warnings, "差异基因表")
    tables = [
        _table(
            "de_genes", "Differential genes per cell type (Welch t + BH)",
            [_col("cell_type", "Cell type"), _col("gene", "Gene"),
             _col("group1", "Group 1"), _col("group2", "Group 2"),
             _col("mean_in", "Mean 1", "number"), _col("mean_out", "Mean 2", "number"),
             _col("logfc", "logFC", "number"), _col("t_stat", "t", "number"),
             _col("pval", "p", "number"), _col("pval_adj", "FDR", "number"),
             _col("pct_in", "pct in 1", "number"), _col("pct_out", "pct in 2", "number"),
             _col("direction", "Direction")],
            rows,
            note="logFC = group1 - group2 的 log 标准化均值差；p 值为双侧 Welch t。",
        ),
        _table(
            "de_summary", "Differential expression summary per cell type",
            [_col("cell_type", "Cell type"), _col("n_in", "n group1", "number"),
             _col("n_out", "n group2", "number"), _col("n_tested", "Genes tested", "number"),
             _col("n_sig", "FDR<0.05", "number"), _col("n_up", "Up", "number"),
             _col("n_down", "Down", "number"), _col("top_gene", "Top gene"),
             _col("top_pval_adj", "Top FDR", "number")],
            _de_summary_rows(data),
        ),
    ]
    panels = []
    for cell_type in data["cell_types"][:6]:
        record = dict(data["per_type"][cell_type])
        record["genes"] = data["genes"]
        panels.append((cell_type, record))
    if not panels:
        raise EngineError("prerequisite", "没有可绘制差异结果的细胞类型")
    columns = 2 if len(panels) > 1 else 1
    rows_count = int(math.ceil(len(panels) / columns))
    figure, axes_grid = _new_figure_grid(rows_count, columns, 4.2 * columns + 1.0, 3.5 * rows_count)
    axes_array = axes_grid.ravel()
    for axis, (cell_type, record) in zip(axes_array, panels):
        _draw_volcano(
            axis, record, _fig_label(cell_type, "cell type"), 0.05, 0.25, 8,
        )
    for axis in axes_array[len(panels):]:
        axis.set_visible(False)
    figure.suptitle("%s vs %s: differential genes per cell type" % (
        _fig_label(data["group1"], "group1"), _fig_label(data["group2"], "group2")), fontsize=11)
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "de_group", "de_volcano",
                            "Differential genes per cell type")]
    bullets = [
        "分组：%s vs %s（列 %s）；检验细胞类型：%s" % (
            data["group1"], data["group2"], data["group_column"], "、".join(data["cell_types"][:10])),
        "逐细胞类型 Welch t + BH FDR；表格含每类 top %d 差异基因与汇总" % top_n,
        "火山图展示前 %d 个细胞类型（按细胞数顺序）" % min(len(panels), 6),
    ]
    metrics = {
        "n_cell_types": len(data["cell_types"]),
        "n_skipped": len(data["skipped"]),
        "n_genes_tested": len(data["genes"]),
        "group1": data["group1"],
        "group2": data["group2"],
    }
    return _result("de_group", bullets, metrics, tables, figures, warnings=warnings)


def _run_volcano(ctx, params, outdir):
    warnings = []
    data = _compute_de(ctx, params, warnings)
    padj_cutoff = float(params.get("padj_cutoff", 0.05) or 0.05)
    logfc_cutoff = float(params.get("logfc_cutoff", 0.25) or 0.0)
    top_labels = int(params.get("top_labels", 10) or 10)
    requested = str(params.get("cell_type") or "").strip()
    obs_columns = {str(column) for column in _ctx_obs(ctx).columns}
    if requested in obs_columns:  # 列名（如默认 cell_type / clusters）表示自动选全部
        requested = ""
    cell_types = data["cell_types"]
    if requested:
        if requested not in data["per_type"]:
            raise EngineError(
                "bad_request", "细胞类型 %s 没有结果（可用：%s）" % (requested, "、".join(cell_types[:20])))
        cell_types = [requested]
    selected = cell_types[:6]
    panels = []
    significant_rows = []
    for cell_type in cell_types:
        record = dict(data["per_type"][cell_type])
        record["genes"] = data["genes"]
        mask = (record["padj"] < padj_cutoff) & (np.abs(record["logfc"]) > logfc_cutoff)
        order = np.argsort(record["padj"], kind="stable")
        for index in order:
            index = int(index)
            if not mask[index]:
                continue
            significant_rows.append({
                "cell_type": cell_type,
                "gene": data["genes"][index],
                "logfc": _num(record["logfc"][index]),
                "pval": _num(record["pvalue"][index]),
                "pval_adj": _num(record["padj"][index]),
                "pct_in": _num(record["pct_in"][index]),
                "pct_out": _num(record["pct_out"][index]),
                "direction": "up" if record["logfc"][index] > 0 else "down",
            })
    for cell_type in selected:
        record = dict(data["per_type"][cell_type])
        record["genes"] = data["genes"]
        panels.append((cell_type, record))
    significant_rows.sort(key=lambda row: (row["pval_adj"] if row["pval_adj"] is not None else 1.0,
                                           row["cell_type"], row["gene"]))
    significant_rows = _cap_rows(significant_rows, 400, warnings, "显著基因表")
    tables = [
        _table(
            "significant_genes", "Significant genes",
            [_col("cell_type", "Cell type"), _col("gene", "Gene"),
             _col("logfc", "logFC", "number"), _col("pval", "p", "number"),
             _col("pval_adj", "FDR", "number"), _col("pct_in", "pct in 1", "number"),
             _col("pct_out", "pct in 2", "number"), _col("direction", "Direction")],
            significant_rows,
            note="阈值：FDR < %g 且 |logFC| > %g。" % (padj_cutoff, logfc_cutoff),
        ),
        _table(
            "de_summary", "Differential expression summary per cell type",
            [_col("cell_type", "Cell type"), _col("n_in", "n group1", "number"),
             _col("n_out", "n group2", "number"), _col("n_sig", "FDR<%.2g" % padj_cutoff, "number"),
             _col("n_up", "Up", "number"), _col("n_down", "Down", "number"),
             _col("top_gene", "Top gene")],
            [
                {"cell_type": row["cell_type"], "n_in": row["n_in"], "n_out": row["n_out"],
                 "n_sig": row["n_sig"], "n_up": row["n_up"], "n_down": row["n_down"],
                 "top_gene": row["top_gene"]}
                for row in _de_summary_rows(data, padj_cutoff=padj_cutoff)
            ],
        ),
    ]
    if not panels:
        raise EngineError("prerequisite", "没有可绘制的细胞类型差异结果")
    columns = 2 if len(panels) > 1 else 1
    rows_count = int(math.ceil(len(panels) / columns))
    figure, axes_grid = _new_figure_grid(rows_count, columns, 4.3 * columns + 1.0, 3.6 * rows_count)
    axes_array = axes_grid.ravel()
    for axis, (cell_type, record) in zip(axes_array, panels):
        _draw_volcano(axis, record, _fig_label(cell_type, "cell type"),
                      padj_cutoff, logfc_cutoff, top_labels)
    for axis in axes_array[len(panels):]:
        axis.set_visible(False)
    figure.suptitle("Volcano: %s vs %s" % (
        _fig_label(data["group1"], "group1"), _fig_label(data["group2"], "group2")), fontsize=11)
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "volcano", "volcano", "Volcano plot")]
    bullets = [
        "分组 %s vs %s；显著基因（FDR<%g，|logFC|>%g）共 %d 个" % (
            data["group1"], data["group2"], padj_cutoff, logfc_cutoff, len(significant_rows)),
        "火山图面板：%s" % "、".join(selected),
        "结果复用 de_group 的 Welch t + BH 统计（同参数自动读缓存）",
    ]
    metrics = {
        "n_significant": len(significant_rows),
        "padj_cutoff": _num(padj_cutoff),
        "logfc_cutoff": _num(logfc_cutoff),
        "group1": data["group1"],
        "group2": data["group2"],
    }
    return _result("volcano", bullets, metrics, tables, figures, warnings=warnings)


# ===========================================================================
# 8. 可视化：proportion_test / proportion_sankey / featureplot
# ===========================================================================


def _limit_levels(levels, masks, limit, warnings, what="分组"):
    limit = int(max(1, limit))
    if len(levels) <= limit:
        return list(levels)
    ordered = sorted(levels, key=lambda label: -int(masks[label].sum()))[:limit]
    keep = set(ordered)
    warnings.append("%s数 %d 超过展示上限 %d，只展示细胞数最多的 %d 个" % (
        what, len(levels), limit, limit))
    return [label for label in levels if label in keep]


def _run_proportion_test(ctx, params, outdir):
    warnings = []
    group_column = str(params.get("group") or "group")
    group_values, group_levels, group_masks = _group_labels(ctx, group_column, warnings, role="分组列")
    cell_type_column, cell_type_series = _obs_column(
        ctx, params.get("cell_type") or "cell_type", role="细胞类型列")
    cell_type_values = np.asarray([str(value) for value in cell_type_series.tolist()], dtype=object)
    cell_type_levels = _ordered_labels(cell_type_values.tolist())
    if len(cell_type_levels) < 2:
        raise EngineError("prerequisite", "细胞类型列 %s 只有一个类别，无法比较比例" % cell_type_column)
    min_cells = int(params.get("min_cells", 10) or 10)
    group_levels = _limit_levels(group_levels, group_masks, 12, warnings, "分组")
    cell_type_levels = [
        level for level in cell_type_levels
        if int(np.sum(cell_type_values == level)) >= min_cells
    ]
    dropped = int(len(_ordered_labels(cell_type_values.tolist()))) - len(cell_type_levels)
    if dropped > 0:
        warnings.append("%d 个细胞类型总细胞数 < %d，已从比例检验中剔除" % (dropped, min_cells))
    if len(cell_type_levels) < 2:
        raise EngineError("prerequisite", "过滤后细胞类型不足 2 个，无法比较比例")
    counts = np.zeros((len(group_levels), len(cell_type_levels)), dtype=np.int64)
    for i, group in enumerate(group_levels):
        for j, cell_type in enumerate(cell_type_levels):
            counts[i, j] = int(np.sum(group_masks[group] & (cell_type_values == cell_type)))
    row_totals = counts.sum(axis=1).astype(np.float64)
    col_totals = counts.sum(axis=0).astype(np.float64)
    grand_total = float(counts.sum())
    expected = np.outer(row_totals, col_totals) / max(grand_total, 1.0)
    if grand_total <= 0:
        raise EngineError("prerequisite", "分组 × 细胞类型列联表为空")
    try:
        chi2, chi_p, dof, expected_full = _stats.chi2_contingency(counts + 0.0, correction=False)
    except Exception:
        chi2, chi_p, dof = 0.0, 1.0, 0
    cramers_v = math.sqrt(chi2 / max(grand_total * (min(counts.shape) - 1), 1e-12)) if chi2 >= 0 else 0.0
    if chi_p < 0.05:
        warnings.append("整体卡方检验显著（p=%.3g），比例分布在组间存在差异" % chi_p)
    low_expected = int(np.sum(expected_full < 5)) if expected_full is not None else 0
    if low_expected:
        warnings.append("有 %d 个格子的期望频数 <5，卡方近似可能不稳（可增大样本）" % low_expected)

    count_rows = []
    for i, group in enumerate(group_levels):
        row = {"group": group, "total_cells": int(row_totals[i])}
        for j, cell_type in enumerate(cell_type_levels):
            row["n_%d" % (j + 1)] = int(counts[i, j])
        count_rows.append(row)
    count_columns = [_col("group", "Group"), _col("total_cells", "Total cells", "number")]
    for j, cell_type in enumerate(cell_type_levels):
        count_columns.append(_col("n_%d" % (j + 1), _fig_label(cell_type, "type %d" % (j + 1)), "number"))

    proportion_rows = []
    pvalues = []
    for j, cell_type in enumerate(cell_type_levels):
        in_type = (cell_type_values == cell_type)
        table = np.zeros((len(group_levels), 2), dtype=np.int64)
        for i, group in enumerate(group_levels):
            inside = int(np.sum(group_masks[group] & in_type))
            table[i] = [inside, int(np.sum(group_masks[group])) - inside]
        try:
            _, type_p, _, _ = _stats.chi2_contingency(table, correction=False)
        except Exception:
            type_p = 1.0
        pvalues.append(_safe_p(type_p))
        standardized = (counts[:, j] - expected[:, j]) / np.sqrt(np.maximum(expected[:, j], 1e-9))
        best_group_index = int(np.argmax(standardized))
        observed = float(counts[best_group_index, j])
        total_in_group = float(row_totals[best_group_index])
        observed_pct = observed / max(total_in_group, 1.0)
        expected_pct = float(col_totals[j]) / max(grand_total, 1.0)
        ratio = _safe_div(observed / max(total_in_group, 1.0), expected_pct, 0.0)
        ci_low, ci_high = _beta_ci(observed, total_in_group)
        direction = "enriched in %s" % group_levels[best_group_index] if standardized[best_group_index] > 0 \
            else "depleted in %s" % group_levels[best_group_index]
        proportion_rows.append({
            "cell_type": cell_type,
            "total_cells": int(col_totals[j]),
            "pct_overall": _num(100.0 * expected_pct),
            "enriched_group": group_levels[best_group_index],
            "pct_in_group": _num(100.0 * observed_pct),
            "log2_obs_exp": _num(math.log2(max(ratio, 1e-9)) if ratio > 0 else -9.0),
            "pct_ci_low": _num(100.0 * ci_low),
            "pct_ci_high": _num(100.0 * ci_high),
            "chi2_p": _num(type_p),
            "pval_adj": None,
            "direction": direction,
        })
    adjusted = _bh_fdr(np.asarray(pvalues, dtype=np.float64))
    for row, value in zip(proportion_rows, adjusted):
        row["pval_adj"] = _num(value)
    proportion_rows.sort(key=lambda row: (row["pval_adj"] if row["pval_adj"] is not None else 1.0,
                                          row["cell_type"]))
    proportion_columns = [
        _col("cell_type", "Cell type"), _col("total_cells", "Cells", "number"),
        _col("pct_overall", "Overall %", "number"), _col("enriched_group", "Enriched group"),
        _col("pct_in_group", "% in group", "number"), _col("log2_obs_exp", "log2(obs/exp)", "number"),
        _col("pct_ci_low", "95% CI low", "number"), _col("pct_ci_high", "95% CI high", "number"),
        _col("chi2_p", "Chi-square p", "number"), _col("pval_adj", "FDR", "number"),
        _col("direction", "Direction"),
    ]
    tables = [
        _table("cell_counts", "Cell counts per group and cell type", count_columns, count_rows),
        _table(
            "proportion_enrichment", "Proportion enrichment per cell type",
            proportion_columns, proportion_rows,
            note="每类做 2×G 卡方检验 + BH；富集组为标准化残差最大的组；区间为 Beta(1,1) 后验 95% 区间。",
        ),
    ]
    # 堆叠比例柱状图
    palette = _palette(len(cell_type_levels))
    figure, axis = _new_figure(max(7.0, 1.4 * len(group_levels) + 3.6), 4.8)
    bottoms = np.zeros(len(group_levels))
    for j, cell_type in enumerate(cell_type_levels):
        proportions = 100.0 * counts[:, j] / np.maximum(row_totals, 1.0)
        axis.bar(
            np.arange(len(group_levels)), proportions, bottom=bottoms, width=0.62,
            color=palette[j], edgecolor="white", linewidth=0.5,
            label=_shorten(_fig_label(cell_type, "type %d" % (j + 1)), 26),
        )
        bottoms += proportions
    axis.set_xticks(np.arange(len(group_levels)))
    axis.set_xticklabels([_fig_label(group, "group %d" % (i + 1)) for i, group in enumerate(group_levels)])
    axis.set_ylabel("Percent of cells (%)")
    axis.set_xlabel("Group")
    axis.set_title("Cell type composition by group (chi-square p=%.3g)" % chi_p)
    axis.set_ylim(0, 100)
    axis.legend(frameon=False, bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=7)
    _style_axis(axis)
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "proportion_test", "proportion_bar",
                            "Cell type composition by group")]
    bullets = [
        "列联表：%d 组 × %d 个细胞类型，整体卡方 p=%.3g（Cramér's V=%.3f）" % (
            len(group_levels), len(cell_type_levels), chi_p, cramers_v),
        "比例显著差异（FDR<0.05）的细胞类型：%d / %d" % (
            int(np.sum(adjusted < 0.05)), len(adjusted)),
        "每组占比与富集方向见 proportion_enrichment 表（含 95%% 后验区间）",
    ]
    metrics = {
        "n_groups": len(group_levels),
        "n_cell_types": len(cell_type_levels),
        "chi2_p": _num(chi_p),
        "cramers_v": _num(cramers_v),
        "n_significant": int(np.sum(adjusted < 0.05)),
    }
    return _result("proportion_test", bullets, metrics, tables, figures, warnings=warnings)


def _alluvial_positions(totals, gap=0.015, top=1.0):
    """把节点总量映射成 [bottom, top] 堆叠区间（归一化到 0..1）。"""
    total = float(sum(totals.values()))
    gap_total = gap * max(0, len(totals) - 1)
    scale = (top - gap_total) / total if total > 0 else 0.0
    positions = {}
    cursor = top
    for label, value in totals.items():
        height = float(value) * scale
        positions[label] = (cursor - height, cursor)
        cursor -= height + gap
    return positions, scale


def _ribbon_patch(x0, x1, y0_top, y0_bottom, y1_top, y1_bottom, color, alpha=0.55):
    xm = 0.5 * (x0 + x1)
    vertices = [
        (x0, y0_top), (xm, y0_top), (xm, y1_top), (x1, y1_top),
        (x1, y1_bottom), (xm, y1_bottom), (xm, y0_bottom), (x0, y0_bottom),
        (x0, y0_top),
    ]
    codes = [
        _MplPath.MOVETO, _MplPath.CURVE4, _MplPath.CURVE4, _MplPath.CURVE4,
        _MplPath.LINETO, _MplPath.CURVE4, _MplPath.CURVE4, _MplPath.CURVE4,
        _MplPath.CLOSEPOLY,
    ]
    patch = _PathPatch(_MplPath(vertices, codes), facecolor=color, edgecolor="none", alpha=alpha)
    return patch


def _run_proportion_sankey(ctx, params, outdir):
    warnings = []
    group_column = str(params.get("group") or "group")
    group_values, group_levels, group_masks = _group_labels(ctx, group_column, warnings, role="分组列")
    cell_type_column, cell_type_series = _obs_column(
        ctx, params.get("cell_type") or "cell_type", role="细胞类型列")
    cell_type_values = np.asarray([str(value) for value in cell_type_series.tolist()], dtype=object)
    cell_type_levels = _ordered_labels(cell_type_values.tolist())
    if len(cell_type_levels) < 2:
        raise EngineError("prerequisite", "细胞类型列 %s 至少需要 2 个类别" % cell_type_column)
    min_cells = int(params.get("min_cells", 5) or 5)
    group_levels = _limit_levels(group_levels, group_masks, 8, warnings, "分组")
    cell_type_levels = [
        level for level in cell_type_levels
        if int(np.sum(cell_type_values == level)) >= min_cells
    ]
    if len(cell_type_levels) < 2:
        raise EngineError("prerequisite", "过滤后细胞类型不足 2 个，无法绘制桑基图")
    if len(cell_type_levels) > 12:
        keep = sorted(cell_type_levels, key=lambda level: -int(np.sum(cell_type_values == level)))[:12]
        keep_set = set(keep)
        warnings.append("细胞类型过多，只展示细胞数最多的 12 个")
        cell_type_levels = [level for level in cell_type_levels if level in keep_set]
    flows = {}
    for group in group_levels:
        for cell_type in cell_type_levels:
            count = int(np.sum(group_masks[group] & (cell_type_values == cell_type)))
            if count >= min_cells:
                flows[(group, cell_type)] = count
    if not flows:
        raise EngineError("prerequisite", "没有达到最小流向细胞数（min_cells=%d）的组-细胞类型组合" % min_cells)
    left_totals = {group: 0 for group in group_levels}
    right_totals = {cell_type: 0 for cell_type in cell_type_levels}
    for (group, cell_type), count in flows.items():
        left_totals[group] += count
        right_totals[cell_type] += count
    left_totals = {key: value for key, value in left_totals.items() if value > 0}
    right_totals = {key: value for key, value in right_totals.items() if value > 0}
    left_positions, scale = _alluvial_positions(left_totals)
    right_positions, _ = _alluvial_positions(right_totals)
    left_cursor = {label: left_positions[label][1] for label in left_totals}
    right_cursor = {label: right_positions[label][1] for label in right_totals}
    group_colors = dict(zip(list(left_totals.keys()), _palette(len(left_totals))))
    cell_colors = dict(zip(list(right_totals.keys()), _palette(len(right_totals))))
    figure, axis = _new_figure(9.4, 5.6)
    for group in left_totals:
        for cell_type in right_totals:
            count = flows.get((group, cell_type))
            if not count:
                continue
            height = count * scale
            left_top = left_cursor[group]
            right_top = right_cursor[cell_type]
            axis.add_patch(_ribbon_patch(
                0.12, 0.88, left_top, left_top - height, right_top, right_top - height,
                group_colors[group], alpha=0.5,
            ))
            left_cursor[group] -= height
            right_cursor[cell_type] -= height
    bar_width = 0.035
    for group, (bottom, top) in left_positions.items():
        axis.add_patch(plt.Rectangle((0.12 - bar_width, bottom), bar_width, top - bottom,
                                     facecolor=group_colors[group], edgecolor="#333333", linewidth=0.5))
        axis.text(0.12 - bar_width - 0.015, 0.5 * (bottom + top),
                  "%s (n=%d)" % (_shorten(_fig_label(group, "group"), 22), left_totals[group]),
                  ha="right", va="center", fontsize=8)
    for cell_type, (bottom, top) in right_positions.items():
        axis.add_patch(plt.Rectangle((0.88, bottom), bar_width, top - bottom,
                                     facecolor=cell_colors[cell_type], edgecolor="#333333", linewidth=0.5))
        axis.text(0.88 + bar_width + 0.015, 0.5 * (bottom + top),
                  "%s (n=%d)" % (_shorten(_fig_label(cell_type, "type"), 24), right_totals[cell_type]),
                  ha="left", va="center", fontsize=8)
    axis.set_xlim(-0.28, 1.28)
    axis.set_ylim(-0.02, 1.02)
    axis.axis("off")
    axis.set_title("Cell type proportions: %s → cell type" % _fig_label(group_column, "group"), fontsize=11)
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "proportion_sankey", "alluvial_flow",
                            "Group to cell type alluvial flow")]
    flow_rows = [
        {"group": group, "cell_type": cell_type, "n_cells": count,
         "pct_of_group": _num(100.0 * count / max(left_totals.get(group, 0), 1)),
         "pct_of_type": _num(100.0 * count / max(right_totals.get(cell_type, 0), 1))}
        for (group, cell_type), count in sorted(
            flows.items(), key=lambda item: (item[0][0], item[0][1]))
    ]
    tables = [_table(
        "flows", "Group to cell type flows",
        [_col("group", "Group"), _col("cell_type", "Cell type"),
         _col("n_cells", "Cells", "number"), _col("pct_of_group", "% of group", "number"),
         _col("pct_of_type", "% of cell type", "number")],
        flow_rows,
    )]
    bullets = [
        "两级桑基图：%s → %s（%d 条流向）" % (group_column, cell_type_column, len(flow_rows)),
        "左侧节点=组（%d 个），右侧节点=细胞类型（%d 个），带宽与细胞数成正比" % (
            len(left_totals), len(right_totals)),
        "小于 min_cells=%d 的流向不绘制（保持图面可读）" % min_cells,
    ]
    metrics = {
        "n_groups": len(left_totals),
        "n_cell_types": len(right_totals),
        "n_flows": len(flow_rows),
        "n_cells": int(sum(flows.values())),
    }
    return _result("proportion_sankey", bullets, metrics, tables, figures, warnings=warnings)


def _parse_gene_list(text, index):
    """解析逗号/空格分隔的基因列表，返回 (positions, found, missing)。"""
    names = [item for item in re.split(r"[,;\s]+", str(text or "")) if item]
    return _gene_positions(index, names)


def _run_featureplot(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    index = _gene_index(ctx)
    max_genes = int(params.get("max_genes", 6) or 6)
    max_genes = min(max(max_genes, 1), 6)
    positions, found, missing = _parse_gene_list(params.get("genes"), index)
    if not positions:
        raise EngineError("bad_request", "基因列表为空或全部不在数据中：%s" % params.get("genes"))
    if len(positions) > max_genes:
        warnings.append("基因数超过上限 %d，只取前 %d 个" % (max_genes, max_genes))
        positions = positions[:max_genes]
        found = found[:max_genes]
    if missing:
        warnings.append("以下基因不在数据中已忽略：%s" % "、".join(missing[:10]))
    kind = str(params.get("kind") or "embedding").lower()
    if kind not in ("embedding", "vln", "dot", "heatmap"):
        raise EngineError("bad_request", "未知图类型 %s（可用 embedding/vln/dot/heatmap）" % kind)
    group_column = str(params.get("groupby") or "clusters")
    try:
        group_values, group_levels, group_masks = _group_labels(ctx, group_column, warnings, role="分组列")
    except EngineError:
        if group_column != "clusters":
            raise
        group_values, group_levels, group_masks = _cluster_labels(ctx), None, None
        group_levels = _ordered_labels(list(group_values))
        group_masks = {label: group_values == label for label in group_levels}
    group_levels = _limit_levels(group_levels, group_masks, 20, warnings, "分组")
    n_obs = int(matrix.shape[0])
    expression = {}
    for position, gene in zip(positions, found):
        expression[gene] = _gene_vector(matrix, position)
    mean_matrix = np.zeros((len(found), len(group_levels)), dtype=np.float64)
    pct_matrix = np.zeros((len(found), len(group_levels)), dtype=np.float64)
    table_rows = []
    for gene_index, gene in enumerate(found):
        values = expression[gene]
        for column_index, level in enumerate(group_levels):
            mask = group_masks[level]
            group_values_expr = values[mask]
            mean_matrix[gene_index, column_index] = float(np.mean(group_values_expr)) if group_values_expr.size else 0.0
            pct_matrix[gene_index, column_index] = float(np.mean(group_values_expr > 0)) if group_values_expr.size else 0.0
            table_rows.append({
                "gene": gene,
                "group": level,
                "n_cells": int(mask.sum()),
                "mean_expr": _num(np.mean(group_values_expr) if group_values_expr.size else 0.0),
                "pct_expr": _num(100.0 * pct_matrix[gene_index, column_index]),
            })
    figures = []
    if kind == "embedding":
        obsm = _ctx_get(ctx, "obsm") or {}
        coordinates = None
        layout = "X_emb"
        if isinstance(obsm, dict) and obsm.get("X_emb") is not None:
            coordinates = np.asarray(obsm["X_emb"], dtype=np.float64)
        elif isinstance(obsm, dict) and obsm.get("X_pca") is not None:
            coordinates = np.asarray(obsm["X_pca"], dtype=np.float64)[:, :2]
            layout = "X_pca"
            warnings.append("没有 X_emb 嵌入，特征图改用 PCA 前两维")
        if coordinates is None:
            raise EngineError("prerequisite", "需要先完成 embed 步骤（obsm 里没有 X_emb/X_pca）")
        if coordinates.shape[0] != n_obs:
            raise EngineError("internal", "嵌入坐标行数与细胞数不一致")
        sample = np.arange(n_obs)
        if n_obs > 6000:
            sample = np.sort(_rng().choice(n_obs, size=6000, replace=False))
            warnings.append("嵌入图点数过多，已随机下采样到 6000 个细胞")
        columns = min(3, len(found))
        rows_count = int(math.ceil(len(found) / columns))
        figure, axes_grid = _new_figure_grid(rows_count, columns, 4.1 * columns + 0.5, 3.7 * rows_count)
        axes_array = axes_grid.ravel()
        for axis, gene in zip(axes_array, found):
            values = expression[gene][sample]
            scatter = axis.scatter(
                coordinates[sample, 0], coordinates[sample, 1],
                c=values, cmap="viridis", s=8, linewidths=0, alpha=0.9,
            )
            axis.set_title(gene, fontsize=10)
            axis.set_xlabel("%s 1" % layout)
            axis.set_ylabel("%s 2" % layout)
            axis.set_xticks([])
            axis.set_yticks([])
            figure.colorbar(scatter, ax=axis, fraction=0.045, pad=0.01)
        for axis in axes_array[len(found):]:
            axis.set_visible(False)
        figure.suptitle("Feature plots (log-normalized expression)", fontsize=11)
        figure.tight_layout()
        figures.append(_save_figure(figure, ctx, outdir, "featureplot", "feature_embedding",
                                    "Feature plots on embedding"))
    elif kind == "vln":
        columns = min(3, len(found))
        rows_count = int(math.ceil(len(found) / columns))
        figure, axes_grid = _new_figure_grid(rows_count, columns, 4.0 * columns + 0.5, 3.3 * rows_count)
        axes_array = axes_grid.ravel()
        for axis, gene in zip(axes_array, found):
            values = expression[gene]
            data = []
            labels_drawn = []
            for level in group_levels:
                mask = group_masks[level]
                selected = values[mask]
                if selected.size > 200:
                    selected = selected[_rng().choice(selected.size, size=200, replace=False)]
                data.append(selected)
                labels_drawn.append(_fig_label(level, "g%d" % len(labels_drawn)))
            if data:
                parts = axis.violinplot(data, showmedians=True, showextrema=False, widths=0.8)
                for body in parts["bodies"]:
                    body.set_facecolor("#4e79a7")
                    body.set_edgecolor("#444444")
                    body.set_alpha(0.7)
                if "cmedians" in parts:
                    parts["cmedians"].set_color("#222222")
            axis.set_xticks(np.arange(1, len(labels_drawn) + 1))
            axis.set_xticklabels(labels_drawn, rotation=45, ha="right", fontsize=7)
            axis.set_title(gene, fontsize=10)
            axis.set_ylabel("Expression")
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
        for axis in axes_array[len(found):]:
            axis.set_visible(False)
        figure.suptitle("Expression by %s (violin)" % _fig_label(group_column, "group"), fontsize=11)
        figure.tight_layout()
        figures.append(_save_figure(figure, ctx, outdir, "featureplot", "feature_violin",
                                    "Gene expression violin plots"))
    elif kind == "dot":
        zvalues = mean_matrix.copy()
        for row in range(zvalues.shape[0]):
            sd = float(np.std(zvalues[row]))
            zvalues[row] = (zvalues[row] - np.mean(zvalues[row])) / (sd if sd > 1e-12 else 1.0)
        figure, axis = _new_figure(max(6.5, 0.55 * len(group_levels) + 3.0), 0.55 * len(found) + 2.4)
        xgrid, ygrid = np.meshgrid(np.arange(len(group_levels)), np.arange(len(found)))
        sizes = 20.0 + 220.0 * pct_matrix
        scatter = axis.scatter(
            xgrid.ravel(), ygrid.ravel(), s=sizes.ravel(), c=zvalues.ravel(),
            cmap="RdBu_r", vmin=-2, vmax=2, edgecolors="#555555", linewidths=0.3,
        )
        axis.set_xticks(np.arange(len(group_levels)))
        axis.set_xticklabels([_fig_label(level, "g%d" % index) for index, level in enumerate(group_levels)],
                             rotation=45, ha="right", fontsize=8)
        axis.set_yticks(np.arange(len(found)))
        axis.set_yticklabels(found, fontsize=8)
        axis.set_title("Dot plot: color=mean z-score, size=pct expressing")
        figure.colorbar(scatter, ax=axis, fraction=0.03, pad=0.02, label="z-score")
        figure.tight_layout()
        figures.append(_save_figure(figure, ctx, outdir, "featureplot", "feature_dot",
                                    "Feature dot plot"))
    else:  # heatmap
        heat = mean_matrix.copy()
        for row in range(heat.shape[0]):
            sd = float(np.std(heat[row]))
            heat[row] = (heat[row] - np.mean(heat[row])) / (sd if sd > 1e-12 else 1.0)
        figure, axis = _new_figure(max(6.5, 0.5 * len(group_levels) + 3.0), 0.55 * len(found) + 2.4)
        image = axis.imshow(heat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        axis.set_xticks(np.arange(len(group_levels)))
        axis.set_xticklabels([_fig_label(level, "g%d" % index) for index, level in enumerate(group_levels)],
                             rotation=45, ha="right", fontsize=8)
        axis.set_yticks(np.arange(len(found)))
        axis.set_yticklabels(found, fontsize=8)
        axis.set_title("Mean expression heatmap (z-score across groups)")
        figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="z-score")
        figure.tight_layout()
        figures.append(_save_figure(figure, ctx, outdir, "featureplot", "feature_heatmap",
                                    "Feature heatmap"))
    tables = [_table(
        "expression_by_group", "Expression by group",
        [_col("gene", "Gene"), _col("group", "Group"), _col("n_cells", "Cells", "number"),
         _col("mean_expr", "Mean expression", "number"), _col("pct_expr", "% expressing", "number")],
        _cap_rows(table_rows, 400, warnings, "基因表达表"),
    )]
    bullets = [
        "基因：%s（数据中命中 %d 个）" % ("、".join(found), len(found)),
        "图类型=%s；分组列=%s（%d 组）" % (kind, group_column, len(group_levels)),
        "表达式为 log 标准化值（缺失基因已忽略）",
    ]
    metrics = {
        "n_genes": len(found),
        "kind": kind,
        "n_groups": len(group_levels),
        "n_cells": n_obs,
    }
    return _result("featureplot", bullets, metrics, tables, figures, warnings=warnings)


# ===========================================================================
# 9. 通路分析：pathway_ora / gsea_prerank
# ===========================================================================


def _cluster_marker_foreground(ctx, matrix, cluster_param, top_per_cluster=5, warnings=None):
    """从各簇 one-vs-rest Welch t 取正向 top 基因作为 ORA 前景（无外部 marker 缓存时的默认）。"""
    labels = _cluster_labels(ctx)
    levels = _ordered_labels(labels.tolist())
    requested = str(cluster_param or "").strip()
    if requested:
        if requested not in levels:
            raise EngineError(
                "prerequisite", "簇 %s 不存在（现有簇：%s）" % (requested, "、".join(levels[:20])))
        levels = [requested]
    gene_names = _var_names(ctx)
    selected = []
    seen = set()
    for level in levels:
        mask = labels == level
        rows_in = np.where(mask)[0]
        rows_out = np.where(~mask)[0]
        if rows_in.size < 5 or rows_out.size < 5:
            continue
        mean_in, var_in, n_in = _group_stats(matrix, rows_in)
        mean_out, var_out, n_out = _group_stats(matrix, rows_out)
        t_stat, _, _ = _welch_t(mean_in, var_in, n_in, mean_out, var_out, n_out)
        order = np.argsort(-np.asarray(t_stat), kind="stable")
        added = 0
        for index in order:
            index = int(index)
            if t_stat[index] <= 0:
                break
            if index in seen:
                continue
            seen.add(index)
            selected.append(gene_names[index])
            added += 1
            if added >= int(top_per_cluster):
                break
    return selected


def _ora_rows(foreground, sets, background, min_set_size):
    """超几何 ORA：返回 (rows, tested_count, n_significant)。"""
    background_set = set(background)
    background_size = len(background_set)
    query = []
    query_seen = set()
    for gene in foreground:
        text = str(gene)
        if text in background_set and text not in query_seen:
            query_seen.add(text)
            query.append(text)
    query_size = len(query)
    if query_size == 0:
        raise EngineError("prerequisite", "前景基因与数据集的基因名没有交集，无法做 ORA")
    pvalues = []
    pending = []
    for term, genes in sets.items():
        member_positions = []
        seen = set()
        for gene in genes:
            text = str(gene)
            if text in background_set and text not in seen:
                seen.add(text)
                member_positions.append(text)
        set_size = len(member_positions)
        if set_size < int(min_set_size):
            continue
        hits = [gene for gene in member_positions if gene in query_seen]
        expected = query_size * set_size / max(background_size, 1)
        try:
            pvalue = float(_stats.hypergeom.sf(len(hits) - 1, background_size, set_size, query_size))
        except Exception:
            pvalue = 1.0
        pvalue = _safe_p(pvalue)
        pvalues.append(pvalue)
        pending.append((term, member_positions, hits, expected, pvalue))
    if not pending:
        raise EngineError(
            "prerequisite",
            "没有基因集满足最小大小 %d（内置集可能与该数据基因符号不匹配）" % int(min_set_size),
        )
    adjusted = _bh_fdr(np.asarray(pvalues, dtype=np.float64))
    rows = []
    for (term, member_positions, hits, expected, pvalue), padj in zip(pending, adjusted):
        if not hits:
            continue
        fold = _safe_div(len(hits), expected, 0.0)
        rows.append({
            "term": term,
            "n_genes_in_data": len(member_positions),
            "n_hits": len(hits),
            "expected": _num(expected),
            "fold_enrichment": _num(fold),
            "pval": _num(pvalue),
            "pval_adj": _num(padj),
            "hits": ", ".join(hits[:15]),
        })
    rows.sort(key=lambda row: (row["pval_adj"] if row["pval_adj"] is not None else 1.0,
                               -(row["fold_enrichment"] or 0.0), row["term"]))
    n_significant = int(np.sum(np.asarray(adjusted) < 0.05))
    return rows, len(pending), n_significant


def _run_pathway_ora(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    index = _gene_index(ctx)
    source = str(params.get("source") or "builtin").strip().lower()
    gene_text = str(params.get("gene_list") or "").strip()
    background = _var_names(ctx)
    foreground = []
    foreground_source = "用户提供的基因列表"
    if gene_text:
        positions, found, missing = _parse_gene_list(gene_text, index)
        if missing:
            warnings.append("%d 个基因不在数据中（如 %s）" % (len(missing), "、".join(missing[:5])))
        foreground = list(found)
    if not foreground:
        foreground = _cluster_marker_foreground(ctx, matrix, params.get("cluster"), 5, warnings)
        foreground_source = "各簇 one-vs-rest top 标记基因"
    if not foreground:
        raise EngineError(
            "prerequisite",
            "没有前景基因：请提供 gene_list，或先完成聚类以生成簇标记基因",
        )
    max_input = int(params.get("top_n", 15) or 15)
    max_input = min(max(max_input, 1), 100)
    if source == "enrichr":
        library = str(params.get("enrichr_library") or "KEGG_2021_Human")
        hits = _enrichr_ora(foreground, library)
        rows = []
        for entry in hits:
            rows.append({
                "term": entry["term"],
                "rank": entry["rank"],
                "n_hits": len(entry["hits"]),
                "combined_score": _num(entry["combined_score"]),
                "pval": _num(entry["pvalue"]),
                "pval_adj": _num(entry["pvalue_adj"]),
                "hits": ", ".join(entry["hits"][:15]),
            })
        rows.sort(key=lambda row: (row["pval_adj"] if row["pval_adj"] is not None else 1.0,
                                   row["rank"]))
        rows = rows[:max(max_input, 15)]
        tables = [_table(
            "enrichment", "Enrichr ORA results (%s)" % library,
            [_col("term", "Term"), _col("rank", "Rank", "number"),
             _col("n_hits", "Overlap", "number"), _col("combined_score", "Combined score", "number"),
             _col("pval", "p", "number"), _col("pval_adj", "FDR", "number"),
             _col("hits", "Overlap genes")],
            rows,
        )]
        top = rows[:15]
        figure, axis = _new_figure(8.2, max(3.2, 0.34 * len(top) + 1.8))
        yvalues = np.arange(len(top))[::-1]
        heights = [_minus_log10(row["pval_adj"]) for row in top]
        axis.barh(yvalues, heights, color="#4e79a7", height=0.68)
        axis.set_yticks(yvalues)
        axis.set_yticklabels([_shorten(row["term"], 38) for row in top], fontsize=7.5)
        axis.set_xlabel("-log10(FDR)")
        axis.set_title("Enrichr ORA (%s)" % library)
        _style_axis(axis, grid_axis="x")
        figure.tight_layout()
        figures = [_save_figure(figure, ctx, outdir, "pathway_ora", "enrichr_bar",
                                "Enrichr ORA bar plot")]
        bullets = [
            "Enrichr 在线 ORA：%s 库，前景基因 %d 个" % (library, len(foreground)),
            "返回 %d 条结果，FDR<0.05 共 %d 条" % (
                len(rows), sum(1 for row in rows if (row["pval_adj"] or 1.0) < 0.05)),
            "在线服务需要联网；失败会给出可读错误并可回退到内置 GMT",
        ]
        return _result("pathway_ora", bullets, {
            "source": "enrichr", "library": library, "n_input_genes": len(foreground),
            "n_results": len(rows),
        }, tables, figures, warnings=warnings)

    sets, source_label = _resolve_gene_sets(params, warnings)
    min_set_size = int(params.get("min_set_size", 3) or 3)
    rows, n_tested, n_significant = _ora_rows(foreground, sets, background, min_set_size)
    rows = rows[: max(max_input, 5)]
    tables = [
        _table(
            "enrichment", "Over-representation (hypergeometric + BH)",
            [_col("term", "Term"), _col("n_genes_in_data", "Set genes in data", "number"),
             _col("n_hits", "Overlap", "number"), _col("expected", "Expected", "number"),
             _col("fold_enrichment", "Fold enrichment", "number"),
             _col("pval", "p", "number"), _col("pval_adj", "FDR", "number"),
             _col("hits", "Overlap genes")],
            rows,
            note="背景 = 数据集全部基因；基因集大小按数据中可匹配的基因数计算。",
        ),
    ]
    top = rows[:15]
    if top:
        figure, axis = _new_figure(8.4, max(3.2, 0.34 * len(top) + 1.8))
        yvalues = np.arange(len(top))[::-1]
        heights = [_minus_log10(row["pval_adj"]) for row in top]
        colors = ["#e15759" if (_num_to_float(row["pval_adj"]) < 0.05) else "#bab0ab" for row in top]
        axis.barh(yvalues, heights, color=colors, height=0.68)
        axis.set_yticks(yvalues)
        axis.set_yticklabels([_shorten(row["term"], 38) for row in top], fontsize=7.5)
        for yvalue, row in zip(yvalues, top):
            axis.text(_minus_log10(row["pval_adj"]) + 0.02, yvalue,
                      "%d/%d" % (row["n_hits"], row["n_genes_in_data"]),
                      va="center", fontsize=6.5, color="#333333")
        axis.set_xlabel("-log10(FDR)")
        axis.set_title("Over-representation analysis (%s)" % source_label)
        _style_axis(axis, grid_axis="x")
        figure.tight_layout()
        figures = [_save_figure(figure, ctx, outdir, "pathway_ora", "ora_bar",
                                "ORA enrichment bar plot")]
    else:
        figures = []
        warnings.append("没有任何基因集与前景基因重叠，未生成条形图")
    bullets = [
        "前景基因 %d 个（来源：%s），其中 %d 个在数据集中" % (
            len(foreground), foreground_source,
            len(set(foreground) & set(background))),
        "检验 %d 个基因集（最小集合大小 %d），FDR<0.05 显著 %d 条" % (
            n_tested, min_set_size, n_significant),
        "来源：%s；超几何检验 + BH FDR" % source_label,
    ]
    metrics = {
        "source": source,
        "n_input_genes": len(foreground),
        "n_genes_in_background": len(background),
        "n_sets_tested": n_tested,
        "n_significant": n_significant,
    }
    return _result("pathway_ora", bullets, metrics, tables, figures, warnings=warnings)


def _num_to_float(value):
    try:
        number = float(value)
    except Exception:
        return 1.0
    return number if math.isfinite(number) else 1.0


def _gene_mean_variance(matrix):
    n_obs = int(matrix.shape[0])
    total = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    total_sq = np.asarray(matrix.multiply(matrix).sum(axis=0)).ravel().astype(np.float64)
    mean = total / max(n_obs, 1)
    variance = (total_sq - total * total / max(n_obs, 1)) / max(n_obs - 1, 1)
    return mean, np.maximum(variance, 0.0)


def _ranking_metric(ctx, matrix, params, warnings):
    """GSEA 预排序指标：PC1 相关 / 归一化离散度 / 自定义。返回 (gene_idx, metric, label)。"""
    gene_names = _var_names(ctx)
    ranking = str(params.get("ranking") or "pca1").strip().lower()
    n_vars = int(matrix.shape[1])
    if ranking == "custom":
        text = str(params.get("gene_scores") or "").strip()
        if not text:
            raise EngineError("bad_request", "ranking=custom 时必须提供 gene_scores（GENE:score）")
        index = _gene_index(ctx)
        entries = {}
        for item in re.split(r"[,;\n]+", text):
            item = item.strip()
            if not item:
                continue
            parts = re.split(r"[:=\t]+", item, maxsplit=1)
            if len(parts) != 2:
                continue
            gene = parts[0].strip()
            try:
                value = float(parts[1])
            except Exception:
                continue
            if not math.isfinite(value):
                continue
            positions, _, _ = _gene_positions(index, [gene])
            if positions:
                entries[int(positions[0])] = value
        if len(entries) < 20:
            raise EngineError("bad_request", "custom 排序至少需要 20 个可识别基因（当前 %d 个）" % len(entries))
        gene_idx = np.array(sorted(entries.keys()), dtype=np.int64)
        metric = np.array([entries[int(index)] for index in gene_idx], dtype=np.float64)
        return gene_idx, metric, "custom gene scores"
    if ranking == "pca1":
        obsm = _ctx_get(ctx, "obsm") or {}
        pca = obsm.get("X_pca") if isinstance(obsm, dict) else None
        n_obs = int(matrix.shape[0])
        if pca is not None and np.asarray(pca).shape[0] == n_obs and np.asarray(pca).shape[1] >= 1:
            pca = np.asarray(pca, dtype=np.float64)
            score = pca[:, 0] - float(np.mean(pca[:, 0]))
            sd_pc = float(np.std(pca[:, 0]))
            if sd_pc <= 0:
                raise EngineError("prerequisite", "PC1 方差为 0，无法用于 GSEA 排序")
            mean, variance = _gene_mean_variance(matrix)
            sd_gene = np.sqrt(variance)
            covariance = np.asarray(matrix.T.dot(score)).ravel() / max(n_obs - 1, 1)
            denominator = sd_gene * sd_pc
            metric = np.zeros_like(covariance)
            good = denominator > 1e-12
            metric[good] = covariance[good] / denominator[good]
            gene_idx = np.arange(n_vars, dtype=np.int64)
            return gene_idx, metric, "Pearson correlation with PC1"
        warnings.append("obsm 没有 X_pca，GSEA 排序改用归一化离散度")
        ranking = "dispersion"
    mean, variance = _gene_mean_variance(matrix)
    log_variance = np.log(np.maximum(variance, 1e-12))
    sd = float(np.std(log_variance))
    if sd <= 0:
        raise EngineError("prerequisite", "基因表达离散度没有变异，无法用于 GSEA 排序")
    metric = (log_variance - float(np.mean(log_variance))) / sd
    gene_idx = np.arange(n_vars, dtype=np.int64)
    return gene_idx, metric, "normalized log dispersion"


def _es_from_positions(positions, weights_by_gene, n_genes):
    """加权 KS 富集分数（observed/置换共用）。返回 (ES, sign, peak_index)。"""
    positions = np.asarray(positions, dtype=np.int64)
    n_set = int(positions.size)
    if n_set == 0 or n_set >= int(n_genes):
        return 0.0, 0, 0
    weights = np.asarray(weights_by_gene[positions], dtype=np.float64)
    total_weight = float(weights.sum())
    if total_weight <= 0:
        weights = np.ones(n_set, dtype=np.float64)
        total_weight = float(n_set)
    increments = weights / total_weight
    decrement = 1.0 / float(int(n_genes) - n_set)
    cumulative = np.cumsum(increments)
    ranks = np.arange(1, n_set + 1, dtype=np.float64)
    running_at = cumulative - ((positions + 1.0) - ranks) * decrement
    running_before = np.empty(n_set, dtype=np.float64)
    running_before[0] = -(float(positions[0])) * decrement
    if n_set > 1:
        running_before[1:] = cumulative[:-1] - (
            positions[1:].astype(np.float64) - np.arange(1, n_set, dtype=np.float64)
        ) * decrement
    positive = float(max(0.0, float(np.max(running_at))))
    negative = float(min(0.0, float(np.min(running_before))))
    if positive >= -negative:
        return positive, 1, int(np.argmax(running_at))
    return negative, -1, int(np.argmin(running_before))


def _es_curve(positions, weights_by_gene, n_genes):
    positions = np.asarray(positions, dtype=np.int64)
    n_set = int(positions.size)
    if n_set == 0 or n_set >= int(n_genes):
        return [0.0, 1.0], [0.0, 0.0]
    weights = np.asarray(weights_by_gene[positions], dtype=np.float64)
    total_weight = float(weights.sum())
    if total_weight <= 0:
        weights = np.ones(n_set)
        total_weight = float(n_set)
    increments = weights / total_weight
    decrement = 1.0 / float(int(n_genes) - n_set)
    xs = [0.0]
    ys = [0.0]
    last = 0
    running = 0.0
    for index, position in enumerate(positions):
        running -= float(int(position) - last) * decrement
        xs.append(float(position) / n_genes)
        ys.append(running)
        running += float(increments[index])
        xs.append(float(position) / n_genes)
        ys.append(running)
        last = int(position) + 1
    running -= float(int(n_genes) - last) * decrement
    xs.append(1.0)
    ys.append(running)
    return xs, ys


def _run_gsea_prerank(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    gene_idx, metric, ranking_label = _ranking_metric(ctx, matrix, params, warnings)
    n_genes = int(gene_idx.size)
    if n_genes < 100:
        raise EngineError("prerequisite", "可用于 GSEA 的基因只有 %d 个（至少 100）" % n_genes)
    ranking = str(params.get("ranking") or "pca1").strip().lower()
    weight = float(params.get("weight", 1.0) or 0.0)
    weight = min(max(weight, 0.0), 2.0)
    n_perm = int(params.get("n_perm", 1000) or 1000)
    n_perm = min(max(n_perm, 10), 10000)
    top_n = int(params.get("top_n", 20) or 20)
    min_set_size = int(params.get("min_set_size", 5) or 5)
    max_set_size = int(params.get("max_set_size", 500) or 500)
    order = np.argsort(-metric, kind="stable")
    abs_weights = np.abs(metric[order]) ** weight
    gene_names = _var_names(ctx)
    ranking_genes = [gene_names[int(gene_idx[int(position)])] for position in order]
    rank_position = np.full(len(gene_names), -1, dtype=np.int64)
    rank_position[gene_idx[order]] = np.arange(n_genes, dtype=np.int64)
    sets_params = dict(params)
    if not sets_params.get("gmt_file") and not sets_params.get("source"):
        sets_params["source"] = "builtin"
    elif sets_params.get("gmt_file") and not sets_params.get("source"):
        sets_params["source"] = "gmt"
    sets, source_label = _resolve_gene_sets(sets_params, warnings)
    tested = []
    for term, genes in sets.items():
        member_positions, _, _ = _gene_positions(_gene_index(ctx), genes)
        rank_positions = np.sort(np.asarray(
            [int(rank_position[position]) for position in member_positions
             if rank_position[position] >= 0], dtype=np.int64))
        if rank_positions.size < min_set_size or rank_positions.size > max_set_size:
            continue
        tested.append((term, rank_positions))
    if not tested:
        raise EngineError(
            "prerequisite",
            "没有基因集落在 [%d, %d] 的基因数范围内（可能与数据基因符号不匹配）"
            % (min_set_size, max_set_size),
        )
    if len(tested) > 100:
        warnings.append("可检验基因集 %d 个超过上限 100，按集合大小取前 100 个" % len(tested))
        tested = sorted(tested, key=lambda item: (-item[1].size, item[0]))[:100]
    observed = []
    for term, positions in tested:
        es, sign, peak = _es_from_positions(positions, abs_weights, n_genes)
        observed.append((es, sign, peak))
    nulls = np.zeros((len(tested), n_perm), dtype=np.float64)
    permutation_rng = _rng()
    for permutation in range(n_perm):
        permuted = abs_weights[permutation_rng.permutation(n_genes)]
        for index, (_, positions) in enumerate(tested):
            es, _, _ = _es_from_positions(positions, permuted, n_genes)
            nulls[index, permutation] = es
    records = []
    pvalues = []
    for index, (term, positions) in enumerate(tested):
        es, sign, peak = observed[index]
        null_values = nulls[index]
        same_sign = null_values[null_values > 0] if sign >= 0 else null_values[null_values < 0]
        if same_sign.size == 0:
            normalized = 0.0
            pvalue = 1.0
        else:
            scale = float(np.mean(np.abs(same_sign)))
            normalized = es / scale if scale > 0 else 0.0
            pvalue = (1.0 + float(np.sum(np.abs(same_sign) >= abs(es)))) / (1.0 + same_sign.size)
        pvalue = _safe_p(pvalue)
        if sign >= 0:
            leading = positions[: peak + 1]
        else:
            leading = positions[peak:]
        leading_names = [ranking_genes[int(position)] for position in leading[:10]]
        peak_fraction = float(positions[peak]) / max(n_genes - 1, 1) if positions.size else 0.0
        records.append({
            "term": term,
            "size": int(positions.size),
            "es": es,
            "nes": normalized,
            "pval": pvalue,
            "direction": "up" if es > 0 else "down",
            "leading_edge": ", ".join(leading_names[:5]),
            "rank_at_peak": _num(100.0 * peak_fraction),
            "_positions": positions,
            "_peak": int(peak),
        })
        pvalues.append(pvalue)
    adjusted = _bh_fdr(np.asarray(pvalues, dtype=np.float64))
    for record, padj in zip(records, adjusted):
        record["pval_adj"] = _safe_p(padj)
    records.sort(key=lambda item: (-abs(item["nes"]), item["term"]))
    n_significant = int(np.sum(np.asarray(adjusted) < 0.05))
    rows = [
        {"term": record["term"], "size": record["size"],
         "es": _num(record["es"]), "nes": _num(record["nes"]),
         "pval": _num(record["pval"]), "pval_adj": _num(record["pval_adj"]),
         "direction": record["direction"], "leading_edge": record["leading_edge"],
         "rank_at_peak": record["rank_at_peak"]}
        for record in records[: max(top_n, 5)]
    ]
    tables = [
        _table(
            "gsea_results", "Pre-ranked GSEA results (weighted KS + permutations)",
            [_col("term", "Term"), _col("size", "Set size", "number"),
             _col("es", "ES", "number"), _col("nes", "NES", "number"),
             _col("pval", "p", "number"), _col("pval_adj", "FDR", "number"),
             _col("direction", "Direction"), _col("leading_edge", "Leading edge"),
             _col("rank_at_peak", "Peak rank %", "number")],
            rows,
            note="排序依据：%s；置换 %d 次，权重指数 %.1f。" % (ranking_label, n_perm, weight),
        ),
        _table(
            "ranking_top", "Ranking metric extremes",
            [_col("gene", "Gene"), _col("rank", "Rank", "number"), _col("metric", "Metric", "number")],
            _cap_rows(
                [{"gene": ranking_genes[index], "rank": index + 1,
                  "metric": _num(float(metric[order[index]]))}
                 for index in list(range(min(10, n_genes)))
                 + list(range(max(0, n_genes - 10), n_genes))],
                40, warnings, "排序基因表"),
        ),
    ]
    top_records = records[:6]
    figures = []
    if top_records:
        columns = 2 if len(top_records) > 1 else 1
        rows_count = int(math.ceil(len(top_records) / columns))
        figure, axes_grid = _new_figure_grid(rows_count, columns, 4.6 * columns + 0.8, 3.2 * rows_count)
        axes_array = axes_grid.ravel()
        for axis, record in zip(axes_array, top_records):
            xs, ys = _es_curve(record["_positions"], abs_weights, n_genes)
            color = "#e15759" if record["es"] > 0 else "#4e79a7"
            axis.plot(xs, ys, color=color, linewidth=1.3)
            axis.axhline(0.0, color="#888888", linewidth=0.6)
            if record["_positions"].size:
                peak_position = int(record["_positions"][record["_peak"]])
                axis.axvline(float(peak_position) / max(n_genes - 1, 1),
                             color="#555555", linestyle="--", linewidth=0.7)
            axis.set_xlabel("Gene rank (fraction)")
            axis.set_ylabel("Running enrichment score")
            axis.set_title("%s\nNES=%.2f, FDR=%.2g" % (
                _shorten(_fig_label(record["term"], "set"), 34),
                record["nes"], _num_to_float(record["pval_adj"])), fontsize=8)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
        for axis in axes_array[len(top_records):]:
            axis.set_visible(False)
        figure.suptitle("GSEA enrichment curves (top %d by |NES|)" % len(top_records), fontsize=11)
        figure.tight_layout()
        figures.append(_save_figure(figure, ctx, outdir, "gsea_prerank", "gsea_curves",
                                    "GSEA enrichment curves"))
    bullets = [
        "排序依据：%s（%d 个基因）；权重指数 %.1f；置换 %d 次" % (ranking_label, n_genes, weight, n_perm),
        "检验 %d 个基因集，FDR<0.05 显著 %d 条；NES 按同号置换零分布标准化" % (len(tested), n_significant),
        "富集曲线图展示 |NES| 最大的 %d 条通路（前 6 条）" % min(len(top_records), 6),
    ]
    metrics = {
        "ranking": ranking,
        "n_genes_ranked": n_genes,
        "n_sets_tested": len(tested),
        "n_significant": n_significant,
        "n_perm": n_perm,
        "weight": _num(weight),
    }
    return _result("gsea_prerank", bullets, metrics, tables, figures, warnings=warnings)


# ===========================================================================
# 10. 基因集合分析：gsva_like / auc_score / cnmf_programs / wgcna_lite / geneset_clinical
# ===========================================================================


def _usable_gene_sets(ctx, params, warnings, min_set_size=5):
    """解析内置/GMT 基因集并映射到数据基因位置；返回 (usable, source_label)。"""
    sets_params = dict(params)
    if params.get("gmt_file") and not params.get("source"):
        sets_params["source"] = "gmt"
    elif not params.get("source"):
        sets_params["source"] = "builtin"
    sets, source_label = _resolve_gene_sets(sets_params, warnings)
    index = _gene_index(ctx)
    usable = []
    for term, genes in sets.items():
        positions, found, _ = _gene_positions(index, genes)
        if len(positions) < int(min_set_size):
            continue
        usable.append((term, np.asarray(positions, dtype=np.int64)))
    if not usable:
        raise EngineError(
            "prerequisite",
            "没有基因集在数据中命中至少 %d 个基因；请检查基因符号或提供自己的 .gmt" % int(min_set_size),
        )
    return usable, source_label


def _rank_ks_scores(matrix, usable, n_vars):
    """每细胞对每个基因集计算秩加权 KS 富集分数（ssGSEA 式简化）。"""
    n_obs = int(matrix.shape[0])
    n_sets = len(usable)
    scores = np.zeros((n_obs, n_sets), dtype=np.float64)
    weight_by_rank = (n_vars - np.arange(n_vars, dtype=np.float64))
    ranks = np.empty(n_vars, dtype=np.int64)
    for cell in range(n_obs):
        row = np.asarray(matrix[cell].todense()).ravel().astype(np.float64)
        order = np.argsort(-row, kind="stable")
        ranks[order] = np.arange(n_vars, dtype=np.int64)
        for set_index, (_, positions) in enumerate(usable):
            set_ranks = np.sort(ranks[positions])
            es, _, _ = _es_from_positions(set_ranks, weight_by_rank, n_vars)
            scores[cell, set_index] = es
    return scores


def _select_sets_by_cluster_variance(ctx, scores, usable, rows, n_sets, warnings):
    """按（簇间均分方差）优先挑出展示用基因集；无簇时按命中基因数。"""
    if len(usable) <= n_sets:
        return list(range(len(usable)))
    clusters = None
    try:
        clusters = _cluster_labels(ctx)
    except EngineError:
        clusters = None
    if clusters is not None and clusters.shape[0] == int(_ctx_get(ctx, "n_obs", clusters.shape[0])):
        sub_clusters = clusters[rows]
        levels = _ordered_labels(sub_clusters.tolist())
        priority = np.zeros(len(usable), dtype=np.float64)
        for set_index in range(len(usable)):
            means = []
            for level in levels:
                mask = sub_clusters == level
                if mask.any():
                    means.append(float(np.mean(scores[mask, set_index])))
            priority[set_index] = float(np.var(means)) if means else 0.0
        order = np.argsort(-priority, kind="stable")[: int(n_sets)]
        warnings.append("按簇间差异挑出 %d 个最可区分的基因集（共 %d 个可用）" % (n_sets, len(usable)))
        return sorted(int(index) for index in order)
    return list(range(int(n_sets)))


def _set_cluster_heatmap(ctx, scores, rows, usable, selected, warnings):
    """返回 (heat 矩阵, 行标签, 列标签)；有簇时按簇均值，否则抽样细胞。"""
    clusters = None
    try:
        clusters = _cluster_labels(ctx)
    except EngineError:
        clusters = None
    set_names = [usable[index][0] for index in selected]
    if clusters is not None and clusters.shape[0] == int(_ctx_get(ctx, "n_obs", clusters.shape[0])):
        sub_clusters = clusters[rows]
        levels = _ordered_labels(sub_clusters.tolist())
        masks = {level: sub_clusters == level for level in levels}
        levels = _limit_levels(levels, masks, 20, warnings, "簇")
        heat = np.zeros((len(selected), len(levels)), dtype=np.float64)
        for row_index, set_index in enumerate(selected):
            for column_index, level in enumerate(levels):
                mask = masks[level]
                heat[row_index, column_index] = float(np.mean(scores[mask, set_index])) if mask.any() else 0.0
        return heat, set_names, [_fig_label(level, "C%d" % index) for index, level in enumerate(levels)], "cluster"
    n_cells = int(scores.shape[0])
    sample = np.arange(n_cells)
    if n_cells > 300:
        sample = np.sort(_rng().choice(n_cells, size=300, replace=False))
    heat = scores[sample][:, selected].T
    return heat, set_names, ["cell %d" % index for index in sample[:0]], "cells"


def _heatmap_figure(heat, row_labels, column_labels, title, cbar_label, symmetric=False):
    n_rows, n_columns = heat.shape
    width = max(6.2, 0.42 * n_columns + 3.4)
    height = max(2.8, 0.36 * n_rows + 2.2)
    figure, axis = _new_figure(width, height)
    if symmetric:
        limit = max(float(np.max(np.abs(heat))) if heat.size else 1.0, 1e-9)
        image = axis.imshow(heat, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    else:
        image = axis.imshow(heat, aspect="auto", cmap="viridis")
    axis.set_yticks(np.arange(n_rows))
    axis.set_yticklabels([_shorten(label, 34) for label in row_labels], fontsize=7.5)
    axis.set_xticks(np.arange(n_columns))
    axis.set_xticklabels([_shorten(label, 22) for label in column_labels],
                         rotation=45, ha="right", fontsize=7.5)
    axis.set_title(title)
    figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label=cbar_label)
    figure.tight_layout()
    return figure


def _run_gsva_like(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    n_vars = int(matrix.shape[1])
    usable, source_label = _usable_gene_sets(ctx, params, warnings, min_set_size=5)
    max_cells = int(params.get("max_cells", 5000) or 5000)
    max_cells = min(max(max_cells, 100), 20000)
    matrix_sub, rows, _ = _limit_cells(ctx, matrix, max_cells, warnings, note="GSVA ")
    scores = _rank_ks_scores(matrix_sub, usable, n_vars)
    n_sets = int(params.get("n_sets", 12) or 12)
    n_sets = min(max(n_sets, 1), 40)
    selected = _select_sets_by_cluster_variance(ctx, scores, usable, rows, n_sets, warnings)
    prefix = str(params.get("obs_prefix") or "gsva")
    n_obs = _n_obs(ctx)
    used = set()
    obs_columns = []
    obs_values = {}
    summary_rows = []
    for set_index in selected:
        term, positions = usable[set_index]
        key = _obs_key(prefix, term, used)
        full_values = [None] * n_obs
        for local_index, original_row in enumerate(rows):
            full_values[int(original_row)] = float(scores[local_index, set_index])
        entry, values_map = _obs_pair(ctx, key, "numeric", full_values)
        obs_columns.append(entry)
        obs_values.update(values_map)
        set_scores = scores[:, set_index]
        summary_rows.append({
            "gene_set": term,
            "n_genes_in_data": int(positions.size),
            "obs_key": key,
            "mean_score": _num(np.mean(set_scores)),
            "sd_score": _num(np.std(set_scores)),
        })
    heat, set_names, column_labels, heat_mode = _set_cluster_heatmap(
        ctx, scores, rows, usable, selected, warnings)
    figure = _heatmap_figure(
        heat, set_names, column_labels,
        "Pathway activity (ssGSEA-style rank score) by %s" % ("cluster" if heat_mode == "cluster" else "cells"),
        "mean enrichment score", symmetric=True)
    figures = [_save_figure(figure, ctx, outdir, "gsva_like", "pathway_heatmap",
                            "Pathway activity heatmap")]
    tables = [_table(
        "set_summary", "Gene set activity summary",
        [_col("gene_set", "Gene set"), _col("n_genes_in_data", "Genes in data", "number"),
         _col("obs_key", "obs column"), _col("mean_score", "Mean score", "number"),
         _col("sd_score", "SD", "number")],
        summary_rows,
    )]
    if heat_mode == "cluster":
        wide_rows = []
        for row_index, term in enumerate(set_names):
            row = {"gene_set": term}
            for column_index, label in enumerate(column_labels):
                row["c_%d" % (column_index + 1)] = _num(heat[row_index, column_index])
            wide_rows.append(row)
        columns = [_col("gene_set", "Gene set")]
        for column_index, label in enumerate(column_labels):
            columns.append(_col("c_%d" % (column_index + 1), label, "number"))
        tables.append(_table("cluster_means", "Mean activity by cluster", columns, wide_rows))
    bullets = [
        "对 %d 个细胞的 %d 个基因集计算秩加权 KS 富集分数（ssGSEA 式简化实现）" % (
            int(scores.shape[0]), len(usable)),
        "选中的 %d 个基因集写入 obs（前缀 %s）；>= 上限的细胞未计算处为缺失" % (len(selected), prefix),
        "来源：%s；热图为簇 × 通路平均分数" % source_label,
    ]
    metrics = {
        "n_cells_scored": int(scores.shape[0]),
        "n_cells_total": n_obs,
        "n_gene_sets_available": len(usable),
        "n_gene_sets_used": len(selected),
        "obs_prefix": prefix,
    }
    return _result("gsva_like", bullets, metrics, tables, figures,
                   obs_columns=obs_columns, obs_values=obs_values, warnings=warnings)


def _rank_auc_scores(matrix, usable, n_vars):
    """AUCell 式 rank-AUC：AUC = 1 - (Σrank_asc - n(n+1)/2) / (n(G-n))。"""
    from scipy.stats import rankdata

    n_obs = int(matrix.shape[0])
    n_sets = len(usable)
    scores = np.zeros((n_obs, n_sets), dtype=np.float64)
    for cell in range(n_obs):
        row = np.asarray(matrix[cell].todense()).ravel().astype(np.float64)
        ranks = rankdata(row, method="average")
        for set_index, (_, positions) in enumerate(usable):
            member_ranks = ranks[positions]
            n_set = int(member_ranks.size)
            denominator = n_set * (n_vars - n_set)
            if denominator <= 0:
                scores[cell, set_index] = 1.0
                continue
            auc = 1.0 - (float(member_ranks.sum()) - n_set * (n_set + 1) / 2.0) / denominator
            scores[cell, set_index] = min(max(auc, 0.0), 1.0)
    return scores


def _run_auc_score(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    n_vars = int(matrix.shape[1])
    usable, source_label = _usable_gene_sets(ctx, params, warnings, min_set_size=5)
    max_cells = int(params.get("max_cells", 5000) or 5000)
    max_cells = min(max(max_cells, 100), 20000)
    matrix_sub, rows, _ = _limit_cells(ctx, matrix, max_cells, warnings, note="AUC ")
    scores = _rank_auc_scores(matrix_sub, usable, n_vars)
    n_sets = int(params.get("n_sets", 12) or 12)
    n_sets = min(max(n_sets, 1), 40)
    selected = _select_sets_by_cluster_variance(ctx, scores, usable, rows, n_sets, warnings)
    prefix = str(params.get("obs_prefix") or "auc")
    threshold_param = float(params.get("threshold", 0.0) or 0.0)
    n_obs = _n_obs(ctx)
    used = set()
    obs_columns = []
    obs_values = {}
    summary_rows = []
    thresholds = {}
    for set_index in selected:
        term, positions = usable[set_index]
        key = _obs_key(prefix, term, used)
        full_values = [None] * n_obs
        for local_index, original_row in enumerate(rows):
            full_values[int(original_row)] = float(scores[local_index, set_index])
        entry, values_map = _obs_pair(ctx, key, "numeric", full_values)
        obs_columns.append(entry)
        obs_values.update(values_map)
        set_scores = scores[:, set_index]
        threshold = threshold_param if threshold_param > 0 else float(np.mean(set_scores) + np.std(set_scores))
        thresholds[set_index] = threshold
        summary_rows.append({
            "gene_set": term,
            "n_genes_in_data": int(positions.size),
            "obs_key": key,
            "mean_auc": _num(np.mean(set_scores)),
            "sd_auc": _num(np.std(set_scores)),
            "threshold": _num(threshold),
            "pct_active": _num(100.0 * float(np.mean(set_scores >= threshold))),
        })
    heat, set_names, column_labels, heat_mode = _set_cluster_heatmap(
        ctx, scores, rows, usable, selected, warnings)
    figure = _heatmap_figure(
        heat, set_names, column_labels,
        "Gene set activity (AUCell-style rank AUC) by %s" % (
            "cluster" if heat_mode == "cluster" else "cells"),
        "mean AUC", symmetric=False)
    figures = [_save_figure(figure, ctx, outdir, "auc_score", "auc_heatmap",
                            "Gene set AUC heatmap")]
    tables = [_table(
        "set_summary", "Gene set AUC summary",
        [_col("gene_set", "Gene set"), _col("n_genes_in_data", "Genes in data", "number"),
         _col("obs_key", "obs column"), _col("mean_auc", "Mean AUC", "number"),
         _col("sd_auc", "SD", "number"), _col("threshold", "Threshold", "number"),
         _col("pct_active", "% active", "number")],
        summary_rows,
    )]
    if heat_mode == "cluster":
        wide_rows = []
        for row_index, term in enumerate(set_names):
            row = {"gene_set": term}
            for column_index, label in enumerate(column_labels):
                row["c_%d" % (column_index + 1)] = _num(heat[row_index, column_index])
            wide_rows.append(row)
        columns = [_col("gene_set", "Gene set")]
        for column_index, label in enumerate(column_labels):
            columns.append(_col("c_%d" % (column_index + 1), label, "number"))
        tables.append(_table("cluster_auc", "Mean AUC by cluster", columns, wide_rows))
    bullets = [
        "对 %d 个细胞的 %d 个基因集计算 AUCell 式 rank-AUC（秩和公式实现）" % (
            int(scores.shape[0]), len(usable)),
        "选中的 %d 个基因集写入 obs（前缀 %s）；阈值=%s" % (
            len(selected), prefix, "均值+1SD" if threshold_param <= 0 else "%.3f" % threshold_param),
        "来源：%s；热图为簇 × 基因集平均 AUC" % source_label,
    ]
    metrics = {
        "n_cells_scored": int(scores.shape[0]),
        "n_cells_total": n_obs,
        "n_gene_sets_available": len(usable),
        "n_gene_sets_used": len(selected),
        "obs_prefix": prefix,
    }
    return _result("auc_score", bullets, metrics, tables, figures,
                   obs_columns=obs_columns, obs_values=obs_values, warnings=warnings)


def _run_cnmf_programs(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    k = int(params.get("k", 6) or 6)
    k = min(max(k, 2), 20)
    n_top_genes = int(params.get("n_top_genes", 1000) or 1000)
    n_top_genes = min(max(n_top_genes, 50), 4000)
    max_iter = int(params.get("max_iter", 300) or 300)
    n_top = int(params.get("top_genes", 15) or 15)
    prefix = str(params.get("obs_prefix") or "nmf")
    n_obs = _n_obs(ctx)
    cell_limit = 10000
    matrix_sub, rows, _ = _limit_cells(ctx, matrix, cell_limit, warnings, note="NMF ")
    matrix_sub, gene_positions = _limit_genes(
        ctx, matrix_sub, np.arange(matrix_sub.shape[0]), n_top_genes, warnings, note="NMF ")
    n_cells, n_genes = int(matrix_sub.shape[0]), int(matrix_sub.shape[1])
    if n_genes < 20:
        raise EngineError("prerequisite", "可用基因只有 %d 个，太少无法做 NMF" % n_genes)
    k = min(k, max(2, n_cells - 1), max(2, n_genes - 1))
    dense = np.asarray(matrix_sub.toarray(), dtype=np.float64)
    dense = np.maximum(dense, 0.0)
    try:
        import warnings as _warnings

        from sklearn.decomposition import NMF
        from sklearn.exceptions import ConvergenceWarning

        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore", ConvergenceWarning)
            model = NMF(n_components=k, init="nndsvda", random_state=0,
                        max_iter=max_iter, tol=1e-4)
            usage = np.asarray(model.fit_transform(dense), dtype=np.float64)
            programs = np.asarray(model.components_, dtype=np.float64)
        try:
            iterations_used = int(getattr(model, "n_iter_", 0))
        except Exception:
            iterations_used = 0
        if iterations_used >= int(max_iter):
            warnings.append(
                "NMF 在 max_iter=%d 内未完全收敛；可提高 max_iter 或降低 k" % int(max_iter))
    except Exception as error:
        raise EngineError("internal", "NMF 分解失败：%s" % error)
    energy = usage.sum(axis=0)
    order = np.argsort(-energy, kind="stable")
    usage = usage[:, order]
    programs = programs[order]
    program_names = ["P%d" % (index + 1) for index in range(k)]
    gene_names = _var_names(ctx)
    selected_genes = [gene_names[int(gene_positions[index])] for index in range(n_genes)]

    obs_columns = []
    obs_values = {}
    used = set()
    obs_keys = []
    for program_index in range(k):
        key = _obs_key(prefix, "program_%d" % (program_index + 1), used)
        obs_keys.append(key)
        full_values = [None] * n_obs
        for local_index, original_row in enumerate(rows):
            full_values[int(original_row)] = float(usage[local_index, program_index])
        entry, values_map = _obs_pair(ctx, key, "numeric", full_values)
        obs_columns.append(entry)
        obs_values.update(values_map)

    gene_rows = []
    program_summary = []
    for program_index in range(k):
        weights = programs[program_index]
        total = float(weights.sum())
        top_order = np.argsort(-weights, kind="stable")[: int(max(n_top, 1))]
        for rank, index in enumerate(top_order, start=1):
            gene_rows.append({
                "program": program_names[program_index],
                "rank": rank,
                "gene": selected_genes[int(index)],
                "weight": _num(weights[int(index)]),
                "weight_share": _num(_safe_div(weights[int(index)], total, 0.0)),
            })
        top_genes_text = ", ".join(selected_genes[int(index)] for index in top_order[:5])
        program_summary.append({
            "program": program_names[program_index],
            "obs_key": obs_keys[program_index],
            "total_usage": _num(float(energy[order][program_index])),
            "mean_usage": _num(float(np.mean(usage[:, program_index]))),
            "top_genes": top_genes_text,
        })
    tables = [
        _table(
            "program_genes", "Top genes per expression program",
            [_col("program", "Program"), _col("rank", "Rank", "number"),
             _col("gene", "Gene"), _col("weight", "Weight", "number"),
             _col("weight_share", "Weight share", "number")],
            _cap_rows(gene_rows, 400, warnings, "程序基因表"),
        ),
        _table(
            "program_summary", "Program summary",
            [_col("program", "Program"), _col("obs_key", "obs column"),
             _col("total_usage", "Total usage", "number"),
             _col("mean_usage", "Mean usage", "number"), _col("top_genes", "Top genes")],
            program_summary,
        ),
    ]
    clusters = None
    try:
        clusters = _cluster_labels(ctx)
    except EngineError:
        clusters = None
    if clusters is not None and clusters.shape[0] == n_obs:
        sub_clusters = clusters[rows]
        levels = _ordered_labels(sub_clusters.tolist())
        masks = {level: sub_clusters == level for level in levels}
        levels = _limit_levels(levels, masks, 20, warnings, "簇")
        cluster_rows = []
        for program_index in range(k):
            row = {"program": program_names[program_index]}
            dominant_level = ""
            dominant_value = -np.inf
            for column_index, level in enumerate(levels):
                mask = masks[level]
                value = float(np.mean(usage[mask, program_index])) if mask.any() else 0.0
                row["c_%d" % (column_index + 1)] = _num(value)
                if value > dominant_value:
                    dominant_value = value
                    dominant_level = level
            row["top_cluster"] = dominant_level
            cluster_rows.append(row)
        columns = [_col("program", "Program")]
        for column_index, level in enumerate(levels):
            columns.append(_col("c_%d" % (column_index + 1), _fig_label(level, "C%d" % column_index),
                                "number"))
        columns.append(_col("top_cluster", "Top cluster"))
        tables.append(_table("usage_by_cluster", "Mean program usage by cluster", columns, cluster_rows))
    # 细胞用法矩阵热图（按主导程序排序）
    n_plot = min(n_cells, 2000)
    if n_cells > n_plot:
        plot_rows = np.sort(_rng().choice(n_cells, size=n_plot, replace=False))
    else:
        plot_rows = np.arange(n_cells)
    dominant = np.argmax(usage[plot_rows], axis=1)
    plot_order = plot_rows[np.lexsort((usage[plot_rows].sum(axis=1), dominant))]
    heat = usage[plot_order].T
    normalized = heat / np.maximum(heat.max(axis=1, keepdims=True), 1e-9)
    figure, axis = _new_figure(8.6, 4.6)
    image = axis.imshow(normalized, aspect="auto", cmap="magma", vmin=0.0, vmax=1.0)
    axis.set_yticks(np.arange(k))
    axis.set_yticklabels(program_names, fontsize=8)
    axis.set_xlabel("Cells (sorted by dominant program)")
    axis.set_title("NMF expression program usage")
    figure.colorbar(image, ax=axis, fraction=0.026, pad=0.01, label="usage (row-normalized)")
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "cnmf_programs", "program_usage_heatmap",
                            "NMF program usage heatmap")]
    bullets = [
        "NMF（sklearn，random_state=0）在 %d 个 HVG 上分解出 %d 个表达程序" % (n_genes, k),
        "程序按总使用量降序命名 P1..Pk；程序基因表每程序 top %d 个基因" % n_top,
        "细胞用法写入 obs（%s）；热图按主导程序排序" % "、".join(obs_keys[:4]) + ("…" if k > 4 else ""),
    ]
    metrics = {
        "k": k,
        "n_top_genes": n_genes,
        "n_cells": n_cells,
        "n_cells_total": n_obs,
        "max_iter": max_iter,
    }
    return _result("cnmf_programs", bullets, metrics, tables, figures,
                   obs_columns=obs_columns, obs_values=obs_values, warnings=warnings)


def _module_eigengene(standardized, seed=0):
    """模块特征基因：标准化表达的第一主成分得分，符号与平均表达对齐。"""
    centered = standardized - standardized.mean(axis=0, keepdims=True)
    try:
        u, s, _ = np.linalg.svd(centered, full_matrices=False)
        eigengene = u[:, 0] * s[0]
    except Exception:
        eigengene = centered.mean(axis=1)
    if np.corrcoef(eigengene, centered.mean(axis=1))[0, 1] < 0:
        eigengene = -eigengene
    return eigengene


def _run_wgcna_lite(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    n_obs = _n_obs(ctx)
    n_top_genes = int(params.get("n_top_genes", 800) or 800)
    n_top_genes = min(max(n_top_genes, 100), 2000)
    min_corr = float(params.get("min_corr", 0.5) or 0.5)
    min_corr = min(max(min_corr, -0.5), 0.95)
    min_module_size = int(params.get("min_module_size", 20) or 20)
    n_hub = int(params.get("n_hub", 5) or 5)
    cell_limit = 5000
    matrix_sub, rows, _ = _limit_cells(ctx, matrix, cell_limit, warnings, note="WGCNA ")
    matrix_sub, gene_positions = _limit_genes(
        ctx, matrix_sub, np.arange(matrix_sub.shape[0]), n_top_genes, warnings, note="WGCNA ")
    n_cells, n_genes = int(matrix_sub.shape[0]), int(matrix_sub.shape[1])
    if n_genes < 30:
        raise EngineError("prerequisite", "可用基因只有 %d 个，太少无法做共表达模块分析" % n_genes)
    dense = np.asarray(matrix_sub.toarray(), dtype=np.float64)
    mean = dense.mean(axis=0)
    sd = dense.std(axis=0)
    sd[sd <= 0] = 1.0
    standardized = (dense - mean) / sd
    correlation = np.corrcoef(standardized, rowvar=False)
    correlation = np.nan_to_num(correlation, nan=0.0, posinf=1.0, neginf=-1.0)
    np.fill_diagonal(correlation, 1.0)
    distance = 1.0 - correlation
    distance = (distance + distance.T) / 2.0
    np.clip(distance, 0.0, 2.0, out=distance)
    np.fill_diagonal(distance, 0.0)
    try:
        condensed = _squareform(distance, checks=False)
        linkage_matrix = _hierarchy.linkage(condensed, method="average")
        height = max(0.05, min(1.99, 1.0 - min_corr))
        raw_labels = _hierarchy.fcluster(linkage_matrix, t=height, criterion="distance")
    except Exception as error:
        raise EngineError("internal", "共表达模块聚类失败：%s" % error)
    module_members = {}
    for gene_index, label in enumerate(raw_labels.tolist()):
        module_members.setdefault(int(label), []).append(int(gene_index))
    big = [(label, members) for label, members in module_members.items()
           if len(members) >= min_module_size]
    big.sort(key=lambda item: (-len(item[1]), item[0]))
    if not big:
        biggest = max(module_members.items(), key=lambda item: len(item[1]))
        big = [biggest]
        warnings.append("没有模块达到最小模块大小 %d，已把最大模块作为 M1 输出" % min_module_size)
    if len(big) == 1:
        warnings.append("只得到 1 个共表达模块；可降低 min_corr 以细分模块")
    gene_names = _var_names(ctx)
    selected_genes = [gene_names[int(gene_positions[index])] for index in range(n_genes)]
    module_eigengenes = {}
    hub_rows = []
    module_rows = []
    for module_index, (label, members) in enumerate(big, start=1):
        module_name = "M%d" % module_index
        module_genes = np.asarray(members, dtype=np.int64)
        eigengene = _module_eigengene(standardized[:, module_genes])
        module_eigengenes[module_name] = (module_genes, eigengene)
        kme = np.array([
            float(np.corrcoef(standardized[:, gene], eigengene)[0, 1])
            if np.std(standardized[:, gene]) > 0 and np.std(eigengene) > 0 else 0.0
            for gene in module_genes
        ])
        kme = np.nan_to_num(kme, nan=0.0)
        hub_order = np.argsort(-kme, kind="stable")[: int(max(n_hub, 1))]
        for rank, position in enumerate(hub_order, start=1):
            hub_rows.append({
                "module": module_name,
                "rank": rank,
                "gene": selected_genes[int(module_genes[int(position)])],
                "kME": _num(kme[int(position)]),
                "module_size": int(module_genes.size),
            })
        module_rows.append({
            "module": module_name,
            "n_genes": int(module_genes.size),
            "mean_kME": _num(float(np.mean(kme))),
            "top_hub_genes": ", ".join(
                selected_genes[int(module_genes[int(position)])]
                for position in hub_order[:3]),
        })
    tables = [
        _table(
            "modules", "Co-expression modules",
            [_col("module", "Module"), _col("n_genes", "Genes", "number"),
             _col("mean_kME", "Mean kME", "number"), _col("top_hub_genes", "Top hub genes")],
            module_rows,
        ),
        _table(
            "hub_genes", "Hub genes per module",
            [_col("module", "Module"), _col("rank", "Rank", "number"),
             _col("gene", "Gene"), _col("kME", "kME", "number"),
             _col("module_size", "Module size", "number")],
            _cap_rows(hub_rows, 400, warnings, "hub 基因表"),
        ),
    ]
    # obs：模块特征基因（最多 8 个模块）
    obs_columns = []
    obs_values = {}
    used = set()
    obs_keys = []
    n_modules_obs = min(len(big), 8)
    for module_index in range(n_modules_obs):
        module_name = "M%d" % (module_index + 1)
        module_genes, _ = module_eigengenes[module_name]
        original_positions = np.asarray(
            [int(gene_positions[int(index)]) for index in module_genes], dtype=np.int64)
        module_dense = np.asarray(matrix[:, original_positions].toarray(), dtype=np.float64)
        module_sd = module_dense.std(axis=0)
        module_sd[module_sd <= 0] = 1.0
        eigengene_all = _module_eigengene((module_dense - module_dense.mean(axis=0)) / module_sd)
        key = _obs_key(prefix="wgcna_ME" if False else "wgcna",
                       term="ME_%d" % (module_index + 1), used=used)
        obs_keys.append(key)
        full_values = [float(value) for value in eigengene_all]
        entry, values_map = _obs_pair(ctx, key, "numeric", full_values)
        obs_columns.append(entry)
        obs_values.update(values_map)
    clusters = None
    try:
        clusters = _cluster_labels(ctx)
    except EngineError:
        clusters = None
    if clusters is not None and clusters.shape[0] == n_obs:
        levels = _ordered_labels(clusters.tolist())
        masks = {level: clusters == level for level in levels}
        levels = _limit_levels(levels, masks, 20, warnings, "簇")
        heat = np.zeros((n_modules_obs, len(levels)), dtype=np.float64)
        for module_index in range(n_modules_obs):
            module_genes, eigengene = module_eigengenes["M%d" % (module_index + 1)]
            original_positions = np.asarray(
                [int(gene_positions[int(index)]) for index in module_genes], dtype=np.int64)
            module_dense = np.asarray(matrix[:, original_positions].toarray(), dtype=np.float64)
            module_sd = module_dense.std(axis=0)
            module_sd[module_sd <= 0] = 1.0
            eigengene_all = _module_eigengene((module_dense - module_dense.mean(axis=0)) / module_sd)
            for column_index, level in enumerate(levels):
                mask = masks[level]
                heat[module_index, column_index] = float(np.mean(eigengene_all[mask])) if mask.any() else 0.0
        row_labels = ["M%d (n=%d)" % (index + 1, len(big[index][1])) for index in range(n_modules_obs)]
        figure = _heatmap_figure(
            heat, row_labels,
            [_fig_label(level, "C%d" % index) for index, level in enumerate(levels)],
            "Module eigengene mean expression by cluster (z-score)",
            "z-score", symmetric=True)
        figures = [_save_figure(figure, ctx, outdir, "wgcna_lite", "module_cluster_heatmap",
                                "Module-cluster eigengene heatmap")]
        wide_rows = []
        for row_index, label in enumerate(row_labels):
            row = {"module": label}
            for column_index, level in enumerate(levels):
                row["c_%d" % (column_index + 1)] = _num(heat[row_index, column_index])
            wide_rows.append(row)
        columns = [_col("module", "Module")]
        for column_index, level in enumerate(levels):
            columns.append(_col("c_%d" % (column_index + 1), _fig_label(level, "C%d" % column_index),
                                "number"))
        tables.append(_table("module_cluster", "Module eigengene by cluster", columns, wide_rows))
    else:
        figure, axis = _new_figure(7.2, 3.8)
        sizes = [len(members) for _, members in big]
        yvalues = np.arange(len(big))
        axis.barh(yvalues, sizes, color=_palette(len(big)), height=0.65)
        axis.set_yticks(yvalues)
        axis.set_yticklabels(
            ["M%d (n=%d)" % (index + 1, len(big[index][1])) for index in range(len(big))], fontsize=8)
        axis.set_xlabel("Genes in module")
        axis.set_title("WGCNA-lite co-expression modules")
        _style_axis(axis, grid_axis="x")
        figure.tight_layout()
        figures = [_save_figure(figure, ctx, outdir, "wgcna_lite", "module_sizes",
                                "Co-expression module sizes")]
    bullets = [
        "在 %d 个 HVG 上计算相关矩阵 → 层次聚类（min_corr=%.2f）得到 %d 个模块" % (
            n_genes, min_corr, len(big)),
        "模块按大小命名 M1..M%d；每模块输出 top %d 枢纽基因（kME）" % (len(big), n_hub),
        "模块特征基因写入 obs（%s）%s" % (
            "、".join(obs_keys), "（最多 8 个）" if len(big) > 8 else ""),
    ]
    metrics = {
        "n_genes": n_genes,
        "n_cells": n_cells,
        "n_cells_total": n_obs,
        "n_modules": len(big),
        "min_corr": _num(min_corr),
        "min_module_size": min_module_size,
    }
    return _result("wgcna_lite", bullets, metrics, tables, figures,
                   obs_columns=obs_columns, obs_values=obs_values, warnings=warnings)


def _is_internal_obs_column(column):
    text = str(column)
    if text.startswith("_"):
        return True
    for prefix in ("phase", "S_score", "G2M_score", "subcluster", "gsva_", "auc_",
                   "nmf_", "wgcna_", "module_"):
        if text.startswith(prefix):
            return True
    return False


def _run_geneset_clinical(ctx, params, outdir):
    warnings = []
    matrix = _expression(ctx, warnings)
    usable, source_label = _usable_gene_sets(ctx, params, warnings, min_set_size=3)
    n_sets = int(params.get("n_sets", 8) or 8)
    n_sets = min(max(n_sets, 1), 30)
    usable.sort(key=lambda item: (-item[1].size, item[0]))
    usable = usable[:n_sets]
    max_levels = int(params.get("max_levels", 12) or 12)
    obs = _ctx_obs(ctx)
    requested = str(params.get("columns") or "").strip()
    if requested:
        candidates = [item for item in re.split(r"[,;\s]+", requested) if item]
        missing = [column for column in candidates if column not in obs.columns]
        if missing:
            raise EngineError(
                "bad_request", "obs 里没有这些元数据列：%s（可用：%s）" % (
                    "、".join(missing), "、".join(str(column) for column in obs.columns[:20])))
    else:
        candidates = [str(column) for column in obs.columns if not _is_internal_obs_column(column)]
    metadata = []
    for column in candidates:
        if column not in obs.columns:
            continue
        series = obs[column]
        if pd.api.types.is_bool_dtype(series):
            series = series.astype("object")
        try:
            unique = int(series.nunique(dropna=True))
        except Exception:
            continue
        if unique < 2:
            continue
        numeric = pd.api.types.is_numeric_dtype(series) and unique > max_levels
        if numeric:
            metadata.append((column, "numeric", series))
        elif unique <= max_levels:
            metadata.append((column, "categorical", series))
        else:
            continue
    if not metadata:
        raise EngineError(
            "prerequisite",
            "obs 里没有可用的元数据列（分类列 ≤%d 个水平，或数值列）；请先 annotate 或提供 columns 参数"
            % max_levels,
        )
    set_scores = {}
    for term, positions in usable:
        values = np.asarray(matrix[:, positions].mean(axis=1)).ravel().astype(np.float64)
        set_scores[term] = _zscore(values)
    association_rows = []
    pvalues = []
    for term, _positions in usable:
        scores = set_scores[term]
        for column, column_type, series in metadata:
            if column_type == "numeric":
                values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
                mask = np.isfinite(values)
                if int(mask.sum()) < 10 or int(np.unique(values[mask]).size) < 3:
                    continue
                try:
                    rho, pvalue = _stats.spearmanr(scores[mask], values[mask])
                except Exception:
                    continue
                if not math.isfinite(float(rho)):
                    continue
                pvalue = _safe_p(pvalue)
                pvalues.append(pvalue)
                association_rows.append({
                    "gene_set": term,
                    "variable": column,
                    "variable_type": "numeric",
                    "test": "Spearman",
                    "statistic": _num(rho),
                    "pval": _num(pvalue),
                    "pval_adj": None,
                    "effect": _num(rho * rho),
                    "n": int(mask.sum()),
                    "detail": "rho=%.3f" % rho,
                })
            else:
                values = series.astype("object")
                texts = np.asarray([str(value) for value in values.tolist()], dtype=object)
                groups = []
                labels = []
                for level in _ordered_labels(texts.tolist()):
                    mask = texts == level
                    count = int(mask.sum())
                    if count < 3:
                        continue
                    groups.append(scores[mask])
                    labels.append(level)
                if len(groups) < 2:
                    continue
                try:
                    f_stat, pvalue = _stats.f_oneway(*groups)
                except Exception:
                    continue
                if not math.isfinite(float(f_stat)):
                    continue
                grand = float(np.mean(scores))
                ss_between = float(sum(group.size * (float(np.mean(group)) - grand) ** 2 for group in groups))
                ss_total = float(np.sum((scores - grand) ** 2))
                eta_squared = _safe_div(ss_between, ss_total, 0.0)
                pvalue = _safe_p(pvalue)
                try:
                    _, kruskal_p = _stats.kruskal(*groups)
                except Exception:
                    kruskal_p = 1.0
                best_group = labels[int(np.argmax([float(np.mean(group)) for group in groups]))]
                pvalues.append(pvalue)
                association_rows.append({
                    "gene_set": term,
                    "variable": column,
                    "variable_type": "categorical",
                    "test": "ANOVA",
                    "statistic": _num(f_stat),
                    "pval": _num(pvalue),
                    "pval_adj": None,
                    "effect": _num(eta_squared),
                    "n": int(sum(group.size for group in groups)),
                    "detail": "highest: %s; Kruskal p=%.3g" % (best_group, _safe_p(kruskal_p)),
                })
    if not association_rows:
        raise EngineError("prerequisite", "没有任何基因集 × 元数据列组合可以检验（样本量或水平数不足）")
    adjusted = _bh_fdr(np.asarray(pvalues, dtype=np.float64))
    for row, padj in zip(association_rows, adjusted):
        row["pval_adj"] = _num(padj)
    association_rows.sort(
        key=lambda row: (row["pval_adj"] if row["pval_adj"] is not None else 1.0,
                         row["gene_set"], row["variable"]))
    tables = [_table(
        "associations", "Gene set score vs metadata associations",
        [_col("gene_set", "Gene set"), _col("variable", "Variable"),
         _col("variable_type", "Variable type"), _col("test", "Test"),
         _col("statistic", "Statistic", "number"), _col("pval", "p", "number"),
         _col("pval_adj", "FDR", "number"), _col("effect", "Effect size", "number"),
         _col("n", "N", "number"), _col("detail", "Detail")],
        _cap_rows(association_rows, 400, warnings, "关联表"),
        note="分类变量用单因素 ANOVA（效应量 eta²），数值变量用 Spearman 相关；BH 校正。",
    )]
    top = association_rows[:15]
    figure, axis = _new_figure(8.6, max(3.2, 0.36 * len(top) + 1.8))
    yvalues = np.arange(len(top))[::-1]
    heights = [_minus_log10(row["pval_adj"]) for row in top]
    colors = ["#f28e2b" if row["variable_type"] == "categorical" else "#4e79a7" for row in top]
    axis.barh(yvalues, heights, color=colors, height=0.66)
    axis.set_yticks(yvalues)
    axis.set_yticklabels(
        ["%s | %s" % (_shorten(row["gene_set"], 26), _shorten(row["variable"], 16)) for row in top],
        fontsize=7.5)
    axis.set_xlabel("-log10(FDR)")
    axis.set_title("Gene set score associations (ANOVA / Spearman)")
    _style_axis(axis, grid_axis="x")
    figure.tight_layout()
    figures = [_save_figure(figure, ctx, outdir, "geneset_clinical", "association_bar",
                            "Gene set vs metadata associations")]
    bullet_pairs = "%d" % len(association_rows)
    bullets = [
        "基因集 %d 个 × 元数据列 %d 个，共 %s 个关联检验（BH 校正）" % (
            len(usable), len(metadata), bullet_pairs),
        "分类变量用 ANOVA（eta²），数值变量用 Spearman（rho）；FDR<0.05 的关联 %d 个" % (
            int(np.sum(adjusted < 0.05))),
        "基因集分数 = 集合内基因 log 标准化表达均值（跨细胞 z-score）；来源：%s" % source_label,
    ]
    metrics = {
        "n_gene_sets": len(usable),
        "n_metadata_columns": len(metadata),
        "n_tests": len(association_rows),
        "n_significant": int(np.sum(adjusted < 0.05)),
    }
    return _result("geneset_clinical", bullets, metrics, tables, figures, warnings=warnings)


# ===========================================================================
# 11. 参数校验与分发入口
# ===========================================================================

_RUNNERS = {
    "cell_cycle": _run_cell_cycle,
    "subcluster": _run_subcluster,
    "proportion_test": _run_proportion_test,
    "de_group": _run_de_group,
    "volcano": _run_volcano,
    "proportion_sankey": _run_proportion_sankey,
    "featureplot": _run_featureplot,
    "pathway_ora": _run_pathway_ora,
    "gsea_prerank": _run_gsea_prerank,
    "gsva_like": _run_gsva_like,
    "auc_score": _run_auc_score,
    "cnmf_programs": _run_cnmf_programs,
    "wgcna_lite": _run_wgcna_lite,
    "geneset_clinical": _run_geneset_clinical,
}


def _validate_params(analysis_id, params):
    spec = _spec(analysis_id)
    params = params if isinstance(params, dict) else {}
    cleaned = {}
    for entry in spec.get("params") or []:
        key = str(entry["key"])
        kind = str(entry.get("type") or "string")
        default = entry.get("default")
        value = params.get(key, default)
        if value is None:
            value = default
        if kind == "number":
            try:
                if isinstance(value, bool):
                    raise ValueError
                number = float(value)
            except Exception:
                raise EngineError(
                    "bad_request", "参数 %s 需要数值，收到 %r" % (key, value))
            if not math.isfinite(number):
                raise EngineError("bad_request", "参数 %s 需要有限数值" % key)
            if isinstance(default, int) and not isinstance(default, bool):
                cleaned[key] = int(round(number))
            else:
                cleaned[key] = number
        elif kind == "boolean":
            cleaned[key] = bool(value)
        elif kind == "select":
            options = [str(option) for option in (entry.get("options") or [])]
            text = str(value)
            if options and text not in options:
                raise EngineError(
                    "bad_request", "参数 %s 取值 %r 不在允许范围：%s" % (key, text, "、".join(options)))
            cleaned[key] = text
        else:
            cleaned[key] = "" if value is None else str(value)
    return cleaned


def run(analysis_id, ctx, params, outdir):
    """模块统一入口（契约 §9.2）：run(analysis_id, ctx, params, outdir) -> dict。"""
    analysis_id = str(analysis_id)
    runner = _RUNNERS.get(analysis_id)
    if runner is None:
        raise EngineError(
            "not_found",
            "未知分析：%s（sc-analysis-a 支持：%s）" % (analysis_id, ", ".join(_RUNNERS)),
        )
    if params is None:
        params = {}
    if not isinstance(params, dict):
        raise EngineError("bad_request", "分析参数必须是 JSON 对象")
    cleaned = _validate_params(analysis_id, params)
    try:
        directory = Path(str(outdir)).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
    except Exception as error:
        raise EngineError("internal", "无法创建图输出目录：%s" % error)
    return runner(ctx, cleaned, directory)


__all__ = ["ANALYSES", "EngineError", "run"]
