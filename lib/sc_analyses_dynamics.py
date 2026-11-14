#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench — sc-analysis-b 高级分析模块（整合 / 轨迹 / 异质性 / 通讯 / 调控 / 虚拟敲除）

契约
----
- CONTRACT.md §9.2 模块 API：``ANALYSES`` 元数据 + ``run(analysis_id, ctx, params, outdir)``；
- 返回键固定：``id/title/method/summary/tables/figures/obs_columns/warnings``；
- 需要写回 obs 的分析额外返回顶层 ``obs_values``：``{列名: [与 ctx.obs.index 等长的值]}``，
  ``obs_columns`` 里同时声明 ``{"key":..., "kind":"numeric|categorical"}``；
- 缺前置数据抛 ``EngineError("prerequisite", ...)``，参数错误抛 ``EngineError("bad_request", ...)``；
- 图写 ``outdir``，白底 dpi≈140，图上文字全英文；结果确定性（固定种子，与调用次数无关）；
- 细胞 ≤ 20000、基因 ≤ 4000，超限随机下采样并在 warnings 写明。

诚实声明
--------
本模块里的高级分析全部是**本插件内的重实现 / 近似**，不是官方工具：
``integration``（ComBat 式位置-尺度 + BBKNN 式批次平衡 kNN，无经验贝叶斯/精确 BBKNN）、
``trajectory_dpt``（扩散图 + DPT 简化实现，非 monocle）、
``trajectory_paga``（簇间边计数/期望比 + 最大生成树，非官方 PAGA）、
``entropy_stemness``（表达熵 + 邻居平滑启发式）、
``cytotrace_lite``（CytoTRACE 思路的简化实现）、
``cnv_score`` / ``cnv_tree`` / ``ith_score``（inferCNV 代理分数：参照均值/标准差 z 分数 +
移动窗口平滑；无基因顺序时是按方差排序的代理分数，不能当真实 CNV 区段解读）、
``lr_communication``（CellChat/NATMI-lite：表达乘积 + 置换 p，不复刻 CellChat 数据库/复合物）、
``regulon_activity``（SCENIC-lite：内置 TF 列表 + 共表达靶基因 + AUCell 式排序分数，无 motif 数据库）、
``virtual_knockout`` / ``virtual_overexpress`` / ``virtual_screen``（scTenifoldKnk-lite：
共表达相关 GRN + 稀疏化 + 影响传播，非 PC regression/张量分解）。
"""

from __future__ import annotations

import importlib.util
import math
import os
import re
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# matplotlib 必须在任何 pyplot 导入前设定后端（引擎也会设，这里兜底以便独立测试）
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "dsh-singlecell-workbench-mpl")
)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

matplotlib.rcParams.update(
    {
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
    }
)

# ---------------------------------------------------------------------------
# 内置生物学先验（同目录 biopriors.py；按文件路径加载，兼容 registry 的动态导入）
# ---------------------------------------------------------------------------


def _load_biopriors():
    cached = sys.modules.get("_dsh_sc_biopriors")
    if cached is not None:
        return cached
    try:  # 正常情况：lib/ 在 sys.path 上
        import biopriors as module  # type: ignore

        sys.modules.setdefault("_dsh_sc_biopriors", module)
        return module
    except Exception:
        pass
    path = Path(__file__).resolve().parent / "biopriors.py"
    if not path.exists():  # pragma: no cover - 防御
        raise EngineError("unsupported", "缺少先验文件：%s" % path)
    spec = importlib.util.spec_from_file_location("_dsh_sc_biopriors", str(path))
    module = importlib.util.module_from_spec(spec)
    sys.modules["_dsh_sc_biopriors"] = module
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

MAX_CELLS = 20000
MAX_GENES = 4000
DEFAULT_SEED = 0

_LR_PAIRS = []  # 延迟填充，避免 import 顺序问题
_TF_LIST = []
_PRIORS_LOADED = False


def _priors():
    """惰性加载 biopriors，返回 (lr_pairs, tf_list)。"""
    global _PRIORS_LOADED, _LR_PAIRS, _TF_LIST
    if not _PRIORS_LOADED:
        module = _load_biopriors()
        _LR_PAIRS = list(module.lr_pairs())
        _TF_LIST = list(module.tf_list())
        _PRIORS_LOADED = True
    return _LR_PAIRS, _TF_LIST


class EngineError(Exception):
    """可预期的业务错误：code ∈ CONTRACT §2.1 的错误码集合。"""

    CODES = {"bad_request", "not_found", "prerequisite", "unsupported", "internal"}

    def __init__(self, code, message, details=None):
        self.code = code if code in self.CODES else "internal"
        self.message = str(message)
        self.details = details or None
        super().__init__(self.message)


# ---------------------------------------------------------------------------
# 通用小工具
# ---------------------------------------------------------------------------


def _num(value, digits=None):
    """把 numpy 数值转成 JSON 安全数值；非有限值返回 None。"""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    if digits is not None:
        return round(number, int(digits))
    return number


def _int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _slug(text):
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(text or "dataset")).strip("-._")
    return cleaned or "dataset"


def _clamp_int(params, key, default, minimum, maximum):
    value = params.get(key, default)
    if value is None or value == "":
        value = default
    try:
        number = int(float(value))
    except (TypeError, ValueError):
        raise EngineError("bad_request", "参数 %s 需要整数，收到 %r" % (key, value))
    return max(int(minimum), min(int(maximum), number))


def _clamp_float(params, key, default, minimum, maximum):
    value = params.get(key, default)
    if value is None or value == "":
        value = default
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise EngineError("bad_request", "参数 %s 需要数值，收到 %r" % (key, value))
    if not math.isfinite(number):
        raise EngineError("bad_request", "参数 %s 需要有限数值" % key)
    return max(float(minimum), min(float(maximum), number))


def _rng(seed=DEFAULT_SEED):
    """分析内部固定种子（不依赖 ctx['rng']，保证同名参数结果一致）。"""
    return np.random.default_rng(int(seed))


def _sorted_categories(labels):
    """给簇名排序：纯数字按数值排，其余按字符串排（稳定、确定性）。"""
    unique = []
    seen = set()
    for label in labels:
        text = str(label)
        if text not in seen:
            seen.add(text)
            unique.append(text)

    def key(text):
        try:
            return (0, float(text), "")
        except (TypeError, ValueError):
            return (1, 0.0, text)

    return sorted(unique, key=key)


_BASE_COLORS = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2", "#EECA3B",
    "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#1F77B4", "#FF7F0E",
    "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
    "#BCBD22", "#17BECF",
]


def _palette(count):
    """确定性调色板（前 20 色固定，更多时用 turbo 均匀采样）。"""
    count = max(1, int(count))
    if count <= len(_BASE_COLORS):
        return list(_BASE_COLORS[:count])
    cmap = plt.get_cmap("turbo")
    return [cmap(value) for value in np.linspace(0.05, 0.95, count)]


def _rank01(values):
    """把数值转成 [0,1] 的秩（并列取平均，确定性）。"""
    array = np.asarray(values, dtype=np.float64)
    if array.size <= 1:
        return np.zeros(array.shape, dtype=np.float64)
    order = np.argsort(array, kind="stable")
    ranks = np.empty(array.size, dtype=np.float64)
    ranks[order] = np.arange(array.size, dtype=np.float64)
    # 简单并列平均（扫描排序后的等值段）
    sorted_values = array[order]
    start = 0
    for position in range(1, array.size + 1):
        if position == array.size or sorted_values[position] != sorted_values[start]:
            if position - start > 1:
                ranks[order[start:position]] = (start + position - 1) / 2.0
            start = position
    return (ranks - ranks.min()) / (ranks.max() - ranks.min() + 1e-12)


def _zscore(values):
    array = np.asarray(values, dtype=np.float64)
    std = float(array.std())
    if not math.isfinite(std) or std <= 1e-12:
        return np.zeros(array.shape, dtype=np.float64)
    return (array - array.mean()) / std


def _png_size(path):
    try:
        with open(path, "rb") as handle:
            head = handle.read(24)
        if head[:8] == b"\x89PNG\r\n\x1a\n" and len(head) >= 24:
            return struct.unpack(">II", head[16:24])
    except Exception:
        pass
    return (0, 0)


def _jsonable(value):
    """递归转成 json.dumps 可接受的结构（NaN/Inf → None，numpy → Python）。"""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, np.ndarray):
        return [_jsonable(item) for item in value.tolist()]
    if isinstance(value, (np.str_,)):
        return str(value)
    if value is None or isinstance(value, (int, str)):
        return value
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


# ---------------------------------------------------------------------------
# ctx / 矩阵视图
# ---------------------------------------------------------------------------


def _log_normalize(matrix, target_sum=1e4):
    """counts → log1p(每细胞归一化到 target_sum)。"""
    matrix = sp.csr_matrix(matrix, dtype=np.float32)
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    scale = np.zeros(totals.shape, dtype=np.float32)
    nonzero = totals > 0
    scale[nonzero] = (float(target_sum) / totals[nonzero]).astype(np.float32)
    normalized = (sp.diags(scale) @ matrix).tocsr().astype(np.float32)
    normalized.data = np.log1p(normalized.data)
    return normalized


def _matrix_lognorm(ctx):
    matrix = ctx.get("X")
    if matrix is None:
        raise EngineError("prerequisite", "当前数据集没有表达矩阵，请先导入数据")
    if sp.issparse(matrix):
        matrix = matrix.tocsr().astype(np.float32, copy=False)
    else:
        matrix = sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    kind = str(ctx.get("X_kind") or "counts").lower()
    if kind != "lognorm":
        matrix = _log_normalize(matrix)
    return matrix


def _counts(ctx):
    matrix = ctx.get("counts")
    if matrix is None:
        matrix = ctx.get("X")
    if matrix is None:
        raise EngineError("prerequisite", "当前数据集没有表达矩阵，请先导入数据")
    if sp.issparse(matrix):
        return matrix.tocsr().astype(np.float32, copy=False)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float32))


def _view(ctx, warnings=None, max_cells=MAX_CELLS, seed=DEFAULT_SEED):
    """构造分析视图：可选细胞下采样，统一 lognorm 矩阵 + obs/obsm/簇标签。"""
    if not isinstance(ctx, dict):
        raise EngineError("internal", "分析 ctx 必须是对象")
    matrix = _matrix_lognorm(ctx)
    n_obs = int(ctx.get("n_obs") or matrix.shape[0])
    if matrix.shape[0] != n_obs:
        n_obs = int(matrix.shape[0])
    var_names = [str(gene) for gene in (ctx.get("var_names") or [])]
    if len(var_names) != matrix.shape[1]:
        var_names = ["gene%d" % (index + 1) for index in range(matrix.shape[1])]
    idx = np.arange(n_obs)
    if max_cells and n_obs > int(max_cells):
        rng = _rng(seed)
        idx = np.sort(rng.choice(n_obs, size=int(max_cells), replace=False))
        if warnings is not None:
            warnings.append(
                "细胞数 %d 超过分析上限 %d，已随机下采样至 %d（固定 seed=%d，结果可复现）"
                % (n_obs, int(max_cells), int(max_cells), int(seed))
            )
    obs = ctx.get("obs")
    if isinstance(obs, pd.DataFrame):
        obs = obs.iloc[idx].copy()
    else:  # pragma: no cover - 防御：没有 obs 时造一个空表
        obs = pd.DataFrame(index=pd.Index(["cell%d" % (i + 1) for i in idx]))
    obsm = ctx.get("obsm") or {}
    pca = obsm.get("X_pca")
    pca = np.asarray(pca, dtype=np.float32)[idx] if pca is not None else None
    emb = obsm.get("X_emb")
    emb = np.asarray(emb, dtype=np.float32)[idx] if emb is not None else None

    clusters = None
    raw = ctx.get("clusters")
    if isinstance(raw, dict) and raw.get("labels") is not None:
        labels = [str(value) for value in raw["labels"]]
        if len(labels) == n_obs:
            clusters = np.asarray(labels, dtype=object)[idx]
    if clusters is None and "clusters" in obs.columns:
        clusters = obs["clusters"].astype(str).to_numpy(dtype=object)

    return {
        "X": matrix[idx] if idx.size != n_obs else matrix,
        "obs": obs,
        "var_names": var_names,
        "pca": pca,
        "emb": emb,
        "clusters": clusters,
        "idx": idx,
        "n_total": n_obs,
        "n_vars": matrix.shape[1],
        "ctx": ctx,
    }


def _hvg_indices(matrix, n_genes):
    """按方差取高变基因索引（升序返回，确定性）。"""
    total = int(matrix.shape[1])
    n_genes = max(1, min(int(n_genes), total))
    if n_genes >= total:
        return np.arange(total, dtype=np.int64)
    mean = np.asarray(matrix.mean(axis=0)).ravel().astype(np.float64)
    mean_square = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel().astype(np.float64)
    variance = np.maximum(mean_square - mean * mean, 0.0)
    order = np.argsort(-variance, kind="stable")[:n_genes]
    return np.sort(order.astype(np.int64))


def _pca(view, n_pcs=30):
    """取 PCA：优先 ctx.obsm['X_pca']，缺失时用 HVG 现场计算（确定性）。"""
    provided = view.get("pca")
    if provided is not None and provided.ndim == 2 and provided.shape[1] >= 2:
        wide = np.asarray(provided, dtype=np.float32)
        return wide[:, : int(min(n_pcs, wide.shape[1]))]
    matrix = view["X"]
    genes = _hvg_indices(matrix, min(2000, matrix.shape[1]))
    sub = matrix[:, genes]
    n_cells = sub.shape[0]
    components = int(min(n_pcs, len(genes) - 1, n_cells - 1))
    if components < 2:
        raise EngineError("prerequisite", "细胞/基因数量太少，无法计算 PCA，请先完成 pca 步骤")
    try:
        from sklearn.decomposition import PCA, TruncatedSVD
    except ImportError as error:  # pragma: no cover - 环境兜底
        raise EngineError("unsupported", "缺少 scikit-learn，无法计算 PCA：%s" % error)
    if n_cells * len(genes) <= 40_000_000:
        model = PCA(n_components=components, svd_solver="randomized", random_state=DEFAULT_SEED)
        return np.asarray(model.fit_transform(sub.toarray()), dtype=np.float32)
    model = TruncatedSVD(n_components=components, random_state=DEFAULT_SEED)
    return np.asarray(model.fit_transform(sub), dtype=np.float32)


def _embedding(view):
    """取已有二维嵌入；没有则用 PCA 前两维。返回 (coords, label)。"""
    emb = view.get("emb")
    if emb is not None and emb.ndim == 2 and emb.shape[1] >= 2:
        return np.asarray(emb, dtype=np.float32)[:, :2], "existing embedding"
    return _pca(view, n_pcs=2)[:, :2], "PCA (first 2 PCs)"


def _knn(points, k):
    """k 近邻（不含自身），返回 (indices, distances)，确定性 brute 距离。"""
    points = np.ascontiguousarray(points, dtype=np.float32)
    n_cells = points.shape[0]
    k = int(max(1, min(int(k), n_cells - 1)))
    from sklearn.neighbors import NearestNeighbors

    model = NearestNeighbors(n_neighbors=k + 1, algorithm="brute", metric="euclidean", n_jobs=1)
    model.fit(points)
    distances, indices = model.kneighbors(points)
    self_mask = indices != np.arange(n_cells)[:, None]
    order = np.argsort(~self_mask, axis=1, kind="stable")[:, :k]
    keep = np.take_along_axis(indices, order, axis=1)
    keep_distances = np.take_along_axis(distances, order, axis=1).astype(np.float32)
    return keep.astype(np.int64), keep_distances


def _clusters(view, required=False, what="该分析"):
    labels = view.get("clusters")
    if labels is None:
        if required:
            raise EngineError("prerequisite", "需要先完成 cluster 步骤，才能做%s" % what)
        return None
    return np.asarray([str(value) for value in labels], dtype=object)


def _group_means(matrix, labels, categories):
    """簇 × 基因 平均表达（稠密小矩阵）与每簇细胞数。"""
    categories = list(categories)
    position = {category: index for index, category in enumerate(categories)}
    codes = np.asarray([position[str(label)] for label in labels], dtype=np.int64)
    rows = np.arange(len(codes), dtype=np.int64)
    onehot = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (codes, rows)),
        shape=(len(categories), len(codes)),
    )
    sums = np.asarray((onehot @ matrix).toarray(), dtype=np.float64)
    counts = np.bincount(codes, minlength=len(categories)).astype(np.float64)
    means = sums / np.maximum(counts[:, None], 1.0)
    return means, counts, codes


# ---------------------------------------------------------------------------
# 出图 / 表格 / obs_values
# ---------------------------------------------------------------------------


def _dataset_slug(ctx):
    name = ctx.get("name") or ctx.get("dataset_id") or "dataset"
    return _slug(name)


def _save_figure(figure, ctx, outdir, analysis_id, slug, title=None):
    """保存 PNG（白底 dpi=140），返回契约 figure 条目。"""
    out_dir = Path(str(outdir)).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    filename = "%s_%s_%s.png" % (_dataset_slug(ctx), analysis_id, _slug(slug))
    path = out_dir / filename
    figure.savefig(str(path), dpi=140, facecolor="white", bbox_inches="tight")
    width, height = _png_size(path)
    plt.close(figure)
    return {
        "id": _slug(slug),
        "title": str(title or slug),
        "file": str(path.resolve()),
        "width": int(width),
        "height": int(height),
    }


def _new_figure(width=5.6, height=4.4):
    figure, axis = plt.subplots(figsize=(float(width), float(height)), dpi=140)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    return figure, axis


def _table(table_id, title, columns, rows, note=None):
    table = {"id": str(table_id), "title": str(title), "columns": list(columns), "rows": list(rows)}
    if note:
        table["note"] = str(note)
    return table


def _column(key, label, kind="string"):
    return {"key": str(key), "label": str(label), "type": str(kind)}


def _full_numeric(values, idx, n_total):
    """把下采样后的数值铺回全长（缺失 = None）。"""
    output = [None] * int(n_total)
    for position, value in zip(np.asarray(idx, dtype=np.int64), np.asarray(values).ravel()):
        number = float(value)
        output[int(position)] = number if math.isfinite(number) else None
    return output


def _full_categorical(values, idx, n_total):
    output = [None] * int(n_total)
    for position, value in zip(np.asarray(idx, dtype=np.int64), np.asarray(values, dtype=object)):
        output[int(position)] = None if value is None else str(value)
    return output


def _cluster_palette(labels, categories):
    colors = _palette(len(categories))
    mapping = {category: colors[index] for index, category in enumerate(categories)}
    return [mapping[str(label)] for label in labels]


def _scatter_embedding(axis, coords, colors, size=6.0, alpha=0.75):
    axis.scatter(
        coords[:, 0], coords[:, 1], c=colors, s=size, alpha=alpha, linewidths=0, rasterized=True
    )


def _annotate_no_data(axis, text="No data available"):
    axis.text(0.5, 0.5, str(text), ha="center", va="center", transform=axis.transAxes, fontsize=9)
    axis.set_xticks([])
    axis.set_yticks([])

# ###########################################################################
# 1. integration — ComBat 式位置-尺度校正 + BBKNN 式批次平衡 kNN
# ###########################################################################


def _batch_mixing_from_neighbors(neighbors, batch_codes, n_batches):
    """根据邻居索引算批次混合指标（逐细胞同批比例 + 归一化熵）。"""
    neighbor_codes = batch_codes[neighbors]
    per_cell_same = np.mean(neighbor_codes == batch_codes[:, None], axis=1)
    n_cells, k = neighbors.shape
    counts = np.zeros((n_cells, n_batches), dtype=np.float64)
    rows = np.repeat(np.arange(n_cells), k)
    np.add.at(counts, (rows, neighbor_codes.ravel()), 1.0)
    probabilities = counts / float(k)
    entropy = -np.sum(
        np.where(probabilities > 0, probabilities * np.log(probabilities + 1e-12), 0.0), axis=1
    )
    normalized = entropy / math.log(max(n_batches, 2))
    return float(per_cell_same.mean()), per_cell_same, normalized


def _combat_location_scale(points, batch_codes, n_batches):
    """ComBat 式批次位置-尺度校正（无经验贝叶斯收缩），逐维对齐到全局均值/标准差。"""
    corrected = np.array(points, dtype=np.float64, copy=True)
    global_mean = corrected.mean(axis=0)
    global_std = corrected.std(axis=0)
    global_std[global_std < 1e-8] = 1.0
    for code in range(n_batches):
        mask = batch_codes == code
        if mask.sum() < 2:
            continue
        batch_mean = corrected[mask].mean(axis=0)
        batch_std = corrected[mask].std(axis=0)
        batch_std[batch_std < 1e-8] = 1.0
        corrected[mask] = (corrected[mask] - batch_mean) / batch_std * global_std + global_mean
    return corrected


def _bbknn_neighbors(points, batch_codes, k_total):
    """BBKNN 式批次平衡 kNN：从全局候选池中按批次配额取邻居。"""
    n_cells = points.shape[0]
    n_batches = int(batch_codes.max()) + 1 if n_cells else 0
    k_total = int(min(max(2, k_total), n_cells - 1))
    per_batch = max(1, int(math.floor(k_total / float(max(n_batches, 1)))))
    pool = int(min(n_cells - 1, max(k_total * 3, k_total + 10)))
    candidates, _distances = _knn(points, pool)
    neighbors = np.empty((n_cells, k_total), dtype=np.int64)
    for cell in range(n_cells):
        chosen = []
        remaining = [per_batch] * n_batches
        chosen_set = set()
        for candidate in candidates[cell]:
            code = int(batch_codes[candidate])
            if remaining[code] > 0:
                chosen.append(int(candidate))
                chosen_set.add(int(candidate))
                remaining[code] -= 1
                if len(chosen) >= k_total:
                    break
        if len(chosen) < k_total:
            for candidate in candidates[cell]:
                value = int(candidate)
                if value not in chosen_set:
                    chosen.append(value)
                    chosen_set.add(value)
                    if len(chosen) >= k_total:
                        break
        neighbors[cell] = np.asarray(chosen[:k_total], dtype=np.int64)
    return neighbors


def _run_integration(view, params, outdir, warnings):
    obs = view["obs"]
    batch_key = str(params.get("batch") or "batch").strip()
    if not batch_key:
        batch_key = "batch"
    if batch_key not in obs.columns:
        available = [str(column) for column in obs.columns if not str(column).startswith("_")]
        raise EngineError(
            "prerequisite",
            "integration 需要 obs 里的批次列 %r（可用列：%s）；请重新导入并带上 batch 列，"
            "或用 params.batch 指定其它分类列" % (batch_key, "、".join(available[:10]) or "无"),
        )
    batch_values = obs[batch_key].astype(str).to_numpy(dtype=object)
    categories = _sorted_categories(batch_values)
    n_batches = len(categories)
    if n_batches < 2:
        raise EngineError(
            "bad_request",
            "批次列 %r 只有 %d 个取值，无法做批次校正（至少需要 2 个批次）" % (batch_key, n_batches),
        )
    position = {category: index for index, category in enumerate(categories)}
    batch_codes = np.asarray([position[str(value)] for value in batch_values], dtype=np.int64)
    n_neighbors = _clamp_int(params, "n_neighbors", 15, 3, 100)
    n_pcs = _clamp_int(params, "n_pcs", 30, 2, 100)
    if n_batches > n_neighbors:
        warnings.append(
            "批次数 %d 大于近邻数 k=%d：BBKNN 式配额每个批次最多 1 个邻居，批次平衡受限，建议增大 n_neighbors"
            % (n_batches, n_neighbors)
        )

    points = _pca(view, n_pcs=n_pcs)
    corrected = _combat_location_scale(points, batch_codes, n_batches)

    # before：批次校正前 PCA 上的普通 kNN
    before_neighbors, _ = _knn(points, n_neighbors)
    before_same, before_same_cells, before_entropy = _batch_mixing_from_neighbors(
        before_neighbors, batch_codes, n_batches
    )
    # after：BBKNN 式批次平衡 kNN（每个批次等配额）
    after_neighbors = _bbknn_neighbors(corrected, batch_codes, n_neighbors)
    after_same, after_same_cells, after_entropy = _batch_mixing_from_neighbors(
        after_neighbors, batch_codes, n_batches
    )

    before_embedding = points[:, :2]
    after_embedding = corrected[:, :2]

    rows = []
    for code, category in enumerate(categories):
        mask = batch_codes == code
        rows.append(
            {
                "batch": category,
                "n_cells": int(mask.sum()),
                "same_batch_frac_before": _num(np.mean(before_same_cells[mask]), 4),
                "same_batch_frac_after": _num(np.mean(after_same_cells[mask]), 4),
                "batch_entropy_before": _num(np.mean(before_entropy[mask]), 4),
                "batch_entropy_after": _num(np.mean(after_entropy[mask]), 4),
            }
        )
    table = _table(
        "batch_mixing",
        "Batch mixing metrics before / after correction",
        [
            _column("batch", "Batch"),
            _column("n_cells", "Cells", "number"),
            _column("same_batch_frac_before", "Same-batch neighbor frac (before)", "number"),
            _column("same_batch_frac_after", "Same-batch neighbor frac (after)", "number"),
            _column("batch_entropy_before", "Batch entropy (before)", "number"),
            _column("batch_entropy_after", "Batch entropy (after)", "number"),
        ],
        rows,
        note="熵按 log(批次数) 归一化；越高表示邻居越跨批次。BBKNN 式配额邻居列表用于 after。",
    )

    figure, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), dpi=140, squeeze=False)
    colors = _palette(n_batches)
    batch_colors = [colors[code] for code in batch_codes]
    for axis, coords, title, same, entropy in (
        (
            axes[0][0],
            before_embedding,
            "Before correction (PCA, colored by batch)",
            before_same,
            before_entropy,
        ),
        (
            axes[0][1],
            after_embedding,
            "After ComBat-style correction (colored by batch)",
            after_same,
            after_entropy,
        ),
    ):
        axis.scatter(
            coords[:, 0],
            coords[:, 1],
            c=batch_colors,
            s=5,
            alpha=0.7,
            linewidths=0,
            rasterized=True,
        )
        axis.set_xlabel("PC1")
        axis.set_ylabel("PC2")
        axis.set_title(
            "%s\nsame-batch=%.2f, batch entropy=%.2f"
            % (title, float(np.mean(same)), float(np.mean(entropy))),
            fontsize=8.5,
        )
        handles = [
            plt.Line2D([], [], marker="o", linestyle="", color=colors[index], label=str(category))
            for index, category in enumerate(categories)
        ]
        axis.legend(handles=handles, frameon=False, loc="best", markerscale=1.4)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "integration", "batch_embedding", "Batch effect before/after"
    )

    # 校正后 PC1/PC2 + 混合熵写回 obs（供嵌入图着色）
    idx = view["idx"]
    n_total = view["n_total"]
    obs_values = {
        "integration_pc1": _full_numeric(after_embedding[:, 0], idx, n_total),
        "integration_pc2": _full_numeric(after_embedding[:, 1], idx, n_total),
        "integration_batch_entropy_before": _full_numeric(before_entropy, idx, n_total),
        "integration_batch_entropy_after": _full_numeric(after_entropy, idx, n_total),
    }
    mean_before_entropy = float(np.mean(before_entropy))
    mean_after_entropy = float(np.mean(after_entropy))
    effective_batches_after = float(math.exp(mean_after_entropy * math.log(max(n_batches, 2))))
    bullets = [
        "批次列：obs.%s（%d 个批次，%d 个细胞）" % (batch_key, n_batches, len(batch_values)),
        "ComBat 式位置-尺度校正（本插件实现，无经验贝叶斯收缩）：同一批次的每个 PC 维对齐到全局均值/标准差",
        "BBKNN 式批次平衡 kNN（本插件实现）：每个细胞按批次等配额挑选邻居，候选池=全局 kNN",
        "批次混合：同批邻居比例 %.3f → %.3f；归一化批次熵 %.3f → %.3f（1=完全混合）"
        % (np.mean(before_same), np.mean(after_same), mean_before_entropy, mean_after_entropy),
        "校正后 PC1/PC2 与批次熵已写入 obs（integration_pc1/pc2、integration_batch_entropy_*），可直接用于嵌入图着色",
    ]
    summary = {
        "bullets": bullets,
        "metrics": {
            "batch_key": batch_key,
            "n_batches": n_batches,
            "same_batch_frac_before": _num(np.mean(before_same), 4),
            "same_batch_frac_after": _num(np.mean(after_same), 4),
            "batch_entropy_before": _num(mean_before_entropy, 4),
            "batch_entropy_after": _num(mean_after_entropy, 4),
            "effective_batches_after": _num(effective_batches_after, 3),
            "n_pcs": int(points.shape[1]),
            "n_neighbors": int(n_neighbors),
        },
    }
    return {
        "summary": summary,
        "tables": [table],
        "figures": [figure_entry],
        "obs_columns": [
            {"key": "integration_pc1", "kind": "numeric"},
            {"key": "integration_pc2", "kind": "numeric"},
            {"key": "integration_batch_entropy_before", "kind": "numeric"},
            {"key": "integration_batch_entropy_after", "kind": "numeric"},
        ],
        "obs_values": obs_values,
    }


# ###########################################################################
# 2. trajectory_dpt — 扩散图 + DPT 伪时间
# ###########################################################################


def _diffusion_components(points, n_neighbors=15, n_components=4):
    """扩散图：自适应高斯核 kNN 图 → 对称归一化 → 前若干特征向量。

    返回 (psi, eigenvalues)：psi 为扩散分量（psi[:,0] 近似常数），
    eigenvalues 降序；psi[:,l] = v_l / sqrt(degree)。
    """
    n_cells = points.shape[0]
    k = int(min(max(3, n_neighbors), max(3, n_cells - 2)))
    neighbors, distances = _knn(points, k)
    rows = np.repeat(np.arange(n_cells), k)
    cols = neighbors.ravel()
    local_scale = np.maximum(distances[:, -1], 1e-6)
    variance = distances.ravel() ** 2
    denominator = local_scale[rows] * local_scale[cols] + 1e-12
    weights = np.exp(-variance / denominator).astype(np.float64)
    kernel = sp.csr_matrix((weights, (rows, cols)), shape=(n_cells, n_cells))
    kernel = kernel.maximum(kernel.T)
    degree = np.asarray(kernel.sum(axis=1)).ravel()
    degree = np.maximum(degree, 1e-12)
    inv_sqrt = 1.0 / np.sqrt(degree)
    symmetric = sp.diags(inv_sqrt) @ kernel @ sp.diags(inv_sqrt)
    symmetric = symmetric.tocsr()
    total = int(min(n_components + 1, n_cells - 2))
    total = max(2, total)
    eigenvalues = None
    eigenvectors = None
    try:
        from scipy.sparse.linalg import eigsh

        # ARPACK 默认随机 v0 → 结果可能不可复现；这里固定起始向量保证确定性
        start_vector = np.random.default_rng(DEFAULT_SEED).normal(size=n_cells)
        eigenvalues, eigenvectors = eigsh(symmetric, k=total, which="LA", v0=start_vector)
    except Exception:  # pragma: no cover - 小矩阵/收敛失败时回退稠密分解
        dense = np.asarray(symmetric.toarray(), dtype=np.float64)
        values, vectors = np.linalg.eigh(dense)
        eigenvalues = values[-total:]
        eigenvectors = vectors[:, -total:]
    order = np.argsort(-np.asarray(eigenvalues))
    eigenvalues = np.asarray(eigenvalues)[order]
    eigenvectors = np.asarray(eigenvectors)[:, order]
    psi = eigenvectors / np.sqrt(degree)[:, None]
    return psi, eigenvalues


def _run_trajectory_dpt(view, params, outdir, warnings):
    n_neighbors = _clamp_int(params, "n_neighbors", 15, 3, 100)
    n_dcs = _clamp_int(params, "n_dcs", 4, 2, 10)
    points = _pca(view, n_pcs=30)
    n_cells = points.shape[0]
    if n_cells < 10:
        raise EngineError("prerequisite", "细胞数太少（<10），无法构建扩散图，请先完成 pca/neighbors")
    psi, eigenvalues = _diffusion_components(points, n_neighbors=n_neighbors, n_components=n_dcs)
    n_components = min(int(psi.shape[1]) - 1, int(n_dcs))
    if n_components < 1:
        raise EngineError("prerequisite", "扩散图分量不足，无法计算伪时间")

    # 根细胞：params.root_cell 指定（obs 索引 / 位置），否则取第一扩散分量极值
    obs_index = [str(value) for value in view["obs"].index]
    root_text = str(params.get("root_cell") or "").strip()
    root_mode = str(params.get("root_mode") or "auto").strip().lower()
    if root_text:
        if root_text in obs_index:
            root = int(obs_index.index(root_text))
        else:
            try:
                root = int(root_text)
            except ValueError:
                root = -1
            if root < 0 or root >= n_cells:
                raise EngineError(
                    "bad_request",
                    "找不到根细胞 %r（需要 obs 索引名，或 0-based 位置 0..%d）" % (root_text, n_cells - 1),
                )
        root_source = "params.root_cell"
    elif root_mode == "max_genes":
        detected = np.asarray((view["X"] > 0).sum(axis=1)).ravel()
        root = int(np.argmax(detected))
        root_source = "auto: max detected genes"
    else:
        component = psi[:, 1]
        root = int(np.argmax(component))
        root_source = "auto: extreme of first diffusion component"

    weights = []
    for component in range(1, n_components + 1):
        value = float(eigenvalues[component]) if component < len(eigenvalues) else 0.0
        value = max(min(value, 0.999999), 0.0)
        weights.append(value / max(1.0 - value, 1e-6))
    weights = np.asarray(weights, dtype=np.float64)
    diffusion = psi[:, 1 : n_components + 1]
    delta = diffusion - diffusion[root : root + 1, :]
    distances = np.sqrt(np.sum((delta * weights[None, :]) ** 2, axis=1))
    if distances.max() <= 1e-12:
        warnings.append("根细胞到所有细胞的扩散距离为 0：扩散图可能连通性异常，伪时间退化为常数")
        pseudotime = np.zeros(n_cells, dtype=np.float64)
    else:
        pseudotime = distances / distances.max()
    pseudotime = np.clip(pseudotime, 0.0, 1.0)

    clusters = _clusters(view, required=False)
    # 分支：在加权扩散分量上做 KMeans（近似分支，不是 monocle 的分支点检测）
    n_branches = 3
    if clusters is not None:
        n_branches = int(min(max(2, len(_sorted_categories(clusters)) // 2), 4))
    n_branches = int(min(max(2, n_branches), max(2, n_cells // 10)))
    branch_labels = None
    try:
        from sklearn.cluster import KMeans

        features = diffusion * weights[None, :]
        model = KMeans(n_clusters=n_branches, random_state=DEFAULT_SEED, n_init=10)
        codes = model.fit_predict(features)
        branch_labels = np.asarray(["branch%d" % (int(code) + 1) for code in codes], dtype=object)
    except Exception as error:  # pragma: no cover - 防御
        warnings.append("分支检测（KMeans）失败，已跳过分支：%s" % error)

    # 簇级伪时间排序表
    tables = []
    if clusters is not None:
        categories = _sorted_categories(clusters)
        rows = []
        for category in categories:
            mask = clusters == category
            rows.append(
                {
                    "cluster": category,
                    "n_cells": int(mask.sum()),
                    "mean_pseudotime": _num(np.mean(pseudotime[mask]), 4),
                    "median_pseudotime": _num(np.median(pseudotime[mask]), 4),
                    "min_pseudotime": _num(np.min(pseudotime[mask]), 4),
                    "max_pseudotime": _num(np.max(pseudotime[mask]), 4),
                }
            )
        rows.sort(key=lambda row: (row["mean_pseudotime"] is None, row["mean_pseudotime"] or 0.0))
        for order_index, row in enumerate(rows):
            row["order"] = order_index + 1
        for row in rows:
            row["order"] = row.pop("order")
        tables.append(
            _table(
                "pseudotime_by_cluster",
                "Pseudotime by cluster (sorted)",
                [
                    _column("order", "Order", "number"),
                    _column("cluster", "Cluster"),
                    _column("n_cells", "Cells", "number"),
                    _column("mean_pseudotime", "Mean pseudotime", "number"),
                    _column("median_pseudotime", "Median pseudotime", "number"),
                    _column("min_pseudotime", "Min", "number"),
                    _column("max_pseudotime", "Max", "number"),
                ],
                rows,
                note="伪时间 = 根细胞在加权扩散空间中的欧氏距离，归一化到 [0,1]。",
            )
        )
    if branch_labels is not None:
        branch_categories = _sorted_categories(branch_labels)
        branch_rows = []
        for category in branch_categories:
            mask = branch_labels == category
            top_clusters = ""
            if clusters is not None:
                values, counts = np.unique(clusters[mask], return_counts=True)
                top = values[np.argsort(-counts)][:3]
                top_clusters = ", ".join(str(value) for value in top)
            branch_rows.append(
                {
                    "branch": category,
                    "n_cells": int(mask.sum()),
                    "mean_pseudotime": _num(np.mean(pseudotime[mask]), 4),
                    "top_clusters": top_clusters,
                }
            )
        tables.append(
            _table(
                "branches",
                "Diffusion branches (KMeans on diffusion components)",
                [
                    _column("branch", "Branch"),
                    _column("n_cells", "Cells", "number"),
                    _column("mean_pseudotime", "Mean pseudotime", "number"),
                    _column("top_clusters", "Top clusters"),
                ],
                branch_rows,
                note="分支为扩散分量上的 KMeans 分群，属于近似分支划分，不等价于 monocle 分支点检验。",
            )
        )

    # 图：嵌入按伪时间着色 + 簇质心轨迹线
    embedding, embedding_label = _embedding(view)
    figure, axis = _new_figure(6.4, 5.2)
    scatter = axis.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=pseudotime,
        cmap="viridis",
        s=7,
        alpha=0.85,
        linewidths=0,
        rasterized=True,
    )
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.03)
    colorbar.set_label("DPT pseudotime")
    axis.scatter(
        [embedding[root, 0]],
        [embedding[root, 1]],
        marker="*",
        s=170,
        c="#E45756",
        edgecolors="white",
        linewidths=0.8,
        zorder=6,
        label="root cell",
    )
    if clusters is not None:
        categories = _sorted_categories(clusters)
        centroids = []
        for category in categories:
            mask = clusters == category
            centroids.append(
                (
                    float(np.mean(embedding[mask, 0])),
                    float(np.mean(embedding[mask, 1])),
                    float(np.mean(pseudotime[mask])),
                    category,
                )
            )
        centroids.sort(key=lambda item: item[2])
        if len(centroids) > 1:
            coordinates = np.asarray([[item[0], item[1]] for item in centroids])
            axis.plot(
                coordinates[:, 0],
                coordinates[:, 1],
                color="#333333",
                linewidth=1.2,
                alpha=0.85,
                zorder=5,
            )
            for index in range(len(centroids) - 1):
                start = coordinates[index]
                end = coordinates[index + 1]
                axis.annotate(
                    "",
                    xy=(end[0], end[1]),
                    xytext=(start[0], start[1]),
                    arrowprops=dict(arrowstyle="-|>", color="#333333", lw=0.9),
                    zorder=5,
                )
        for item in centroids:
            axis.text(item[0], item[1], str(item[3]), fontsize=7.5, ha="center", va="center", zorder=7,
                      bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#999999", alpha=0.75, lw=0.5))
    axis.set_xlabel("%s 1" % embedding_label)
    axis.set_ylabel("%s 2" % embedding_label)
    axis.set_title("Diffusion pseudotime (%d cells, root=%s)" % (n_cells, root_source), fontsize=9)
    axis.legend(frameon=False, loc="best")
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "trajectory_dpt", "pseudotime", "Diffusion pseudotime"
    )

    idx = view["idx"]
    n_total = view["n_total"]
    obs_values = {"dpt_pseudotime": _full_numeric(pseudotime, idx, n_total)}
    obs_columns = [{"key": "dpt_pseudotime", "kind": "numeric"}]
    if branch_labels is not None:
        obs_values["dpt_branch"] = _full_categorical(branch_labels, idx, n_total)
        obs_columns.append({"key": "dpt_branch", "kind": "categorical"})

    bullets = [
        "扩散图：%d 个细胞、k=%d 的自适应高斯核图，取 %d 个扩散分量（本插件实现，非官方 DPT/monocle）"
        % (n_cells, n_neighbors, n_components),
        "根细胞：%s（位置 %d）" % (root_source, root),
        "伪时间范围 [0,1]；根细胞=%0.2f，最远细胞=%0.2f" % (float(pseudotime[root]), float(pseudotime.max())),
    ]
    if branch_labels is not None:
        bullets.append("分支：在加权扩散分量上做 KMeans（%d 类）得到近似分支标签 dpt_branch" % n_branches)
    bullets.append("伪时间已写入 obs 列 dpt_pseudotime" + ("，分支写入 dpt_branch" if branch_labels is not None else ""))
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_cells": int(n_cells),
            "n_neighbors": int(n_neighbors),
            "n_diffusion_components": int(n_components),
            "root_index": int(root),
            "root_source": root_source,
            "pseudotime_min": _num(float(pseudotime.min()), 4),
            "pseudotime_max": _num(float(pseudotime.max()), 4),
            "n_branches": int(n_branches) if branch_labels is not None else None,
        },
    }
    return {
        "summary": summary,
        "tables": tables,
        "figures": [figure_entry],
        "obs_columns": obs_columns,
        "obs_values": obs_values,
    }


# ###########################################################################
# 3. trajectory_paga — 簇间连通度图 + MST/ELG + 簇级伪时间
# ###########################################################################


def _separate_points(points, min_distance, iterations=80):
    """把重合的二维点确定性推开（用于簇级图布局，避免节点互相遮挡）。"""
    positions = np.array(points, dtype=np.float64, copy=True)
    n_points = positions.shape[0]
    if n_points <= 1 or min_distance <= 0:
        return positions
    for _ in range(int(iterations)):
        moved = False
        for left in range(n_points):
            for right in range(left + 1, n_points):
                delta = positions[right] - positions[left]
                distance = float(np.linalg.norm(delta))
                if distance >= min_distance:
                    continue
                if distance < 1e-9:
                    angle = 2.0 * math.pi * ((left + 1) * 37 + (right + 1) * 11) / 97.0
                    direction = np.asarray([math.cos(angle), math.sin(angle)])
                    distance = 1e-9
                else:
                    direction = delta / distance
                push = (min_distance - distance) / 2.0
                positions[left] -= direction * push
                positions[right] += direction * push
                moved = True
        if not moved:
            break
    return positions


def _maximum_spanning_tree(connectivity):
    """Prim 最大生成树（k×k 稠密矩阵，确定性 tie-break）。"""
    n = connectivity.shape[0]
    visited = np.zeros(n, dtype=bool)
    visited[0] = True
    edges = []
    for _ in range(n - 1):
        best = None
        for left in range(n):
            if not visited[left]:
                continue
            for right in range(n):
                if visited[right]:
                    continue
                score = float(connectivity[left, right])
                key = (-score, left, right)
                if best is None or key < best[0]:
                    best = (key, left, right, score)
        if best is None:
            break
        _key, left, right, score = best
        visited[right] = True
        edges.append((left, right, score))
    return edges


def _dijkstra_on_tree(n, edges):
    """MST 上的最短路径距离（边权 = 1/max(connectivity, floor)）。"""
    adjacency = [[] for _ in range(n)]
    for left, right, score in edges:
        length = 1.0 / max(float(score), 0.05)
        adjacency[left].append((right, length))
        adjacency[right].append((left, length))
    distances = np.full(n, np.inf)
    distances[0] = 0.0
    visited = np.zeros(n, dtype=bool)
    for _ in range(n):
        candidates = np.where(~visited, distances, np.inf)
        node = int(np.argmin(candidates))
        if not math.isfinite(float(candidates[node])):
            break
        visited[node] = True
        for neighbor, length in adjacency[node]:
            if distances[node] + length < distances[neighbor]:
                distances[neighbor] = distances[node] + length
    finite = np.isfinite(distances)
    if finite.any() and distances[finite].max() > 0:
        distances = np.where(finite, distances, distances[finite].max())
        distances = distances / distances.max()
    else:
        distances = np.zeros(n, dtype=np.float64)
    return distances


def _run_trajectory_paga(view, params, outdir, warnings):
    clusters = _clusters(view, required=True, what=" PAGA-like 簇级轨迹")
    categories = _sorted_categories(clusters)
    k_clusters = len(categories)
    if k_clusters < 2:
        raise EngineError("bad_request", "只有 1 个簇，无法构建簇间轨迹图（PAGA-like）")
    n_neighbors = _clamp_int(params, "n_neighbors", 15, 3, 100)
    points = _pca(view, n_pcs=20)
    n_cells = points.shape[0]
    neighbors, _distances = _knn(points, n_neighbors)
    position = {category: index for index, category in enumerate(categories)}
    codes = np.asarray([position[str(label)] for label in clusters], dtype=np.int64)
    neighbor_codes = codes[neighbors]

    counts_matrix = np.zeros((k_clusters, k_clusters), dtype=np.float64)
    rows = np.repeat(np.arange(n_cells), neighbors.shape[1])
    np.add.at(
        counts_matrix,
        (codes[rows], neighbor_codes.ravel()),
        1.0,
    )
    observed = 0.5 * (counts_matrix + counts_matrix.T)
    np.fill_diagonal(observed, 0.0)
    cluster_sizes = np.bincount(codes, minlength=k_clusters).astype(np.float64)
    expected = np.zeros((k_clusters, k_clusters), dtype=np.float64)
    for left in range(k_clusters):
        for right in range(k_clusters):
            if left != right:
                expected[left, right] = (
                    n_neighbors * cluster_sizes[left] * cluster_sizes[right] / max(n_cells - 1, 1)
                )
    connectivity = np.where(expected > 1e-9, observed / np.maximum(expected, 1e-9), 0.0)
    np.fill_diagonal(connectivity, 0.0)

    # 簇在 PCA 空间完全分离时，k=15 的 kNN 图可能没有任何跨簇边 → connectivity 全 0。
    # 此时退化为"簇质心距离"近似：connectivity = exp(-距离/中位距离)，并在 warnings 说明。
    centroid_pca = np.zeros((k_clusters, 2), dtype=np.float64)
    for index in range(k_clusters):
        mask = codes == index
        centroid_pca[index] = np.mean(points[mask, :2], axis=0)
    connectivity_source = "knn_edges"
    off_diagonal = ~np.eye(k_clusters, dtype=bool)
    positive_fraction = float(np.mean(connectivity[off_diagonal] > 0)) if k_clusters > 1 else 0.0
    if positive_fraction < 0.25 or float(np.max(connectivity)) <= 1e-9:
        pairwise_distance = np.linalg.norm(
            centroid_pca[:, None, :] - centroid_pca[None, :, :], axis=2
        )
        positive = pairwise_distance[pairwise_distance > 1e-9]
        scale = float(np.median(positive)) if positive.size else 1.0
        connectivity = np.exp(-pairwise_distance / max(scale, 1e-9))
        np.fill_diagonal(connectivity, 0.0)
        connectivity_source = "centroid_distance_fallback"
        warnings.append(
            "k=%d 时只有 %.0f%% 的簇对存在跨簇 kNN 边（簇在 PCA 空间高度分离）：PAGA connectivity 证据不足，"
            "已退化为簇质心距离近似（connectivity=exp(-距离/中位距离)）；connectivity 表仍保留观测边数=0"
            % (n_neighbors, positive_fraction * 100.0)
        )

    mst_edges = _maximum_spanning_tree(connectivity)
    mst_pairs = {tuple(sorted((left, right))) for left, right, _score in mst_edges}
    threshold = float(np.max(connectivity)) * 0.5 if np.max(connectivity) > 0 else 0.0
    elg_pairs = {
        (left, right)
        for left in range(k_clusters)
        for right in range(left + 1, k_clusters)
        if connectivity[left, right] >= max(threshold, 1e-9)
    }

    # 根簇：params.root_cluster 或自动选"检测基因数最多"的簇（最不分化）
    detected = np.asarray((view["X"] > 0).sum(axis=1)).ravel()
    cluster_median_genes = np.asarray(
        [float(np.median(detected[codes == index])) for index in range(k_clusters)]
    )
    root_text = str(params.get("root_cluster") or "").strip()
    if root_text:
        if root_text not in position:
            raise EngineError(
                "bad_request",
                "找不到根簇 %r（可用簇：%s）" % (root_text, "、".join(categories)),
            )
        root = position[root_text]
        root_source = "params.root_cluster"
    else:
        root = int(np.argmax(cluster_median_genes))
        root_source = "auto: cluster with most detected genes"
    # 重排簇顺序，让根簇成为 MST 起点（保证伪时间从根开始）
    reorder = [root] + [index for index in range(k_clusters) if index != root]
    remap = {old: new for new, old in enumerate(reorder)}
    ordered_connectivity = connectivity[np.ix_(reorder, reorder)]
    ordered_edges = _maximum_spanning_tree(ordered_connectivity)
    ordered_distances = _dijkstra_on_tree(k_clusters, ordered_edges)
    paga_pseudotime = np.zeros(k_clusters, dtype=np.float64)
    for new, old in enumerate(reorder):
        paga_pseudotime[old] = ordered_distances[new]

    # 边表
    edge_rows = []
    for left in range(k_clusters):
        for right in range(left + 1, k_clusters):
            value = float(connectivity[left, right])
            if value <= 0 and (left, right) not in mst_pairs:
                continue
            pair = tuple(sorted((left, right)))
            edge_rows.append(
                {
                    "source": categories[left],
                    "target": categories[right],
                    "n_edges": _num(observed[left, right], 1),
                    "expected_edges": _num(expected[left, right], 2),
                    "connectivity": _num(value, 4),
                    "mst": bool(pair in mst_pairs),
                    "elg": bool(pair in elg_pairs),
                }
            )
    edge_rows.sort(key=lambda row: (-(row["connectivity"] or 0.0), row["source"], row["target"]))
    edge_table = _table(
        "cluster_connectivity",
        "Cluster-cluster connectivity (observed / expected kNN edges)",
        [
            _column("source", "Source"),
            _column("target", "Target"),
            _column("n_edges", "Observed edges", "number"),
            _column("expected_edges", "Expected edges", "number"),
            _column("connectivity", "Connectivity", "number"),
            _column("mst", "In MST"),
            _column("elg", "In ELG"),
        ],
        edge_rows,
        note=(
            "connectivity = 观察到的簇间 kNN 边 / 随机期望；MST=最大生成树，ELG=connectivity ≥ 50%% 最大值的边。"
            + (
                " 本次簇间无 kNN 边，connectivity 为簇质心距离的指数近似（见 warnings）。"
                if connectivity_source != "knn_edges"
                else ""
            )
        ),
    )
    cluster_rows = []
    for index, category in enumerate(categories):
        degree = 0.0
        for left, right, _score in mst_edges:
            if index in (left, right):
                degree += 1.0
        cluster_rows.append(
            {
                "cluster": category,
                "n_cells": int(cluster_sizes[index]),
                "paga_pseudotime": _num(paga_pseudotime[index], 4),
                "median_detected_genes": _num(cluster_median_genes[index], 1),
                "mst_degree": int(degree),
                "root": bool(index == root),
            }
        )
    cluster_rows.sort(key=lambda row: (row["paga_pseudotime"] or 0.0, row["cluster"]))
    for order_index, row in enumerate(cluster_rows):
        row["order"] = order_index + 1
    cluster_table = _table(
        "cluster_pseudotime",
        "Cluster-level pseudotime (MST distance from root)",
        [
            _column("order", "Order", "number"),
            _column("cluster", "Cluster"),
            _column("n_cells", "Cells", "number"),
            _column("paga_pseudotime", "PAGA-like pseudotime", "number"),
            _column("median_detected_genes", "Median detected genes", "number"),
            _column("mst_degree", "MST degree", "number"),
            _column("root", "Root"),
        ],
        cluster_rows,
        note="簇级伪时间 = 根簇在最大生成树上的最短路径（边权 = 1/connectivity），归一化到 [0,1]。",
    )

    # 图：簇间网络（节点位置 = 簇质心前两维 PCA；退化时用环形布局）
    centroids = np.array(centroid_pca, dtype=np.float64, copy=True)
    spread = float(np.std(centroids, axis=0).max()) if k_clusters > 1 else 0.0
    if spread < 1e-6:
        angles = np.linspace(0, 2 * math.pi, k_clusters, endpoint=False)
        centroids = np.column_stack([np.cos(angles), np.sin(angles)])
        spread = 1.0
    else:
        # 簇质心可能重合（重复/亚群簇）：确定性推开，避免节点互相遮挡
        centroids = _separate_points(centroids, min_distance=0.22 * float(spread))

    figure, axis = _new_figure(6.6, 5.6)
    maximum = max(float(np.max(connectivity)), 1e-9)
    zero_connectivity_mst = 0
    for left in range(k_clusters):
        for right in range(left + 1, k_clusters):
            value = float(connectivity[left, right])
            pair = tuple(sorted((left, right)))
            is_mst = pair in mst_pairs
            if value <= 0 and not is_mst:
                continue
            if is_mst and value <= 0:
                zero_connectivity_mst += 1
            axis.plot(
                [centroids[left, 0], centroids[right, 0]],
                [centroids[left, 1], centroids[right, 1]],
                color="#333333" if is_mst else "#BBBBBB",
                linewidth=(0.9 if (is_mst and value <= 0) else 0.6 + 2.6 * (value / maximum))
                if is_mst
                else 0.3 + 1.0 * (value / maximum),
                linestyle=":" if (is_mst and value <= 0) else "-",
                alpha=0.9 if is_mst else 0.45,
                zorder=1 if not is_mst else 2,
            )
    scatter = axis.scatter(
        centroids[:, 0],
        centroids[:, 1],
        s=120.0 + 900.0 * (cluster_sizes / max(cluster_sizes.max(), 1.0)),
        c=paga_pseudotime,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.6,
        zorder=3,
    )
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.03)
    colorbar.set_label("PAGA-like pseudotime")
    for index, category in enumerate(categories):
        axis.text(
            centroids[index, 0],
            centroids[index, 1],
            str(category),
            fontsize=8,
            ha="center",
            va="center",
            zorder=4,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7),
        )
    axis.set_xlabel("PC1 (cluster centroid)")
    axis.set_ylabel("PC2 (cluster centroid)")
    axis.set_title(
        "PAGA-like cluster graph: %d clusters, %d MST edges, root=%s"
        % (k_clusters, len(mst_edges), categories[root]),
        fontsize=9,
    )
    axis.set_xticks([])
    axis.set_yticks([])
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "trajectory_paga", "cluster_graph", "PAGA-like cluster graph"
    )

    # 每细胞也写入其簇的伪时间
    idx = view["idx"]
    n_total = view["n_total"]
    cell_pseudotime = np.asarray([paga_pseudotime[position[str(label)]] for label in clusters])
    obs_values = {"paga_pseudotime": _full_numeric(cell_pseudotime, idx, n_total)}

    bullets = [
        "%d 个簇、%d 条 kNN 簇间边；connectivity = 观察边数 / 随机期望（本插件实现，非官方 PAGA）"
        % (k_clusters, len(edge_rows)),
        "最大生成树 %d 条边；ELG 阈值 %.2f（最大 connectivity 的 50%%）" % (len(mst_edges), threshold),
        "根簇 %r（%s），簇级伪时间从根簇沿 MST 计算" % (categories[root], root_source),
        "每细胞的簇级伪时间已写入 obs 列 paga_pseudotime",
    ]
    if connectivity_source != "knn_edges":
        bullets.append(
            "注意：簇间无 kNN 边，connectivity 使用簇质心距离指数近似（不是官方 PAGA 的边计数/期望比）"
        )
    elif zero_connectivity_mst:
        warnings.append(
            "%d/%d 条 MST 边 connectivity=0（k=%d 时跨簇 kNN 边稀疏）：图中以虚线示意这些连接，"
            "其方向/距离证据较弱" % (zero_connectivity_mst, len(mst_edges), n_neighbors)
        )
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_clusters": int(k_clusters),
            "n_edges": int(len(edge_rows)),
            "n_mst_edges": int(len(mst_edges)),
            "n_elg_edges": int(len(elg_pairs)),
            "root_cluster": str(categories[root]),
            "root_source": root_source,
            "connectivity_source": connectivity_source,
            "max_connectivity": _num(maximum, 4),
        },
    }
    return {
        "summary": summary,
        "tables": [edge_table, cluster_table],
        "figures": [figure_entry],
        "obs_columns": [{"key": "paga_pseudotime", "kind": "numeric"}],
        "obs_values": obs_values,
    }


# ###########################################################################
# 4. entropy_stemness — 表达熵 + 邻居平滑
# ###########################################################################


def _expression_entropy(matrix_genes):
    """每细胞在给定基因集上的归一化 Shannon 熵（稀疏行内计算）。"""
    matrix = matrix_genes.tocsr()
    n_cells = matrix.shape[0]
    row_sums = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    row_lengths = np.diff(matrix.indptr)
    if matrix.nnz == 0:
        return np.zeros(n_cells, dtype=np.float64)
    row_sum_repeat = np.repeat(row_sums, row_lengths)
    probabilities = matrix.data.astype(np.float64) / np.maximum(row_sum_repeat, 1e-12)
    contribution = np.where(probabilities > 0, probabilities * np.log(probabilities + 1e-300), 0.0)
    entropy = -np.bincount(
        np.repeat(np.arange(n_cells), row_lengths), weights=contribution, minlength=n_cells
    )
    # 基因数越多熵上限越大：按 log(有效基因数) 归一化（用平均检出基因数近似有效维度）
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel().astype(np.float64)
    normalizer = np.log(np.maximum(detected, 2))
    return np.clip(entropy / np.maximum(normalizer, 1e-9), 0.0, 1.0)


def _run_entropy_stemness(view, params, outdir, warnings):
    n_genes = _clamp_int(params, "n_genes", 2000, 100, MAX_GENES)
    n_neighbors = _clamp_int(params, "n_neighbors", 15, 3, 100)
    genes = _hvg_indices(view["X"], n_genes)
    expression = _expression_entropy(view["X"][:, genes])
    points = _pca(view, n_pcs=20)
    neighbors, _distances = _knn(points, n_neighbors)
    smoothed = expression[neighbors].mean(axis=1)
    combined = 0.5 * _zscore(expression) + 0.5 * _zscore(smoothed)
    score = _rank01(combined)

    clusters = _clusters(view, required=False)
    tables = []
    if clusters is not None:
        categories = _sorted_categories(clusters)
        rows = []
        for category in categories:
            mask = clusters == category
            rows.append(
                {
                    "cluster": category,
                    "n_cells": int(mask.sum()),
                    "mean_score": _num(np.mean(score[mask]), 4),
                    "median_score": _num(np.median(score[mask]), 4),
                    "mean_expression_entropy": _num(np.mean(expression[mask]), 4),
                }
            )
        rows.sort(key=lambda row: -(row["mean_score"] or 0.0))
        tables.append(
            _table(
                "entropy_by_cluster",
                "Entropy stemness score by cluster",
                [
                    _column("cluster", "Cluster"),
                    _column("n_cells", "Cells", "number"),
                    _column("mean_score", "Mean entropy score", "number"),
                    _column("median_score", "Median score", "number"),
                    _column("mean_expression_entropy", "Mean expression entropy", "number"),
                ],
                rows,
                note="score=0.5*z(表达熵)+0.5*z(邻居平滑熵)，再秩变换到 [0,1]；高=表达更弥散、干性样。",
            )
        )
    else:
        tables.append(
            _table(
                "entropy_summary",
                "Entropy score summary",
                [_column("metric", "Metric"), _column("value", "Value", "number")],
                [
                    {"metric": "n_cells", "value": int(len(score))},
                    {"metric": "mean_score", "value": _num(np.mean(score), 4)},
                    {"metric": "mean_expression_entropy", "value": _num(np.mean(expression), 4)},
                ],
            )
        )

    embedding, embedding_label = _embedding(view)
    figure, axis = _new_figure(5.9, 4.9)
    scatter = axis.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=score,
        cmap="magma",
        s=7,
        alpha=0.85,
        linewidths=0,
        rasterized=True,
    )
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.03)
    colorbar.set_label("entropy stemness score")
    axis.set_xlabel("%s 1" % embedding_label)
    axis.set_ylabel("%s 2" % embedding_label)
    axis.set_title("Entropy-based stemness score", fontsize=9)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "entropy_stemness", "entropy_score", "Entropy stemness score"
    )

    idx = view["idx"]
    n_total = view["n_total"]
    bullets = [
        "表达熵：每细胞在 %d 个高变基因上的归一化 Shannon 熵（稀疏行内计算）" % len(genes),
        "邻居平滑：k=%d 的 PCA kNN 均值（本插件启发式，非已发表干性打分）" % n_neighbors,
        "综合打分 = 0.5*z(表达熵) + 0.5*z(平滑熵)，秩变换到 [0,1]；高=表达更弥散（干性样）",
        "打分已写入 obs 列 entropy_stemness_score，原始表达熵写入 expression_entropy",
    ]
    if clusters is not None:
        top = max(tables[0]["rows"], key=lambda row: row["mean_score"] or 0.0) if tables[0]["rows"] else None
        if top:
            bullets.append("熵分数最高的簇：%s（均值 %.3f）" % (top["cluster"], top["mean_score"]))
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_cells": int(len(score)),
            "n_genes": int(len(genes)),
            "n_neighbors": int(n_neighbors),
            "mean_expression_entropy": _num(np.mean(expression), 4),
            "mean_score": _num(np.mean(score), 4),
        },
    }
    return {
        "summary": summary,
        "tables": tables,
        "figures": [figure_entry],
        "obs_columns": [
            {"key": "entropy_stemness_score", "kind": "numeric"},
            {"key": "expression_entropy", "kind": "numeric"},
        ],
        "obs_values": {
            "entropy_stemness_score": _full_numeric(score, idx, n_total),
            "expression_entropy": _full_numeric(expression, idx, n_total),
        },
    }


# ###########################################################################
# 5. cytotrace_lite — 基因数-表达相关性平滑 + 近邻池化
# ###########################################################################


def _column_ranks(dense):
    """按列取平均秩（scipy.stats.rankdata 的轻量包装）。"""
    from scipy.stats import rankdata

    return rankdata(dense, axis=0, method="average")


def _run_cytotrace_lite(view, params, outdir, warnings):
    n_genes = _clamp_int(params, "n_genes", 2000, 100, MAX_GENES)
    n_neighbors = _clamp_int(params, "n_neighbors", 15, 3, 100)
    top_signature = _clamp_int(params, "signature_size", 100, 10, 500)
    matrix = view["X"]
    n_cells = matrix.shape[0]
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel().astype(np.float64)
    points = _pca(view, n_pcs=20)
    neighbors, _distances = _knn(points, n_neighbors)
    smoothed_detected = 0.5 * detected + 0.5 * detected[neighbors].mean(axis=1)

    genes = _hvg_indices(matrix, n_genes)
    ranks_smoothed = _column_ranks(smoothed_detected.reshape(-1, 1)).ravel()
    correlations = np.zeros(len(genes), dtype=np.float64)
    chunk = 500
    for start in range(0, len(genes), chunk):
        index = genes[start : start + chunk]
        block = matrix[:, index].toarray().astype(np.float64)
        block_ranks = _column_ranks(block)
        centered = block_ranks - block_ranks.mean(axis=0, keepdims=True)
        smoothed_centered = ranks_smoothed - ranks_smoothed.mean()
        denominator = np.sqrt(
            np.sum(centered ** 2, axis=0) * np.sum(smoothed_centered ** 2)
        )
        correlations[start : start + chunk] = np.where(
            denominator > 1e-12,
            centered.T @ smoothed_centered / np.maximum(denominator, 1e-12),
            0.0,
        )
    positive_positions = np.where(correlations > 0.02)[0]
    positive_positions = positive_positions[
        np.argsort(-correlations[positive_positions], kind="stable")
    ]
    signature_positions = positive_positions[: int(min(top_signature, positive_positions.size))]
    if signature_positions.size < 5:
        raise EngineError(
            "prerequisite",
            "细胞数量或表达结构不足：找不到与检出基因数正相关的基因，无法计算 CytoTRACE-like 分数",
        )
    signature_genes = genes[signature_positions]
    weights = np.maximum(correlations[signature_positions], 0.0)
    weights = weights / max(weights.sum(), 1e-12)

    signature_matrix = matrix[:, signature_genes].toarray().astype(np.float64)
    mean = signature_matrix.mean(axis=0)
    std = signature_matrix.std(axis=0)
    std[std < 1e-9] = 1.0
    standardized = (signature_matrix - mean) / std
    raw_score = standardized @ weights
    pooled_score = 0.5 * raw_score + 0.5 * raw_score[neighbors].mean(axis=1)
    score = _rank01(pooled_score)

    from scipy.stats import spearmanr

    gene_count_correlation = float(spearmanr(score, detected).statistic) if n_cells > 2 else 0.0
    if not math.isfinite(gene_count_correlation):
        gene_count_correlation = 0.0

    clusters = _clusters(view, required=False)
    tables = []
    if clusters is not None:
        categories = _sorted_categories(clusters)
        rows = []
        for category in categories:
            mask = clusters == category
            rows.append(
                {
                    "cluster": category,
                    "n_cells": int(mask.sum()),
                    "mean_differentiation_potential": _num(np.mean(score[mask]), 4),
                    "median_differentiation_potential": _num(np.median(score[mask]), 4),
                    "mean_detected_genes": _num(np.mean(detected[mask]), 1),
                }
            )
        rows.sort(key=lambda row: -(row["mean_differentiation_potential"] or 0.0))
        tables.append(
            _table(
                "cytotrace_by_cluster",
                "CytoTRACE-like differentiation potential by cluster",
                [
                    _column("cluster", "Cluster"),
                    _column("n_cells", "Cells", "number"),
                    _column("mean_differentiation_potential", "Mean potential", "number"),
                    _column("median_differentiation_potential", "Median potential", "number"),
                    _column("mean_detected_genes", "Mean detected genes", "number"),
                ],
                rows,
                note="分数越高 = 分化潜能越高（越干/越未分化）；为 CytoTRACE 思路的简化实现。",
            )
        )
    else:
        tables.append(
            _table(
                "cytotrace_summary",
                "CytoTRACE-like summary",
                [_column("metric", "Metric"), _column("value", "Value", "number")],
                [
                    {"metric": "n_cells", "value": int(n_cells)},
                    {"metric": "mean_potential", "value": _num(np.mean(score), 4)},
                    {"metric": "signature_genes", "value": int(len(signature_genes))},
                    {"metric": "spearman_score_vs_detected_genes", "value": _num(gene_count_correlation, 4)},
                ],
            )
        )

    embedding, embedding_label = _embedding(view)
    figure, axis = _new_figure(5.9, 4.9)
    scatter = axis.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=score,
        cmap="coolwarm",
        s=7,
        alpha=0.85,
        linewidths=0,
        rasterized=True,
    )
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.03)
    colorbar.set_label("CytoTRACE-like potential")
    axis.set_xlabel("%s 1" % embedding_label)
    axis.set_ylabel("%s 2" % embedding_label)
    axis.set_title("CytoTRACE-like differentiation potential", fontsize=9)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure,
        view["ctx"],
        outdir,
        "cytotrace_lite",
        "differentiation_potential",
        "CytoTRACE-like differentiation potential",
    )

    idx = view["idx"]
    n_total = view["n_total"]
    top_signature_names = [view["var_names"][int(gene)] for gene in signature_genes[:8]]
    bullets = [
        "思路：检出基因数 → kNN 平滑 → 与每个基因算 Spearman 相关 → 取正相关 top %d 基因做加权表达签名（CytoTRACE 近似）"
        % len(signature_genes),
        "签名基因示例：%s" % "、".join(top_signature_names),
        "分数 = 0.5*签名表达 + 0.5*邻居平滑，再秩变换到 [0,1]；高=分化潜能高（干性样）",
        "与检出基因数的 Spearman 相关 = %.3f（CytoTRACE 的核心假设是分化程度与检出基因数正相关）"
        % gene_count_correlation,
        "分数已写入 obs 列 cytotrace_lite_score；检出基因数写入 cytotrace_n_genes",
    ]
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_cells": int(n_cells),
            "n_genes_used": int(len(genes)),
            "signature_genes": int(len(signature_genes)),
            "n_neighbors": int(n_neighbors),
            "spearman_score_vs_detected_genes": _num(gene_count_correlation, 4),
            "mean_potential": _num(np.mean(score), 4),
        },
    }
    return {
        "summary": summary,
        "tables": tables,
        "figures": [figure_entry],
        "obs_columns": [
            {"key": "cytotrace_lite_score", "kind": "numeric"},
            {"key": "cytotrace_n_genes", "kind": "numeric"},
        ],
        "obs_values": {
            "cytotrace_lite_score": _full_numeric(score, idx, n_total),
            "cytotrace_n_genes": _full_numeric(detected, idx, n_total),
        },
    }

# ###########################################################################
# 6. cnv_score / 7. cnv_tree / 8. ith_score — inferCNV 代理实现
# ###########################################################################

_IMMUNE_REFERENCE_WORDS = (
    "t cell", "cd4", "cd8", "nk", "b cell", "plasma", "monocyte", "macrophage",
    "dendritic", "mast", "neutrophil", "platelet", "erythrocyte", "endothelial",
    "fibroblast", "low quality", "lymphocyte", "immune", "normal",
)


def _auto_cnv_reference(obs, warnings):
    """自动挑选免疫/非肿瘤样参照细胞（cell_type/annotation 等列的关键词启发式）。

    找不到可信参照时退回全部细胞（均值参照），并在 warnings 里写明。
    """
    columns = [str(column) for column in obs.columns]
    for column in ("cell_type", "celltype", "cell_type_major", "annotation", "celltype_major"):
        if column not in columns:
            continue
        values = obs[column].astype(str).to_numpy(dtype=object)
        immune = np.asarray(
            [
                any(word in str(value).lower() for word in _IMMUNE_REFERENCE_WORDS)
                for value in values
            ],
            dtype=bool,
        )
        fraction = float(immune.mean()) if len(immune) else 0.0
        if 0.05 <= fraction <= 0.95:
            warnings.append(
                "未显式指定参考细胞：自动把 obs.%s 中的免疫/非肿瘤样细胞（%d/%d）作为 CNV 参照；"
                "如需指定请用 params.reference='列名:取值'" % (column, int(immune.sum()), len(immune))
            )
            return immune, "auto: obs.%s immune-like cells (%d cells)" % (column, int(immune.sum()))
    warnings.append("未显式指定参考细胞：使用全部细胞作为 CNV 参照（均值/标准差），结果只反映相对偏离")
    return np.ones(len(obs), dtype=bool), "all cells"


def _reference_mask(view, params, warnings):
    text = str(params.get("reference") or "").strip()
    obs = view["obs"]
    if text:
        if ":" not in text:
            raise EngineError(
                "bad_request",
                "reference 参数格式应为 'obs列名:取值'（如 'cell_type:Monocyte' 或 'clusters:0'）",
            )
        column, value = text.split(":", 1)
        column = column.strip()
        value = value.strip()
        if column not in obs.columns:
            available = [str(item) for item in obs.columns if not str(item).startswith("_")]
            raise EngineError(
                "bad_request",
                "reference 指定的列 %r 不在 obs 里（可用列：%s）" % (column, "、".join(available[:10])),
            )
        mask = obs[column].astype(str).str.strip().to_numpy(dtype=object) == value
        if int(mask.sum()) < 10:
            raise EngineError(
                "bad_request",
                "reference='%s' 只匹配到 %d 个细胞（至少需要 10 个）" % (text, int(mask.sum())),
            )
        return np.asarray(mask, dtype=bool), "params.reference (%d cells)" % int(mask.sum())
    return _auto_cnv_reference(obs, warnings)


def _parse_gene_order_file(path):
    """解析基因位置文件：两列 gene, chr:start 或三列 gene, chr, start[ end]。

    返回 {gene: (chrom_rank, start)}。
    """
    target = Path(str(path)).expanduser()
    if not target.exists():
        raise EngineError("bad_request", "gene_order 文件不存在：%s" % target)
    if target.is_dir():
        raise EngineError("bad_request", "gene_order 需要是文件，收到目录：%s" % target)
    try:
        text = target.read_text(encoding="utf-8", errors="replace")
    except Exception as error:  # pragma: no cover - 权限/编码兜底
        raise EngineError("bad_request", "无法读取 gene_order 文件 %s：%s" % (target, error))
    chrom_rank = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    order = {}
    for line_number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [part for part in re.split(r"[,\t ]+", line) if part]
        if len(parts) < 2:
            continue
        gene = parts[0].strip()
        location = parts[1].strip()
        if ":" in location:
            chrom_text, start_text = location.split(":", 1)
        elif len(parts) >= 3:
            chrom_text, start_text = location, parts[2]
        else:
            chrom_text, start_text = location, "0"
        chrom_text = re.sub(r"^chr", "", chrom_text.strip(), flags=re.IGNORECASE).upper()
        if chrom_text in chrom_rank:
            rank = chrom_rank[chrom_text]
        else:
            try:
                rank = int(float(chrom_text))
            except (TypeError, ValueError):
                rank = 999
        try:
            start = float(start_text)
        except (TypeError, ValueError):
            start = float("inf")
        order[gene] = (rank, start)
    if not order:
        raise EngineError("bad_request", "gene_order 文件没有解析到任何基因位置：%s" % target)
    return order


def _cnv_gene_selection(view, params, warnings):
    """选 CNV 分析用的基因：有位置文件→按染色体排序；否则→按方差排序的代理顺序。"""
    matrix = view["X"]
    var_names = view["var_names"]
    max_genes = _clamp_int(params, "max_genes", MAX_GENES, 100, MAX_GENES)
    path = str(params.get("gene_order") or "").strip()
    if path:
        order = _parse_gene_order_file(path)
        matched = [(index, gene) for index, gene in enumerate(var_names) if gene in order]
        if len(matched) < 10:
            raise EngineError(
                "bad_request",
                "gene_order 文件与表达矩阵基因名交集过少（%d 个）：请确认文件里的基因符号与 "
                "var_names 一致（人源 HGNC symbol）" % len(matched),
            )
        matched.sort(key=lambda item: (order[item[1]][0], order[item[1]][1], item[1]))
        if len(matched) > max_genes:
            step = len(matched) / float(max_genes)
            matched = [matched[int(position * step)] for position in range(max_genes)]
            warnings.append(
                "gene_order 命中的基因 %d 个超过上限 %d，已按染色体位置等间隔抽取 %d 个"
                % (len(order), max_genes, len(matched))
            )
        gene_indices = np.asarray([item[0] for item in matched], dtype=np.int64)
        gene_labels = [item[1] for item in matched]
        ranks = [order[item[1]][0] for item in matched]
        rank_names = {
            23: "chrX", 24: "chrY", 25: "chrM",
        }
        boundaries = []
        segment_start = 0
        for position in range(1, len(ranks) + 1):
            if position == len(ranks) or ranks[position] != ranks[segment_start]:
                label = rank_names.get(ranks[segment_start], "chr%s" % ranks[segment_start])
                boundaries.append((label, segment_start, position))
                segment_start = position
        return {
            "gene_indices": gene_indices,
            "gene_labels": gene_labels,
            "ordered": True,
            "boundaries": boundaries,
            "gene_order_path": str(path),
        }
    # 无位置文件：按方差排序的代理分数（不能当真实 CNV 区段解读）
    total = int(matrix.shape[1])
    mean = np.asarray(matrix.mean(axis=0)).ravel().astype(np.float64)
    mean_square = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel().astype(np.float64)
    variance = np.maximum(mean_square - mean * mean, 0.0)
    order = np.argsort(-variance, kind="stable")[: int(min(max_genes, total))]
    warnings.append(
        "未提供 gene_order 基因位置文件：CNV 分数是**按方差排序的代理分数**，"
        "基因顺序不代表染色体坐标，不能解读为具体染色体区段/CNV 事件；"
        "若要做 inferCNV 式区段分析，请提供 gene/chr:start 文件（params.gene_order）"
    )
    return {
        "gene_indices": np.asarray(order, dtype=np.int64),
        "gene_labels": [var_names[int(index)] for index in order],
        "ordered": False,
        "boundaries": [],
        "gene_order_path": "",
    }


def _cnv_core(view, params, warnings, sample_rows=None):
    """inferCNV-like CNV 代理分数核心计算。

    - 参照细胞：params.reference（'列:值'）> 自动免疫样关键词 > 全部细胞；
    - 逐基因 z 分数 = (表达 - 参照均值) / 参照标准差；
    - 沿基因顺序做移动窗口平滑；每细胞 CNV 分数 = 平滑后 |z| 的均值；
    - 返回簇均值谱（热图）、每细胞分数、可选样本细胞的平滑谱（ITH 用）。
    """
    matrix = view["X"]
    n_cells = matrix.shape[0]
    selection = _cnv_gene_selection(view, params, warnings)
    gene_indices = selection["gene_indices"]
    gene_labels = selection["gene_labels"]
    n_genes = len(gene_indices)
    reference_mask, reference_note = _reference_mask(view, params, warnings)
    reference_count = int(reference_mask.sum())
    if reference_count < 10:
        raise EngineError("prerequisite", "CNV 参照细胞不足（%d 个，至少需要 10 个）" % reference_count)
    reference_matrix = matrix[reference_mask][:, gene_indices]
    reference_mean = np.asarray(reference_matrix.mean(axis=0)).ravel().astype(np.float64)
    reference_square = np.asarray(reference_matrix.multiply(reference_matrix).mean(axis=0)).ravel()
    reference_std = np.sqrt(np.maximum(reference_square - reference_mean ** 2, 0.0))
    reference_std[reference_std < 1e-6] = 1.0

    window = _clamp_int(params, "window", 11, 3, 201)
    if window % 2 == 0:
        window += 1
    try:
        from scipy.ndimage import uniform_filter1d
    except ImportError:  # pragma: no cover - scipy 一定存在，防御
        uniform_filter1d = None

    clusters = _clusters(view, required=False)
    cluster_categories = _sorted_categories(clusters) if clusters is not None else None
    if cluster_categories is not None:
        position = {category: index for index, category in enumerate(cluster_categories)}
        cluster_codes = np.asarray([position[str(label)] for label in clusters], dtype=np.int64)
        cluster_sums = np.zeros((len(cluster_categories), n_genes), dtype=np.float64)
        cluster_counts = np.bincount(cluster_codes, minlength=len(cluster_categories)).astype(np.float64)
    else:
        cluster_codes = None
        cluster_sums = None
        cluster_counts = None

    scores = np.zeros(n_cells, dtype=np.float64)
    sample_positions = None
    sample_profiles = None
    if sample_rows is not None and len(sample_rows):
        sample_rows = np.asarray(sample_rows, dtype=np.int64)
        sample_positions = {int(row): position for position, row in enumerate(sample_rows)}
        sample_profiles = np.zeros((len(sample_rows), n_genes), dtype=np.float32)

    chunk = max(1, min(2000, int(4_000_000 / max(n_genes, 1))))
    for start in range(0, n_cells, chunk):
        end = min(n_cells, start + chunk)
        block = matrix[start:end][:, gene_indices].toarray().astype(np.float64)
        z_scores = (block - reference_mean[None, :]) / reference_std[None, :]
        if uniform_filter1d is not None:
            smoothed = uniform_filter1d(z_scores, size=window, axis=1, mode="nearest")
        else:  # pragma: no cover - 防御性回退
            smoothed = z_scores
        # 逐细胞中位数居中：去掉细胞整体表达水平/复杂度造成的全基因组偏移，
        # 让热图更接近"相对参照的区段性偏离"（inferCNV 也做逐细胞基线处理）
        smoothed = smoothed - np.median(smoothed, axis=1, keepdims=True)
        scores[start:end] = np.mean(np.abs(smoothed), axis=1)
        if cluster_sums is not None:
            np.add.at(cluster_sums, cluster_codes[start:end], smoothed)
        if sample_profiles is not None:
            for row in range(start, end):
                target = sample_positions.get(int(row))
                if target is not None:
                    sample_profiles[target] = smoothed[row - start].astype(np.float32)

    # 簇均值谱可能因 np.add.at 与浮点误差为 0，做一次防御性归一
    cluster_profiles = None
    if cluster_sums is not None:
        cluster_profiles = cluster_sums / np.maximum(cluster_counts[:, None], 1.0)
    return {
        **selection,
        "cell_scores": scores,
        "cluster_profiles": cluster_profiles,
        "cluster_categories": cluster_categories,
        "cluster_counts": cluster_counts,
        "sample_profiles": sample_profiles,
        "sample_rows": sample_rows,
        "window": int(window),
        "reference_count": int(reference_count),
        "reference_note": reference_note,
        "n_genes_used": int(n_genes),
    }


def _cnv_heatmap_ticks(selection, n_genes):
    boundaries = selection.get("boundaries") or []
    if not boundaries:
        return np.array([]), []
    positions = []
    labels = []
    for label, start, end in boundaries:
        positions.append((start + end - 1) / 2.0)
        labels.append(str(label))
    # 染色体太多时只标前 12 个，避免重叠
    if len(positions) > 12:
        step = int(math.ceil(len(positions) / 12.0))
        positions = positions[::step]
        labels = labels[::step]
    return np.asarray(positions), labels


def _run_cnv_score(view, params, outdir, warnings):
    clusters = _clusters(view, required=False)
    sample_rows = None
    if clusters is None:
        # 无簇时：抽样部分细胞展示单细胞 CNV 谱
        rng = _rng(DEFAULT_SEED)
        n_cells = view["X"].shape[0]
        count = int(min(200, n_cells))
        sample_rows = (
            np.sort(rng.choice(n_cells, size=count, replace=False))
            if n_cells > count
            else np.arange(n_cells)
        )
    core = _cnv_core(view, params, warnings, sample_rows=sample_rows)
    scores = core["cell_scores"]
    gene_labels = core["gene_labels"]
    n_genes = core["n_genes_used"]
    selection = {"boundaries": core["boundaries"]}
    positions, tick_labels = _cnv_heatmap_ticks(selection, n_genes)

    tables = []
    if core["cluster_profiles"] is not None:
        rows = []
        for index, category in enumerate(core["cluster_categories"]):
            mask = clusters == category
            rows.append(
                {
                    "cluster": category,
                    "n_cells": int(core["cluster_counts"][index]),
                    "mean_cnv_score": _num(np.mean(scores[mask]), 4),
                    "median_cnv_score": _num(np.median(scores[mask]), 4),
                    "sd_cnv_score": _num(np.std(scores[mask]), 4),
                }
            )
        rows.sort(key=lambda row: -(row["mean_cnv_score"] or 0.0))
        tables.append(
            _table(
                "cnv_by_cluster",
                "CNV proxy score by cluster",
                [
                    _column("cluster", "Cluster"),
                    _column("n_cells", "Cells", "number"),
                    _column("mean_cnv_score", "Mean CNV score", "number"),
                    _column("median_cnv_score", "Median CNV score", "number"),
                    _column("sd_cnv_score", "SD", "number"),
                ],
                rows,
                note="CNV 分数 = 沿基因顺序移动窗口平滑后 |z| 的均值；参照细胞见 summary。"
                + ("" if core["ordered"] else " 无基因位置文件时为按方差排序的代理分数。"),
            )
        )
    else:
        tables.append(
            _table(
                "cnv_summary",
                "CNV proxy summary",
                [_column("metric", "Metric"), _column("value", "Value", "number")],
                [
                    {"metric": "n_cells", "value": int(len(scores))},
                    {"metric": "mean_cnv_score", "value": _num(np.mean(scores), 4)},
                    {"metric": "n_genes_used", "value": int(n_genes)},
                    {"metric": "window", "value": int(core["window"])},
                ],
            )
        )

    figure, axis = _new_figure(8.2, 4.6)
    if core["cluster_profiles"] is not None:
        heat = core["cluster_profiles"]
        image = axis.imshow(
            heat, aspect="auto", cmap="RdBu_r", interpolation="nearest",
            vmin=-float(np.max(np.abs(heat)) + 1e-9), vmax=float(np.max(np.abs(heat)) + 1e-9),
        )
        axis.set_yticks(np.arange(len(core["cluster_categories"])))
        axis.set_yticklabels([str(item) for item in core["cluster_categories"]], fontsize=7.5)
        axis.set_ylabel("Cluster")
        title_prefix = "Cluster CNV profiles"
    else:
        # 无簇：按 CNV 分数排序展示抽样细胞的平滑谱
        heat = np.asarray(core["sample_profiles"], dtype=np.float64)
        heat = heat - heat.mean(axis=1, keepdims=True)
        heat = heat / (heat.std(axis=1, keepdims=True) + 1e-9)
        image = axis.imshow(
            heat, aspect="auto", cmap="RdBu_r", interpolation="nearest",
            vmin=-float(np.max(np.abs(heat)) + 1e-9), vmax=float(np.max(np.abs(heat)) + 1e-9),
        )
        axis.set_ylabel("Cells (sampled)")
        axis.set_yticks([])
        title_prefix = "Single-cell CNV profiles (sampled)"
    if len(positions):
        axis.set_xticks(positions)
        axis.set_xticklabels(tick_labels, fontsize=7)
        axis.set_xlabel("Genomic position (chromosome order from gene_order file)")
    else:
        axis.set_xticks([])
        axis.set_xlabel(
            "Genes ordered by variance (PROXY: no chromosome order; not a real CNV profile)"
            if not core["ordered"]
            else "Genes in provided order"
        )
    colorbar = figure.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
    colorbar.set_label("Smoothed z-score")
    axis.set_title("%s, window=%d genes" % (title_prefix, core["window"]), fontsize=9)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "cnv_score", "cnv_heatmap", "CNV proxy heatmap"
    )

    idx = view["idx"]
    n_total = view["n_total"]
    bullets = [
        "inferCNV-like 代理实现：%d 个基因、移动窗口 %d 平滑，参照细胞 %d 个（%s）"
        % (n_genes, core["window"], core["reference_count"], core["reference_note"]),
        "每细胞 CNV 分数 = 平滑 z 分数绝对值的均值；分数高低只表示相对参照的偏离程度",
    ]
    if core["ordered"]:
        bullets.append(
            "已按 gene_order 文件的染色体位置排序（%d 条染色体区段），热图 x 轴为真实基因组顺序"
            % len(core["boundaries"])
        )
    else:
        bullets.append(
            "**没有 gene_order 文件：这是按方差排序的代理分数（proxy），基因顺序不对应染色体坐标，不能当作真实 CNV 区段解读**"
        )
    bullets.append("CNV 分数已写入 obs 列 cnv_score")
    method_note = "inferCNV-like（本插件代理实现）"
    summary = {
        "bullets": bullets,
        "metrics": {
            "method_scope": method_note,
            "ordered_by_chromosome": bool(core["ordered"]),
            "n_chromosome_segments": int(len(core["boundaries"])) if core["ordered"] else 0,
            "n_genes_used": int(n_genes),
            "window": int(core["window"]),
            "reference_cells": int(core["reference_count"]),
            "reference_note": core["reference_note"],
            "mean_cnv_score": _num(np.mean(scores), 4),
            "max_cnv_score": _num(np.max(scores), 4),
        },
    }
    return {
        "summary": summary,
        "tables": tables,
        "figures": [figure_entry],
        "obs_columns": [{"key": "cnv_score", "kind": "numeric"}],
        "obs_values": {"cnv_score": _full_numeric(scores, idx, n_total)},
    }


def _linkage_newick(linkage_matrix, labels):
    """把 scipy linkage 矩阵转成 Newick 字符串（分支长度 = 高度差）。"""
    count = len(labels)

    def safe(text):
        return re.sub(r"[\s()\[\]{},;:'\"]+", "_", str(text))

    nodes = {
        index: {"text": safe(labels[index]), "height": 0.0} for index in range(count)
    }
    for row_index, row in enumerate(np.asarray(linkage_matrix)):
        left_id = int(row[0])
        right_id = int(row[1])
        height = float(row[2])
        left = nodes.pop(left_id)
        right = nodes.pop(right_id)
        text = "(%s:%.6f,%s:%.6f)" % (
            left["text"],
            max(height - left["height"], 0.0),
            right["text"],
            max(height - right["height"], 0.0),
        )
        nodes[count + row_index] = {"text": text, "height": height}
    root = nodes[max(nodes.keys())]
    return root["text"] + ";"


def _run_cnv_tree(view, params, outdir, warnings):
    clusters = _clusters(view, required=True, what=" CNV 聚类树")
    categories = _sorted_categories(clusters)
    if len(categories) < 2:
        raise EngineError("bad_request", "只有 1 个簇，无法构建 CNV 距离树（至少需要 2 个簇）")
    core = _cnv_core(view, params, warnings)
    if core["cluster_profiles"] is None:
        raise EngineError("prerequisite", "缺少簇标签，无法构建 CNV 聚类树")
    profiles = np.asarray(core["cluster_profiles"], dtype=np.float64)
    # 去掉全零行（参照均值谱），防止距离矩阵退化
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import pdist
    except ImportError as error:  # pragma: no cover
        raise EngineError("unsupported", "缺少 scipy，无法做层次聚类：%s" % error)
    distances = pdist(profiles, metric="euclidean")
    linkage_matrix = linkage(distances, method="average", optimal_ordering=True)
    newick = _linkage_newick(linkage_matrix, categories)

    merge_rows = []
    for step, row in enumerate(np.asarray(linkage_matrix)):
        left_id, right_id, height, size = (int(row[0]), int(row[1]), float(row[2]), int(row[3]))

        def name_of(node_id):
            if node_id < len(categories):
                return str(categories[node_id])
            return "node%d" % node_id

        merge_rows.append(
            {
                "step": step + 1,
                "left": name_of(left_id),
                "right": name_of(right_id),
                "distance": _num(height, 5),
                "cluster_size": size,
                "new_node": "node%d" % (len(categories) + step),
            }
        )
    merge_table = _table(
        "cnv_tree_merges",
        "Hierarchical clustering merges (average linkage, Euclidean on CNV profiles)",
        [
            _column("step", "Step", "number"),
            _column("left", "Left"),
            _column("right", "Right"),
            _column("distance", "Distance", "number"),
            _column("cluster_size", "Merged size", "number"),
            _column("new_node", "New node"),
        ],
        merge_rows,
        note="叶节点 = 簇名；newick 字符串见 summary.metrics.newick。",
    )

    figure, axis = _new_figure(7.2, 3.9)
    dendrogram(
        linkage_matrix,
        labels=[str(item) for item in categories],
        ax=axis,
        leaf_font_size=8,
        color_threshold=0.0,
        above_threshold_color="#333333",
    )
    axis.set_ylabel("Euclidean distance (CNV profile)")
    axis.set_title("CNV-profile cluster tree (inferCNV-like proxy)", fontsize=9)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "cnv_tree", "cnv_dendrogram", "CNV cluster tree"
    )
    bullets = [
        "对 %d 个簇的 CNV 平滑谱（%d 个基因）做欧氏距离 + average-linkage 层次聚类（scipy）"
        % (len(categories), core["n_genes_used"]),
        "树是无根的表型聚类，用于展示簇间 CNV 相似性，不是真正的克隆演化树（无突变/等位基因信息）",
    ]
    if not core["ordered"]:
        bullets.append("注意：本次 CNV 无基因位置文件，谱为按方差排序的代理顺序")
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_clusters": int(len(categories)),
            "n_genes_used": int(core["n_genes_used"]),
            "linkage": "average",
            "metric": "euclidean",
            "newick": newick,
            "max_merge_distance": _num(float(np.max(linkage_matrix[:, 2])), 5),
        },
    }
    return {
        "summary": summary,
        "tables": [merge_table],
        "figures": [figure_entry],
        "obs_columns": [],
    }


def _run_ith_score(view, params, outdir, warnings):
    clusters = _clusters(view, required=True, what=" ITH（瘤内异质性）打分")
    categories = _sorted_categories(clusters)
    if len(categories) < 2:
        raise EngineError("bad_request", "只有 1 个簇，无法计算瘤内异质性（ITH）分数")
    per_cluster = _clamp_int(params, "max_cells_per_cluster", 120, 10, 500)
    rng = _rng(DEFAULT_SEED)
    sample_rows = []
    cluster_sample_indices = {}
    for category in categories:
        indices = np.where(clusters == category)[0]
        if indices.size > per_cluster:
            chosen = np.sort(rng.choice(indices, size=per_cluster, replace=False))
        else:
            chosen = indices
        cluster_sample_indices[category] = chosen
        sample_rows.extend(int(value) for value in chosen)
    sample_rows = np.asarray(sorted(set(sample_rows)), dtype=np.int64)
    core = _cnv_core(view, params, warnings, sample_rows=sample_rows)
    profiles = core["sample_profiles"]
    if profiles is None or profiles.shape[0] == 0:
        raise EngineError("internal", "ITH 计算失败：没有可用的 CNV 平滑谱")
    row_lookup = {int(row): position for position, row in enumerate(core["sample_rows"])}
    n_genes = core["n_genes_used"]

    rows = []
    cell_scores = np.zeros(view["X"].shape[0], dtype=np.float64)
    for category in categories:
        indices = cluster_sample_indices[category]
        positions = [row_lookup[int(index)] for index in indices if int(index) in row_lookup]
        matrix = profiles[positions].astype(np.float64)
        if len(positions) > 1:
            from scipy.spatial.distance import pdist

            pairwise = pdist(matrix, metric="euclidean")
            mean_pairwise = float(np.mean(pairwise))
            median_pairwise = float(np.median(pairwise))
        else:
            mean_pairwise = 0.0
            median_pairwise = 0.0
        normalized = mean_pairwise / math.sqrt(max(n_genes, 1))
        cell_scores[clusters == category] = normalized
        rows.append(
            {
                "cluster": category,
                "n_cells": int(np.sum(clusters == category)),
                "n_sampled": int(len(positions)),
                "mean_pairwise_cnv_distance": _num(mean_pairwise, 5),
                "median_pairwise_cnv_distance": _num(median_pairwise, 5),
                "ith_score": _num(normalized, 5),
            }
        )
    rows.sort(key=lambda row: -(row["ith_score"] or 0.0))
    table = _table(
        "ith_by_cluster",
        "Intra-tumor heterogeneity by cluster (mean pairwise CNV-profile distance)",
        [
            _column("cluster", "Cluster"),
            _column("n_cells", "Cells", "number"),
            _column("n_sampled", "Sampled", "number"),
            _column("mean_pairwise_cnv_distance", "Mean pairwise distance", "number"),
            _column("median_pairwise_cnv_distance", "Median pairwise distance", "number"),
            _column("ith_score", "ITH score (distance / sqrt(genes))", "number"),
        ],
        rows,
        note="簇内平均成对欧氏距离（基于 CNV 平滑谱，每簇最多抽样 %d 个细胞）。" % per_cluster,
    )

    figure, axis = _new_figure(6.2, 3.8)
    ordered = list(reversed(rows))
    y_positions = np.arange(len(ordered))
    values = [row["ith_score"] or 0.0 for row in ordered]
    colors = _palette(len(ordered))
    axis.barh(y_positions, values, color=colors, height=0.65, edgecolor="black", linewidth=0.4)
    axis.set_yticks(y_positions)
    axis.set_yticklabels([str(row["cluster"]) for row in ordered], fontsize=8)
    axis.set_xlabel("ITH score (mean pairwise CNV distance / sqrt(n_genes))")
    axis.set_ylabel("Cluster")
    axis.set_title("Intra-tumor heterogeneity by cluster", fontsize=9)
    top = rows[0] if rows else None
    if top:
        axis.text(
            0.98,
            0.04,
            "top cluster: %s (%.3f)" % (top["cluster"], top["ith_score"] or 0.0),
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
        )
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "ith_score", "ith_bar", "ITH score by cluster"
    )

    idx = view["idx"]
    n_total = view["n_total"]
    bullets = [
        "ITH 分数 = 簇内平均成对 CNV 平滑谱欧氏距离 / sqrt(基因数)，每簇最多抽样 %d 个细胞" % per_cluster,
        "分数越高表示该簇内部 CNV 组成越异质（亚克隆混合程度越高）；这是表型代理，不是克隆推断",
    ]
    if not core["ordered"]:
        bullets.append("注意：本次 CNV 无基因位置文件，谱为按方差排序的代理顺序")
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_clusters": int(len(categories)),
            "n_genes_used": int(n_genes),
            "max_cells_per_cluster": int(per_cluster),
            "top_cluster": top["cluster"] if top else None,
            "top_ith_score": top["ith_score"] if top else None,
        },
    }
    return {
        "summary": summary,
        "tables": [table],
        "figures": [figure_entry],
        "obs_columns": [{"key": "ith_score", "kind": "numeric"}],
        "obs_values": {"ith_score": _full_numeric(cell_scores, idx, n_total)},
    }

# ###########################################################################
# 9. lr_communication — CellChat/NATMI-lite 配体-受体通讯
# ###########################################################################


def _bh_qvalues(pvalues):
    """Benjamini–Hochberg FDR（输入 list[float]，返回同序 list[float]）。"""
    values = np.asarray(pvalues, dtype=np.float64)
    count = values.size
    if count == 0:
        return []
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    qvalues = np.empty(count, dtype=np.float64)
    previous = 1.0
    for position in range(count - 1, -1, -1):
        rank = position + 1
        adjusted = min(previous, ranked[position] * count / rank)
        qvalues[position] = adjusted
        previous = adjusted
    output = np.empty(count, dtype=np.float64)
    output[order] = np.clip(qvalues, 0.0, 1.0)
    return output.tolist()


def _run_lr_communication(view, params, outdir, warnings):
    clusters = _clusters(view, required=True, what="细胞通讯分析（lr_communication）")
    categories = _sorted_categories(clusters)
    k_clusters = len(categories)
    if k_clusters < 2:
        raise EngineError(
            "bad_request",
            "细胞通讯分析至少需要 2 个簇/细胞类型（当前只有 %d 个）" % k_clusters,
        )
    n_perm = _clamp_int(params, "n_perm", 50, 0, 200)
    top_n = _clamp_int(params, "top_n", 40, 5, 200)
    matrix = view["X"]
    n_cells = matrix.shape[0]
    var_names = view["var_names"]
    gene_index = {gene: index for index, gene in enumerate(var_names)}
    pairs = [
        (ligand, receptor, family)
        for ligand, receptor, family in _priors()[0]
        if ligand in gene_index and receptor in gene_index
    ]
    if len(pairs) < 3:
        raise EngineError(
            "prerequisite",
            "表达矩阵中只找到 %d 个内置配体-受体对（%d 对人的 LR 对）：lr_communication 需要人源基因符号；"
            "请确认 var_names 是 HGNC symbol（如 CCL5/CCR5、IL6/IL6R）" % (len(pairs), len(_priors()[0])),
        )

    lr_genes = []
    seen = set()
    for ligand, receptor, _family in pairs:
        for gene in (ligand, receptor):
            if gene not in seen:
                seen.add(gene)
                lr_genes.append(gene)
    lr_columns = np.asarray([gene_index[gene] for gene in lr_genes], dtype=np.int64)
    lr_matrix = matrix[:, lr_columns].tocsr()
    means, counts, codes = _group_means(lr_matrix, clusters, categories)
    gene_position = {gene: index for index, gene in enumerate(lr_genes)}
    small_clusters = [categories[index] for index in range(k_clusters) if counts[index] < 5]
    if small_clusters:
        warnings.append(
            "以下簇细胞数 <5，通讯强度估计不稳定：%s" % "、".join(str(item) for item in small_clusters)
        )

    observed = {}
    for pair_index, (ligand, receptor, _family) in enumerate(pairs):
        ligand_values = means[:, gene_position[ligand]]
        receptor_values = means[:, gene_position[receptor]]
        observed[pair_index] = np.outer(ligand_values, receptor_values)

    # 置换零分布：打乱簇标签后重算簇均值，构造 LR 强度的全球零分布
    permuted_means = []
    if n_perm > 0:
        rng = _rng(DEFAULT_SEED)
        rows = np.arange(n_cells, dtype=np.int64)
        for _ in range(n_perm):
            permuted = rng.permutation(codes)
            onehot = sp.csr_matrix(
                (np.ones(n_cells, dtype=np.float32), (permuted, rows)),
                shape=(k_clusters, n_cells),
            )
            sums = np.asarray((onehot @ lr_matrix).toarray(), dtype=np.float64)
            permuted_means.append(sums / np.maximum(counts[:, None], 1.0))

    records = []
    if n_perm == 0:
        warnings.append("n_perm=0：未做置换检验，p 值全部为 1（建议至少 20 次置换）")
    for pair_index, (ligand, receptor, family) in enumerate(pairs):
        strength = observed[pair_index]
        if n_perm > 0:
            null_stack = np.asarray(
                [
                    np.outer(block[:, gene_position[ligand]], block[:, gene_position[receptor]])
                    for block in permuted_means
                ]
            )
            pvalues = (1.0 + np.sum(null_stack >= strength[None, :, :], axis=0)) / (n_perm + 1.0)
        else:
            pvalues = np.ones_like(strength)
        for source in range(k_clusters):
            for target in range(k_clusters):
                value = float(strength[source, target])
                if value <= 0:
                    continue
                records.append(
                    {
                        "pair_index": pair_index,
                        "ligand": ligand,
                        "receptor": receptor,
                        "family": family,
                        "source": categories[source],
                        "target": categories[target],
                        "strength": value,
                        "p_value": float(pvalues[source, target]),
                        "source_index": source,
                        "target_index": target,
                    }
                )
    if not records:
        raise EngineError("internal", "没有计算到任何非零通讯强度（LR 基因在簇均值上全为 0）")
    qvalues = _bh_qvalues([record["p_value"] for record in records])
    for record, qvalue in zip(records, qvalues):
        record["q_value"] = float(qvalue)
        record["significant"] = bool(qvalue < 0.05 and record["strength"] > 0)
    significant_records = [record for record in records if record["significant"]]
    if not significant_records:
        warnings.append("置换检验后没有 q<0.05 的显著通讯；表格仍按强度排序展示，解释需谨慎")

    records.sort(key=lambda record: (-record["strength"], record["ligand"], record["receptor"]))
    interaction_rows = []
    for record in records[:top_n]:
        interaction_rows.append(
            {
                "source": record["source"],
                "target": record["target"],
                "ligand": record["ligand"],
                "receptor": record["receptor"],
                "family": record["family"],
                "strength": _num(record["strength"], 5),
                "p_value": _num(record["p_value"], 4),
                "q_value": _num(record["q_value"], 4),
                "significant": "yes" if record["significant"] else "no",
            }
        )
    interaction_table = _table(
        "lr_interactions",
        "Top ligand-receptor interactions",
        [
            _column("source", "Source"),
            _column("target", "Target"),
            _column("ligand", "Ligand"),
            _column("receptor", "Receptor"),
            _column("family", "Family"),
            _column("strength", "Strength", "number"),
            _column("p_value", "p (permutation)", "number"),
            _column("q_value", "q (BH)", "number"),
            _column("significant", "Significant"),
        ],
        interaction_rows,
        note="strength = 配体簇均值 × 受体簇均值；p 来自打乱簇标签的置换零分布（全球零），q 为 BH 校正。",
    )

    # 细胞类型对聚合强度
    total_strength = np.zeros((k_clusters, k_clusters), dtype=np.float64)
    significant_count = np.zeros((k_clusters, k_clusters), dtype=np.int64)
    for record in records:
        total_strength[record["source_index"], record["target_index"]] += record["strength"]
        if record["significant"]:
            significant_count[record["source_index"], record["target_index"]] += 1
    aggregate_rows = []
    for source in range(k_clusters):
        for target in range(k_clusters):
            if source == target and significant_count[source, target] == 0 and total_strength[source, target] <= 0:
                continue
            aggregate_rows.append(
                {
                    "source": categories[source],
                    "target": categories[target],
                    "total_strength": _num(total_strength[source, target], 5),
                    "n_significant_lr": int(significant_count[source, target]),
                }
            )
    aggregate_rows.sort(key=lambda row: -(row["total_strength"] or 0.0))
    aggregate_table = _table(
        "cell_cell_strength",
        "Aggregated communication strength between cell groups",
        [
            _column("source", "Source"),
            _column("target", "Target"),
            _column("total_strength", "Total strength", "number"),
            _column("n_significant_lr", "Significant LR pairs", "number"),
        ],
        aggregate_rows,
        note="每个 source→target 方向对所有 LR 对的强度求和。",
    )

    # 网络图：节点 = 细胞类型（环形布局），边 = 聚合通讯强度
    angles = np.linspace(0.0, 2.0 * math.pi, k_clusters, endpoint=False)
    node_positions = np.column_stack([np.cos(angles), np.sin(angles)])
    figure, axis = _new_figure(6.6, 6.4)
    significant_pairs = np.argwhere(significant_count > 0)
    if significant_pairs.size:
        draw_pairs = [(int(left), int(right)) for left, right in significant_pairs]
    else:
        flat = np.dstack(np.unravel_index(np.argsort(-total_strength, axis=None), total_strength.shape))[0]
        draw_pairs = [(int(left), int(right)) for left, right in flat[: min(30, len(flat))] if total_strength[left, right] > 0]
    maximum_pair = max(
        [float(total_strength[left, right]) for left, right in draw_pairs] or [1.0]
    )
    for left, right in draw_pairs:
        value = float(total_strength[left, right])
        width = 0.6 + 4.0 * math.sqrt(max(value, 0.0) / max(maximum_pair, 1e-9))
        if left == right:
            # 自分泌：画一个小环
            position = node_positions[left]
            circle = plt.Circle(
                (position[0] * 1.18, position[1] * 1.18),
                0.06 + 0.05 * width,
                fill=False,
                color="#F58518",
                linewidth=width,
                alpha=0.8,
            )
            axis.add_patch(circle)
            continue
        axis.annotate(
            "",
            xy=(node_positions[right, 0], node_positions[right, 1]),
            xytext=(node_positions[left, 0], node_positions[left, 1]),
            arrowprops=dict(
                arrowstyle="-|>",
                color="#4C78A8",
                lw=width,
                alpha=0.75,
                shrinkA=14,
                shrinkB=14,
                connectionstyle="arc3,rad=0.08",
            ),
            zorder=2,
        )
    node_sizes = 220.0 + 1200.0 * (counts / max(counts.max(), 1.0))
    node_colors = _palette(k_clusters)
    axis.scatter(
        node_positions[:, 0],
        node_positions[:, 1],
        s=node_sizes,
        c=node_colors,
        edgecolors="black",
        linewidths=0.7,
        zorder=3,
    )
    for index, category in enumerate(categories):
        axis.text(
            node_positions[index, 0] * 1.25,
            node_positions[index, 1] * 1.25,
            str(category),
            fontsize=8,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#BBBBBB", lw=0.5, alpha=0.9),
            zorder=4,
        )
    axis.set_xlim(-1.55, 1.55)
    axis.set_ylim(-1.55, 1.55)
    axis.set_aspect("equal")
    axis.axis("off")
    if significant_records:
        title = "Cell-cell communication (CellChat/NATMI-lite)\n%d significant LR interactions, %d groups" % (
            len(significant_records),
            k_clusters,
        )
    else:
        title = (
            "Cell-cell communication (CellChat/NATMI-lite)\n"
            "no significant LR after permutation; edges = top aggregated strengths, %d groups" % k_clusters
        )
    axis.set_title(title, fontsize=9)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "lr_communication", "communication_network", "Communication network"
    )

    # 每细胞的簇级通讯分数（出度+入度，显著 LR 聚合）写回 obs，供嵌入图着色
    outgoing = significant_count.sum(axis=1).astype(np.float64)
    incoming = significant_count.sum(axis=0).astype(np.float64)
    cluster_score = outgoing + incoming
    cell_scores = np.asarray([cluster_score[categories.index(str(label))] for label in clusters])
    idx = view["idx"]
    n_total = view["n_total"]

    bullets = [
        "内置 %d 个配体-受体对（该数据集命中 %d 个），通讯强度 = 配体簇均值 × 受体簇均值（CellChat/NATMI-lite）"
        % (len(_priors()[0]), len(pairs)),
        "置换检验：打乱簇标签 %d 次（全球零分布），BH 校正后 q<0.05 视为显著；本实现不复刻 CellChat 数据库与复合物"
        % n_perm,
        "显著 LR 交互 %d 条；最强：%s→%s %s-%s（strength=%.3f）"
        % (
            len(significant_records),
            records[0]["source"],
            records[0]["target"],
            records[0]["ligand"],
            records[0]["receptor"],
            records[0]["strength"],
        ),
        "簇级通讯分数（出+入度显著 LR 数）已写入 obs 列 lr_communication_score",
    ]
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_lr_pairs_total": int(len(_priors()[0])),
            "n_lr_pairs_present": int(len(pairs)),
            "n_permutations": int(n_perm),
            "n_interactions": int(len(records)),
            "n_significant": int(len(significant_records)),
            "n_groups": int(k_clusters),
        },
    }
    return {
        "summary": summary,
        "tables": [interaction_table, aggregate_table],
        "figures": [figure_entry],
        "obs_columns": [{"key": "lr_communication_score", "kind": "numeric"}],
        "obs_values": {"lr_communication_score": _full_numeric(cell_scores, idx, n_total)},
    }


# ###########################################################################
# 10. regulon_activity — SCENIC-lite 转录调控
# ###########################################################################


def _sanitize_obs_key(prefix, name):
    return prefix + re.sub(r"[^A-Za-z0-9_]+", "_", str(name)).strip("_")


def _run_regulon_activity(view, params, outdir, warnings):
    n_hvg = _clamp_int(params, "n_genes", 1500, 100, MAX_GENES)
    n_targets = _clamp_int(params, "n_targets", 20, 3, 100)
    n_regulons = _clamp_int(params, "n_regulons", 30, 3, 200)
    top_columns = _clamp_int(params, "top_regulons", 8, 1, 30)
    matrix = view["X"]
    n_cells = matrix.shape[0]
    var_names = view["var_names"]
    gene_index = {gene: index for index, gene in enumerate(var_names)}
    tf_all = _priors()[1]
    tfs_present = [tf for tf in tf_all if tf in gene_index]
    if len(tfs_present) < 3:
        raise EngineError(
            "prerequisite",
            "表达矩阵中只找到 %d 个内置人源 TF（内置列表 %d 个）：regulon_activity 需要人源基因符号；"
            "请确认 var_names 是 HGNC symbol（如 STAT1、IRF1、CEBPA）" % (len(tfs_present), len(tf_all)),
        )

    hvg = _hvg_indices(matrix, n_hvg)
    tf_columns = np.asarray([gene_index[tf] for tf in tfs_present], dtype=np.int64)
    extra = np.setdiff1d(tf_columns, hvg)
    gene_columns = np.concatenate([hvg, extra])
    if gene_columns.size > MAX_GENES:
        cut = int(MAX_GENES - extra.size)
        if cut >= 100:
            gene_columns = np.concatenate([hvg[:cut], extra])
        else:
            gene_columns = gene_columns[:MAX_GENES]
        warnings.append("regulon 基因集超过 %d 上限，已截断到 %d 个基因" % (MAX_GENES, gene_columns.size))
    gene_labels = [var_names[int(index)] for index in gene_columns]
    label_position = {label: position for position, label in enumerate(gene_labels)}
    tfs_used = [tf for tf in tfs_present if tf in label_position]
    if len(tfs_used) < 3:
        raise EngineError("internal", "regulon 构建失败：可用 TF 太少（%d）" % len(tfs_used))

    # TF × 基因 共表达相关（genes 含所有 TF，保证 TF 自己也在列里）
    tf_dense = matrix[:, [gene_index[tf] for tf in tfs_used]].toarray().astype(np.float64)
    tf_mean = tf_dense.mean(axis=0)
    tf_std = tf_dense.std(axis=0)
    tf_std[tf_std < 1e-9] = 1.0
    z_tf = (tf_dense - tf_mean[None, :]) / tf_std[None, :]
    correlations = np.zeros((len(tfs_used), len(gene_columns)), dtype=np.float64)
    block_size = 400
    for start in range(0, len(gene_columns), block_size):
        columns = gene_columns[start : start + block_size]
        block = matrix[:, columns].toarray().astype(np.float64)
        mean = block.mean(axis=0)
        std = block.std(axis=0)
        std[std < 1e-9] = 1.0
        z_block = (block - mean[None, :]) / std[None, :]
        correlations[:, start : start + block_size] = (z_tf.T @ z_block) / max(n_cells - 1, 1)

    regulons = []
    for row, tf in enumerate(tfs_used):
        values = correlations[row].copy()
        own = label_position[tf]
        values[own] = -np.inf
        order = np.argsort(-values, kind="stable")
        selected = [position for position in order[:n_targets] if values[position] > 0.02]
        if len(selected) < 3:
            continue
        selected = np.asarray(selected, dtype=np.int64)
        score = float(np.mean(values[selected]))
        regulons.append(
            {
                "tf": tf,
                "target_positions": selected,
                "target_labels": [gene_labels[int(position)] for position in selected],
                "target_correlations": [float(values[int(position)]) for position in selected],
                "score": score,
            }
        )
    if not regulons:
        raise EngineError(
            "prerequisite",
            "没有构建出任何 regulon：TF 与高变基因的共表达相关过弱（可能是数据过度稀疏或基因名不匹配）",
        )
    regulons.sort(key=lambda item: -item["score"])
    regulons = regulons[:n_regulons]

    # AUCell 式活性：每细胞在背景基因集内的秩（1=最高），regulon 基因平均恢复率
    activity = np.zeros((n_cells, len(regulons)), dtype=np.float32)
    target_positions = [item["target_positions"] for item in regulons]
    background = len(gene_columns)
    chunk = 1000
    for start in range(0, n_cells, chunk):
        end = min(n_cells, start + chunk)
        block = matrix[start:end][:, gene_columns].toarray().astype(np.float32)
        order = np.argsort(-block, axis=1)
        ranks = np.empty_like(order)
        row_index = np.arange(block.shape[0])[:, None]
        ranks[row_index, order] = np.arange(background, dtype=np.int64)[None, :]
        for regulon_index, positions in enumerate(target_positions):
            recovery = 1.0 - ranks[:, positions].astype(np.float64) / float(background)
            activity[start:end, regulon_index] = recovery.mean(axis=1)

    variances = activity.var(axis=0)
    top_indices = np.argsort(-variances, kind="stable")[:top_columns]

    clusters = _clusters(view, required=False)
    table_rows = []
    heat_matrix = None
    heat_columns = None
    if clusters is not None:
        categories = _sorted_categories(clusters)
        position = {category: index for index, category in enumerate(categories)}
        codes = np.asarray([position[str(label)] for label in clusters], dtype=np.int64)
        cluster_means = np.zeros((len(regulons), len(categories)), dtype=np.float64)
        for regulon_index in range(len(regulons)):
            for cluster_index in range(len(categories)):
                mask = codes == cluster_index
                cluster_means[regulon_index, cluster_index] = (
                    float(np.mean(activity[mask, regulon_index])) if mask.any() else 0.0
                )
        heat_matrix = cluster_means
        heat_columns = categories
        table_categories = categories[:8]
        for regulon_index, regulon in enumerate(regulons):
            row = {
                "tf": regulon["tf"],
                "family": str(_priors_family(regulon["tf"])),
                "n_targets": int(len(regulon["target_labels"])),
                "mean_activity": _num(float(np.mean(activity[:, regulon_index])), 4),
                "sd_activity": _num(float(np.std(activity[:, regulon_index])), 4),
                "variance": _num(float(variances[regulon_index]), 6),
            }
            for cluster_index, category in enumerate(table_categories):
                row["cluster_%s" % category] = _num(cluster_means[regulon_index, cluster_index], 4)
            table_rows.append(row)
        columns = [
            _column("tf", "TF"),
            _column("family", "Family"),
            _column("n_targets", "Targets", "number"),
            _column("mean_activity", "Mean activity", "number"),
            _column("sd_activity", "SD", "number"),
            _column("variance", "Variance", "number"),
        ] + [
            _column("cluster_%s" % category, "Activity in %s" % category, "number")
            for category in table_categories
        ]
        activity_table = _table(
            "regulon_activity",
            "Regulon activity (AUCell-like) by cluster",
            columns,
            table_rows,
            note="活性 = regulon 基因在背景基因集内秩恢复率的均值（0-1），不是官方 AUCell 的 AUC。"
            " 簇多于 8 个时表格只展示前 8 个簇。",
        )
    else:
        for regulon_index, regulon in enumerate(regulons):
            table_rows.append(
                {
                    "tf": regulon["tf"],
                    "family": str(_priors_family(regulon["tf"])),
                    "n_targets": int(len(regulon["target_labels"])),
                    "mean_activity": _num(float(np.mean(activity[:, regulon_index])), 4),
                    "sd_activity": _num(float(np.std(activity[:, regulon_index])), 4),
                    "variance": _num(float(variances[regulon_index]), 6),
                }
            )
        activity_table = _table(
            "regulon_activity",
            "Regulon activity (AUCell-like)",
            [
                _column("tf", "TF"),
                _column("family", "Family"),
                _column("n_targets", "Targets", "number"),
                _column("mean_activity", "Mean activity", "number"),
                _column("sd_activity", "SD", "number"),
                _column("variance", "Variance", "number"),
            ],
            table_rows,
            note="活性 = regulon 基因在背景基因集内秩恢复率的均值（0-1）。",
        )

    targets_rows = []
    for regulon in regulons:
        for target, corr in list(zip(regulon["target_labels"], regulon["target_correlations"]))[:3]:
            targets_rows.append(
                {
                    "tf": regulon["tf"],
                    "target": str(target),
                    "correlation": _num(corr, 4),
                }
            )
    targets_table = _table(
        "regulon_targets",
        "Top co-expression targets per regulon",
        [
            _column("tf", "TF"),
            _column("target", "Target gene"),
            _column("correlation", "Pearson r (co-expression)", "number"),
        ],
        targets_rows,
        note="靶基因仅由共表达相关选出，未使用 motif/cisTarget 数据库（SCENIC-lite 的关键简化）。",
    )

    # 热图
    figure, axis = _new_figure(7.6, 4.6)
    if heat_matrix is not None:
        image = axis.imshow(heat_matrix, aspect="auto", cmap="viridis", interpolation="nearest")
        axis.set_yticks(np.arange(len(regulons)))
        axis.set_yticklabels([regulon["tf"] for regulon in regulons], fontsize=7.5)
        axis.set_xticks(np.arange(len(heat_columns)))
        axis.set_xticklabels([str(item) for item in heat_columns], rotation=45, ha="right", fontsize=7.5)
        axis.set_xlabel("Cluster")
    else:
        rng = _rng(DEFAULT_SEED)
        order = np.argsort(-activity.mean(axis=1))[: min(300, n_cells)]
        if order.size > 300:
            order = order[:: max(1, order.size // 300)][:300]
        image = axis.imshow(activity[order].T, aspect="auto", cmap="viridis", interpolation="nearest")
        axis.set_yticklabels([])
        axis.set_yticks(np.arange(len(regulons)))
        axis.set_yticklabels([regulon["tf"] for regulon in regulons], fontsize=7.5)
        axis.set_xlabel("Cells ordered by mean regulon activity (sampled)")
    colorbar = figure.colorbar(image, ax=axis, fraction=0.035, pad=0.02)
    colorbar.set_label("Regulon activity (0-1)")
    axis.set_ylabel("Regulon (TF)")
    axis.set_title(
        "SCENIC-lite regulon activity (%d regulons, %d background genes)"
        % (len(regulons), background),
        fontsize=9,
    )
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "regulon_activity", "regulon_heatmap", "Regulon activity heatmap"
    )

    idx = view["idx"]
    n_total = view["n_total"]
    obs_values = {}
    obs_columns = []
    for regulon_index in top_indices:
        regulon = regulons[int(regulon_index)]
        key = _sanitize_obs_key("regulon_", regulon["tf"])
        if key in obs_values:
            continue
        obs_values[key] = _full_numeric(activity[:, regulon_index], idx, n_total)
        obs_columns.append({"key": key, "kind": "numeric"})
    obs_values["regulon_mean_activity"] = _full_numeric(activity.mean(axis=1), idx, n_total)
    obs_columns.append({"key": "regulon_mean_activity", "kind": "numeric"})

    top_row = table_rows[0]
    bullets = [
        "内置 %d 个人源 TF（该数据集命中 %d 个），每个 TF 取共表达 |r| 最高的 %d 个正相关基因作为靶基因"
        % (len(tf_all), len(tfs_present), n_targets),
        "构建 %d 个 regulon（按平均靶基因相关排序），活性为 AUCell 式秩恢复分数（背景 %d 个基因）"
        % (len(regulons), background),
        "最强 regulon：%s（靶基因 %d 个，平均相关 %.3f）"
        % (top_row["tf"], top_row["n_targets"], top_row.get("mean_activity") or 0.0),
        "SCENIC-lite：**只用表达共表达，没有 motif/cisTarget 数据库**，不等于官方 SCENIC/pySCENIC 结果",
        "前 %d 个 regulon 活性已写入 obs（%s 等）" % (len(obs_columns) - 1, ", ".join(
            item["key"] for item in obs_columns[:3]
        )),
    ]
    summary = {
        "bullets": bullets,
        "metrics": {
            "n_tfs_in_data": int(len(tfs_present)),
            "n_regulons": int(len(regulons)),
            "n_genes_background": int(background),
            "n_targets_per_regulon": int(n_targets),
            "regulon_columns": [item["key"] for item in obs_columns],
        },
    }
    return {
        "summary": summary,
        "tables": [activity_table, targets_table],
        "figures": [figure_entry],
        "obs_columns": obs_columns,
        "obs_values": obs_values,
    }


def _priors_family(tf):
    """TF 家族（来自 biopriors.TF_FAMILY；缺失返回 "other"）。"""
    module = _load_biopriors()
    return module.tf_family(tf) or "other"


# ###########################################################################
# 11/12/13. virtual_knockout / virtual_overexpress / virtual_screen
# ###########################################################################


def _grn_matrices(view, params, warnings, required_genes=None):
    """构建共表达 GRN（scTenifoldKnk-lite 用）。

    - 取方差最高的 n_genes 个基因（默认 1000，上限 3000），并强制包含 required_genes；
    - 用稀疏均值/标准差计算基因间 Pearson 相关（cells 分块累加 Gram 矩阵）；
    - 每个基因保留 top_k 强相关边并对称化（带符号，行按绝对值归一化）。
    """
    n_genes = _clamp_int(params, "n_genes", 1000, 100, 3000)
    top_k = _clamp_int(params, "top_k", 10, 3, 50)
    matrix = view["X"]
    n_cells = matrix.shape[0]
    var_names = view["var_names"]
    gene_index = {gene: index for index, gene in enumerate(var_names)}
    selected = list(_hvg_indices(matrix, n_genes))
    if required_genes:
        for gene in required_genes:
            column = gene_index.get(str(gene))
            if column is None:
                raise EngineError(
                    "bad_request",
                    "目标基因 %r 不在表达矩阵的基因名里（virtual_* 需要基因符号精确匹配）" % str(gene),
                )
            if column not in selected:
                selected.append(column)
    selected = np.asarray(sorted(set(int(value) for value in selected)), dtype=np.int64)
    if selected.size < 50:
        raise EngineError(
            "prerequisite",
            "可用基因太少（%d），无法构建稳定的共表达 GRN（至少需要 50 个）" % selected.size,
        )
    labels = [var_names[int(index)] for index in selected]
    label_position = {label: position for position, label in enumerate(labels)}

    # 均值/标准差（稀疏）
    sub_matrix = matrix[:, selected]
    mean = np.asarray(sub_matrix.mean(axis=0)).ravel().astype(np.float64)
    mean_square = np.asarray(sub_matrix.multiply(sub_matrix).mean(axis=0)).ravel().astype(np.float64)
    variance = np.maximum(mean_square - mean * mean, 0.0)
    std = np.sqrt(variance)
    std[std < 1e-9] = 1.0

    # 相关矩阵：分块 cells，累加 (x-μ)(y-μ)
    gene_count = len(selected)
    gram = np.zeros((gene_count, gene_count), dtype=np.float64)
    chunk = 800
    for start in range(0, n_cells, chunk):
        end = min(n_cells, start + chunk)
        block = matrix[start:end][:, selected].toarray().astype(np.float64)
        centered = block - mean[None, :]
        gram += centered.T @ centered
    correlation = gram / ((max(n_cells - 1, 1)) * std[:, None] * std[None, :])
    np.fill_diagonal(correlation, 0.0)
    correlation = np.clip(correlation, -1.0, 1.0)
    if not np.all(np.isfinite(correlation)):
        correlation = np.nan_to_num(correlation, nan=0.0, posinf=0.0, neginf=0.0)

    # top-k 稀疏化 + 对称化（保留绝对值更大的符号）
    weights = np.zeros_like(correlation)
    absolute = np.abs(correlation)
    k = int(min(top_k, gene_count - 1))
    for row in range(gene_count):
        values = absolute[row]
        candidate = np.argpartition(-values, k)[:k]
        weights[row, candidate] = correlation[row, candidate]
    flipped = np.where(np.abs(weights) >= np.abs(weights.T), weights, weights.T)
    weights = flipped
    np.fill_diagonal(weights, 0.0)
    row_norm = np.sum(np.abs(weights), axis=1, keepdims=True)
    weights = weights / np.maximum(row_norm, 1e-12)
    return {
        "labels": labels,
        "label_position": label_position,
        "correlation": correlation,
        "weights": weights,
        "n_genes": gene_count,
        "mean": mean,
        "std": std,
        "top_k": int(k),
        "n_cells": int(n_cells),
    }


def _propagate_effect(weights, direct_effect, alpha, steps):
    """沿 GRN 传播扰动效应：total = direct + Σ α^t (Wᵀ)^t direct。"""
    total = np.array(direct_effect, dtype=np.float64, copy=True)
    current = np.array(direct_effect, dtype=np.float64, copy=True)
    for _ in range(int(steps)):
        current = alpha * (weights.T @ current)
        total += current
    return total


def _perturbation_analysis(view, params, outdir, warnings, mode):
    """KO / 过表达共享实现：direct effect（一阶线性）+ GRN 传播。"""
    target = str(params.get("target") or params.get("gene") or "").strip()
    if not target:
        raise EngineError(
            "bad_request",
            "%s 需要 params.target 指定目标基因（如 'STAT1' 或 'CEBPA'）"
            % ("virtual_knockout" if mode == "knockout" else "virtual_overexpress"),
        )
    grn = _grn_matrices(view, params, warnings, required_genes=[target])
    if target not in grn["label_position"]:
        raise EngineError("bad_request", "目标基因 %r 不在分析基因集里" % target)
    position = grn["label_position"][target]
    alpha = _clamp_float(params, "alpha", 0.5, 0.05, 0.95)
    steps = _clamp_int(params, "steps", 4, 1, 12)
    n_report = _clamp_int(params, "top_n", 30, 5, 100)
    mean_value = float(grn["mean"][position])
    std_value = float(grn["std"][position])
    if std_value <= 1e-9:
        raise EngineError("bad_request", "目标基因 %r 在所有细胞中表达恒定，无法做扰动分析" % target)
    direction = -1.0 if mode == "knockout" else 1.0
    # 一阶线性响应：目标基因变化 direction × mean → 其他基因 Δ ≈ r_ij × (Δ目标/sd目标) × sd_j
    direct_effect = (
        direction
        * (mean_value / std_value)
        * grn["correlation"][position, :]
        * grn["std"]
    )
    total_effect = _propagate_effect(grn["weights"], direct_effect, alpha, steps)

    rows = []
    for index, label in enumerate(grn["labels"]):
        if index == position:
            continue
        effect = float(total_effect[index])
        rows.append(
            {
                "gene": str(label),
                "corr_with_target": _num(grn["correlation"][position, index], 4),
                "direct_effect": _num(direct_effect[index], 5),
                "total_effect": _num(effect, 5),
                "abs_effect": abs(effect),
                "direction": "up" if effect >= 0 else "down",
            }
        )
    rows.sort(key=lambda row: -row["abs_effect"])
    reported = rows[:n_report]
    for rank, row in enumerate(reported, start=1):
        row["rank"] = rank
        row.pop("abs_effect")
    up_count = int(sum(1 for row in rows if row["total_effect"] and row["total_effect"] > 0))
    down_count = int(sum(1 for row in rows if row["total_effect"] and row["total_effect"] < 0))

    table = _table(
        "affected_genes",
        "Most affected genes after in-silico %s of %s"
        % ("knockout" if mode == "knockout" else "overexpression", target),
        [
            _column("rank", "Rank", "number"),
            _column("gene", "Gene"),
            _column("corr_with_target", "Correlation with target", "number"),
            _column("direct_effect", "Direct effect (first order)", "number"),
            _column("total_effect", "Propagated effect", "number"),
            _column("direction", "Direction"),
        ],
        reported,
        note="direct effect = 一阶线性近似（r × Δ目标/sd目标 × sd_j）；total effect = direct + GRN 传播（α=%s, steps=%d）。"
        % (alpha, steps),
    )

    figure, axis = _new_figure(6.2, 4.8)
    plot_rows = list(reversed(reported[:15]))
    y_positions = np.arange(len(plot_rows))
    values = [row["total_effect"] or 0.0 for row in plot_rows]
    colors = ["#D62728" if value > 0 else "#4C78A8" for value in values]
    axis.barh(y_positions, values, color=colors, height=0.68, edgecolor="black", linewidth=0.35)
    axis.set_yticks(y_positions)
    axis.set_yticklabels([row["gene"] for row in plot_rows], fontsize=8)
    axis.axvline(0, color="#333333", linewidth=0.7)
    axis.set_xlabel("Predicted expression change (log-normalized scale)")
    axis.set_title(
        "In-silico %s of %s: top affected genes"
        % ("knockout" if mode == "knockout" else "overexpression", target),
        fontsize=9,
    )
    figure.tight_layout()
    figure_entry = _save_figure(
        figure,
        view["ctx"],
        outdir,
        "virtual_knockout" if mode == "knockout" else "virtual_overexpress",
        "affected_genes",
        "Affected genes",
    )

    action = "敲除（表达置 0）" if mode == "knockout" else "过表达（+1 倍均值）"
    bullets = [
        "scTenifoldKnk-lite：对 %d 个高变基因构建共表达相关 GRN（top_k=%d 稀疏化、行列归一化）"
        % (grn["n_genes"], grn["top_k"]),
        "对目标基因 %s 做 in-silico %s，一阶线性效应沿 GRN 传播 α=%.2f、steps=%d（本插件近似实现，非官方 scTenifoldKnk）"
        % (target, action, alpha, steps),
        "受影响基因：上调 %d 个、下调 %d 个（共计入 %d 个基因）" % (up_count, down_count, len(rows)),
        "目标基因 %s 平均表达 %.3f（log-normalized），标准差 %.3f"
        % (target, mean_value, std_value),
        "**这是基于共表达的预测，不是实验验证的敲除结果**；GRN 为相关/稀疏化近似，未使用染色质先验与因果推断",
    ]
    summary = {
        "bullets": bullets,
        "metrics": {
            "mode": "knockout" if mode == "knockout" else "overexpression",
            "target_gene": target,
            "n_genes": int(grn["n_genes"]),
            "top_k": int(grn["top_k"]),
            "alpha": float(alpha),
            "steps": int(steps),
            "n_up": up_count,
            "n_down": down_count,
            "target_mean_expression": _num(mean_value, 4),
        },
    }
    return {
        "summary": summary,
        "tables": [table],
        "figures": [figure_entry],
        "obs_columns": [],
    }


def _run_virtual_knockout(view, params, outdir, warnings):
    return _perturbation_analysis(view, params, outdir, warnings, mode="knockout")


def _run_virtual_overexpress(view, params, outdir, warnings):
    return _perturbation_analysis(view, params, outdir, warnings, mode="overexpress")


def _cluster_signature_genes(view, cluster_id, top_n=50):
    """简单簇 marker：簇均值 - 其余细胞均值，取差值最大的正基因。"""
    clusters = _clusters(view, required=False)
    matrix = view["X"]
    var_names = view["var_names"]
    if clusters is None:
        raise EngineError(
            "prerequisite",
            "virtual_screen 默认签名来自某簇标记基因，需要先完成 cluster；"
            "或者用 params.signature 显式给一组基因（逗号/空格分隔）",
        )
    categories = _sorted_categories(clusters)
    if cluster_id is None or str(cluster_id).strip() == "":
        cluster_id = categories[0]
    cluster_id = str(cluster_id)
    if cluster_id not in categories:
        raise EngineError(
            "bad_request", "找不到簇 %r（可用簇：%s）" % (cluster_id, "、".join(categories))
        )
    mask = clusters == cluster_id
    if int(mask.sum()) < 10:
        raise EngineError("bad_request", "簇 %r 只有 %d 个细胞，太少，无法取签名" % (cluster_id, int(mask.sum())))
    in_mean = np.asarray(matrix[mask].mean(axis=0)).ravel()
    out_mean = np.asarray(matrix[~mask].mean(axis=0)).ravel() if (~mask).any() else np.zeros_like(in_mean)
    difference = in_mean - out_mean
    order = np.argsort(-difference, kind="stable")
    selected = [int(index) for index in order[: int(top_n)] if difference[index] > 0]
    if len(selected) < 5:
        raise EngineError("bad_request", "簇 %r 找不到足够的上调标记基因，无法作为 screen 签名" % cluster_id)
    return selected, [var_names[position] for position in selected], cluster_id


def _run_virtual_screen(view, params, outdir, warnings):
    n_candidates = _clamp_int(params, "n_candidates", 30, 3, 100)
    n_genes = _clamp_int(params, "n_genes", 1000, 100, 3000)
    alpha = _clamp_float(params, "alpha", 0.5, 0.05, 0.95)
    steps = _clamp_int(params, "steps", 4, 1, 12)
    signature_text = str(params.get("signature") or "").strip()
    cluster_id = str(params.get("cluster") or "").strip()
    matrix = view["X"]
    var_names = view["var_names"]

    if signature_text:
        requested = [gene for gene in re.split(r"[,\s;]+", signature_text) if gene]
        gene_index = {gene: index for index, gene in enumerate(var_names)}
        signature_columns = [gene_index[gene] for gene in requested if gene in gene_index]
        missing = [gene for gene in requested if gene not in gene_index]
        if missing:
            warnings.append("signature 中以下基因不在表达矩阵里，已忽略：%s" % "、".join(missing[:10]))
        if len(signature_columns) < 3:
            raise EngineError(
                "bad_request",
                "params.signature 解析后只剩 %d 个有效基因（至少需要 3 个）" % len(signature_columns),
            )
        signature_source = "params.signature (%d genes)" % len(signature_columns)
        signature_labels = [var_names[index] for index in signature_columns]
    else:
        signature_columns, signature_labels, cluster_id = _cluster_signature_genes(
            view, cluster_id or None, top_n=_clamp_int(params, "signature_size", 50, 5, 200)
        )
        signature_source = "top markers of cluster %r (%d genes)" % (cluster_id, len(signature_columns))

    # 候选扰动基因：默认取 GRN 基因里内置 TF；命中太少则退化为高变基因
    candidates_text = str(params.get("candidates") or "").strip()
    gene_index = {gene: index for index, gene in enumerate(var_names)}
    if candidates_text:
        requested = [gene for gene in re.split(r"[,\s;]+", candidates_text) if gene]
        candidates = [gene for gene in requested if gene in gene_index]
        if len(candidates) < 1:
            raise EngineError("bad_request", "params.candidates 里没有基因存在于表达矩阵")
        candidate_source = "params.candidates (%d genes)" % len(candidates)
    else:
        tf_all = _priors()[1]
        candidates = [tf for tf in tf_all if tf in gene_index]
        if len(candidates) < 3:
            warnings.append(
                "数据集里内置 TF 命中 %d 个，virtual_screen 候选退化为方差最高的高变基因" % len(candidates)
            )
            candidates = [var_names[int(index)] for index in _hvg_indices(matrix, 60)]
            candidate_source = "fallback: top variable genes (%d genes)" % len(candidates)
        else:
            candidate_source = "builtin TF list (%d candidates)" % len(candidates)
    candidates = candidates[: max(3, min(len(candidates), n_candidates * 2))]

    grn = _grn_matrices(view, params, warnings, required_genes=[signature_labels[0]])
    signature_positions = [grn["label_position"][label] for label in signature_labels if label in grn["label_position"]]
    if len(signature_positions) < 3:
        raise EngineError(
            "internal",
            "签名基因只有 %d 个落在 GRN 基因集里，无法计算 screen 分数" % len(signature_positions),
        )
    candidates = [gene for gene in candidates if gene in grn["label_position"]][:n_candidates]
    if not candidates:
        raise EngineError("bad_request", "候选基因都不在 GRN 基因集里，无法做虚拟筛选")

    signature_positions = np.asarray(signature_positions, dtype=np.int64)
    rows = []
    for gene in candidates:
        position = grn["label_position"][gene]
        if grn["std"][position] <= 1e-9:
            continue
        direct_effect = (
            -(grn["mean"][position] / grn["std"][position]) * grn["correlation"][position, :] * grn["std"]
        )
        total_effect = _propagate_effect(grn["weights"], direct_effect, alpha, steps)
        signature_effect = float(np.mean(total_effect[signature_positions]))
        down_in_signature = int(np.sum(total_effect[signature_positions] < 0))
        rows.append(
            {
                "candidate": str(gene),
                "signature_score": signature_effect,
                "abs_signature_effect": abs(signature_effect),
                "n_signature_down": down_in_signature,
                "n_signature_total": int(len(signature_positions)),
                "n_genes_shifted": int(np.sum(np.abs(total_effect) > 0.05)),
            }
        )
    if not rows:
        raise EngineError("internal", "virtual_screen 没有可评估的候选基因")
    rows.sort(key=lambda row: (row["signature_score"], row["candidate"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        row.pop("abs_signature_effect")
        row["signature_score"] = _num(row["signature_score"], 5)
    reported = rows
    table = _table(
        "screen_results",
        "Virtual screen: predicted effect of each candidate perturbation on the target signature",
        [
            _column("rank", "Rank", "number"),
            _column("candidate", "Candidate gene"),
            _column("signature_score", "Signature effect (KO)", "number"),
            _column("n_signature_down", "Signature genes down", "number"),
            _column("n_signature_total", "Signature genes", "number"),
            _column("n_genes_shifted", "Genes shifted (|Δ|>0.05)", "number"),
        ],
        reported,
        note="负分 = 敲除该基因会压低目标签名（推定为该签名的正向调控因子）；分数为一阶效应 + GRN 传播的近似。",
    )

    figure, axis = _new_figure(6.2, 4.8)
    plot_rows = list(reversed(rows[:15]))
    y_positions = np.arange(len(plot_rows))
    values = [row["signature_score"] or 0.0 for row in plot_rows]
    colors = ["#4C78A8" if value < 0 else "#D62728" for value in values]
    axis.barh(y_positions, values, color=colors, height=0.68, edgecolor="black", linewidth=0.35)
    axis.set_yticks(y_positions)
    axis.set_yticklabels([row["candidate"] for row in plot_rows], fontsize=8)
    axis.axvline(0, color="#333333", linewidth=0.7)
    axis.set_xlabel("Mean predicted signature change after knockout")
    axis.set_ylabel("Candidate gene")
    axis.set_title("Virtual screen (scTenifoldKnk-lite)", fontsize=9)
    figure.tight_layout()
    figure_entry = _save_figure(
        figure, view["ctx"], outdir, "virtual_screen", "screen_ranking", "Virtual screen ranking"
    )

    top = rows[0]
    bullets = [
        "目标签名：%s；来源：%s" % ("、".join(signature_labels[:8]) + ("…" if len(signature_labels) > 8 else ""), signature_source),
        "候选扰动基因：%d 个（%s），每个做 in-silico 敲除并沿 GRN 传播（α=%.2f、steps=%d）"
        % (len(rows), candidate_source, alpha, steps),
        "按对签名的平均预测变化排序；负分 = 敲除后签名下调（推定为正调控因子）",
        "最显著候选：%s（signature score %.4f，%d/%d 个签名基因下调）"
        % (top["candidate"], top["signature_score"] or 0.0, top["n_signature_down"], top["n_signature_total"]),
        "**这是共表达 GRN 上的预测性筛选，不是实验验证**；与官方 scTenifoldKnk 的口径不同（无张量分解/PC regression）",
    ]
    summary = {
        "bullets": bullets,
        "metrics": {
            "signature_size": int(len(signature_positions)),
            "signature_source": signature_source,
            "n_candidates": int(len(rows)),
            "candidate_source": candidate_source,
            "alpha": float(alpha),
            "steps": int(steps),
            "top_candidate": top["candidate"],
            "top_signature_score": top["signature_score"],
        },
    }
    return {
        "summary": summary,
        "tables": [table],
        "figures": [figure_entry],
        "obs_columns": [],
    }

# ###########################################################################
# ANALYSES 元数据（与 CONTRACT §9.3 sc-analysis-b 表逐项对应）
# 注意：title/method 会被 run() 原样写进结果，UI/Agent 直接展示这些字段。
# ###########################################################################

_LR_METHOD = (
    "CellChat/NATMI-lite（本插件实现）：内置人源配体-受体对 + 簇均值表达乘积打分 + "
    "打乱簇标签置换 p + BH 校正；不含 CellChat 官方数据库/复合物/trimean 细节"
)

ANALYSES = [
    {
        "id": "integration",
        "title": "多样本整合 / 去除批次效应",
        "group": "基本分析",
        "method": (
            "ComBat 式位置-尺度批次校正（本插件实现，无经验贝叶斯收缩）+ "
            "BBKNN 式批次平衡 kNN（按批次配额取邻居）；近似实现，结果不与官方 ComBat/BBKNN 逐值一致"
        ),
        "params": [
            {"key": "batch", "label": "批次列（obs 分类列）", "type": "column", "default": "batch",
             "hint": "obs 里的批次/样本列，例如 batch、sample、donor"},
            {"key": "n_neighbors", "label": "近邻数 k", "type": "number", "default": 15, "hint": "3-100"},
            {"key": "n_pcs", "label": "PCA 维数", "type": "number", "default": 30, "hint": "2-100"},
        ],
        "needs": ["groups:batch"],
    },
    {
        "id": "trajectory_dpt",
        "title": "扩散伪时间（monocle-like）",
        "group": "发育轨迹分析",
        "method": (
            "扩散图（自适应高斯核 + 对称归一化）+ 加权扩散分量距离的 DPT 简化实现；"
            "分支为扩散分量上的 KMeans 近似。本插件实现，非 monocle/官方 DPT"
        ),
        "params": [
            {"key": "root_cell", "label": "根细胞（obs 索引名，可留空自动）", "type": "text", "default": "",
             "hint": "留空 = 取第一扩散分量极值；也可填细胞条码或 0-based 位置"},
            {"key": "root_mode", "label": "自动选根方式", "type": "select", "default": "auto",
             "options": ["auto", "max_genes"],
             "hint": "auto=第一扩散分量极值；max_genes=检出基因数最多的细胞"},
            {"key": "n_neighbors", "label": "近邻数 k", "type": "number", "default": 15, "hint": "3-100"},
            {"key": "n_dcs", "label": "扩散分量数", "type": "number", "default": 4, "hint": "2-10"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "trajectory_paga",
        "title": "簇级轨迹图（PAGA-like）",
        "group": "发育轨迹分析",
        "method": (
            "PAGA-like（本插件实现）：簇间 kNN 边计数/随机期望得 connectivity，"
            "最大生成树 + ELG 边 + MST 距离簇级伪时间；非官方 PAGA"
        ),
        "params": [
            {"key": "root_cluster", "label": "根簇（簇名，可留空自动）", "type": "text", "default": "",
             "hint": "留空 = 检出基因数中位数最高的簇（最不分化）"},
            {"key": "n_neighbors", "label": "近邻数 k", "type": "number", "default": 15, "hint": "3-100"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "entropy_stemness",
        "title": "基于熵的干性分析",
        "group": "发育轨迹分析",
        "method": (
            "表达熵（高变基因上的归一化 Shannon 熵）+ kNN 邻居平滑的启发式打分；"
            "本插件启发式实现，不是已发表的干性评分方法"
        ),
        "params": [
            {"key": "n_genes", "label": "高变基因数", "type": "number", "default": 2000, "hint": "100-4000"},
            {"key": "n_neighbors", "label": "平滑近邻数 k", "type": "number", "default": 15, "hint": "3-100"},
        ],
        "needs": [],
    },
    {
        "id": "cytotrace_lite",
        "title": "CytoTRACE-like 分化潜能",
        "group": "发育轨迹分析",
        "method": (
            "CytoTRACE 思路的简化实现（本插件）：检出基因数 → kNN 平滑 → 与基因表达做 Spearman 相关 → "
            "取正相关 top 基因加权表达签名 → 邻居池化 → 秩变换；不是官方 CytoTRACE"
        ),
        "params": [
            {"key": "n_genes", "label": "参与相关的基因数", "type": "number", "default": 2000, "hint": "100-4000"},
            {"key": "n_neighbors", "label": "平滑近邻数 k", "type": "number", "default": 15, "hint": "3-100"},
            {"key": "signature_size", "label": "签名基因数", "type": "number", "default": 100, "hint": "10-500"},
        ],
        "needs": [],
    },
    {
        "id": "cnv_score",
        "title": "CNV 信号/ITH（inferCNV-like，需基因顺序）",
        "group": "肿瘤异质性分析",
        "method": (
            "inferCNV-like 代理实现（本插件）：参照细胞逐基因 z 分数 + 逐细胞中位数居中 + 基因顺序移动窗口平滑 + "
            "每细胞 |z| 均值；提供 gene_order 时按染色体位置排序，否则为按方差排序的**代理分数**，不能当真实 CNV 区段"
        ),
        "params": [
            {"key": "gene_order", "label": "基因位置文件（可选；gene  chr:start 或 gene chr start end）",
             "type": "file", "default": "", "hint": "留空 = 按方差排序的代理分数"},
            {"key": "reference", "label": "参照细胞（obs列:取值，可留空）", "type": "text", "default": "",
             "hint": "如 cell_type:Monocyte；留空=自动找免疫样细胞，找不到则全部细胞"},
            {"key": "window", "label": "平滑窗口（基因数，奇数）", "type": "number", "default": 11, "hint": "3-201"},
            {"key": "max_genes", "label": "最多使用基因数", "type": "number", "default": 4000, "hint": "100-4000"},
        ],
        "needs": ["gene_order"],
    },
    {
        "id": "cnv_tree",
        "title": "CNVscore/进化树",
        "group": "肿瘤异质性分析",
        "method": (
            "CNV 代理谱的簇级欧氏距离 + average-linkage 层次聚类（scipy.cluster.hierarchy），"
            "输出 newick/合并表/树图；是表型相似性树，不是真实克隆演化树"
        ),
        "params": [
            {"key": "gene_order", "label": "基因位置文件（可选）", "type": "file", "default": ""},
            {"key": "reference", "label": "参照细胞（obs列:取值，可留空）", "type": "text", "default": ""},
            {"key": "window", "label": "平滑窗口（基因数）", "type": "number", "default": 11, "hint": "3-201"},
            {"key": "max_genes", "label": "最多使用基因数", "type": "number", "default": 4000, "hint": "100-4000"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "ith_score",
        "title": "ITH_score 肿瘤异质性",
        "group": "肿瘤异质性分析",
        "method": (
            "基于 CNV 代理谱的簇内平均成对欧氏距离（每簇抽样上限可调），"
            "归一化后作为 ITH 代理分数；非克隆推断/非官方 ITH 工具"
        ),
        "params": [
            {"key": "gene_order", "label": "基因位置文件（可选）", "type": "file", "default": ""},
            {"key": "reference", "label": "参照细胞（obs列:取值，可留空）", "type": "text", "default": ""},
            {"key": "window", "label": "平滑窗口（基因数）", "type": "number", "default": 11, "hint": "3-201"},
            {"key": "max_genes", "label": "最多使用基因数", "type": "number", "default": 4000, "hint": "100-4000"},
            {"key": "max_cells_per_cluster", "label": "每簇抽样细胞数", "type": "number", "default": 120,
             "hint": "10-500，越大越慢"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "lr_communication",
        "title": "细胞通讯（CellChat/NATMI-lite）",
        "group": "细胞通讯分析",
        "method": _LR_METHOD,
        "params": [
            {"key": "n_perm", "label": "置换次数", "type": "number", "default": 50, "hint": "0-200；0=不做检验"},
            {"key": "top_n", "label": "展示交互条数", "type": "number", "default": 40, "hint": "5-200"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "regulon_activity",
        "title": "转录调控（SCENIC-lite）",
        "group": "转录调控分析",
        "method": (
            "SCENIC-lite（本插件）：内置 %d 个人源 TF + 共表达相关选靶基因 + AUCell 式秩恢复活性；"
            "**无 motif/cisTarget 数据库**，不等于官方 SCENIC/pySCENIC" % (len(_priors()[1]),)
        ),
        "params": [
            {"key": "n_genes", "label": "背景基因数（高变基因）", "type": "number", "default": 1500, "hint": "100-4000"},
            {"key": "n_targets", "label": "每个 regulon 靶基因数", "type": "number", "default": 20, "hint": "3-100"},
            {"key": "n_regulons", "label": "构建 regulon 数", "type": "number", "default": 30, "hint": "3-200"},
            {"key": "top_regulons", "label": "写回 obs 的 regulon 数", "type": "number", "default": 8, "hint": "1-30"},
        ],
        "needs": [],
    },
    {
        "id": "virtual_knockout",
        "title": "虚拟敲除（scTenifoldKnk-lite）",
        "group": "虚拟敲除",
        "method": (
            "scTenifoldKnk-lite（本插件）：共表达相关 GRN top-k 稀疏化 + 一阶线性扰动 + GRN 传播；"
            "非 PC regression/张量分解实现，结果为计算预测而非实验验证"
        ),
        "params": [
            {"key": "target", "label": "目标基因（必填）", "type": "text", "default": "", "hint": "如 STAT1、CEBPA"},
            {"key": "n_genes", "label": "GRN 基因数", "type": "number", "default": 1000, "hint": "100-3000"},
            {"key": "top_k", "label": "每个基因保留边数", "type": "number", "default": 10, "hint": "3-50"},
            {"key": "alpha", "label": "传播衰减 α", "type": "number", "default": 0.5, "hint": "0.05-0.95"},
            {"key": "steps", "label": "传播步数", "type": "number", "default": 4, "hint": "1-12"},
            {"key": "top_n", "label": "展示受影响基因数", "type": "number", "default": 30, "hint": "5-100"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "virtual_overexpress",
        "title": "虚拟过表达",
        "group": "虚拟敲除",
        "method": (
            "scTenifoldKnk-lite（本插件）：与虚拟敲除同一 GRN，方向取反（目标基因 +1 倍均值表达）；"
            "结果为计算预测而非实验验证"
        ),
        "params": [
            {"key": "target", "label": "目标基因（必填）", "type": "text", "default": "", "hint": "如 STAT1、CEBPA"},
            {"key": "n_genes", "label": "GRN 基因数", "type": "number", "default": 1000, "hint": "100-3000"},
            {"key": "top_k", "label": "每个基因保留边数", "type": "number", "default": 10, "hint": "3-50"},
            {"key": "alpha", "label": "传播衰减 α", "type": "number", "default": 0.5, "hint": "0.05-0.95"},
            {"key": "steps", "label": "传播步数", "type": "number", "default": 4, "hint": "1-12"},
            {"key": "top_n", "label": "展示受影响基因数", "type": "number", "default": 30, "hint": "5-100"},
        ],
        "needs": ["clusters"],
    },
    {
        "id": "virtual_screen",
        "title": "虚拟筛选",
        "group": "虚拟敲除",
        "method": (
            "scTenifoldKnk-lite 批量版（本插件）：对目标签名（默认某簇 top markers）逐一虚拟敲除候选 TF，"
            "按签名平均预测变化排序；共表达 GRN 近似，非实验筛选"
        ),
        "params": [
            {"key": "signature", "label": "目标签名基因（逗号/空格分隔；留空=用某簇 markers）",
             "type": "text", "default": "", "hint": "如 CD3D,CD3E,LCK"},
            {"key": "cluster", "label": "签名来源簇（signature 留空时使用）", "type": "text", "default": "",
             "hint": "簇名，留空=第一个簇"},
            {"key": "candidates", "label": "候选扰动基因（留空=内置 TF 列表）", "type": "text", "default": ""},
            {"key": "n_candidates", "label": "候选基因数", "type": "number", "default": 30, "hint": "3-100"},
            {"key": "signature_size", "label": "簇 marker 签名大小", "type": "number", "default": 50, "hint": "5-200"},
            {"key": "n_genes", "label": "GRN 基因数", "type": "number", "default": 1000, "hint": "100-3000"},
            {"key": "alpha", "label": "传播衰减 α", "type": "number", "default": 0.5, "hint": "0.05-0.95"},
            {"key": "steps", "label": "传播步数", "type": "number", "default": 4, "hint": "1-12"},
        ],
        "needs": ["clusters"],
    },
]

_SPEC_BY_ID = {str(spec["id"]): spec for spec in ANALYSES}

_RUNNERS = {
    "integration": _run_integration,
    "trajectory_dpt": _run_trajectory_dpt,
    "trajectory_paga": _run_trajectory_paga,
    "entropy_stemness": _run_entropy_stemness,
    "cytotrace_lite": _run_cytotrace_lite,
    "cnv_score": _run_cnv_score,
    "cnv_tree": _run_cnv_tree,
    "ith_score": _run_ith_score,
    "lr_communication": _run_lr_communication,
    "regulon_activity": _run_regulon_activity,
    "virtual_knockout": _run_virtual_knockout,
    "virtual_overexpress": _run_virtual_overexpress,
    "virtual_screen": _run_virtual_screen,
}


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------


def _validate_result(payload, ctx, spec):
    """校验固定返回键、obs_values 长度和图片文件（开发期兜底，异常抛 internal）。"""
    n_obs = int(ctx.get("n_obs") or 0)
    summary = payload.get("summary")
    if not isinstance(summary, dict):
        raise EngineError("internal", "%s 返回的 summary 不是对象" % spec["id"])
    summary.setdefault("bullets", [])
    summary.setdefault("metrics", {})
    payload["summary"] = summary
    tables = payload.get("tables")
    if not isinstance(tables, list):
        raise EngineError("internal", "%s 返回的 tables 不是数组" % spec["id"])
    for table in tables:
        if not isinstance(table, dict) or "id" not in table or "title" not in table or "columns" not in table or "rows" not in table:
            raise EngineError("internal", "%s 返回的 table 结构不合法：%r" % (spec["id"], table))
        keys = {str(column.get("key")) for column in table["columns"] if isinstance(column, dict)}
        for row in table["rows"]:
            if not isinstance(row, dict):
                raise EngineError("internal", "%s 表格行不是对象" % spec["id"])
            if not set(row).issubset(keys):
                raise EngineError(
                    "internal",
                    "%s 表格 %s 的行里有未声明的列：%s"
                    % (spec["id"], table["id"], sorted(set(row) - keys)),
                )
    figures = payload.get("figures")
    if not isinstance(figures, list):
        raise EngineError("internal", "%s 返回的 figures 不是数组" % spec["id"])
    for figure in figures:
        if not isinstance(figure, dict) or not figure.get("file"):
            raise EngineError("internal", "%s 返回的 figure 缺少 file 字段" % spec["id"])
        path = Path(str(figure["file"]))
        if not path.exists():
            raise EngineError("internal", "%s 的图片文件不存在：%s" % (spec["id"], path))
        if not figure.get("width") or not figure.get("height"):
            width, height = _png_size(path)
            figure["width"] = int(width)
            figure["height"] = int(height)
        if Path(str(path)).suffix.lower() != ".png":
            raise EngineError("internal", "%s 的图片必须是 PNG：%s" % (spec["id"], path))
    obs_columns = payload.get("obs_columns")
    if not isinstance(obs_columns, list):
        raise EngineError("internal", "%s 返回的 obs_columns 不是数组" % spec["id"])
    declared = {}
    for column in obs_columns:
        if not isinstance(column, dict) or not column.get("key"):
            raise EngineError("internal", "%s 的 obs_columns 条目缺少 key" % spec["id"])
        kind = str(column.get("kind") or "numeric")
        if kind not in ("numeric", "categorical"):
            raise EngineError("internal", "%s 的 obs_columns.kind 非法：%s" % (spec["id"], kind))
        column["kind"] = kind
        declared[str(column["key"])] = kind
    obs_values = payload.get("obs_values")
    if obs_values is not None:
        if not isinstance(obs_values, dict):
            raise EngineError("internal", "%s 返回的 obs_values 不是对象" % spec["id"])
        for key, values in obs_values.items():
            if not isinstance(values, (list, tuple, np.ndarray)):
                raise EngineError("internal", "%s.obs_values[%s] 不是数组" % (spec["id"], key))
            if n_obs and len(values) != n_obs:
                raise EngineError(
                    "internal",
                    "%s.obs_values[%s] 长度 %d != ctx.n_obs %d（必须与 ctx.obs.index 等长对齐）"
                    % (spec["id"], key, len(values), n_obs),
                )
            if key not in declared:
                raise EngineError(
                    "internal",
                    "%s.obs_values[%s] 没有在 obs_columns 里声明" % (spec["id"], key),
                )
            if declared[key] == "categorical":
                cleaned = [None if value is None else str(value) for value in values]
            else:
                cleaned = [
                    None
                    if value is None
                    else (float(value) if math.isfinite(float(value)) else None)
                    for value in values
                ]
            obs_values[key] = cleaned
        payload["obs_values"] = obs_values
    warnings = payload.get("warnings")
    if not isinstance(warnings, list):
        warnings = []
    payload["warnings"] = [str(item) for item in warnings]
    return payload


def run(analysis_id, ctx, params, outdir):
    """执行一个 sc-analysis-b 分析（契约见 CONTRACT.md §9.2）。"""
    spec = _SPEC_BY_ID.get(str(analysis_id))
    if spec is None:
        raise EngineError(
            "bad_request",
            "未知分析：%s（sc-analysis-b 可用：%s）"
            % (analysis_id, ", ".join(spec["id"] for spec in ANALYSES)),
        )
    if params is None:
        params = {}
    if not isinstance(params, dict):
        raise EngineError("bad_request", "params 必须是 JSON 对象")
    if not isinstance(ctx, dict):
        raise EngineError("internal", "分析 ctx 必须是对象")
    if not outdir:
        raise EngineError("internal", "缺少输出目录 outdir，无法保存图片")
    warnings = []
    view = _view(ctx, warnings)
    runner = _RUNNERS[spec["id"]]
    result = runner(view, params, outdir, warnings) or {}
    payload = {
        "id": str(spec["id"]),
        "title": str(spec["title"]),
        "method": str(spec["method"]),
        "summary": result.get("summary") or {"bullets": [], "metrics": {}},
        "tables": result.get("tables") or [],
        "figures": result.get("figures") or [],
        "obs_columns": result.get("obs_columns") or [],
        "warnings": list(warnings) + list(result.get("warnings") or []),
    }
    if "obs_values" in result:
        payload["obs_values"] = result["obs_values"]
    payload = _validate_result(payload, ctx, spec)
    return _jsonable(payload)


__all__ = ["ANALYSES", "run", "EngineError"]
