#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell Python compute worker.

JSON-lines over stdin/stdout（一行一个请求 / 响应），长驻进程。
┌ 协议契约：plugins/dsh-singlecell/CONTRACT.md §2（v1）与 §6（v2 追加）
├ stdout 只允许协议行；日志一律 stderr
├ 单条命令异常只影响该条响应，进程不退出（只有 shutdown / stdin EOF 退出）
└ builtin 引擎只依赖 numpy / scipy / scikit-learn / anndata / matplotlib；
   检测到 scanpy 时 mode="scanpy" 并尽量走 scanpy 路径，失败自动回退 builtin。

启动方式（Host 负责）：
    python3 -u lib/sc_engine.py --workdir <abs_dir>
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

# ===========================================================================
# 0. 环境准备（必须在 numpy / matplotlib 导入之前完成）
# ===========================================================================


def _parse_cli(argv):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--workdir", default=None)
    known, _ = parser.parse_known_args(argv)
    return known


_CLI = _parse_cli(sys.argv[1:])
WORKDIR = Path(
    _CLI.workdir or os.environ.get("DS_SC_WORKDIR") or os.getcwd()
).expanduser().resolve()
FIG_DIR = WORKDIR / "figures"
EXPORT_DIR = WORKDIR / "exports"
MPL_DIR = WORKDIR / "mpl-cache"
for _dir in (WORKDIR, FIG_DIR, EXPORT_DIR, MPL_DIR):
    try:
        _dir.mkdir(parents=True, exist_ok=True)
    except Exception:  # pragma: no cover - 目录不可写时交给后续操作报错
        pass

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

# 不往插件目录写 __pycache__（安装目录可能只读；也让产物更干净）
sys.dont_write_bytecode = True

try:  # 统一 UTF-8，避免中文响应在窄 locale 下炸掉
    sys.stdout.reconfigure(encoding="utf-8", errors="replace", newline="\n")
    sys.stdin.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # pragma: no cover
    pass

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import scipy.sparse.csgraph as _csgraph  # noqa: E402
from scipy import stats as _stats  # noqa: E402

import sklearn  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.manifold import TSNE  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402

try:
    import anndata as ad  # noqa: E402
except Exception:  # pragma: no cover
    ad = None

import matplotlib  # noqa: E402

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402

try:  # scanpy / umap 均为可选依赖（本机未安装）
    import scanpy as _scanpy  # noqa: E402

    SCANPY = _scanpy
except Exception:  # pragma: no cover
    SCANPY = None

try:  # pragma: no cover - 可选
    import umap as _umap  # noqa: E402

    HAS_UMAP = True
except Exception:  # pragma: no cover
    HAS_UMAP = False

try:  # pragma: no cover - 可选
    import leidenalg as _leidenalg  # noqa: E402

    HAS_LEIDENALG = True
except Exception:  # pragma: no cover
    HAS_LEIDENALG = False

ENGINE_MODE = "scanpy" if SCANPY is not None else "builtin"
SEED = 0
PY_VERSION = "%d.%d.%d" % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
STARTED_AT = time.time()
PLUGIN_VERSION = "0.1.0"


def _plugin_version():
    """从 package.json 读插件版本（缺失时用常量兜底）。"""
    package_json = Path(__file__).resolve().parent.parent / "package.json"
    try:
        data = json.loads(package_json.read_text(encoding="utf-8"))
        version = data.get("version")
        if isinstance(version, str) and version:
            return version
    except Exception:
        pass
    return PLUGIN_VERSION


try:
    PLUGIN_VERSION = _plugin_version()
except Exception:  # pragma: no cover
    pass


def _pkg_version(distribution, module=None, attr="__version__"):
    try:
        from importlib.metadata import version as _version

        return _version(distribution)
    except Exception:
        pass
    try:
        module = module or distribution
        mod = __import__(module)
        value = getattr(mod, attr, None)
        if isinstance(value, str) and value:
            return value
    except Exception:
        pass
    return None


PACKAGE_VERSIONS = {
    "numpy": _pkg_version("numpy"),
    "pandas": _pkg_version("pandas"),
    "scipy": _pkg_version("scipy"),
    "sklearn": _pkg_version("scikit-learn", "sklearn"),
    "anndata": _pkg_version("anndata"),
    "scanpy": _pkg_version("scanpy", "scanpy") if SCANPY is not None else None,
    "matplotlib": _pkg_version("matplotlib"),
}

LAYOUTS = ["pca", "tsne"] + (["umap"] if (HAS_UMAP or SCANPY is not None) else [])

STARTUP_WARNINGS = []
if ENGINE_MODE == "builtin":
    STARTUP_WARNINGS.append(
        "未检测到 scanpy：使用内置 numpy/scipy/scikit-learn 引擎（无 UMAP/Leiden，"
        "聚类为纯 Python Louvain 实现）"
    )
STARTUP_WARNINGS.append(f"matplotlib 使用 Agg 后端，配置缓存目录：{MPL_DIR}")

# 图表全局样式：白底、小字号、英文标签（本机无 CJK 字体）
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


# ===========================================================================
# 1. 协议基础工具
# ===========================================================================


class EngineError(Exception):
    """可预期的业务错误：code ∈ CONTRACT §2.1 / §8 的错误码集合。"""

    CODES = {"bad_request", "not_found", "prerequisite", "unsupported", "internal", "ambiguous"}

    def __init__(self, code, message, details=None):
        self.code = code if code in self.CODES else "internal"
        self.message = str(message)
        self.details = details or None
        super().__init__(self.message)


def _log(message):
    try:
        sys.stderr.write("[sc-engine %s] %s\n" % (time.strftime("%H:%M:%S"), message))
        sys.stderr.flush()
    except Exception:  # pragma: no cover
        pass


def _jsonable(value):
    """把 numpy / pandas 标量与容器递归转换成可 JSON 序列化的纯 Python 对象。"""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, np.ndarray):
        return [_jsonable(item) for item in value.tolist()]
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if isinstance(value, (pd.Series, pd.Index)):
        return [_jsonable(item) for item in value.tolist()]
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _emit(payload):
    try:
        text = json.dumps(_jsonable(payload), ensure_ascii=False, allow_nan=False)
    except Exception:  # pragma: no cover - 极端兜底：再洗一次
        text = json.dumps(_jsonable(json.loads(json.dumps(payload, default=str))), ensure_ascii=False)
    sys.stdout.write(text + "\n")
    sys.stdout.flush()


def _now_iso():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def _fmt_num(value):
    """给 Methods / 报告用的数字格式化。"""
    try:
        number = float(value)
    except Exception:
        return str(value)
    if not math.isfinite(number):
        return str(value)
    if abs(number - round(number)) < 1e-9:
        return str(int(round(number)))
    if abs(number) >= 1000 or abs(number) < 0.001:
        return ("%.4g" % number)
    return ("%.4g" % number)


def _as_int(value, name):
    try:
        if isinstance(value, bool):
            raise ValueError
        number = float(value)
    except Exception:
        raise EngineError("bad_request", "参数 %s 需要整数，收到 %r" % (name, value))
    if not math.isfinite(number):
        raise EngineError("bad_request", "参数 %s 需要有限数值" % name)
    return int(round(number))


def _as_float(value, name):
    try:
        if isinstance(value, bool):
            raise ValueError
        number = float(value)
    except Exception:
        raise EngineError("bad_request", "参数 %s 需要数值，收到 %r" % (name, value))
    if not math.isfinite(number):
        raise EngineError("bad_request", "参数 %s 需要有限数值" % name)
    return number


def _as_optional_int(value, name):
    if value is None or value == "":
        return None
    return _as_int(value, name)


def _as_optional_float(value, name):
    if value is None or value == "":
        return None
    return _as_float(value, name)


def _limited_float(value, name, minimum, maximum):
    number = _as_float(value, name)
    if number < minimum or number > maximum:
        raise EngineError(
            "bad_request", "参数 %s 需要在 [%s, %s] 之间，收到 %s" % (name, minimum, maximum, number)
        )
    return number


# ---------------------------------------------------------------------------
# 调色板（数据色：服务端下发，UI 直接使用）
# ---------------------------------------------------------------------------

_TAB20 = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948",
    "#b07aa1", "#ff9da7", "#9c755f", "#bab0ab", "#1f77b4", "#ff7f0e",
    "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
]


def _palette(count):
    """确定性调色板：前 20 个用 Tableau 色，之后用黄金角 HSV。"""
    colors = list(_TAB20[: max(0, count)])
    if count <= len(_TAB20):
        return colors
    import colorsys

    index = len(colors)
    while len(colors) < count:
        hue = (index * 0.618033988749895) % 1.0
        red, green, blue = colorsys.hsv_to_rgb(hue, 0.58, 0.86)
        colors.append("#%02x%02x%02x" % (int(red * 255), int(green * 255), int(blue * 255)))
        index += 1
    return colors


def _gradient_colors(count):
    cmap = matplotlib.colormaps.get_cmap("viridis")
    if count <= 1:
        return [to_hex(cmap(0.5))]
    return [to_hex(cmap(index / (count - 1))) for index in range(count)]


def _png_dimensions(path):
    """只读 PNG 头（IHDR）拿宽高，避免引入 PIL 依赖。"""
    try:
        with open(path, "rb") as handle:
            header = handle.read(24)
        if header[:8] != b"\x89PNG\r\n\x1a\n":
            return None
        width = int.from_bytes(header[16:20], "big")
        height = int.from_bytes(header[20:24], "big")
        return width, height
    except Exception:
        return None


def _bh_fdr(pvalues):
    """Benjamini-Hochberg FDR。"""
    pvalues = np.asarray(pvalues, dtype=np.float64)
    total = pvalues.size
    if total == 0:
        return pvalues
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order]
    ranks = np.arange(1, total + 1, dtype=np.float64)
    adjusted = ranked * total / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    output = np.empty_like(adjusted)
    output[order] = adjusted
    return output


# ===========================================================================
# 2. 纯 Python Louvain（模块度优化 + resolution）
# ===========================================================================


def _adjacency_lists(matrix):
    """CSR → 每个节点的邻居索引 / 权重列表。"""
    count = matrix.shape[0]
    indptr = matrix.indptr
    indices = matrix.indices
    data = matrix.data
    neighbours = [indices[indptr[i]:indptr[i + 1]] for i in range(count)]
    weights = [data[indptr[i]:indptr[i + 1]].astype(np.float64) for i in range(count)]
    return neighbours, weights


def _louvain_local_moving(neighbours, weights, resolution, max_passes=100):
    """一层 Louvain 局部移动。返回每个节点的社区编号（原始、未重排）。"""
    count = len(neighbours)
    degree = np.array([float(w.sum()) for w in weights], dtype=np.float64)
    two_m = float(degree.sum())
    communities = np.arange(count, dtype=np.int64)
    if count <= 1 or two_m <= 0:
        return communities
    sigma_total = degree.copy()
    epsilon = 1e-12
    for _pass in range(max_passes):
        moved = 0
        for node in range(count):
            node_degree = degree[node]
            if node_degree <= 0:
                continue
            current = int(communities[node])
            sigma_total[current] -= node_degree
            # 统计搬到每个邻居社区后的“内部权重”
            weight_to = {}
            for neighbour, weight in zip(neighbours[node], weights[node]):
                community = int(communities[neighbour])
                weight_to[community] = weight_to.get(community, 0.0) + float(weight)
            best_community = current
            best_gain = weight_to.get(current, 0.0) - resolution * sigma_total[current] * node_degree / two_m
            for community, weight in weight_to.items():
                if community == current:
                    continue
                gain = weight - resolution * sigma_total[community] * node_degree / two_m
                if gain > best_gain + epsilon or (
                    abs(gain - best_gain) <= epsilon and community < best_community
                ):
                    best_gain = gain
                    best_community = community
            communities[node] = best_community
            sigma_total[best_community] += node_degree
            if best_community != current:
                moved += 1
        if moved == 0:
            break
    return communities


def _aggregate_louvain_graph(neighbours, weights, communities, community_count):
    """把社区聚合成超节点图（自环 = 社区内部权重和）。"""
    rows = []
    cols = []
    values = []
    for node in range(len(neighbours)):
        source = int(communities[node])
        row_weights = {}
        for neighbour, weight in zip(neighbours[node], weights[node]):
            target = int(communities[neighbour])
            row_weights[target] = row_weights.get(target, 0.0) + float(weight)
        for target, weight in row_weights.items():
            rows.append(source)
            cols.append(target)
            values.append(weight)
    matrix = sp.csr_matrix(
        (np.asarray(values, dtype=np.float64), (rows, cols)),
        shape=(community_count, community_count),
    )
    matrix = ((matrix + matrix.T) * 0.5).tocsr()
    return matrix


def _relabel_by_size(labels):
    """簇编号按大小降序重排为 0,1,2…（同大小按原编号升序，保证确定性）。"""
    labels = np.asarray(labels)
    if labels.size == 0:
        return labels.astype(np.int32)
    values, counts = np.unique(labels, return_counts=True)
    order = sorted(range(values.size), key=lambda i: (-int(counts[i]), int(values[i])))
    mapping = {int(values[index]): new for new, index in enumerate(order)}
    return np.array([mapping[int(label)] for label in labels], dtype=np.int32)


def louvain_clusters(matrix, resolution=1.0, max_levels=20, max_passes=100):
    """在加权无向图上跑纯 Python Louvain，返回按大小降序重排的 int32 标签。"""
    node_count = matrix.shape[0]
    if node_count == 0:
        return np.zeros(0, dtype=np.int32)
    if node_count == 1:
        return np.zeros(1, dtype=np.int32)
    node_map = np.arange(node_count, dtype=np.int64)
    current = matrix.tocsr()
    for _level in range(max_levels):
        current_count = current.shape[0]
        if current_count <= 1:
            break
        neighbours, weights = _adjacency_lists(current)
        raw = _louvain_local_moving(neighbours, weights, float(resolution), max_passes=max_passes)
        values = np.unique(raw)
        community_count = int(values.size)
        remap = np.zeros(current_count, dtype=np.int64)
        remap[values] = np.arange(community_count, dtype=np.int64)
        communities = remap[raw]
        node_map = communities[node_map]
        if community_count >= current_count:
            break
        current = _aggregate_louvain_graph(neighbours, weights, communities, community_count)
    return _relabel_by_size(node_map)


# ===========================================================================
# 3. 数据集对象与状态
# ===========================================================================

DATASETS = {}
_DATASET_COUNTER = 0


def _as_csr_float32(matrix):
    if matrix is None:
        raise EngineError("internal", "矩阵为空")
    if sp.issparse(matrix):
        return sp.csr_matrix(matrix, dtype=np.float32)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float32))


def _integer_like(matrix) -> bool:
    """矩阵是否像原始计数（非负整数）。"""
    try:
        if sp.issparse(matrix):
            if not matrix.nnz:
                return True
            sample = matrix.data[: min(2000, matrix.nnz)]
        else:
            flat = np.asarray(matrix, dtype=np.float64).ravel()
            if not flat.size:
                return True
            sample = flat[: min(2000, flat.size)]
        sample = np.asarray(sample, dtype=np.float64)
        sample = sample[np.isfinite(sample)]
        if sample.size == 0:
            return True
        if np.any(sample < -1e-6):
            return False
        return bool(np.all(np.abs(sample - np.round(sample)) < 1e-6))
    except Exception:  # noqa: BLE001
        return False


class Dataset:
    """一个常驻内存的数据集（含原始矩阵、当前矩阵、obs/var/obsm/uns、历史）。"""

    def __init__(self, dataset_id, name, path, full_x, full_counts, var_names, obs,
                 var=None, warnings=None, source="h5ad"):
        self.id = str(dataset_id)
        self.name = str(name)
        self.path = str(path) if path else None
        self.source = str(source)
        self.full_counts = _as_csr_float32(full_counts)
        self.full_x = _as_csr_float32(full_x if full_x is not None else full_counts)
        if self.full_x.shape != self.full_counts.shape:
            self.full_x = self.full_counts.copy()
        self.full_var_names = [str(gene) for gene in var_names]
        if len(self.full_var_names) != self.full_counts.shape[1]:
            raise EngineError("bad_request", "基因名数量与表达矩阵列数不一致")
        self.full_var = (
            var.copy() if var is not None else pd.DataFrame(index=pd.Index(self.full_var_names, name="gene"))
        )
        self.full_obs = obs.copy()
        if self.full_obs.shape[0] != self.full_counts.shape[0]:
            raise EngineError("bad_request", "obs 行数与表达矩阵行数不一致")

        # 当前（过滤后）矩阵 / 注释
        self.counts = self.full_counts.copy()
        self.X = self.full_x.copy()
        self.var_names = list(self.full_var_names)
        self.var = self.full_var.copy()
        self.obs = self.full_obs.copy()
        self.layers = {"counts": self.counts.copy()}

        self.obsm = {}
        self.uns = {}
        self.adj = None
        self.knn = None
        self.markers = None
        self.layout = None
        self.steps = []
        self.params = {}
        # EXPORT-01：结构化执行记录 —— Methods / 导出脚本 / 复现清单的唯一事实来源
        # 每条记录包含 requested_method / actual_method / 软件版本 / 实际生效参数 /
        # 输入输出矩阵尺度 / 聚类 resolution / 种子 / 回退原因 / 数据版本。
        self.execution_record = []
        self.qc_metrics = None
        self.filter_info = None
        self.version = 0
        # SC-04：任何数据/注释变更都要让结果缓存失效，所以单独维护 revision
        self.revision = 0
        self.history = []
        self.figures = []
        self.resolution_scan = None
        self.warnings = list(warnings or [])
        # 导入数据若自带 cell_type（真实注释），聚类步骤不应把它覆盖成簇号
        self.cell_type_default = "cell_type" not in self.full_obs.columns
        self.source_format = None
        # SC-02：counts 的语义必须显式记录，不能再靠「矩阵名字」或「有就当作原始计数」
        self.counts_source = "counts"
        self.counts_verified = True
        self.matrix_kind = "unknown"
        self.notes = []
        self.meta = {}
        self.analyses_run = []
        self.created_at = _now_iso()

    def touch(self, reason=""):
        """标记数据已变更（分组/注释/过滤/子集…），用于结果缓存失效（SC-04）。"""
        self.revision = int(getattr(self, "revision", 0)) + 1
        if reason:
            try:
                self.uns.setdefault("revisions", []).append({"at": _now_iso(), "reason": str(reason)})
            except Exception:  # noqa: BLE001
                pass
        return self.revision

    # -- 便捷状态 ---------------------------------------------------------
    @property
    def n_obs(self):
        return int(self.obs.shape[0])

    @property
    def n_vars(self):
        return int(len(self.var_names))

    def obs_all(self):
        """全部细胞（未过滤）的 obs，并附上 QC 指标列。"""
        obs = self.full_obs.copy()
        if self.qc_metrics is not None:
            for column in self.qc_metrics.columns:
                obs[column] = self.qc_metrics[column].reindex(obs.index).to_numpy()
        return obs

    def var_index(self):
        return {str(gene): index for index, gene in enumerate(self.var_names)}

    def gene_column(self, gene):
        """取某个基因在当前矩阵里的列（稀疏列 → 稠密一维数组）。"""
        index = self.var_index().get(str(gene))
        if index is None:
            return None
        column = self.X[:, index]
        if sp.issparse(column):
            return np.asarray(column.todense()).ravel()
        return np.asarray(column).ravel()

    def cluster_codes(self):
        """('0','1',…) 有序簇 id + 每细胞 code；没有聚类时返回 (None, None)。"""
        if "clusters" not in self.obs.columns:
            return None, None
        values = [str(value) for value in self.obs["clusters"].tolist()]

        def sort_key(item):
            try:
                return (0, float(item))
            except Exception:
                return (1, str(item))

        categories = sorted(set(values), key=sort_key)
        lookup = {label: index for index, label in enumerate(categories)}
        codes = np.array([lookup[value] for value in values], dtype=np.int64)
        return categories, codes

    def push_history(self, entry):
        self.history.append(entry)
        if len(self.history) > 50:
            self.history = self.history[-50:]


def _next_dataset_id():
    global _DATASET_COUNTER
    _DATASET_COUNTER += 1
    return "ds%d" % _DATASET_COUNTER


def _resolve_dataset(args, required=True):
    """取 args['dataset'] 指向的数据集；省略时回退到最近打开/唯一的数据集。"""
    dataset_id = args.get("dataset")
    if dataset_id is not None:
        key = str(dataset_id)
        dataset = DATASETS.get(key)
        if dataset is None:
            raise EngineError("not_found", "数据集不存在：%s" % key)
        return dataset
    if DATASETS:
        return list(DATASETS.values())[-1]
    if required:
        raise EngineError("bad_request", "还没有载入任何数据集，请先 open 或载入 demo")
    return None


def _obs_columns(dataset):
    columns = []
    for name in dataset.obs.columns:
        text = str(name)
        if text.startswith("_"):
            continue
        series = dataset.obs[name]
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            kind = "numeric"
        else:
            kind = "categorical"
        try:
            n_unique = int(series.nunique(dropna=True))
        except Exception:
            n_unique = 0
        entry = {"key": text, "kind": kind, "n_unique": n_unique}
        if kind == "numeric":
            try:
                numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
                finite = numeric[np.isfinite(numeric)]
                if finite.size:
                    entry["min"] = float(finite.min())
                    entry["max"] = float(finite.max())
            except Exception:
                pass
        columns.append(entry)
    return columns


def _dataset_summary(dataset):
    summary = {
        "id": dataset.id,
        "name": dataset.name,
        "path": dataset.path,
        "n_obs": dataset.n_obs,
        "n_vars": dataset.n_vars,
        "obs_columns": _obs_columns(dataset),
        "steps": list(dataset.steps),
        "current": dataset.steps[-1] if dataset.steps else None,
        "warnings": list(dataset.warnings),
        "history": dataset.history[-20:],
        "source_format": getattr(dataset, "source_format", None) or dataset.source,
        "notes": list(getattr(dataset, "notes", []) or []),
        "analyses_run": list(getattr(dataset, "analyses_run", []) or []),
        "samples": list((getattr(dataset, "uns", {}) or {}).get("samples") or []),
        "groups": _summary_groups(dataset),
        "group_column": (getattr(dataset, "uns", {}) or {}).get("group_column"),
        "import_state": (getattr(dataset, "uns", {}) or {}).get("import_state"),
    }
    return summary


def _summary_groups(dataset):
    uns = getattr(dataset, "uns", {}) or {}
    column = uns.get("group_column")
    if column and column in dataset.obs.columns:
        return sorted({str(item) for item in dataset.obs[column].dropna().unique().tolist()})
    return []


# ===========================================================================
# 4. 数据加载（h5ad / demo）
# ===========================================================================


def _make_sample_module():
    """按绝对路径加载 test/make_sample.py（demo 生成器与 marker 面板的唯一来源）。"""
    path = Path(__file__).resolve().parent.parent / "test" / "make_sample.py"
    if not path.exists():
        raise EngineError("unsupported", "缺少 demo 生成器 test/make_sample.py：%s" % path)
    module_name = "_dsh_sc_make_sample"
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise EngineError("unsupported", "demo 生成器无法加载：%s" % path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as error:  # pragma: no cover
        sys.modules.pop(module_name, None)
        raise EngineError("internal", "demo 生成器加载失败：%s" % error)
    return module


def _load_demo_dataset():
    module = _make_sample_module()
    try:
        data = module.generate_demo()
    except Exception as error:  # pragma: no cover
        raise EngineError("internal", "demo 数据生成失败：%s" % error)
    warnings = []
    return Dataset(
        _next_dataset_id(),
        "demo-pbmc",
        None,
        data["counts"],
        data["counts"],
        data["var_names"],
        data["obs"],
        warnings=warnings,
        source="demo",
    )


def _read_h5_group_with_adata(path):
    """余力支持：普通 .h5 文件里如果有 `adata` 组（h5ad 结构），尽力读出来。

    做法是把该组先复制成 workdir 下的临时 .h5ad 再用 anndata 读，成功后删临时文件。
    """
    try:
        import h5py
    except Exception:
        return None
    try:
        with h5py.File(path, "r") as source:
            if "adata" not in source:
                return None
            group = source["adata"]
            if not isinstance(group, h5py.Group) or "X" not in group or "obs" not in group:
                return None
            temp_dir = WORKDIR / "tmp"
            temp_dir.mkdir(parents=True, exist_ok=True)
            digest = hashlib.md5(str(path).encode("utf-8")).hexdigest()[:12]
            temp_path = temp_dir / ("h5open-%s.h5ad" % digest)
            try:
                with h5py.File(temp_path, "w") as target:
                    for key in group.keys():
                        source.copy(group[key], target, key)
                    for attr_key, attr_value in group.attrs.items():
                        target.attrs[attr_key] = attr_value
                return ad.read_h5ad(str(temp_path))
            finally:
                try:
                    temp_path.unlink()
                except Exception:
                    pass
    except Exception:
        return None


def _load_h5ad_dataset(path, name=None):
    if ad is None:
        raise EngineError("unsupported", "未安装 anndata，无法读取 .h5ad 文件")
    resolved = Path(os.path.expanduser(str(path))).resolve()
    if not resolved.exists():
        raise EngineError("bad_request", "文件不存在：%s" % path)
    if not resolved.is_file():
        raise EngineError("bad_request", "路径不是文件：%s" % path)
    try:
        adata = ad.read_h5ad(str(resolved))
    except Exception as error:
        adata = _read_h5_group_with_adata(resolved)
        if adata is None:
            raise EngineError("bad_request", "无法解析 .h5ad 文件（%s）：%s" % (resolved.name, error))
        if name is None:
            name = "%s::adata" % resolved.name
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise EngineError("bad_request", "h5ad 里没有任何细胞或基因：%s" % resolved.name)

    matrix = adata.X
    if matrix is None:
        raise EngineError("bad_request", "h5ad 里没有 X 表达矩阵：%s" % resolved.name)
    warnings = []
    counts = None
    if "counts" in getattr(adata, "layers", {}):
        counts = adata.layers["counts"]
    elif getattr(adata, "raw", None) is not None:
        counts = adata.raw.X
        warnings.append("未找到 layers['counts']，已使用 adata.raw.X 作为原始计数")
    else:
        counts = matrix
        warnings.append("未找到 layers['counts']，已把 X 当作原始计数（若 X 已标准化，QC 指标可能失真）")

    var_names = [str(gene) for gene in adata.var_names]
    obs = adata.obs.copy()
    if obs.index.empty:
        obs = pd.DataFrame(index=pd.Index(["cell%d" % i for i in range(adata.n_obs)], name="cell"))
    var = adata.var.copy()
    if len(var) != len(var_names):
        var = None
    try:
        matrix_dtype = matrix.dtype
        is_integer_like = np.issubdtype(matrix_dtype, np.integer)
    except Exception:
        is_integer_like = False
    if not is_integer_like and sp.issparse(matrix) and matrix.nnz:
        sample = matrix.data[: min(1000, matrix.nnz)]
        if np.all(np.abs(sample - np.round(sample)) < 1e-6) and float(np.max(sample)) > 0:
            pass
        elif float(np.max(sample)) < 30:
            warnings.append("X 看起来已经过 log 标准化；内置流程的 normalize 步骤会基于 counts 重新计算")
    return Dataset(
        _next_dataset_id(),
        name or resolved.name,
        str(resolved),
        matrix,
        counts,
        var_names,
        obs,
        var=var,
        warnings=warnings,
        source="h5ad",
    )


# ===========================================================================
# 5. 步骤参数（默认值必须与 CONTRACT §2.2 表格逐字一致）
# ===========================================================================

STEP_ORDER = [
    "qc", "filter", "normalize", "hvg", "pca",
    "neighbors", "cluster", "embed", "markers", "annotate",
]

STEP_DEPS = {
    "qc": [],
    "filter": ["qc"],
    "normalize": ["filter"],
    "hvg": ["normalize"],
    "pca": ["hvg"],
    "neighbors": ["pca"],
    "cluster": ["neighbors"],
    "embed": ["pca"],
    "markers": ["cluster"],
    "annotate": ["cluster"],
}

STEP_DEFAULTS = {
    "qc": {"min_genes": 200, "max_genes": None, "min_counts": None, "max_pct_mt": 20.0},
    "filter": {"min_genes": 200, "max_genes": None, "min_counts": None, "max_pct_mt": 20.0, "min_cells": 3},
    "normalize": {"target_sum": 10000.0},
    "hvg": {"n_top_genes": 2000, "flavor": "seurat"},
    "pca": {"n_pcs": 30},
    "neighbors": {"n_neighbors": 15, "metric": "cosine"},
    "cluster": {"resolution": 0.2, "method": "louvain"},
    "embed": {"layout": "umap", "perplexity": 30.0, "n_iter": 1000},
    "markers": {"groupby": "clusters", "top_n": 10, "method": "t-test"},
    "annotate": {"mapping": {}},
}

STEP_TITLES = {
    "qc": "QC 指标",
    "filter": "过滤",
    "normalize": "标准化",
    "hvg": "高变基因",
    "pca": "PCA",
    "neighbors": "邻接图",
    "cluster": "聚类",
    "embed": "二维嵌入",
    "markers": "标记基因",
    "annotate": "细胞类型注释",
}


def _merge_params(step, user_params):
    """合并默认参数并做类型/范围校验，返回 effective params。"""
    if user_params is None:
        user_params = {}
    if not isinstance(user_params, dict):
        raise EngineError("bad_request", "params 必须是 JSON 对象")
    merged = dict(STEP_DEFAULTS[step])
    for key, value in user_params.items():
        if value is None:
            # null 视为“使用默认值”（避免 UI 传 null 时把阈值误当成 0 / 关闭）
            merged[key] = STEP_DEFAULTS[step].get(key)
            continue
        merged[key] = value

    if step in ("qc", "filter"):
        merged["min_genes"] = max(0, _as_int(merged.get("min_genes", 0), "min_genes"))
        merged["max_genes"] = _as_optional_int(merged.get("max_genes"), "max_genes")
        if merged["max_genes"] is not None and merged["max_genes"] < 0:
            raise EngineError("bad_request", "max_genes 不能为负")
        merged["min_counts"] = _as_optional_int(merged.get("min_counts"), "min_counts")
        merged["max_pct_mt"] = _limited_float(
            merged.get("max_pct_mt", 100.0) if merged.get("max_pct_mt") is not None else 100.0,
            "max_pct_mt", 0.0, 100.0,
        )
        if step == "filter":
            merged["min_cells"] = max(0, _as_int(merged.get("min_cells", 3), "min_cells"))
    elif step == "normalize":
        merged["target_sum"] = _limited_float(merged.get("target_sum", 10000.0), "target_sum", 1e-9, 1e9)
    elif step == "hvg":
        merged["n_top_genes"] = _as_int(merged.get("n_top_genes", 2000), "n_top_genes")
        if merged["n_top_genes"] < 1:
            raise EngineError("bad_request", "n_top_genes 必须 ≥ 1")
        merged["flavor"] = str(merged.get("flavor") or "seurat")
    elif step == "pca":
        merged["n_pcs"] = _as_int(merged.get("n_pcs", 30), "n_pcs")
        if merged["n_pcs"] < 1:
            raise EngineError("bad_request", "n_pcs 必须 ≥ 1")
    elif step == "neighbors":
        merged["n_neighbors"] = _as_int(merged.get("n_neighbors", 15), "n_neighbors")
        if merged["n_neighbors"] < 1:
            raise EngineError("bad_request", "n_neighbors 必须 ≥ 1")
        merged["metric"] = str(merged.get("metric") or "cosine")
    elif step == "cluster":
        merged["resolution"] = _limited_float(merged.get("resolution", 0.2), "resolution", 1e-4, 100.0)
        merged["method"] = str(merged.get("method") or "louvain").lower()
    elif step == "embed":
        merged["layout"] = str(merged.get("layout") or "tsne").lower()
        merged["perplexity"] = _limited_float(merged.get("perplexity", 30.0), "perplexity", 2.0, 1000.0)
        merged["n_iter"] = max(250, min(2000, _as_int(merged.get("n_iter", 1000), "n_iter")))
    elif step == "markers":
        merged["groupby"] = str(merged.get("groupby") or "clusters")
        merged["top_n"] = _as_int(merged.get("top_n", 10), "top_n")
        merged["top_n"] = max(1, min(200, merged["top_n"]))
        merged["method"] = str(merged.get("method") or "t-test")
    elif step == "annotate":
        mapping = merged.get("mapping") or {}
        if not isinstance(mapping, dict):
            raise EngineError("bad_request", "annotate 的 mapping 必须是对象：{\"0\":\"T cell\"}")
        merged["mapping"] = {str(key): str(value) for key, value in mapping.items()}
    return merged


def _ordered_steps(step_set):
    return [step for step in STEP_ORDER if step in step_set]


# ===========================================================================
# 6. 各步骤计算实现
# ===========================================================================


def _step_qc(dataset, params):
    counts = dataset.full_counts
    n_obs, n_vars = counts.shape
    n_genes = np.diff(counts.indptr).astype(np.int64)
    n_counts = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    mito_mask = np.array([str(gene).upper().startswith("MT-") for gene in dataset.full_var_names], dtype=bool)
    mito_index = np.flatnonzero(mito_mask)
    if mito_index.size:
        mito_counts = np.asarray(counts[:, mito_index].sum(axis=1)).ravel().astype(np.float64)
    else:
        mito_counts = np.zeros(n_obs, dtype=np.float64)
    safe_total = np.maximum(n_counts, 1e-12)
    pct_counts_mt = np.where(n_counts > 0, 100.0 * mito_counts / safe_total, 0.0)

    qc = pd.DataFrame(
        {"n_genes": n_genes, "n_counts": n_counts, "pct_counts_mt": pct_counts_mt},
        index=dataset.full_obs.index,
    )
    dataset.qc_metrics = qc
    dataset.obs = dataset.full_obs.copy()
    for column in qc.columns:
        dataset.obs[column] = qc[column].to_numpy()

    median_genes = float(np.median(n_genes)) if n_obs else 0.0
    median_counts = float(np.median(n_counts)) if n_obs else 0.0
    median_pct = float(np.median(pct_counts_mt)) if n_obs else 0.0
    summary = {
        "title": STEP_TITLES["qc"],
        "bullets": [
            "%d 个细胞 · %d 个基因" % (n_obs, n_vars),
            "每细胞基因数中位数 %s，总计数中位数 %s" % (_fmt_num(median_genes), _fmt_num(median_counts)),
            "线粒体计数占比中位数 %.2f%%（%d 个 MT- 基因；过滤阈值 max_pct_mt=%s）"
            % (median_pct, int(mito_mask.sum()), _fmt_num(params["max_pct_mt"])),
        ],
        "metrics": {
            "n_obs": int(n_obs),
            "n_vars": int(n_vars),
            "n_mito_genes": int(mito_mask.sum()),
            "median_genes": median_genes,
            "median_counts": median_counts,
            "median_pct_mt": median_pct,
        },
    }
    return summary


def _step_filter(dataset, params):
    if dataset.qc_metrics is None:
        raise EngineError("prerequisite", "filter 需要先完成 qc")
    obs_all = dataset.obs_all()
    n_genes = pd.to_numeric(obs_all["n_genes"], errors="coerce").to_numpy(dtype=np.float64)
    n_counts = pd.to_numeric(obs_all["n_counts"], errors="coerce").to_numpy(dtype=np.float64)
    pct_mt = pd.to_numeric(obs_all["pct_counts_mt"], errors="coerce").to_numpy(dtype=np.float64)
    n_genes = np.nan_to_num(n_genes, nan=0.0)
    n_counts = np.nan_to_num(n_counts, nan=0.0)
    pct_mt = np.nan_to_num(pct_mt, nan=100.0)

    keep = n_genes >= params["min_genes"]
    if params.get("max_genes") is not None:
        keep &= n_genes <= params["max_genes"]
    if params.get("min_counts") is not None:
        keep &= n_counts >= params["min_counts"]
    keep &= pct_mt < params["max_pct_mt"]

    n_obs_before = int(keep.size)
    n_obs_after = int(keep.sum())
    if n_obs_after == 0:
        raise EngineError(
            "bad_request",
            "过滤后没有剩余细胞，请放宽阈值（min_genes=%s, max_genes=%s, min_counts=%s, max_pct_mt=%s）"
            % (params["min_genes"], params.get("max_genes"), params.get("min_counts"), params["max_pct_mt"]),
        )
    filtered_counts = dataset.full_counts[keep]
    gene_n_cells = np.asarray((filtered_counts > 0).sum(axis=0)).ravel()
    gene_keep = gene_n_cells >= params["min_cells"]
    n_vars_before = int(gene_keep.size)
    n_vars_after = int(gene_keep.sum())
    if n_vars_after == 0:
        raise EngineError("bad_request", "过滤后没有剩余基因，请降低 min_cells（当前 %s）" % params["min_cells"])

    gene_index = np.flatnonzero(gene_keep)
    keep_index = np.flatnonzero(keep)
    dataset.counts = filtered_counts[:, gene_index].tocsr()
    dataset.X = dataset.full_x[keep][:, gene_index].tocsr()
    dataset.var_names = [dataset.full_var_names[i] for i in gene_index]
    if len(dataset.full_var) == len(dataset.full_var_names):
        var = dataset.full_var.iloc[gene_index].copy()
        var.index = pd.Index(dataset.var_names, name="gene")
        dataset.var = var
    else:
        dataset.var = pd.DataFrame(index=pd.Index(dataset.var_names, name="gene"))
    dataset.obs = obs_all.iloc[keep_index].copy()
    dataset.layers = {"counts": dataset.counts.copy()}
    dataset.filter_info = {
        "n_obs_before": n_obs_before,
        "n_obs_after": n_obs_after,
        "n_vars_before": n_vars_before,
        "n_vars_after": n_vars_after,
        "params": dict(params),
    }

    summary = {
        "title": STEP_TITLES["filter"],
        "bullets": [
            "保留 %d/%d 个细胞（过滤掉 %d 个）" % (n_obs_after, n_obs_before, n_obs_before - n_obs_after),
            "保留 %d/%d 个基因（%d 个基因在 < %d 个细胞中检出）"
            % (n_vars_after, n_vars_before, n_vars_before - n_vars_after, params["min_cells"]),
        ],
        "metrics": {
            "n_obs_before": n_obs_before,
            "n_obs_after": n_obs_after,
            "n_vars_before": n_vars_before,
            "n_vars_after": n_vars_after,
            "n_cells_removed": n_obs_before - n_obs_after,
            "n_genes_removed": n_vars_before - n_vars_after,
            "min_genes": params["min_genes"],
            "max_genes": params.get("max_genes"),
            "min_counts": params.get("min_counts"),
            "max_pct_mt": params["max_pct_mt"],
            "min_cells": params["min_cells"],
        },
    }
    return summary


def _step_normalize(dataset, params):
    # SC-02：没有经过验证的原始 counts 时，不允许把已标准化的表达再当计数做 log-normalize
    if not getattr(dataset, "counts_verified", True):
        raise EngineError(
            "prerequisite",
            "当前数据没有可验证的原始 counts（来源：%s）。log-normalize 需要真实计数矩阵，"
            "请改用带 counts 层的文件重新导入，或先确认导入的不是已标准化的表达矩阵。"
            % (getattr(dataset, "counts_source", "unknown")),
        )
    counts = dataset.counts
    row_sums = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    factor = np.zeros_like(row_sums)
    nonzero = row_sums > 0
    factor[nonzero] = params["target_sum"] / row_sums[nonzero]
    normalized = counts.multiply(factor[:, None]).tocsr().astype(np.float32)
    if normalized.nnz:
        normalized.data = np.log1p(normalized.data).astype(np.float32)
        normalized.sort_indices()
    dataset.X = normalized
    dataset.uns["normalize"] = {"target_sum": float(params["target_sum"]), "log1p": True}
    return {
        "title": STEP_TITLES["normalize"],
        "bullets": [
            "每细胞缩放到目标总计数 %s，再做 log1p" % _fmt_num(params["target_sum"]),
            "标准化矩阵非零元素 %s 个" % _fmt_num(int(normalized.nnz)),
        ],
        "metrics": {
            "target_sum": float(params["target_sum"]),
            "n_nonzero": int(normalized.nnz),
            "n_cells_zero_total": int((~nonzero).sum()),
        },
    }


def _step_hvg(dataset, params):
    if dataset.X is None:
        raise EngineError("prerequisite", "hvg 需要先完成 normalize")
    matrix = dataset.X.tocsr()
    n_obs = matrix.shape[0]
    n_vars = matrix.shape[1]
    if n_vars == 0:
        raise EngineError("bad_request", "没有基因可用于高变基因选择")
    means = np.asarray(matrix.mean(axis=0)).ravel().astype(np.float64)
    squared = matrix.multiply(matrix)
    mean_squares = np.asarray(squared.mean(axis=0)).ravel().astype(np.float64)
    variances = np.maximum(mean_squares - means ** 2, 0.0)
    if n_obs > 1:
        variances = variances * n_obs / (n_obs - 1)
    dispersions = np.where(means > 0, variances / np.maximum(means, 1e-12), 0.0)

    # 按均值分 20 箱，箱内对离散度做 z-score（Seurat 风格）
    n_bins = 20
    with np.errstate(invalid="ignore"):
        binned = pd.qcut(pd.Series(means), n_bins, labels=False, duplicates="drop")
    binned = np.asarray(binned, dtype=np.float64)
    dispersions_norm = np.zeros(n_vars, dtype=np.float64)
    valid_bins = binned[~np.isnan(binned)]
    for bin_id in np.unique(valid_bins):
        mask = binned == bin_id
        values = dispersions[mask]
        mean_value = float(values.mean()) if values.size else 0.0
        std_value = float(values.std(ddof=0)) if values.size else 0.0
        dispersions_norm[mask] = (values - mean_value) / std_value if std_value > 0 else 0.0
    dispersions_norm[np.isnan(binned)] = 0.0

    valid = np.flatnonzero(means > 0)
    if valid.size == 0:
        raise EngineError("bad_request", "所有基因均无表达，无法选择高变基因")
    order = valid[np.argsort(-dispersions_norm[valid], kind="stable")]
    take = int(min(params["n_top_genes"], order.size))
    mask = np.zeros(n_vars, dtype=bool)
    mask[order[:take]] = True

    dataset.var["highly_variable"] = mask
    dataset.var["means"] = means.astype(np.float64)
    dataset.var["dispersions"] = dispersions.astype(np.float64)
    dataset.var["dispersions_norm"] = dispersions_norm
    dataset.uns["hvg"] = {
        "n_top_genes": int(params["n_top_genes"]),
        "flavor": params["flavor"],
        "requested_flavor": params["flavor"],
        "actual_method": "seurat_binned_dispersion(builtin)",
        "fallback_reason": (
            "builtin engine always uses the Seurat-style mean-binned normalised dispersion; "
            "requested flavor %r was not applied" % params["flavor"]
        ) if params["flavor"] != "seurat" else None,
        "n_selected": take,
    }
    bullets = ["筛选出 %d 个高变基因（n_top_genes=%d，flavor=%s）"
               % (take, params["n_top_genes"], params["flavor"])]
    if take < int(params["n_top_genes"]):
        bullets.append("可用基因只有 %d 个（少于 n_top_genes=%d），已取全部" % (n_vars, params["n_top_genes"]))
    if params["flavor"] != "seurat":
        bullets.append("builtin 引擎统一使用 Seurat 风格均值分箱离散度（请求 flavor=%s）" % params["flavor"])
    return {
        "title": STEP_TITLES["hvg"],
        "bullets": bullets,
        "metrics": {
            "n_top_genes": int(params["n_top_genes"]),
            "n_selected": take,
            "n_vars": int(n_vars),
            "flavor": params["flavor"],
        },
    }


def _seurat_engine(params):
    """流程步骤是否走 Seurat 官方引擎（CONTRACT §10.5）。"""
    return str((params or {}).get("engine") or "").lower() in ("seurat", "seurat6", "r")


def _run_seurat_step(dataset, step, params):
    """把一个流程步骤交给 Seurat 6.0 官方实现，并把结果写回 dataset。"""
    params = dict(params or {})
    # pca 步骤选定的整合方式（harmony/rpca/none）要贯穿后续 neighbors/cluster/embed
    if step != "pca":
        uns = getattr(dataset, "uns", {}) or {}
        if not params.get("integration") and uns.get("seurat_integration"):
            params["integration"] = uns["seurat_integration"]
        if not params.get("reduction") and uns.get("seurat_reduction"):
            params["reduction"] = uns["seurat_reduction"]
    if step == "neighbors":
        # Seurat 路径下 FindNeighbors 在 cluster 内部完成，这里只做状态标记
        dataset.adj = None
        dataset.knn = None
        return {
            "title": STEP_TITLES.get(step, step),
            "bullets": ["Seurat 6.0 模式：邻接图由 cluster 步骤内部的 FindNeighbors 完成"],
            "metrics": {
                "engine": "seurat6",
                "actual_method": "Seurat::FindNeighbors",
                "method_provenance": "Seurat 6.0 官方实现（FindNeighbors/FindClusters 内部完成）",
            },
        }
    external = _load_lib_module("sc_external")
    try:
        result = external.seurat_step(dataset, step, params, str(WORKDIR))
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))
    payload = result.get("payload") or {}
    if step == "pca" and result.get("pca") is not None:
        dataset.obsm["X_pca"] = np.asarray(result["pca"], dtype=np.float32)
        dataset.uns["pca_engine"] = "seurat6"
        metrics_pca = (payload.get("summary") or {}).get("metrics") or {}
        dataset.uns["seurat_integration"] = metrics_pca.get("integration") or params.get("integration") or "none"
        dataset.uns["seurat_reduction"] = metrics_pca.get("reduction") or "pca"
        dataset.uns.pop("pca_variance_ratio", None)
    if step == "cluster":
        labels = result.get("clusters")
        if labels is None:
            raise EngineError("internal", "Seurat 6.0 没有返回聚类结果")
        text = [str(item) for item in labels]
        dataset.obs["clusters"] = pd.Series(text, index=dataset.obs.index, dtype=object)
        metrics_cluster = (payload.get("summary") or {}).get("metrics") or {}
        method = str(metrics_cluster.get("integration") or "none").lower()
        if method in ("harmony", "rpca"):
            dataset.obs["clusters_" + method] = pd.Series(text, index=dataset.obs.index, dtype=object)
        if getattr(dataset, "cell_type_default", True) or "cell_type" not in dataset.obs.columns:
            dataset.obs["cell_type"] = pd.Series(text, index=dataset.obs.index, dtype=object)
            dataset.cell_type_default = True
        dataset.markers = None
        dataset.resolution_scan = None
    if step == "embed" and result.get("embed") is not None:
        embedding = np.asarray(result["embed"], dtype=np.float32)
        if embedding.ndim == 2 and embedding.shape[1] >= 2:
            dataset.obsm["X_emb"] = embedding[:, :2]
            dataset.obs["_emb_x"] = embedding[:, 0]
            dataset.obs["_emb_y"] = embedding[:, 1]
            dataset.layout = "umap(seurat6)"
    summary = payload.get("summary") or {}
    bullets = ["%s · %s" % (STEP_TITLES.get(step, step), payload.get("method") or "Seurat 6.0 官方实现")]
    bullets += [str(item) for item in (summary.get("bullets") or [])]
    metrics = dict(summary.get("metrics") or {})
    metrics["engine"] = "seurat6"
    # EXPORT-01：外部实现要留下 actual_method + provenance（Methods 不得再写内置 Louvain）
    seurat_methods = {
        "pca": "Seurat::RunPCA",
        "cluster": "Seurat::FindClusters",
        "embed": "Seurat::RunUMAP",
    }
    actual_method = seurat_methods.get(step, "Seurat 6.0 官方实现")
    if step == "pca":
        integration = str(metrics.get("integration") or "none")
        actual_method += " (integration=%s)" % integration
    elif step == "cluster":
        integration = str(metrics.get("integration") or "none")
        actual_method += " (Louvain algorithm, integration=%s, resolution=%s)" % (
            integration, _fmt_num(metrics.get("resolution")),
        )
    metrics["actual_method"] = actual_method
    metrics["method_provenance"] = payload.get("method") or "Seurat 6.0 官方实现"
    if not metrics.get("method_provenance"):
        metrics["method_provenance"] = "Seurat 6.0 官方实现（版本未记录）"
    return {"title": STEP_TITLES.get(step, step), "bullets": bullets, "metrics": metrics}


def _step_pca(dataset, params):
    if _seurat_engine(params):
        return _run_seurat_step(dataset, "pca", params)
    if "highly_variable" not in dataset.var.columns:
        raise EngineError("prerequisite", "pca 需要先完成 hvg")
    mask = np.asarray(dataset.var["highly_variable"].to_numpy(), dtype=bool)
    hvg_index = np.flatnonzero(mask)
    if hvg_index.size == 0:
        raise EngineError("bad_request", "高变基因列表为空，无法计算 PCA")
    matrix = dataset.X.tocsr()[:, hvg_index]
    n_obs = matrix.shape[0]
    if n_obs < 2:
        raise EngineError("bad_request", "细胞数太少，无法计算 PCA")
    n_pcs = int(min(params["n_pcs"], n_obs - 1, hvg_index.size))
    n_pcs = max(1, n_pcs)
    solver = "randomized"
    note = None
    try:
        if n_obs * hvg_index.size <= 30_000_000:
            dense = np.asarray(matrix.toarray() if sp.issparse(matrix) else matrix, dtype=np.float32)
            pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=SEED)
            embedding = pca.fit_transform(dense)
        else:
            solver = "covariance_eigh"
            note = "矩阵较大（>3000 万元素），PCA 使用 covariance_eigh 求解器以控制内存"
            pca = PCA(n_components=n_pcs, svd_solver="covariance_eigh")
            embedding = pca.fit_transform(matrix)
    except Exception as error:
        raise EngineError("internal", "PCA 计算失败：%s" % error)
    embedding = np.asarray(embedding, dtype=np.float32)
    ratios = np.asarray(pca.explained_variance_ratio_, dtype=np.float64)
    dataset.obsm["X_pca"] = embedding
    dataset.uns["pca_variance_ratio"] = [float(value) for value in ratios]
    dataset.uns["pca_n_pcs"] = int(n_pcs)
    dataset.uns["pca_solver"] = solver
    cumulative = float(np.sum(ratios))
    bullets = [
        "取 %d 个主成分" % n_pcs,
        "PC1 解释方差 %.1f%%" % (float(ratios[0]) * 100.0 if ratios.size else 0.0),
        "前 %d 个主成分累计解释方差 %.1f%%" % (n_pcs, cumulative * 100.0),
    ]
    if note:
        bullets.append(note)
    return {
        "title": STEP_TITLES["pca"],
        "bullets": bullets,
        "metrics": {
            "n_pcs": n_pcs,
            "pc1_variance_ratio": float(ratios[0]) if ratios.size else 0.0,
            "cumulative_variance_ratio": cumulative,
            "n_hvg": int(hvg_index.size),
            "solver": solver,
        },
    }


def _step_neighbors(dataset, params):
    if _seurat_engine(params):
        return _run_seurat_step(dataset, "neighbors", params)
    if "X_pca" not in dataset.obsm:
        raise EngineError("prerequisite", "neighbors 需要先完成 pca")
    representation = np.asarray(dataset.obsm["X_pca"], dtype=np.float32)
    n_obs = representation.shape[0]
    if n_obs < 3:
        raise EngineError("bad_request", "细胞数太少，无法构建 kNN 图")
    k = int(max(1, min(params["n_neighbors"], n_obs - 1)))
    metric = str(params["metric"])
    try:
        model = NearestNeighbors(n_neighbors=k + 1, metric=metric)
        model.fit(representation)
        distances, indices = model.kneighbors(representation)
    except Exception as error:
        raise EngineError("bad_request", "kNN 构建失败（metric=%s）：%s" % (metric, error))
    indices = np.asarray(indices[:, 1:], dtype=np.int32)
    distances = np.asarray(distances[:, 1:], dtype=np.float32)
    rows = np.repeat(np.arange(n_obs, dtype=np.int64), k)
    cols = indices.ravel().astype(np.int64)
    adjacency = sp.csr_matrix(
        (np.ones(rows.size, dtype=np.float32), (rows, cols)), shape=(n_obs, n_obs)
    )
    adjacency = adjacency.maximum(adjacency.T).tocsr()
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    adjacency.sort_indices()
    edge_count = int(adjacency.nnz // 2)
    n_components = int(_csgraph.connected_components(adjacency, directed=False)[0])
    dataset.adj = adjacency
    dataset.knn = {"indices": indices, "distances": distances, "k": k, "metric": metric}
    dataset.uns["neighbors"] = {"n_neighbors": k, "metric": metric, "n_edges": edge_count,
                                "n_components": n_components}
    bullets = [
        "k=%d（metric=%s）的 kNN 图" % (k, metric),
        "邻接图 %s 条边，%d 个连通分量" % (_fmt_num(edge_count), n_components),
    ]
    if n_components > 1:
        bullets.append("图不连通（%d 个分量），聚类会在分量内分别优化" % n_components)
    return {
        "title": STEP_TITLES["neighbors"],
        "bullets": bullets,
        "metrics": {
            "n_neighbors": k,
            "metric": metric,
            "n_edges": edge_count,
            "n_components": n_components,
        },
    }


def _leiden_via_scanpy(dataset, resolution):
    """scanpy + leidenalg 可用时走真 Leiden；返回 int32 标签（按大小重排）。"""
    if SCANPY is None or not HAS_LEIDENALG:
        raise EngineError("unsupported", "scanpy/leidenalg 不可用")
    n_obs = dataset.n_obs
    matrix = ad.AnnData(X=sp.csr_matrix((n_obs, 1), dtype=np.float32))
    matrix.obs_names = [str(name) for name in dataset.obs.index]
    matrix.obsm["X_pca"] = np.asarray(dataset.obsm["X_pca"], dtype=np.float32)
    matrix.obsp["connectivities"] = dataset.adj
    matrix.uns["neighbors"] = {"connectivities_key": "connectivities", "params": {"method": "umap"}}
    SCANPY.tl.leiden(matrix, resolution=float(resolution), random_state=SEED, key_added="leiden", flavor="leidenalg")
    labels = matrix.obs["leiden"].astype(str).to_numpy()
    codes, _uniques = pd.factorize(labels)
    return _relabel_by_size(codes)


def _step_cluster(dataset, params):
    if _seurat_engine(params):
        return _run_seurat_step(dataset, "cluster", params)
    if dataset.adj is None:
        raise EngineError("prerequisite", "cluster 需要先完成 neighbors")
    resolution = float(params["resolution"])
    requested_method = str(params["method"]).lower()
    method_effective = "louvain"
    note = None
    fallback_reason = None
    labels = None
    if requested_method == "leiden":
        try:
            labels = _leiden_via_scanpy(dataset, resolution)
            method_effective = "leiden(scanpy)"
        except Exception:
            labels = None
            method_effective = "louvain(builtin)"
            note = "leidenalg/scanpy 不可用，已回退内置纯 Python Louvain（resolution 语义一致）"
            fallback_reason = (
                "scanpy/leidenalg unavailable; requested method Leiden fell back to builtin Louvain "
                "(same resolution semantics, but a different algorithm -> results are not interchangeable)"
            )
    if labels is None:
        labels = louvain_clusters(dataset.adj, resolution)
        method_effective = "louvain(builtin)" if requested_method == "louvain" else method_effective
    labels = np.asarray(labels, dtype=np.int32)
    n_clusters = int(labels.max() + 1) if labels.size else 0
    cluster_ids = np.array([str(value) for value in labels], dtype=object)
    dataset.obs["clusters"] = pd.Series(cluster_ids, index=dataset.obs.index, dtype=object)
    if getattr(dataset, "cell_type_default", True) or "cell_type" not in dataset.obs.columns:
        dataset.obs["cell_type"] = pd.Series(cluster_ids, index=dataset.obs.index, dtype=object)
        dataset.cell_type_default = True
    dataset.markers = None
    dataset.resolution_scan = None
    if n_clusters > 30:
        dataset.warnings.append("聚类数 %d 偏多（参考 20–30）：建议把 resolution 调低（当前 %.2f）" % (n_clusters, resolution))
    elif n_clusters < 10:
        dataset.warnings.append("聚类数 %d 偏少（参考 20–30）：建议把 resolution 调高（当前 %.2f）" % (n_clusters, resolution))
    sizes = np.bincount(labels, minlength=n_clusters).astype(int).tolist() if labels.size else []
    dataset.uns["cluster"] = {
        "method": method_effective,
        "requested_method": requested_method,
        "fallback_reason": fallback_reason,
        "resolution": resolution,
        "seed": SEED,
        "input_scale": _matrix_scale(dataset),
        "n_clusters": n_clusters,
    }
    size_text = "，".join("%d:%d" % (index, size) for index, size in enumerate(sizes[:10]))
    bullets = ["resolution=%s → %d 个簇" % (_fmt_num(resolution), n_clusters),
               "簇大小（前 %d）：%s" % (min(10, n_clusters), size_text)]
    if note:
        bullets.append(note)
    return {
        "title": STEP_TITLES["cluster"],
        "bullets": bullets,
        "metrics": {
            "n_clusters": n_clusters,
            "resolution": resolution,
            "method": method_effective,
            "cluster_sizes": {str(index): int(size) for index, size in enumerate(sizes)},
        },
    }


def _tsne_embedding(representation, perplexity, n_iter):
    n_obs = representation.shape[0]
    effective_perplexity = float(min(max(perplexity, 2.0), max(2.0, (n_obs - 1) / 3.0), 50.0))
    effective_perplexity = min(effective_perplexity, max(2.0, n_obs - 1.0))
    model = TSNE(
        n_components=2,
        perplexity=effective_perplexity,
        init="pca",
        random_state=SEED,
        max_iter=int(n_iter),
    )
    embedding = model.fit_transform(np.asarray(representation, dtype=np.float32))
    return np.asarray(embedding, dtype=np.float32), effective_perplexity


def _umap_embedding(representation):
    """UMAP：优先 umap-learn，其次 scanpy；都不可用抛 unsupported。"""
    if HAS_UMAP:
        model = _umap.UMAP(n_components=2, random_state=SEED)
        return np.asarray(model.fit_transform(np.asarray(representation, dtype=np.float32)), dtype=np.float32)
    if SCANPY is not None:
        matrix = ad.AnnData(X=np.asarray(representation, dtype=np.float32))
        SCANPY.pp.neighbors(matrix, n_neighbors=15, use_rep="X")
        SCANPY.tl.umap(matrix, random_state=SEED)
        return np.asarray(matrix.obsm["X_umap"], dtype=np.float32)
    raise EngineError("unsupported", "未安装 umap-learn / scanpy")


def _step_embed(dataset, params):
    if _seurat_engine(params):
        return _run_seurat_step(dataset, "embed", params)
    if "X_pca" not in dataset.obsm:
        raise EngineError("prerequisite", "embed 需要先完成 pca")
    representation = np.asarray(dataset.obsm["X_pca"], dtype=np.float32)
    n_obs = representation.shape[0]
    requested = str(params["layout"]).lower()
    layout = requested
    notes = []
    n_iter = int(params["n_iter"])
    if n_obs > 15000:
        n_iter = min(n_iter, 750)
        notes.append("细胞数 > 15000，t-SNE 迭代降到 %d，可能仍较慢" % n_iter)
    effective_perplexity = None
    embedding = None
    if layout == "umap":
        try:
            embedding = _umap_embedding(representation)
        except Exception:
            layout = "tsne"
            notes.append("UMAP 不可用（未安装 umap-learn/scanpy），已回退 t-SNE")
    if layout == "tsne":
        if n_obs < 5:
            layout = "pca"
            notes.append("细胞数 < 5，无法运行 t-SNE，已回退 PCA 前两维")
        else:
            try:
                embedding, effective_perplexity = _tsne_embedding(
                    representation, float(params["perplexity"]), n_iter
                )
            except Exception as error:
                layout = "pca"
                notes.append("t-SNE 失败（%s），已回退 PCA 前两维" % error)
    if embedding is None or layout == "pca":
        embedding = np.zeros((n_obs, 2), dtype=np.float32)
        dims = int(min(2, representation.shape[1]))
        embedding[:, :dims] = representation[:, :dims]
    embedding = np.nan_to_num(np.asarray(embedding, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    dataset.obsm["X_emb"] = embedding
    dataset.obs["_emb_x"] = embedding[:, 0]
    dataset.obs["_emb_y"] = embedding[:, 1]
    dataset.layout = layout
    fallback_reason = None
    if layout != requested:
        fallback_reason = "; ".join(str(item) for item in notes) or (
            "requested layout %r unavailable; fell back to %r" % (requested, layout)
        )
    dataset.uns["embed"] = {
        "layout": layout,
        "requested_layout": requested,
        "fallback_reason": fallback_reason,
        "perplexity": effective_perplexity,
        "n_iter": n_iter,
    }
    bullets = ["%s 二维嵌入完成（%d 个细胞）" % (layout.upper(), n_obs)]
    if layout == "tsne":
        bullets.append("perplexity=%s，max_iter=%d" % (_fmt_num(effective_perplexity), n_iter))
    elif layout == "pca":
        bullets.append("使用 PCA 前两维作为二维嵌入")
    bullets.extend(notes)
    return {
        "title": STEP_TITLES["embed"],
        "bullets": bullets,
        "metrics": {
            "layout": layout,
            "n_points": int(n_obs),
            "perplexity": effective_perplexity,
            "n_iter": n_iter if layout == "tsne" else None,
        },
    }


def _group_labels_and_order(dataset, groupby):
    """返回 (labels 数组, 有序分组 id 列表, {id: 掩码})。"""
    if groupby not in dataset.obs.columns:
        raise EngineError("not_found", "obs 列不存在：%s（可用：%s）"
                          % (groupby, ", ".join(str(c) for c in dataset.obs.columns)))
    values = [str(value) for value in dataset.obs[groupby].tolist()]
    if groupby == "clusters":
        def sort_key(item):
            try:
                return (0, float(item))
            except Exception:
                return (1, str(item))

        categories = sorted(set(values), key=sort_key)
    else:
        counts = pd.Series(values).value_counts()
        categories = [str(label) for label in counts.index.tolist()]
    masks = {label: np.array([value == label for value in values], dtype=bool) for label in categories}
    return values, categories, masks


def _compute_markers(dataset, groupby, top_n, method):
    """one-vs-rest Welch t 检验 + BH FDR + pct_in/pct_out，返回 groups（已缓存）。"""
    if dataset.X is None:
        raise EngineError("prerequisite", "markers 需要先完成 normalize")
    method = str(method or "t-test")
    if method not in ("t-test", "t_test", "welch", "welch-t"):
        # 其它统计方法不在 builtin 引擎能力内：退化为 Welch t 并在结果里标注
        method_effective = "t-test(builtin)"
    else:
        method_effective = "t-test(builtin)"
    values, categories, masks = _group_labels_and_order(dataset, groupby)
    matrix = dataset.X.tocsr()
    n_obs = matrix.shape[0]
    n_vars = matrix.shape[1]
    gene_names = dataset.var_names

    total_sum = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    total_sq = np.asarray(matrix.multiply(matrix).sum(axis=0)).ravel().astype(np.float64)
    total_nnz = np.asarray(matrix.getnnz(axis=0)).ravel().astype(np.int64)

    groups = []
    for label in categories:
        mask = masks[label]
        n_group = int(mask.sum())
        n_rest = n_obs - n_group
        sub = matrix[mask]
        sum_group = np.asarray(sub.sum(axis=0)).ravel().astype(np.float64)
        sq_group = np.asarray(sub.multiply(sub).sum(axis=0)).ravel().astype(np.float64)
        nnz_group = np.asarray(sub.getnnz(axis=0)).ravel().astype(np.int64)
        sum_rest = total_sum - sum_group
        sq_rest = np.maximum(total_sq - sq_group, 0.0)
        nnz_rest = total_nnz - nnz_group

        mean_group = sum_group / max(n_group, 1)
        mean_rest = sum_rest / max(n_rest, 1)
        var_group = np.maximum(sq_group - sum_group ** 2 / max(n_group, 1), 0.0) / max(n_group - 1, 1)
        var_rest = np.maximum(sq_rest - sum_rest ** 2 / max(n_rest, 1), 0.0) / max(n_rest - 1, 1)

        # Welch t 统计量
        numerator = mean_group - mean_rest
        if n_group >= 2 and n_rest >= 2:
            denom = np.sqrt(var_group / max(n_group, 1) + var_rest / max(n_rest, 1))
            with np.errstate(divide="ignore", invalid="ignore"):
                t_stat = np.where(denom > 0, numerator / np.maximum(denom, 1e-30), 0.0)
                t_stat = np.nan_to_num(t_stat, nan=0.0, posinf=0.0, neginf=0.0)
                a = var_group / max(n_group, 1)
                b = var_rest / max(n_rest, 1)
                df_denom = a ** 2 / max(n_group - 1, 1) + b ** 2 / max(n_rest - 1, 1)
                with np.errstate(divide="ignore", invalid="ignore"):
                    df = np.where(df_denom > 0, (a + b) ** 2 / np.maximum(df_denom, 1e-30), 1.0)
                df = np.nan_to_num(df, nan=1.0, posinf=1.0, neginf=1.0)
            pvalues = 2.0 * _stats.t.sf(np.abs(t_stat), np.maximum(df, 1.0))
            pvalues = np.nan_to_num(pvalues, nan=1.0)
        else:
            t_stat = np.zeros(n_vars, dtype=np.float64)
            pvalues = np.ones(n_vars, dtype=np.float64)
        pvalues_adj = _bh_fdr(pvalues)

        # logFC：基于 expm1 的均值（log 数据 → 原始尺度 → log2 ratio）
        if sub.nnz:
            exp_sum = np.asarray(sub.expm1().sum(axis=0)).ravel().astype(np.float64)
        else:
            exp_sum = np.zeros(n_vars, dtype=np.float64)
        mean_group_raw = exp_sum / max(n_group, 1)
        # 其余细胞的 expm1 均值 = 全矩阵 expm1 列和 - 组内列和
        rest_exp_sum = np.maximum(total_expm1_matrix(dataset) - exp_sum, 0.0)
        mean_rest_raw = rest_exp_sum / max(n_rest, 1)
        logfc = np.log2((mean_group_raw + 1e-9) / (mean_rest_raw + 1e-9))
        logfc = np.nan_to_num(logfc, nan=0.0, posinf=0.0, neginf=0.0)

        pct_in = nnz_group.astype(np.float64) / max(n_group, 1)
        pct_out = nnz_rest.astype(np.float64) / max(n_rest, 1)

        order = np.argsort(-t_stat, kind="stable")[: int(min(top_n, n_vars))]
        genes = []
        for gene_index in order:
            genes.append(
                {
                    "gene": str(gene_names[int(gene_index)]),
                    "score": float(t_stat[gene_index]),
                    "logfc": float(logfc[gene_index]),
                    "pval_adj": float(pvalues_adj[gene_index]),
                    "pct_in": float(pct_in[gene_index]),
                    "pct_out": float(pct_out[gene_index]),
                }
            )
        significant = int(np.sum((pvalues_adj < 0.05) & (logfc > 0)))
        groups.append(
            {
                "id": str(label),
                "label": str(label),
                "n_cells": n_group,
                "n_significant": significant,
                "genes": genes,
            }
        )
    dataset.markers = {
        "groupby": str(groupby),
        "top_n": int(top_n),
        "method": method_effective,
        "groups": groups,
    }
    return groups


_EXPM1_CACHE = {}


def total_expm1_matrix(dataset):
    """当前矩阵逐元素 expm1 的列和（按数据集 id + 版本缓存）。"""
    key = (dataset.id, dataset.version, dataset.X.shape, dataset.X.nnz)
    cached = _EXPM1_CACHE.get(key)
    if cached is not None:
        return cached
    matrix = dataset.X.tocsr()
    if matrix.nnz:
        values = np.asarray(matrix.expm1().sum(axis=0)).ravel().astype(np.float64)
    else:
        values = np.zeros(matrix.shape[1], dtype=np.float64)
    _EXPM1_CACHE.clear()  # 只保留最新版本，避免内存累积
    _EXPM1_CACHE[key] = values
    return values


def _step_markers(dataset, params):
    if "clusters" not in dataset.obs.columns:
        raise EngineError("prerequisite", "markers 需要先完成 cluster")
    requested_method = str(params.get("method") or "t-test")
    groups = _compute_markers(dataset, params["groupby"], params["top_n"], params["method"])
    n_markers = int(sum(len(group["genes"]) for group in groups))
    n_significant = int(sum(group.get("n_significant", 0) for group in groups))
    fallback_reason = None
    if requested_method.lower() not in ("t-test", "t_test", "welch", "welch-t"):
        fallback_reason = (
            "requested marker test %r is not implemented by the builtin engine; "
            "fell back to Welch t-test" % requested_method
        )
    dataset.uns["markers"] = {
        "groupby": params["groupby"],
        "top_n": int(params["top_n"]),
        "requested_method": requested_method,
        "actual_method": "Welch t-test(builtin, one-vs-rest)+Benjamini-Hochberg",
        "fallback_reason": fallback_reason,
        "n_groups": len(groups),
    }
    return {
        "title": STEP_TITLES["markers"],
        "bullets": [
            "%d 个分组 × top %d = %d 个标记基因" % (len(groups), params["top_n"], n_markers),
            "显著正向标记（adj.P<0.05 且 logFC>0）：%d 个" % n_significant,
        ],
        "metrics": {
            "n_markers": n_markers,
            "n_groups": len(groups),
            "n_significant": n_significant,
            "groupby": params["groupby"],
            "method": "t-test(builtin)",
        },
    }


def _step_annotate(dataset, params):
    if "clusters" not in dataset.obs.columns:
        raise EngineError("prerequisite", "annotate 需要先完成 cluster")
    mapping = params["mapping"]
    clusters = [str(value) for value in dataset.obs["clusters"].tolist()]
    labels = [mapping.get(cluster, cluster) for cluster in clusters]
    dataset.obs["cell_type"] = pd.Series(labels, index=dataset.obs.index, dtype=object)
    dataset.cell_type_default = False
    unique_labels = sorted(set(labels))
    bullets = [
        "已为 %d 个簇写入细胞类型（共 %d 种类型）" % (len(mapping), len(unique_labels)),
        "类型：%s" % "、".join(unique_labels[:12]),
    ]
    if not mapping:
        bullets.append("mapping 为空：cell_type 暂时保持簇号，可在 UI / Agent 里补充映射后重跑")
    return {
        "title": STEP_TITLES["annotate"],
        "bullets": bullets,
        "metrics": {
            "n_types": len(unique_labels),
            "n_mapped": len(mapping),
            "cell_types": unique_labels,
        },
    }


STEP_HANDLERS = {
    "qc": _step_qc,
    "filter": _step_filter,
    "normalize": _step_normalize,
    "hvg": _step_hvg,
    "pca": _step_pca,
    "neighbors": _step_neighbors,
    "cluster": _step_cluster,
    "embed": _step_embed,
    "markers": _step_markers,
    "annotate": _step_annotate,
}


# ---------------------------------------------------------------------------
# 状态清理：重跑某一步时清空它自己 + 所有下游状态
# ---------------------------------------------------------------------------


def _clear_step_state(dataset, step):
    if step == "filter":
        dataset.obs = dataset.obs_all()
        dataset.counts = dataset.full_counts.copy()
        dataset.X = dataset.full_x.copy()
        dataset.var_names = list(dataset.full_var_names)
        dataset.var = dataset.full_var.copy()
        dataset.layers = {"counts": dataset.counts.copy()}
        dataset.filter_info = None
    elif step == "normalize":
        dataset.X = dataset.counts.copy()
        dataset.uns.pop("normalize", None)
    elif step == "hvg":
        for column in ("highly_variable", "means", "dispersions", "dispersions_norm"):
            if column in dataset.var.columns:
                dataset.var.drop(columns=column, inplace=True)
        dataset.uns.pop("hvg", None)
    elif step == "pca":
        dataset.obsm.pop("X_pca", None)
        for key in ("pca_variance_ratio", "pca_n_pcs", "pca_solver"):
            dataset.uns.pop(key, None)
    elif step == "neighbors":
        dataset.adj = None
        dataset.knn = None
        dataset.uns.pop("neighbors", None)
    elif step == "cluster":
        columns = ["clusters"]
        if getattr(dataset, "cell_type_default", True):
            columns.append("cell_type")
        for column in columns:
            if column in dataset.obs.columns:
                dataset.obs.drop(columns=column, inplace=True)
        dataset.markers = None
        dataset.resolution_scan = None
        dataset.uns.pop("cluster", None)
        dataset.uns.pop("markers", None)
    elif step == "embed":
        dataset.obsm.pop("X_emb", None)
        for column in ("_emb_x", "_emb_y"):
            if column in dataset.obs.columns:
                dataset.obs.drop(columns=column, inplace=True)
        dataset.layout = None
        dataset.uns.pop("embed", None)
    elif step == "markers":
        dataset.markers = None
        dataset.uns.pop("markers", None)
    elif step == "annotate":
        if "clusters" in dataset.obs.columns and getattr(dataset, "cell_type_default", True):
            dataset.obs["cell_type"] = pd.Series(
                [str(value) for value in dataset.obs["clusters"].tolist()],
                index=dataset.obs.index,
                dtype=object,
            )
    elif step == "qc":
        dataset.obs = dataset.full_obs.copy()
        dataset.qc_metrics = None
        dataset.counts = dataset.full_counts.copy()
        dataset.X = dataset.full_x.copy()
        dataset.var_names = list(dataset.full_var_names)
        dataset.var = dataset.full_var.copy()
        dataset.layers = {"counts": dataset.counts.copy()}
        dataset.filter_info = None


def _downstream_steps(step):
    """按依赖图（而不是线性顺序）找出所有下游步骤。"""
    downstream = set()
    changed = True
    while changed:
        changed = False
        for candidate in STEP_ORDER:
            if candidate == step or candidate in downstream:
                continue
            dependencies = STEP_DEPS.get(candidate) or []
            if any(dependency == step or dependency in downstream for dependency in dependencies):
                downstream.add(candidate)
                changed = True
    return downstream


def _reset_from_step(dataset, step):
    """清空 step 及其所有下游状态（按依赖图；顺序从下游往上游清理）。"""
    targets = {step} | _downstream_steps(step)
    removed = [candidate for candidate in dataset.steps if candidate in targets]
    dataset.steps = [
        candidate for candidate in STEP_ORDER
        if candidate in dataset.steps and candidate not in targets
    ]
    for candidate in sorted(removed, key=lambda item: STEP_ORDER.index(item), reverse=True):
        _clear_step_state(dataset, candidate)
    # EXPORT-01：旧步骤的执行记录必须随状态一起失效，Methods/导出不得描述旧版本
    _prune_execution_records(dataset)


def _summary_line(summary):
    bullets = summary.get("bullets") if isinstance(summary, dict) else None
    if isinstance(bullets, list) and bullets:
        return "；".join(str(item) for item in bullets[:2])[:300]
    if isinstance(summary, dict) and isinstance(summary.get("title"), str):
        return summary["title"]
    return ""


def _engine_info():
    return {
        "mode": ENGINE_MODE,
        "python": PY_VERSION,
        "seed": SEED,
        "versions": dict(PACKAGE_VERSIONS),
    }


# ===========================================================================
# 6b. EXPORT-01 / EXPORT-02：结构化执行记录 + 可重放清单
# ---------------------------------------------------------------------------
# Methods、报告、protocol 与导出脚本都只从 dataset.execution_record 生成；记录
# 在每一步真实成功执行后写入。取不到来源时明确写「未记录」，不凭 UI 名称补写。
# ===========================================================================

_OPTIONAL_VERSIONS = {
    "umap": PACKAGE_VERSIONS.get("umap") if HAS_UMAP else None,
    "leidenalg": _pkg_version("leidenalg") if HAS_LEIDENALG else None,
    "scanpy": PACKAGE_VERSIONS.get("scanpy"),
}

_SHA256_CACHE = {}

REPLAY_MANIFEST_VERSION = 1

# 重放时可以安全重新调用的命令（只覆盖数据/元数据操作，不含外部工具）
REPLAYABLE_COMMANDS = (
    "samples_update",
    "samples_auto_group",
    "obs_set_group",
    "obs_derive",
    "obs_join",
    "samples_remove",
    "markers",
    "analysis",
)

REPLAY_OPEN_ARG_KEYS = (
    "path", "url", "demo", "paths", "merge_dir", "samples",
    "format", "member", "object", "transpose", "name", "delimiter", "gse",
)

TRACKED_OBS_COLUMNS = (
    "sample", "sample_name", "group", "cell_type", "clusters", "batch",
    "donor", "donor_id", "sample_id", "orig.ident", "seurat_clusters", "leiden",
)


def _sha256_file(path):
    """文件校验值（带 (path,size,mtime) 缓存，避免同一导出重复读大文件）。"""
    target = Path(str(path))
    try:
        stat = target.stat()
    except Exception:  # noqa: BLE001
        return None
    key = (str(target), int(stat.st_size), int(stat.st_mtime_ns))
    cached = _SHA256_CACHE.get(key)
    if cached:
        return cached
    digest = hashlib.sha256()
    try:
        with open(target, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except Exception:  # noqa: BLE001
        return None
    value = digest.hexdigest()
    if len(_SHA256_CACHE) > 64:
        _SHA256_CACHE.clear()
    _SHA256_CACHE[key] = value
    return value


def _sha256_json(value):
    text = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _matrix_fingerprint(matrix):
    """确定性矩阵指纹：形状 + 稀疏结构 + float32 值，另附容差比较用的统计量。"""
    if matrix is None:
        return None
    m = matrix if sp.issparse(matrix) else sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    m = m.tocsr(copy=True)
    m.sum_duplicates()
    m.sort_indices()
    data = np.ascontiguousarray(m.data, dtype="<f4")
    digest = hashlib.sha256()
    digest.update(("%d,%d,%d" % (m.shape[0], m.shape[1], int(m.indptr.size))).encode("utf-8"))
    digest.update(np.ascontiguousarray(m.indptr, dtype="<i8").tobytes())
    digest.update(np.ascontiguousarray(m.indices, dtype="<i8").tobytes())
    digest.update(data.tobytes())
    return {
        "shape": [int(m.shape[0]), int(m.shape[1])],
        "nnz": int(m.nnz),
        "sha256": digest.hexdigest(),
        "sum": float(np.sum(data, dtype=np.float64)) if data.size else 0.0,
        "min": float(np.min(data)) if data.size else 0.0,
        "max": float(np.max(data)) if data.size else 0.0,
    }


def _obs_digest(series):
    return _sha256_json([str(value) for value in series.tolist()])


def _record_software():
    """执行记录里的软件/依赖版本（取不到的值保留 None，Methods 写「未记录」）。"""
    return {
        "dsh-singlecell": PLUGIN_VERSION,
        "python": PY_VERSION,
        "engine_mode": ENGINE_MODE,
        "numpy": PACKAGE_VERSIONS.get("numpy"),
        "scipy": PACKAGE_VERSIONS.get("scipy"),
        "pandas": PACKAGE_VERSIONS.get("pandas"),
        "sklearn": PACKAGE_VERSIONS.get("sklearn"),
        "anndata": PACKAGE_VERSIONS.get("anndata"),
        "matplotlib": PACKAGE_VERSIONS.get("matplotlib"),
        "scanpy": PACKAGE_VERSIONS.get("scanpy"),
        "umap": _OPTIONAL_VERSIONS.get("umap"),
        "leidenalg": _OPTIONAL_VERSIONS.get("leidenalg"),
        "seurat6": None,
    }


def _software_fingerprint():
    """决定「逐值比较」还是「容差比较」的环境指纹。"""
    return {
        "plugin_version": PLUGIN_VERSION,
        "engine_mode": ENGINE_MODE,
        "python": "%d.%d" % (sys.version_info[0], sys.version_info[1]),
        "numpy": PACKAGE_VERSIONS.get("numpy"),
        "scipy": PACKAGE_VERSIONS.get("scipy"),
        "pandas": PACKAGE_VERSIONS.get("pandas"),
        "sklearn": PACKAGE_VERSIONS.get("sklearn"),
        "anndata": PACKAGE_VERSIONS.get("anndata"),
        "scanpy": PACKAGE_VERSIONS.get("scanpy"),
        "umap": _OPTIONAL_VERSIONS.get("umap"),
        "leidenalg": _OPTIONAL_VERSIONS.get("leidenalg"),
        "seed": SEED,
    }


def _matrix_scale(dataset):
    """当前 X 的实际尺度（只认已执行的转换记录，不靠步骤名前缀猜测）。"""
    normalize = (getattr(dataset, "uns", {}) or {}).get("normalize") or {}
    if "normalize" in getattr(dataset, "steps", []) and normalize.get("log1p"):
        return "lognorm(log1p, target_sum=%s)" % _fmt_num(normalize.get("target_sum"))
    import_state = (getattr(dataset, "uns", {}) or {}).get("import_state") or {}
    if (import_state.get("details") or {}).get("x_lognormalized"):
        return "imported_lognorm"
    kind = str(getattr(dataset, "matrix_kind", None) or "unknown")
    if kind in ("counts", "integer"):
        return "counts"
    return kind


def _tracked_obs_columns(dataset):
    ordered = [column for column in TRACKED_OBS_COLUMNS if column in dataset.obs.columns]
    extra = set()
    for entry in getattr(dataset, "execution_record", []) or []:
        for column in entry.get("obs_columns") or []:
            extra.add(str(column))
    for column in sorted(extra):
        if column in dataset.obs.columns and column not in ordered:
            ordered.append(column)
    return ordered


def _record_step_execution(dataset, step, user_params, params, summary, before_scale):
    """写入一条结构化执行记录（只在步骤真正成功后调用）。"""
    details = _step_method_details(dataset, step, params, summary)
    record = {
        "seq": len(dataset.execution_record) + 1,
        "at": _now_iso(),
        "kind": "step",
        "step": str(step),
        "requested_method": details.get("requested_method"),
        "actual_method": details.get("actual_method"),
        "method_provenance": details.get("method_provenance"),
        "fallback_reason": details.get("fallback_reason"),
        "params_requested": _jsonable(user_params or {}),
        "params_effective": _jsonable(params or {}),
        # 重放用参数：显式锁定「实际生效」的算法分支（环境差异导致回退时仍然复现实际方法）
        "replay_params": _jsonable(_replay_params_for_step(dataset, step, params)),
        "input_scale": before_scale,
        "output_scale": _matrix_scale(dataset),
        "seed": SEED,
        "software": _record_software(),
        "dataset_version": int(dataset.version),
        "dataset_revision": int(dataset.revision),
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
        "n_obs_after": int(dataset.n_obs),
        "n_vars_after": int(dataset.n_vars),
        "obs_columns": [],
    }
    dataset.execution_record.append(record)
    return record


def _replay_params_for_step(dataset, step, params):
    """把实际生效的算法分支写进重放参数：Leiden 回退就是 Louvain，UMAP 回退就是 t-SNE。"""
    replay = dict(params or {})
    replay.pop("engine", None)  # 外部 Seurat 步骤标为不可重放，重放只走内置实现
    if step == "cluster":
        replay["method"] = "leiden" if str((dataset.uns.get("cluster") or {}).get("method") or "").lower().startswith("leiden") else "louvain"
    elif step == "hvg":
        # builtin 实现始终是 Seurat 风格均值分箱离散度
        replay["flavor"] = "seurat"
    elif step == "embed":
        layout = str((dataset.uns.get("embed") or {}).get("layout") or dataset.layout or replay.get("layout") or "tsne")
        replay["layout"] = layout.split("(")[0]
    elif step == "markers":
        replay["method"] = "t-test"
    return replay


def _step_method_details(dataset, step, params, summary):
    """从真实执行结果里取出 requested / actual / provenance / fallback。"""
    params = params or {}
    metrics = (summary or {}).get("metrics") or {}
    engine_param = str(params.get("engine") or "").lower()
    seurat = engine_param in ("seurat", "seurat6", "r")
    requested = None
    actual = None
    provenance = metrics.get("method_provenance")
    fallback = None
    if step == "qc":
        requested = "builtin"
        actual = "qc_metrics(builtin)"
    elif step == "filter":
        requested = "threshold_filter"
        actual = "threshold_filter(builtin)"
    elif step == "normalize":
        requested = "log1p"
        actual = "library_size_scaling+log1p(builtin, target_sum=%s)" % _fmt_num(params.get("target_sum"))
    elif step == "hvg":
        requested = str(params.get("flavor") or "seurat")
        hvg_meta = dataset.uns.get("hvg") or {}
        actual = "seurat_binned_dispersion(builtin)"
        fallback = hvg_meta.get("fallback_reason")
    elif step == "pca":
        if seurat or str(dataset.uns.get("pca_engine") or "") == "seurat6":
            requested = "seurat6"
            actual = metrics.get("actual_method") or "Seurat::RunPCA"
        else:
            requested = "builtin"
            solver = dataset.uns.get("pca_solver") or metrics.get("solver") or "randomized"
            actual = "PCA(scikit-learn, svd_solver=%s)" % solver
    elif step == "neighbors":
        if seurat:
            requested = "seurat6"
            actual = metrics.get("actual_method") or "Seurat::FindNeighbors"
        else:
            requested = "builtin"
            neighbor_meta = dataset.uns.get("neighbors") or {}
            actual = "knn_graph(scikit-learn NearestNeighbors, metric=%s)" % neighbor_meta.get("metric", params.get("metric"))
    elif step == "cluster":
        cluster_meta = dataset.uns.get("cluster") or {}
        if seurat:
            requested = "seurat6"
            actual = metrics.get("actual_method") or "Seurat::FindClusters"
        else:
            requested = str(params.get("method") or "louvain").lower()
            actual = str(cluster_meta.get("method") or "louvain(builtin)")
            if not actual.startswith("louvain") and not actual.startswith("leiden"):
                actual = "louvain(builtin)"
        fallback = cluster_meta.get("fallback_reason")
    elif step == "embed":
        embed_meta = dataset.uns.get("embed") or {}
        if seurat:
            requested = "seurat6"
            actual = metrics.get("actual_method") or "Seurat::RunUMAP"
        else:
            requested = str(params.get("layout") or "tsne").lower()
            layout = str(embed_meta.get("layout") or dataset.layout or requested)
            if layout == "tsne":
                actual = "t-SNE(scikit-learn, random_state=%d)" % SEED
            elif layout == "umap":
                actual = "UMAP(umap-learn, random_state=%d)" % SEED
            elif layout.startswith("umap"):
                actual = "Seurat::RunUMAP"
            else:
                actual = "PCA(first two components)"
        fallback = embed_meta.get("fallback_reason")
    elif step == "markers":
        method = str(params.get("method") or "t-test")
        requested = method
        actual = "Welch t-test(builtin, one-vs-rest)+Benjamini-Hochberg"
        if method.lower() not in ("t-test", "t_test", "welch", "welch-t"):
            fallback = "requested marker test %r is not implemented by the builtin engine; fell back to Welch t-test" % method
    elif step == "annotate":
        requested = "manual"
        actual = "manual_annotation"
    return {
        "requested_method": requested,
        "actual_method": actual,
        "method_provenance": provenance if provenance else None,
        "fallback_reason": fallback if fallback else metrics.get("fallback_reason"),
    }


def execution_records(dataset):
    """当前数据版本仍有效的执行记录（按时间顺序；已失效步骤的记录被剪掉）。"""
    current_steps = set(getattr(dataset, "steps", []) or [])
    return [
        entry for entry in (getattr(dataset, "execution_record", []) or [])
        if entry.get("kind") != "step" or entry.get("step") in current_steps
    ]


def _record_operation(dataset, kind, command, args, actual_method, *,
                      method_provenance=None, fallback_reason=None, replayable=True,
                      not_replayable_reason=None, obs_columns=None, extra=None):
    """记录一次手动/外部操作（分组、注释、样本移除、外部分析等）。"""
    entry = {
        "seq": len(dataset.execution_record) + 1,
        "at": _now_iso(),
        "kind": str(kind),
        "command": str(command),
        "args": _jsonable(args or {}),
        "requested_method": None,
        "actual_method": actual_method,
        "method_provenance": method_provenance,
        "fallback_reason": fallback_reason,
        "replayable": bool(replayable),
        "not_replayable_reason": not_replayable_reason,
        "obs_columns": [str(column) for column in (obs_columns or [])],
        "seed": SEED,
        "software": _record_software(),
        "dataset_version": int(dataset.version),
        "dataset_revision": int(dataset.revision),
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
        "n_obs_after": int(dataset.n_obs),
        "n_vars_after": int(dataset.n_vars),
    }
    if extra:
        entry.update(_jsonable(extra))
    dataset.execution_record.append(entry)
    return entry


def _prune_execution_records(dataset):
    """数据版本变化后清掉已失效步骤的记录（Methods 不得描述旧版本步骤）。"""
    current_steps = set(getattr(dataset, "steps", []) or [])
    dataset.execution_record = [
        entry for entry in (getattr(dataset, "execution_record", []) or [])
        if entry.get("kind") != "step" or entry.get("step") in current_steps
    ]
    for index, entry in enumerate(dataset.execution_record):
        entry["seq"] = index + 1


def _new_obs_columns(dataset, before_columns):
    return [str(column) for column in dataset.obs.columns if column not in before_columns]


def _input_files(dataset, args=None):
    """记录/校验输入文件清单（本地路径去重后返回绝对路径）。"""
    args = args or {}
    candidates = []
    for key in ("path", "paths"):
        value = args.get(key)
        if isinstance(value, str):
            candidates.append(value)
        elif isinstance(value, (list, tuple)):
            candidates.extend(str(item) for item in value if item)
    samples = args.get("samples")
    if isinstance(samples, list):
        for item in samples:
            if isinstance(item, dict) and item.get("path"):
                candidates.append(str(item["path"]))
    if not candidates and getattr(dataset, "path", None):
        candidates.append(str(dataset.path))
    files = []
    seen = set()
    for candidate in candidates:
        path = Path(os.path.expanduser(str(candidate))).resolve()
        if path.is_dir():
            continue
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        files.append(key)
    return files


def _open_replayable(dataset, args):
    args = args or {}
    if str(getattr(dataset, "source", "")) == "subset":
        return False, "子集数据集需要先重放父数据集，当前导出物只能作为记录"
    if args.get("url") or str(getattr(dataset, "path", "") or "").startswith(("http://", "https://", "ftp://")):
        return False, "URL 输入依赖网络与远端文件，当前导出物只能作为记录"
    files = _input_files(dataset, args)
    missing = [path for path in files if not Path(path).exists()]
    if missing:
        return False, "输入文件不存在，当前导出物只能作为记录：%s" % ", ".join(missing[:3])
    return True, None


def _record_import(dataset, args, source_format=None):
    replay_args = {key: args[key] for key in REPLAY_OPEN_ARG_KEYS if args.get(key) is not None}
    replayable, reason = _open_replayable(dataset, args)
    entry = _record_operation(
        dataset, "open", "open", replay_args,
        "import(%s)" % (source_format or getattr(dataset, "source_format", None) or dataset.source or "auto"),
        replayable=replayable,
        not_replayable_reason=reason,
        extra={
            "source_format": source_format or getattr(dataset, "source_format", None),
            "input_n_obs": int(getattr(dataset, "full_obs", dataset.obs).shape[0]),
            "input_n_vars": int(len(getattr(dataset, "full_var_names", dataset.var_names))),
            "matrix_kind": getattr(dataset, "matrix_kind", None),
            "counts_source": getattr(dataset, "counts_source", None),
            "counts_verified": bool(getattr(dataset, "counts_verified", True)),
        },
    )
    dataset.uns["import"] = dict(entry)
    return entry


def _history_entry(dataset, step, params, ok, ms, summary=None, error=None):
    return {
        "at": _now_iso(),
        "step": str(step),
        "params": _jsonable(params or {}),
        "n_obs": dataset.n_obs,
        "n_vars": dataset.n_vars,
        "ok": bool(ok),
        "ms": int(max(0, round(ms))),
        "engine": _engine_info(),
        "summary": _summary_line(summary) if ok else str(error or "步骤失败"),
    }


def _execute_step(dataset, step, user_params):
    """执行一个步骤：校验依赖 → 清空下游 → 计算 → 更新 steps/params → 写执行记录。"""
    if step not in STEP_ORDER:
        raise EngineError("bad_request", "未知步骤：%s（可用：%s）" % (step, ", ".join(STEP_ORDER)))
    missing = [dependency for dependency in STEP_DEPS[step] if dependency not in dataset.steps]
    if missing:
        raise EngineError(
            "prerequisite",
            "步骤 %s 需要先完成：%s" % (step, ", ".join(missing)),
        )
    params = _merge_params(step, user_params)
    _reset_from_step(dataset, step)
    before_scale = _matrix_scale(dataset)
    summary = STEP_HANDLERS[step](dataset, params)
    dataset.steps = _ordered_steps(set(dataset.steps) | {step})
    dataset.params[step] = params
    dataset.version += 1
    dataset.touch("step:%s" % step)
    # EXPORT-01：只有真正成功执行的步骤才进入结构化记录
    _record_step_execution(dataset, step, user_params, params, summary, before_scale)
    return params, summary


# ===========================================================================
# 7. 二维嵌入数据（UI canvas / embedding 图共用）
# ===========================================================================


def _downsample_indices(n_obs, max_points):
    max_points = int(max_points)
    if n_obs <= max_points:
        return np.arange(n_obs, dtype=np.int64), False
    rng = np.random.default_rng(SEED)
    indices = rng.choice(n_obs, size=max_points, replace=False)
    return np.sort(indices).astype(np.int64), True


def _category_order(labels, prefer_numeric=False):
    unique = sorted(set(str(label) for label in labels))
    if prefer_numeric:
        try:
            return sorted(unique, key=lambda item: (0, float(item)))
        except Exception:
            pass
    counts = pd.Series(unique).value_counts()
    return [str(label) for label in counts.index]


def _categorical_payload(labels_all, indices, legend_title, prefer_numeric=False):
    selected = np.asarray(labels_all, dtype=object)[indices]
    categories = _category_order(selected.tolist(), prefer_numeric=prefer_numeric)
    lookup = {label: index for index, label in enumerate(categories)}
    group_index = np.array([lookup[str(value)] for value in selected.tolist()], dtype=np.int64)
    palette = _palette(len(categories))
    counts = np.bincount(group_index, minlength=len(categories)) if categories else np.zeros(0, dtype=np.int64)
    payload_categories = [
        {
            "id": str(label),
            "label": str(label),
            "color": palette[index],
            "n": int(counts[index]) if index < counts.size else 0,
        }
        for index, label in enumerate(categories)
    ]
    return {
        "continuous": False,
        "legend_title": str(legend_title),
        "categories": payload_categories,
        "group_index": group_index.tolist(),
    }


def _continuous_payload(values_all, indices, legend_title, warnings):
    values = np.asarray(values_all, dtype=np.float64)[indices]
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    if values.size:
        vmin = float(values.min())
        vmax = float(values.max())
    else:
        vmin, vmax = 0.0, 1.0
    if not math.isfinite(vmin):
        vmin = 0.0
    if not math.isfinite(vmax):
        vmax = 1.0
    if vmax - vmin < 1e-12:
        vmax = vmin + 1.0
    n_bins = 8
    edges = np.linspace(vmin, vmax, n_bins + 1)
    bin_index = np.clip(np.digitize(values, edges[1:-1]), 0, n_bins - 1)
    colors = _gradient_colors(n_bins)
    categories = [
        {
            "id": str(index),
            "label": "%.2f–%.2f" % (edges[index], edges[index + 1]),
            "color": colors[index],
            "n": int(np.sum(bin_index == index)),
        }
        for index in range(n_bins)
    ]
    payload = {
        "continuous": True,
        "legend_title": str(legend_title),
        "values": [round(float(value), 4) for value in values],
        "vmin": vmin,
        "vmax": vmax,
        "categories": categories,
        "group_index": bin_index.astype(int).tolist(),
    }
    if warnings:
        payload["warnings"] = list(warnings)
    return payload


def _embedding_payload(dataset, color_by, max_points=12000, include_order=False):
    """构造 /embedding 返回体；coordinate 保留 3 位小数。"""
    if "X_emb" not in dataset.obsm:
        raise EngineError("prerequisite", "尚未运行二维嵌入（embed），请先运行 embed 步骤")
    if not isinstance(color_by, dict):
        color_by = {"kind": "cluster", "key": "clusters"}
    kind = str(color_by.get("kind") or "cluster")
    key = color_by.get("key")
    if key is None:
        key = "clusters" if kind == "cluster" else ""
    key = str(key)
    embedding = np.asarray(dataset.obsm["X_emb"], dtype=np.float64)
    n_obs = int(embedding.shape[0])
    max_points = max(100, min(int(max_points), 200000))
    indices, downsampled = _downsample_indices(n_obs, max_points)
    x = np.round(embedding[indices, 0], 3)
    y = np.round(embedding[indices, 1], 3)
    warnings = []

    if kind == "gene":
        if "normalize" not in dataset.steps:
            raise EngineError("prerequisite", "按基因着色需要先完成 normalize 步骤")
        if not key:
            raise EngineError("bad_request", "按基因着色需要提供 color_by.key（基因名）")
        column = dataset.gene_column(key)
        if column is None:
            warnings.append("基因 %s 不在当前数据集中，已按全 0 处理" % key)
            column = np.zeros(n_obs, dtype=np.float64)
        payload = _continuous_payload(column, indices, key, warnings)
    elif kind == "obs":
        if key not in dataset.obs.columns:
            raise EngineError(
                "not_found",
                "obs 列不存在：%s（可用：%s）"
                % (key, ", ".join(str(column) for column in dataset.obs.columns if not str(column).startswith("_"))),
            )
        series = dataset.obs[key]
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
            payload = _continuous_payload(values, indices, key, warnings)
        else:
            labels = [str(value) for value in series.tolist()]
            payload = _categorical_payload(labels, indices, key, prefer_numeric=False)
    elif kind == "cluster":
        if "clusters" not in dataset.obs.columns:
            raise EngineError("prerequisite", "聚类结果不存在，请先运行 cluster 步骤")
        if not key or key == "clusters":
            key = "clusters"
        if key not in dataset.obs.columns:
            raise EngineError("not_found", "obs 列不存在：%s" % key)
        labels = [str(value) for value in dataset.obs[key].tolist()]
        payload = _categorical_payload(labels, indices, key, prefer_numeric=(key == "clusters"))
    else:
        raise EngineError("bad_request", "未知 color_by.kind：%s（可用 cluster/obs/gene）" % kind)

    result = {
        "layout": str(dataset.layout or "tsne"),
        "n_points": int(indices.size),
        "downsampled": bool(downsampled),
        "x": [float(value) for value in x.tolist()],
        "y": [float(value) for value in y.tolist()],
        "continuous": bool(payload.get("continuous")),
        "legend_title": payload.get("legend_title", key),
        "categories": payload.get("categories", []),
        "group_index": payload.get("group_index", []),
    }
    if result["continuous"]:
        result["values"] = payload.get("values", [])
        result["vmin"] = payload.get("vmin")
        result["vmax"] = payload.get("vmax")
    if payload.get("warnings"):
        result["warnings"] = payload["warnings"]
    if include_order:
        result["_order"] = indices.tolist()
    return result


# ===========================================================================
# 8. 建议下一步（确定性，不调用模型）
# ===========================================================================


def _suggestion(identifier, title, detail, step=None, params=None, action=None):
    entry = {"id": identifier, "title": title, "detail": detail, "params": params or {}}
    if step:
        entry["step"] = step
    if action:
        entry["action"] = action
    return entry


def _suggestions_for(dataset):
    if dataset is None:
        return [
            _suggestion(
                "demo", "载入示例数据",
                "内置 2400 细胞 × 1200 基因的合成 PBMC 数据，可直接跑完整流程",
                action="demo",
            ),
            _suggestion(
                "open", "打开自己的 h5ad",
                "在顶栏输入 .h5ad 绝对路径，或直接上传文件",
                action="open",
            ),
        ]
    steps = set(dataset.steps)
    if "qc" not in steps:
        return [_suggestion("run:qc", "运行 QC", "先算每细胞基因数、总计数与线粒体占比，再决定过滤阈值", step="qc")]
    if "filter" not in steps:
        return [_suggestion(
            "run:filter", "过滤低质量细胞",
            "建议 min_genes=200、max_pct_mt=20、min_cells=3；可先在 QC 图里看分布再调",
            step="filter", params={"min_genes": 200, "max_pct_mt": 20, "min_cells": 3},
        )]
    if "normalize" not in steps:
        return [_suggestion(
            "run:normalize", "标准化", "按每细胞总计数缩放到 10000 后做 log1p",
            step="normalize", params={"target_sum": 10000},
        )]
    if "hvg" not in steps:
        return [_suggestion(
            "run:hvg", "选择高变基因", "按均值分箱的标准化离散度取 top 2000 高变基因",
            step="hvg", params={"n_top_genes": 2000, "flavor": "seurat"},
        )]
    if "pca" not in steps:
        return [_suggestion(
            "run:pca", "计算 PCA", "在高变基因上取 30 个主成分",
            step="pca", params={"n_pcs": 30},
        )]
    if "neighbors" not in steps:
        return [_suggestion(
            "run:neighbors", "构建邻接图", "k=15、metric=cosine 的 kNN 图",
            step="neighbors", params={"n_neighbors": 15, "metric": "cosine"},
        )]
    if "cluster" not in steps:
        return [_suggestion(
            "run:cluster", "运行聚类", "邻接图已就绪，建议 resolution=1.0 的 Louvain 聚类",
            step="cluster", params={"resolution": 1.0},
        )]
    if "embed" not in steps:
        return [_suggestion(
            "run:embed", "计算二维嵌入", "t-SNE（perplexity=30）把细胞投到二维平面上",
            step="embed", params={"layout": "tsne", "perplexity": 30},
        )]
    if "markers" not in steps:
        return [_suggestion(
            "run:markers", "计算标记基因", "每个簇做 one-vs-rest Welch t 检验，取 top 10",
            step="markers", params={"groupby": "clusters", "top_n": 10, "method": "t-test"},
        )]
    suggestions = []
    if "annotate" not in steps:
        suggestions.append(_suggestion(
            "autoannotate:clusters", "自动注释（参考 marker 面板）",
            "用内置 16 类参考 marker 面板给每个簇打分，确认后可用 annotate 写回 cell_type",
        ))
    if dataset.resolution_scan is None:
        suggestions.append(_suggestion(
            "resolution_scan", "做一次分辨率扫描",
            "在既有 kNN 图上比较 0.4–1.6 的簇数与 ARI，验证当前聚类是否稳健",
        ))
    suggestions.append(_suggestion(
        "plot:dotplot", "生成标记基因 dotplot",
        "把各簇 top marker 画成 cluster × gene 点图，用于确认注释",
    ))
    suggestions.append(_suggestion("export:h5ad", "导出 h5ad", "把当前标准化矩阵、簇、嵌入与原始计数一起导出", action="export"))
    return suggestions


# ===========================================================================
# 9. 自检 audit（reviewer 的轻量版）
# ===========================================================================


def _check(check_id, level, title, detail, hint=""):
    return {"id": check_id, "level": level, "title": title, "detail": detail, "hint": hint}


def _audit(dataset):
    checks = []
    steps = set(dataset.steps)
    has_qc = "qc" in steps
    has_filter = "filter" in steps
    has_cluster = "cluster" in steps
    n_obs = dataset.n_obs
    n_vars = dataset.n_vars

    # 1. 先质控再聚类
    if has_cluster and (not has_qc or not has_filter):
        checks.append(_check(
            "qc_before_cluster", "fail", "先质控再聚类",
            "聚类前没有完成 %s" % ("qc" if not has_qc else "filter"),
            "先运行 qc 与 filter，再重新聚类，避免低质量细胞主导簇结构",
        ))
    elif has_qc and has_filter:
        checks.append(_check("qc_before_cluster", "ok", "先质控再聚类", "已完成 qc/filter，聚类在过滤后的数据上运行", ""))
    else:
        checks.append(_check(
            "qc_before_cluster", "warn", "先质控再聚类",
            "尚未完成 QC/过滤", "先运行 qc 与 filter",
        ))

    # 2. 线粒体比例
    if not has_qc:
        checks.append(_check("mito_high", "warn", "线粒体比例", "尚未运行 QC，无法评估 pct_counts_mt", "先运行 qc"))
    else:
        try:
            pct = pd.to_numeric(dataset.obs["pct_counts_mt"], errors="coerce").to_numpy(dtype=np.float64)
            median_pct = float(np.nanmedian(pct)) if pct.size else 0.0
        except Exception:
            median_pct = 0.0
        if median_pct > 20.0:
            checks.append(_check(
                "mito_high", "warn", "线粒体比例",
                "过滤后 pct_counts_mt 中位数 %.1f%%（> 20%%），提示样本质量或阈值偏松" % median_pct,
                "考虑把 max_pct_mt 降到 15 或检查样本制备",
            ))
        else:
            checks.append(_check("mito_high", "ok", "线粒体比例", "过滤后 pct_counts_mt 中位数 %.2f%%" % median_pct, ""))

    # 3. 细胞数
    if n_obs < 200:
        checks.append(_check(
            "cells_low", "fail", "细胞数",
            "过滤后只剩 %d 个细胞（< 200），统计检验不稳定" % n_obs,
            "放宽过滤阈值或合并样本",
        ))
    else:
        checks.append(_check("cells_low", "ok", "细胞数", "过滤后 %d 个细胞 × %d 个基因" % (n_obs, n_vars), ""))

    # 4. 簇数
    categories_cluster, codes = dataset.cluster_codes()
    if not has_cluster or codes is None:
        checks.append(_check("clusters_few", "warn", "簇数", "尚未聚类", "先运行 neighbors 与 cluster"))
        checks.append(_check("cluster_skew", "warn", "簇大小分布", "尚未聚类", "先运行 cluster"))
    else:
        n_clusters = int(len(categories_cluster))
        if n_clusters < 2:
            checks.append(_check(
                "clusters_few", "fail", "簇数", "只有 %d 个簇（< 2），聚类分辨率可能过低" % n_clusters,
                "提高 resolution 或检查邻接图是否连通",
            ))
        else:
            checks.append(_check("clusters_few", "ok", "簇数", "%d 个簇" % n_clusters, ""))
        sizes = np.bincount(codes, minlength=n_clusters).astype(np.float64)
        max_fraction = float(sizes.max() / max(1.0, sizes.sum()))
        if max_fraction > 0.7:
            checks.append(_check(
                "cluster_skew", "warn", "簇大小分布",
                "最大簇占比 %.0f%%（> 70%%），可能存在未被拆分的混合群体" % (max_fraction * 100.0),
                "提高 resolution，或对最大簇做二次聚类",
            ))
        else:
            checks.append(_check(
                "cluster_skew", "ok", "簇大小分布",
                "最大簇占比 %.0f%%，分布均衡" % (max_fraction * 100.0), "",
            ))

    # 5. 标记基因
    if not has_cluster:
        checks.append(_check("markers_missing", "warn", "标记基因", "尚未聚类，无法计算标记基因", "先运行 cluster"))
    elif dataset.markers is None:
        checks.append(_check("markers_missing", "warn", "标记基因", "尚未计算标记基因", "运行 markers 步骤"))
    else:
        groups = dataset.markers.get("groups") or []
        empty_groups = [group["id"] for group in groups if int(group.get("n_significant", 0)) == 0]
        if empty_groups:
            checks.append(_check(
                "markers_missing", "warn", "标记基因",
                "簇 %s 没有显著正向标记（adj.P<0.05 且 logFC>0）" % "、".join(empty_groups[:8]),
                "检查这些簇是否为混合群体，或提高 top_n / 降低过滤强度",
            ))
        else:
            checks.append(_check("markers_missing", "ok", "标记基因", "%d 个簇都有显著标记基因" % len(groups), ""))

    # 6. 二维嵌入
    checks.append(
        _check("embedding_missing", "ok", "二维嵌入", "已生成 %s 二维嵌入" % str(dataset.layout).upper(), "")
        if "embed" in steps
        else _check("embedding_missing", "warn", "二维嵌入", "尚未运行 embed", "运行 embed 生成 t-SNE / UMAP 图")
    )

    # 7. 类型注释
    if "annotate" in steps:
        checks.append(_check("annotate_missing", "ok", "细胞类型注释", "已完成细胞类型注释", ""))
    elif has_cluster:
        checks.append(_check(
            "annotate_missing", "warn", "细胞类型注释", "尚未写回 cell_type（仍是簇号）",
            "先运行 autoannotate 看参考面板建议，再用 annotate 写回",
        ))
    else:
        checks.append(_check("annotate_missing", "warn", "细胞类型注释", "尚未聚类，无法注释", "先运行 cluster"))

    # 8. 双细胞风险（n_genes 极高细胞占比）
    if not has_qc:
        checks.append(_check("doublet_risk", "warn", "双细胞风险", "尚未运行 QC，无法评估 n_genes 分布", "先运行 qc"))
    else:
        try:
            genes_per_cell = pd.to_numeric(dataset.obs["n_genes"], errors="coerce").to_numpy(dtype=np.float64)
            genes_per_cell = genes_per_cell[np.isfinite(genes_per_cell)]
        except Exception:
            genes_per_cell = np.array([], dtype=np.float64)
        if genes_per_cell.size < 10:
            checks.append(_check("doublet_risk", "warn", "双细胞风险", "细胞数太少，无法评估", ""))
        else:
            q1, q3 = np.percentile(genes_per_cell, [25, 75])
            threshold = q3 + 1.5 * (q3 - q1)
            fraction = float(np.mean(genes_per_cell > threshold))
            if fraction > 0.05:
                checks.append(_check(
                    "doublet_risk", "warn", "双细胞风险",
                    "n_genes 极高的细胞占 %.1f%%（> 5%%），可能存在双细胞" % (fraction * 100.0),
                    "考虑做 doublet 检测（如 scDblFinder / scrublet）后再聚类",
                ))
            else:
                checks.append(_check(
                    "doublet_risk", "ok", "双细胞风险",
                    "n_genes 极高（> Q3+1.5·IQR）的细胞仅占 %.2f%%" % (fraction * 100.0), "",
                ))

    # 9. 分辨率稳健性
    if dataset.resolution_scan is None:
        checks.append(_check(
            "resolution_untested", "warn", "分辨率稳健性", "尚未做分辨率扫描",
            "运行 resolution_scan 比较 0.4–1.6 的簇数与 ARI",
        ))
    else:
        results = dataset.resolution_scan.get("results") or []
        checks.append(_check(
            "resolution_untested", "ok", "分辨率稳健性",
            "已扫描分辨率：%s" % ", ".join(_fmt_num(item.get("resolution")) for item in results),
            "",
        ))

    failures = sum(1 for check in checks if check["level"] == "fail")
    warnings = sum(1 for check in checks if check["level"] == "warn")
    score = int(max(0, min(100, 100 - 40 * failures - 10 * warnings)))
    if score >= 90:
        level = "ok"
    elif score >= 60:
        level = "warn"
    else:
        level = "fail"
    return {"score": score, "level": level, "checks": checks}


# ===========================================================================
# 10. 自动注释（内置参考 marker 面板）
# ===========================================================================

PANEL_ID = "builtin-marker-panel-v1"


def _autoannotate(dataset, groupby="clusters"):
    if groupby not in dataset.obs.columns:
        raise EngineError("prerequisite", "自动注释需要先完成聚类（obs 里没有 %s 列）" % groupby)
    if "normalize" not in dataset.steps:
        raise EngineError("prerequisite", "自动注释需要先完成 normalize")
    module = _make_sample_module()
    panel = getattr(module, "BUILTIN_MARKER_PANEL", None)
    if not isinstance(panel, dict) or not panel:
        raise EngineError("unsupported", "内置参考 marker 面板不可用")
    values, categories, masks = _group_labels_and_order(dataset, groupby)
    matrix = dataset.X.tocsr()
    n_vars = matrix.shape[1]
    means_all = np.asarray(matrix.mean(axis=0)).ravel().astype(np.float64)
    mean_squares = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel().astype(np.float64)
    std_all = np.sqrt(np.maximum(mean_squares - means_all ** 2, 0.0))
    var_index = dataset.var_index()

    gene_positions = {}
    gene_order = []
    for class_name in sorted(panel.keys()):
        for gene in panel[class_name]:
            if gene in var_index and gene not in gene_positions:
                gene_positions[gene] = int(var_index[gene])
                gene_order.append(gene)
    if not gene_order:
        raise EngineError("unsupported", "内置面板基因均不在当前数据集中，无法自动注释")
    gene_positions_in_all = {gene: index for index, gene in enumerate(gene_order)}
    column_index = np.array([gene_positions[gene] for gene in gene_order], dtype=np.int64)
    informative = std_all[column_index] > 1e-8
    panel_matrix = matrix[:, column_index]

    class_positions = {}
    for class_name in sorted(panel.keys()):
        positions = [gene_positions_in_all[gene] for gene in panel[class_name] if gene in gene_positions_in_all]
        if len(positions) >= 2:
            class_positions[class_name] = {
                "positions": np.asarray(positions, dtype=np.int64),
                "coverage": len(positions) / float(len(panel[class_name])),
            }
    if not class_positions:
        raise EngineError("unsupported", "当前数据集里没有足够的参考面板基因（每类至少 2 个）")

    suggestions = []
    for label in categories:
        mask = masks[label]
        if not mask.any():
            continue
        cluster_mean = np.asarray(panel_matrix[mask].mean(axis=0)).ravel().astype(np.float64)
        z_scores = np.zeros(cluster_mean.size, dtype=np.float64)
        good = informative
        z_scores[good] = (cluster_mean[good] - means_all[column_index][good]) / std_all[column_index][good]
        class_names = []
        class_scores = []
        for class_name, info in class_positions.items():
            positions = info["positions"]
            usable = good[positions]
            if not usable.any():
                continue
            mean_z = float(z_scores[positions][usable].mean())
            # 覆盖度惩罚：面板基因缺失越多，得分越往 0 收缩
            adjusted = mean_z * (0.6 + 0.4 * float(info["coverage"]))
            class_names.append(class_name)
            class_scores.append(adjusted)
        if not class_names:
            continue
        order = np.argsort(-np.asarray(class_scores), kind="stable")
        top_index = int(order[0])
        second_index = int(order[1]) if order.size > 1 else top_index
        raw_scores = np.asarray(class_scores, dtype=np.float64)
        shifted = raw_scores - raw_scores.max()
        probabilities = np.exp(shifted)
        probabilities = probabilities / probabilities.sum()
        top_label = class_names[top_index]
        evidence_positions = class_positions[top_label]["positions"]
        usable_positions = evidence_positions[good[evidence_positions]]
        if usable_positions.size:
            # evidence = 该簇里表达最高的 5 个面板基因（按簇内平均表达排序）
            evidence_order = usable_positions[np.argsort(-cluster_mean[usable_positions], kind="stable")]
            evidence = [gene_order[int(position)] for position in evidence_order[:5]]
        else:
            evidence = [gene_order[int(position)] for position in evidence_positions[:5]]
        suggestions.append(
            {
                "cluster": str(label),
                "label": str(top_label),
                "score": float(probabilities[top_index]),
                "margin": float(probabilities[top_index] - probabilities[second_index]),
                "z": float(raw_scores[top_index]),
                "evidence": evidence,
                "alternative": str(class_names[second_index]) if order.size > 1 else "",
            }
        )
    return {"panel": PANEL_ID, "groupby": str(groupby), "suggestions": suggestions}


# ===========================================================================
# 11. 分辨率扫描（参数对比，不改当前聚类）
# ===========================================================================


DEFAULT_RESOLUTIONS = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6]


def _resolution_scan(dataset, resolutions):
    if dataset.adj is None:
        raise EngineError("prerequisite", "分辨率扫描需要先完成 neighbors")
    categories_cluster, codes = dataset.cluster_codes()
    if codes is None:
        raise EngineError("prerequisite", "分辨率扫描需要先完成 cluster")
    if not resolutions:
        scan_values = list(DEFAULT_RESOLUTIONS)
    else:
        if not isinstance(resolutions, (list, tuple)):
            raise EngineError("bad_request", "resolutions 需要是数字数组")
        scan_values = []
        for value in list(resolutions)[:8]:
            scan_values.append(_limited_float(value, "resolutions", 1e-4, 100.0))
        if not scan_values:
            raise EngineError("bad_request", "resolutions 不能为空")
    current_labels = [str(value) for value in dataset.obs["clusters"].tolist()]
    current_resolution = float(dataset.params.get("cluster", {}).get("resolution", 1.0))
    n_current = int(len(categories_cluster))
    results = []
    group = "louvain_scan"
    for resolution in scan_values:
        labels = [str(int(value)) for value in louvain_clusters(dataset.adj, float(resolution))]
        n_clusters = int(len(set(labels)))
        ari = float(adjusted_rand_score(current_labels, labels))
        results.append(
            {
                "resolution": float(resolution),
                "n_clusters": n_clusters,
                "ari_vs_current": ari,
            }
        )
    dataset.resolution_scan = {
        "current": {"resolution": current_resolution, "n_clusters": n_current},
        "results": results,
        "scanned_at": _now_iso(),
        "group": group,
    }
    return {
        "current": {"resolution": current_resolution, "n_clusters": n_current},
        "results": results,
    }


# ===========================================================================
# 12. 预设流程（recipes，本地版 skills）
# ===========================================================================

RECIPES = [
    {
        "id": "standard",
        "title": "标准流程",
        "detail": "QC → 过滤 → 标准化 → 高变基因 → PCA → 邻接图 → 聚类 → 二维嵌入 → 标记基因",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {"min_genes": 200, "max_pct_mt": 20, "min_cells": 3}},
            {"step": "normalize", "params": {"target_sum": 10000}},
            {"step": "hvg", "params": {"n_top_genes": 2000, "flavor": "seurat"}},
            {"step": "pca", "params": {"n_pcs": 30}},
            {"step": "neighbors", "params": {"n_neighbors": 15, "metric": "cosine"}},
            {"step": "cluster", "params": {"resolution": 1.0, "method": "louvain"}},
            {"step": "embed", "params": {"layout": "tsne", "perplexity": 30}},
            {"step": "markers", "params": {"groupby": "clusters", "top_n": 10, "method": "t-test"}},
        ],
    },
    {
        "id": "quick",
        "title": "快速流程",
        "detail": "QC → 过滤 → 标准化 → 高变基因 → PCA → 邻接图 → 聚类 → 二维嵌入（跳过标记基因，适合先看结构）",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {"min_genes": 200, "max_pct_mt": 20, "min_cells": 3}},
            {"step": "normalize", "params": {"target_sum": 10000}},
            {"step": "hvg", "params": {"n_top_genes": 2000, "flavor": "seurat"}},
            {"step": "pca", "params": {"n_pcs": 30}},
            {"step": "neighbors", "params": {"n_neighbors": 15, "metric": "cosine"}},
            {"step": "cluster", "params": {"resolution": 1.0, "method": "louvain"}},
            {"step": "embed", "params": {"layout": "tsne", "perplexity": 30}},
        ],
    },
    {
        "id": "qc",
        "title": "质控流程",
        "detail": "QC 指标 + 过滤（之后可在 QC tab 生成 QC 图，再决定是否继续）",
        "steps": [
            {"step": "qc", "params": {}},
            {"step": "filter", "params": {"min_genes": 200, "max_pct_mt": 20, "min_cells": 3}},
        ],
    },
]


_PRO_QC = {"min_genes": 200, "max_genes": 5000, "min_counts": 500, "max_pct_mt": 20, "min_cells": 3}
_PRO_BASE = [
    {"step": "qc", "params": {}},
    {"step": "filter", "params": dict(_PRO_QC)},
    {"step": "normalize", "params": {"target_sum": 10000}},
    {"step": "hvg", "params": {"n_top_genes": 3000, "flavor": "seurat"}},
]

# 对齐用户 pro.R：未整合(PCA, res 0.5) / Harmony(res 0.8) / RPCA(res 0.8)，dims 1:30
SEURAT6_RECIPE = {
    "id": "seurat6",
    "title": "Seurat 6.0 + Harmony（对齐 pro.R）",
    "detail": "QC(nFeature 200-5000, 线粒体<20%, nCount>500) → HVG 3000 → PCA 30 → IntegrateLayers(Harmony) → 聚类 resolution 0.8 → UMAP(30)",
    "steps": _PRO_BASE + [
        {"step": "pca", "params": {"n_pcs": 30, "n_top_genes": 3000, "dims": 30, "integration": "harmony", "engine": "seurat"}},
        {"step": "neighbors", "params": {"n_neighbors": 30, "dims": 30, "integration": "harmony", "engine": "seurat"}},
        {"step": "cluster", "params": {"resolution": 0.8, "dims": 30, "integration": "harmony", "engine": "seurat"}},
        {"step": "embed", "params": {"layout": "umap", "n_neighbors": 30, "dims": 30, "integration": "harmony", "engine": "seurat"}},
    ],
}
PRO_RPCA_RECIPE = {
    "id": "pro-rpca",
    "title": "Seurat 6.0 + RPCA（对齐 pro.R）",
    "detail": "QC 同 pro.R → HVG 3000 → PCA 30 → IntegrateLayers(RPCA) → 聚类 resolution 0.8 → UMAP(30)",
    "steps": _PRO_BASE + [
        {"step": "pca", "params": {"n_pcs": 30, "n_top_genes": 3000, "dims": 30, "integration": "rpca", "engine": "seurat"}},
        {"step": "neighbors", "params": {"n_neighbors": 30, "dims": 30, "integration": "rpca", "engine": "seurat"}},
        {"step": "cluster", "params": {"resolution": 0.8, "dims": 30, "integration": "rpca", "engine": "seurat"}},
        {"step": "embed", "params": {"layout": "umap", "n_neighbors": 30, "dims": 30, "integration": "rpca", "engine": "seurat"}},
    ],
}
PRO_UNINTEGRATED_RECIPE = {
    "id": "pro-unintegrated",
    "title": "Seurat 6.0 未整合（对齐 pro.R）",
    "detail": "QC 同 pro.R → HVG 3000 → PCA 30 → 未整合 PCA → 聚类 resolution 0.5 → UMAP(30)",
    "steps": _PRO_BASE + [
        {"step": "pca", "params": {"n_pcs": 30, "n_top_genes": 3000, "dims": 30, "integration": "none", "engine": "seurat"}},
        {"step": "neighbors", "params": {"n_neighbors": 30, "dims": 30, "integration": "none", "engine": "seurat"}},
        {"step": "cluster", "params": {"resolution": 0.5, "dims": 30, "integration": "none", "engine": "seurat"}},
        {"step": "embed", "params": {"layout": "umap", "n_neighbors": 30, "dims": 30, "integration": "none", "engine": "seurat"}},
    ],
}


def _recipes_payload(recipe_id=None):
    recipes = list(RECIPES) + [SEURAT6_RECIPE, PRO_RPCA_RECIPE, PRO_UNINTEGRATED_RECIPE]
    if recipe_id is None:
        return {"recipes": [dict(recipe) for recipe in recipes]}
    wanted = str(recipe_id)
    for recipe in recipes:
        if recipe["id"] == wanted:
            return dict(recipe)
    raise EngineError(
        "bad_request",
        "未知预设流程：%s（可用：%s）" % (wanted, ", ".join(recipe["id"] for recipe in RECIPES)),
    )


# ===========================================================================
# 13. 出图（7 种，matplotlib / Agg / 白底 / 全英文）
# ===========================================================================

PLOT_KINDS = ["qc_violin", "embedding", "dotplot", "heatmap", "cluster_bar", "pca_variance", "qc_scatter"]
QC_METRICS = ["n_genes", "n_counts", "pct_counts_mt"]


def _require_step(dataset, step, message):
    if step not in dataset.steps:
        raise EngineError("prerequisite", message)


def _new_figure(width, height):
    return plt.subplots(figsize=(float(width), float(height)), dpi=140)


def _ensure_cluster_markers(dataset, top_n):
    top_n = int(max(1, top_n))
    if "clusters" not in dataset.obs.columns:
        raise EngineError("prerequisite", "需要先完成 cluster 步骤")
    cached = dataset.markers
    if (
        isinstance(cached, dict)
        and cached.get("groupby") == "clusters"
        and int(cached.get("top_n", 0)) >= top_n
        and cached.get("groups")
    ):
        return cached["groups"]
    if "normalize" not in dataset.steps:
        raise EngineError("prerequisite", "计算标记基因需要先完成 normalize")
    return _compute_markers(dataset, "clusters", max(top_n, 10), "t-test")


def _cluster_groups_by_id(dataset):
    categories, codes = dataset.cluster_codes()
    if categories is None:
        raise EngineError("prerequisite", "需要先完成 cluster 步骤")
    return categories, codes


def _plot_qc_violin(dataset, params):
    _require_step(dataset, "qc", "生成 QC 图需要先完成 qc 步骤")
    metrics = params.get("metrics") or list(QC_METRICS)
    if not isinstance(metrics, (list, tuple)):
        raise EngineError("bad_request", "qc_violin 的 metrics 需要是数组")
    metrics = [str(metric) for metric in metrics]
    if not metrics:
        raise EngineError("bad_request", "qc_violin 的 metrics 不能为空")
    for metric in metrics:
        if metric not in dataset.obs.columns:
            raise EngineError("bad_request", "QC 指标不存在：%s（可用：%s）" % (metric, ", ".join(QC_METRICS)))
    thresholds = params.get("thresholds") or {}
    if not isinstance(thresholds, dict):
        thresholds = {}
    count = max(1, len(metrics))
    fig, axes = plt.subplots(1, count, figsize=(3.0 * count + 0.6, 3.5), dpi=140, squeeze=False)
    colors = _palette(count)
    for index, metric in enumerate(metrics):
        axis = axes[0][index]
        values = pd.to_numeric(dataset.obs[metric], errors="coerce").to_numpy(dtype=np.float64)
        values = values[np.isfinite(values)]
        if values.size:
            parts = axis.violinplot([values], showmedians=True, showextrema=False, widths=0.75)
            for body in parts["bodies"]:
                body.set_facecolor(colors[index])
                body.set_edgecolor("#555555")
                body.set_alpha(0.75)
            if "cmedians" in parts:
                parts["cmedians"].set_color("#222222")
                parts["cmedians"].set_linewidth(1.2)
        threshold = thresholds.get(metric)
        if threshold is not None:
            try:
                threshold_value = float(threshold)
                axis.axhline(threshold_value, color="#e15759", linestyle="--", linewidth=1.0,
                             label="threshold %s" % _fmt_num(threshold_value))
                axis.legend(loc="best", frameon=False)
            except Exception:
                pass
        axis.set_title(metric)
        axis.set_xticks([])
        axis.set_ylabel("value")
    fig.suptitle("QC metrics", fontsize=11)
    fig.tight_layout()
    return fig


def _plot_embedding(dataset, params):
    _require_step(dataset, "embed", "生成嵌入图需要先完成 embed 步骤")
    color_by = params.get("color_by") or {"kind": "cluster", "key": "clusters"}
    payload = _embedding_payload(dataset, color_by, max_points=8000, include_order=True)
    order = np.asarray(payload.pop("_order", []), dtype=np.int64)
    x = np.asarray(payload["x"], dtype=np.float64)
    y = np.asarray(payload["y"], dtype=np.float64)
    fig, axis = plt.subplots(figsize=(6.2, 5.0), dpi=140)
    if payload.get("continuous"):
        values = np.asarray(payload["values"], dtype=np.float64)
        scatter = axis.scatter(x, y, c=values, cmap="viridis", s=5, linewidths=0, alpha=0.85)
        colorbar = fig.colorbar(scatter, ax=axis, fraction=0.046, pad=0.02)
        colorbar.set_label(payload["legend_title"])
    else:
        categories = payload.get("categories") or []
        group_index = np.asarray(payload.get("group_index", []), dtype=np.int64)
        palette = [category["color"] for category in categories]
        point_colors = [palette[index] if 0 <= index < len(palette) else "#999999" for index in group_index]
        axis.scatter(x, y, c=point_colors, s=5, linewidths=0, alpha=0.85)
        if 0 < len(categories) <= 15:
            handles = [
                matplotlib.lines.Line2D(
                    [0], [0], marker="o", linestyle="", markersize=5,
                    markerfacecolor=category["color"], markeredgewidth=0,
                    label="%s (n=%d)" % (category["label"], category["n"]),
                )
                for category in categories
            ]
            axis.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
                        frameon=False, fontsize=7)
    cluster_categories, cluster_codes = dataset.cluster_codes()
    if cluster_codes is not None and order.size:
        subset_codes = cluster_codes[order]
        for code, label in enumerate(cluster_categories):
            mask = subset_codes == code
            if not mask.any():
                continue
            center_x = float(np.median(x[mask]))
            center_y = float(np.median(y[mask]))
            axis.text(center_x, center_y, str(label), fontsize=8, fontweight="bold",
                      ha="center", va="center",
                      bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                                edgecolor="#888888", linewidth=0.5, alpha=0.8))
    axis.set_xlabel("dimension 1")
    axis.set_ylabel("dimension 2")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_title("Embedding (%s) - colored by %s" % (payload.get("layout", ""), payload.get("legend_title", "")))
    fig.tight_layout()
    return fig


def _marker_gene_matrix(dataset, groups, top_n):
    """cluster × gene 的平均表达 / 表达比例矩阵（供 dotplot / heatmap 使用）。"""
    gene_index = dataset.var_index()
    genes = []
    for group in groups:
        for record in group["genes"][:top_n]:
            gene = str(record["gene"])
            if gene in gene_index and gene not in genes:
                genes.append(gene)
    if not genes:
        raise EngineError("bad_request", "标记基因结果为空，无法出图")
    positions = np.array([int(gene_index[gene]) for gene in genes], dtype=np.int64)
    matrix = dataset.X.tocsr()[:, positions]
    clusters = [str(group["id"]) for group in groups]
    mean_matrix = np.zeros((len(clusters), len(genes)), dtype=np.float64)
    pct_matrix = np.zeros_like(mean_matrix)
    values = [str(value) for value in dataset.obs["clusters"].tolist()]
    for row, label in enumerate(clusters):
        mask = np.array([value == label for value in values], dtype=bool)
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        subset = matrix[mask]
        mean_matrix[row] = np.asarray(subset.mean(axis=0)).ravel()
        pct_matrix[row] = np.asarray(subset.getnnz(axis=0)).ravel() / float(n_cells)
    return genes, clusters, mean_matrix, pct_matrix


def _selected_gene_list(dataset, params):
    """params.genes 指定的基因列表（支持数组 / 逗号/空格分隔字符串），校验并去重。"""
    raw = params.get("genes") if isinstance(params, dict) else None
    if isinstance(raw, str):
        raw = [item for item in raw.replace(",", " ").replace(";", " ").replace("\n", " ").split() if item]
    if not isinstance(raw, (list, tuple)) or not raw:
        return []
    index = dataset.var_index()
    genes, missing = [], []
    for item in raw:
        gene = str(item).strip()
        if not gene:
            continue
        if gene in index:
            if gene not in genes:
                genes.append(gene)
        elif gene not in missing:
            missing.append(gene)
    if missing:
        dataset.warnings.append("出图忽略不在表达矩阵中的基因：%s" % ", ".join(missing[:10]))
    return genes[:120]


def _parse_gene_groups(params):
    """解析 params.gene_groups = [{"label": "T cell", "genes": ["CD3D", ...]}, ...]。"""
    raw = params.get("gene_groups") if isinstance(params, dict) else None
    if not isinstance(raw, (list, tuple)) or not raw:
        return []
    groups = []
    for item in raw:
        if not isinstance(item, dict):
            continue
        label = str(item.get("label") or item.get("name") or "").strip()
        genes = item.get("genes")
        if isinstance(genes, str):
            genes = [value for value in genes.replace(",", " ").replace(";", " ").split() if value]
        if not isinstance(genes, (list, tuple)):
            continue
        clean = [str(gene).strip() for gene in genes if str(gene).strip()]
        if clean:
            groups.append({"label": label, "genes": clean})
    return groups


def _flatten_gene_groups(dataset, groups):
    """把 gene_groups 展平成可用基因列表 + 分组边界；不在矩阵里的基因丢弃。"""
    index = dataset.var_index()
    ordered, bounds, labels, missing = [], [], [], []
    for group in groups:
        genes = []
        for gene in group["genes"]:
            if gene in index and gene not in ordered:
                genes.append(gene)
            elif gene not in index and gene not in missing:
                missing.append(gene)
        if not genes:
            continue
        start = len(ordered)
        ordered.extend(genes)
        bounds.append((start, len(ordered)))
        labels.append(group["label"] or ("Group %d" % len(labels)))
    if missing:
        dataset.warnings.append("点图/热图忽略不在表达矩阵中的基因：%s" % ", ".join(missing[:10]))
    return ordered, bounds, labels


def _expression_by_cluster(dataset, genes, clusters=None):
    """选中基因 × 簇 的平均表达 / 表达比例矩阵（用户手动选 marker 出图用）。"""
    index = dataset.var_index()
    positions = np.array([int(index[gene]) for gene in genes], dtype=np.int64)
    matrix = dataset.X.tocsr()[:, positions]
    values = [str(value) for value in dataset.obs["clusters"].tolist()]
    if clusters is None:
        clusters = sorted(set(values), key=lambda item: (len(item), item))
    else:
        clusters = [str(item) for item in clusters]
    mean_matrix = np.zeros((len(clusters), len(genes)), dtype=np.float64)
    pct_matrix = np.zeros_like(mean_matrix)
    for row, label in enumerate(clusters):
        mask = np.array([value == label for value in values], dtype=bool)
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        subset = matrix[mask]
        mean_matrix[row] = np.asarray(subset.mean(axis=0)).ravel()
        pct_matrix[row] = np.asarray(subset.getnnz(axis=0)).ravel() / float(n_cells)
    return list(genes), clusters, mean_matrix, pct_matrix


def _wrap_group_label(label, width=10):
    text = str(label or "").strip()
    if not text:
        return ""
    words, lines, current = text.replace("/", " / ").split(), [], ""
    for word in words:
        candidate = (current + " " + word).strip()
        if len(candidate) <= width or not current:
            current = candidate
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    return "\n".join(lines)


def _plot_dotplot(dataset, params):
    """对齐用户「点图plus版.r」：基因按细胞类型分面（列）｜ cluster 为行｜
    点大小 = 表达比例，颜色 = 蓝-白-红（按基因 z-score）。"""
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    _require_step(dataset, "cluster", "生成 dotplot 需要先完成 cluster 步骤")
    gene_groups = _parse_gene_groups(params)
    bounds, groups = [], []
    if gene_groups:
        genes, bounds, labels = _flatten_gene_groups(dataset, gene_groups)
        if not genes:
            gene_groups = []
    if not gene_groups:
        selected = _selected_gene_list(dataset, params)
        if selected:
            genes = selected
        else:
            top_n = max(1, min(20, _as_int(params.get("top_n", 5), "top_n")))
            marker_groups = _ensure_cluster_markers(dataset, top_n)
            genes, _clusters, _mean, _pct = _marker_gene_matrix(dataset, marker_groups, top_n)
        bounds = [(0, len(genes))]
        labels = [""]
    genes, clusters, mean_matrix, pct_matrix = _expression_by_cluster(dataset, genes)

    # 按基因在 cluster 间 z-score（= Seurat DotPlot 的 avg.exp.scaled）
    centered = mean_matrix - mean_matrix.mean(axis=0, keepdims=True)
    scale = mean_matrix.std(axis=0, keepdims=True)
    scale[scale < 1e-9] = 1.0
    scaled = np.clip(centered / scale, -2.5, 2.5)

    # 基因位置：组内连续、组间留空隙
    positions = np.zeros(len(genes), dtype=float)
    cursor = 0.0
    spans = []
    for gi, (start, stop) in enumerate(bounds):
        stop = min(stop, len(genes))
        if stop <= start:
            continue
        width = stop - start
        positions[start:stop] = np.arange(width, dtype=float) + cursor
        spans.append((labels[gi], start, stop, cursor, cursor + width - 1))
        cursor += width + 0.9
    n_cols = cursor if cursor > 0 else float(len(genes))

    n_rows = len(clusters)
    width_in = max(6.0, 0.28 * max(1.0, n_cols) + 3.2)
    height_in = max(3.2, 0.30 * n_rows + 1.6)
    fig, axis = plt.subplots(figsize=(width_in, height_in), dpi=150)
    cmap = LinearSegmentedColormap.from_list("dsh_dotplot", ["#3366CC", "#FFFFFF", "#DC3912"])
    xx, yy = np.meshgrid(positions, np.arange(n_rows))
    scatter = axis.scatter(
        xx.ravel(), yy.ravel(),
        s=(pct_matrix.ravel() * 170.0) + 5.0,
        c=scaled.ravel(), cmap=cmap, vmin=-2.5, vmax=2.5,
        linewidths=0.25, edgecolors="#666666",
    )
    colorbar = fig.colorbar(scatter, ax=axis, fraction=0.025, pad=0.015)
    colorbar.set_label("average expression level", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)

    # 细胞类型分面分隔线 + 粗体标题（自动换行）
    for label, start, stop, _left, _right in spans:
        if start > 0:
            axis.axvline(float(positions[start]) - 0.55, color="#999999", linewidth=0.8, linestyle="--", alpha=0.8)
        if label:
            center = (float(positions[start]) + float(positions[stop - 1])) / 2.0
            axis.text(center, n_rows - 0.45, _wrap_group_label(label),
                      ha="center", va="bottom", fontsize=9, fontweight="bold",
                      linespacing=1.15, clip_on=False)

    axis.set_xticks(positions)
    axis.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    axis.set_yticks(np.arange(n_rows))
    axis.set_yticklabels(clusters, fontsize=9)
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.set_xlim(float(positions.min()) - 0.8, float(positions.max()) + 0.8)
    axis.set_ylim(-0.6, n_rows - 0.4)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(length=0)

    # 表达比例大小图例
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777",
               markeredgecolor="#666666", markeredgewidth=0.25, markersize=(p * 170.0 + 5.0) ** 0.5,
               label="%d%%" % int(p * 100))
        for p in (0.25, 0.5, 0.75, 1.0)
    ]
    axis.legend(handles=handles, title="Expressing proportions (%)", loc="upper left",
                bbox_to_anchor=(1.12, 0.55), frameon=False, fontsize=7, title_fontsize=8)
    fig.tight_layout()
    return fig


def _plot_heatmap(dataset, params):
    _require_step(dataset, "cluster", "生成 heatmap 需要先完成 cluster 步骤")
    groups_raw = _parse_gene_groups(params)
    selected = [] if groups_raw else _selected_gene_list(dataset, params)
    if groups_raw:
        ordered, _bounds_h, _labels_h = _flatten_gene_groups(dataset, groups_raw)
        if ordered:
            genes, clusters, mean_matrix, _pct = _expression_by_cluster(dataset, ordered)
            title = "Marker genes by cell type (heatmap)"
        else:
            genes = None
    else:
        genes = None
    if genes is not None:
        pass
    elif selected:
        genes, clusters, mean_matrix, _pct = _expression_by_cluster(dataset, selected)
        title = "Selected marker genes heatmap (%d)" % len(genes)
    else:
        top_n = max(1, min(20, _as_int(params.get("top_n", 5), "top_n")))
        groups = _ensure_cluster_markers(dataset, top_n)
        genes, clusters, mean_matrix, _pct = _marker_gene_matrix(dataset, groups, top_n)
        title = "Marker gene heatmap (top %d per cluster)" % top_n
    centered = mean_matrix - mean_matrix.mean(axis=0, keepdims=True)
    scale = mean_matrix.std(axis=0, keepdims=True)
    scale[scale < 1e-9] = 1.0
    z_matrix = centered / scale
    z_matrix = np.clip(z_matrix, -2.5, 2.5)
    width = max(4.2, 0.28 * len(genes) + 2.4)
    height = max(3.0, 0.30 * len(clusters) + 1.6)
    fig, axis = plt.subplots(figsize=(width, height), dpi=140)
    image = axis.imshow(z_matrix.T, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    colorbar = fig.colorbar(image, ax=axis, fraction=0.03, pad=0.02)
    colorbar.set_label("z-score of mean expression")
    axis.set_xticks(np.arange(len(clusters)))
    axis.set_xticklabels(clusters, fontsize=8)
    axis.set_yticks(np.arange(len(genes)))
    axis.set_yticklabels(genes, fontsize=7)
    axis.set_xlabel("cluster")
    axis.set_ylabel("marker gene")
    axis.set_title(title)
    fig.tight_layout()
    return fig


def _plot_cluster_bar(dataset, params):
    categories, codes = _cluster_groups_by_id(dataset)
    sizes = np.bincount(codes, minlength=len(categories))
    palette = _palette(len(categories))
    fig, axis = plt.subplots(figsize=(max(4.2, 0.55 * len(categories) + 2.2), 3.6), dpi=140)
    positions = np.arange(len(categories))
    axis.bar(positions, sizes, color=palette, edgecolor="#444444", linewidth=0.6)
    for index, size in enumerate(sizes):
        axis.text(index, size, str(int(size)), ha="center", va="bottom", fontsize=7)
    axis.set_xticks(positions)
    axis.set_xticklabels(categories, fontsize=8, rotation=0 if len(categories) <= 12 else 45)
    axis.set_xlabel("cluster")
    axis.set_ylabel("number of cells")
    axis.set_title("Cells per cluster (n=%d)" % int(sizes.sum()))
    axis.margins(y=0.12)
    fig.tight_layout()
    return fig


def _plot_pca_variance(dataset, params):
    if "X_pca" not in dataset.obsm:
        raise EngineError("prerequisite", "生成 PCA 碎石图需要先完成 pca 步骤")
    ratios = np.asarray(dataset.uns.get("pca_variance_ratio") or [], dtype=np.float64)
    if ratios.size == 0:
        raise EngineError("prerequisite", "缺少 PCA 解释方差信息，请重跑 pca")
    positions = np.arange(1, ratios.size + 1)
    fig, axis = plt.subplots(figsize=(6.0, 3.8), dpi=140)
    axis.bar(positions, ratios * 100.0, color="#4e79a7", edgecolor="#2c4a66", linewidth=0.4)
    axis.plot(positions, ratios * 100.0, color="#e15759", linewidth=1.0, marker="o", markersize=3)
    axis.set_xlabel("principal component")
    axis.set_ylabel("explained variance (%)")
    axis.set_title("PCA explained variance (PC1 = %.1f%%)" % (float(ratios[0]) * 100.0))
    axis.set_xticks(positions[:: max(1, ratios.size // 15)])
    fig.tight_layout()
    return fig


def _plot_qc_scatter(dataset, params):
    _require_step(dataset, "qc", "生成 QC 散点图需要先完成 qc 步骤")
    x_metric = str(params.get("x_metric") or "n_genes")
    y_metric = str(params.get("y_metric") or "n_counts")
    for metric in (x_metric, y_metric):
        if metric not in dataset.obs.columns:
            raise EngineError("bad_request", "QC 指标不存在：%s" % metric)
    color_by_pct_mt = params.get("color_by_pct_mt", True)
    color_by_pct_mt = bool(color_by_pct_mt) if isinstance(color_by_pct_mt, (bool, int)) else True
    x_values = pd.to_numeric(dataset.obs[x_metric], errors="coerce").to_numpy(dtype=np.float64)
    y_values = pd.to_numeric(dataset.obs[y_metric], errors="coerce").to_numpy(dtype=np.float64)
    fig, axis = plt.subplots(figsize=(5.8, 4.4), dpi=140)
    if color_by_pct_mt and "pct_counts_mt" in dataset.obs.columns:
        pct_values = pd.to_numeric(dataset.obs["pct_counts_mt"], errors="coerce").to_numpy(dtype=np.float64)
        scatter = axis.scatter(x_values, y_values, c=pct_values, cmap="magma_r", s=5, linewidths=0, alpha=0.8)
        colorbar = fig.colorbar(scatter, ax=axis, fraction=0.046, pad=0.02)
        colorbar.set_label("pct_counts_mt (%)")
    else:
        axis.scatter(x_values, y_values, c="#4e79a7", s=5, linewidths=0, alpha=0.6)
    axis.set_xlabel(x_metric)
    axis.set_ylabel(y_metric)
    axis.set_title("QC scatter: %s vs %s" % (y_metric, x_metric))
    fig.tight_layout()
    return fig


PLOT_BUILDERS = {
    "qc_violin": _plot_qc_violin,
    "embedding": _plot_embedding,
    "dotplot": _plot_dotplot,
    "heatmap": _plot_heatmap,
    "cluster_bar": _plot_cluster_bar,
    "pca_variance": _plot_pca_variance,
    "qc_scatter": _plot_qc_scatter,
}


def _figure_hash(dataset, kind, params):
    payload = json.dumps(
        {"kind": kind, "params": _jsonable(params), "version": dataset.version},
        sort_keys=True, ensure_ascii=False,
    )
    return hashlib.md5(payload.encode("utf-8")).hexdigest()[:12]


def _canonical_plot_params(kind, params):
    """把缺省参数补全，保证 UI / report / 不同调用方式的同参数图复用同一文件。"""
    params = dict(params or {})
    if kind == "embedding":
        color_by = params.get("color_by")
        if not isinstance(color_by, dict) or not color_by:
            params["color_by"] = {"kind": "cluster", "key": "clusters"}
    elif kind == "qc_violin":
        if not params.get("metrics"):
            params["metrics"] = list(QC_METRICS)
        if not isinstance(params.get("thresholds"), dict):
            params.pop("thresholds", None)
    elif kind in ("dotplot", "heatmap"):
        params["top_n"] = params.get("top_n", 5)
    elif kind == "qc_scatter":
        params.setdefault("x_metric", "n_genes")
        params.setdefault("y_metric", "n_counts")
        params.setdefault("color_by_pct_mt", True)
    return params


def _make_plot(dataset, kind, params):
    if kind not in PLOT_BUILDERS:
        raise EngineError("bad_request", "未知图类型：%s（可用：%s）" % (kind, ", ".join(PLOT_KINDS)))
    if params is None:
        params = {}
    if not isinstance(params, dict):
        raise EngineError("bad_request", "plot params 必须是 JSON 对象")
    params = _canonical_plot_params(kind, _jsonable(params))
    digest = _figure_hash(dataset, kind, params)
    path = FIG_DIR / ("fig-%s-%s-%s.png" % (dataset.id, kind, digest))
    if path.exists() and path.stat().st_size > 0:
        dimensions = _png_dimensions(path)
        if dimensions:
            return {
                "file": str(path),
                "width": int(dimensions[0]),
                "height": int(dimensions[1]),
                "bytes": int(path.stat().st_size),
                "kind": kind,
                "reused": True,
            }
    figure = PLOT_BUILDERS[kind](dataset, params)
    try:
        figure.savefig(str(path), dpi=140, facecolor="white")
    finally:
        plt.close(figure)
    dimensions = _png_dimensions(path) or (0, 0)
    entry = {
        "file": str(path),
        "width": int(dimensions[0]),
        "height": int(dimensions[1]),
        "bytes": int(path.stat().st_size),
        "kind": kind,
        "reused": False,
    }
    if not any(item.get("file") == entry["file"] for item in dataset.figures):
        dataset.figures.append(
            {
                "kind": kind,
                "file": str(path),
                "params": params,
                "width": entry["width"],
                "height": entry["height"],
                "bytes": entry["bytes"],
            }
        )
        if len(dataset.figures) > 60:
            dataset.figures = dataset.figures[-60:]
    return entry


# ===========================================================================
# 14. 导出（h5ad / csv / json + v2 的 protocol / methods / report / script）
# ===========================================================================

EXPORT_EXTENSIONS = {
    "h5ad": ".h5ad",
    "markers_csv": ".csv",
    "obs_csv": ".csv",
    "summary_json": ".json",
    "protocol": ".json",
    "methods": ".md",
    "report": ".md",
    "script": ".py",
}
EXPORT_WHATS = list(EXPORT_EXTENSIONS.keys())


def _unique_names(index):
    names = [str(value) for value in index.tolist()]
    seen = set()
    output = []
    counters = {}
    for name in names:
        if name not in seen:
            seen.add(name)
            output.append(name)
            continue
        counter = counters.get(name, 1)
        candidate = "%s-%d" % (name, counter)
        while candidate in seen:
            counter += 1
            candidate = "%s-%d" % (name, counter)
        counters[name] = counter + 1
        seen.add(candidate)
        output.append(candidate)
    return output


def _unique_index(index, default_name):
    return pd.Index(_unique_names(index), name=index.name or default_name)


def _sanitize_obs_for_h5ad(obs):
    obs = obs.copy()
    obs.index = _unique_index(obs.index, "cell")
    for column in list(obs.columns):
        series = obs[column]
        if series.dtype == object:
            try:
                obs[column] = series.astype(str)
            except Exception:
                pass
    return obs


def _default_export_path(dataset, what):
    extension = EXPORT_EXTENSIONS[what]
    stamp = time.strftime("%Y%m%d-%H%M%S") + "-%03d" % (int(time.time() * 1000) % 1000)
    return EXPORT_DIR / ("%s-%s-%s%s" % (dataset.id, what, stamp, extension))


def _resolve_export_path(dataset, what, user_path):
    if user_path is None or str(user_path).strip() == "":
        path = _default_export_path(dataset, what)
    else:
        raw = os.path.expanduser(str(user_path).strip())
        path = Path(raw)
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        if path.is_dir():
            path = path / ("%s-%s" % (dataset.id, what))
        if path.suffix == "":
            path = path.with_suffix(EXPORT_EXTENSIONS[what])
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as error:
            raise EngineError("bad_request", "无法创建导出目录（%s）：%s" % (path.parent, error))
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except Exception as error:
        raise EngineError("bad_request", "无法创建导出目录（%s）：%s" % (path.parent, error))
    return path


def _export_h5ad(dataset, path):
    if ad is None:
        raise EngineError("unsupported", "未安装 anndata，无法导出 h5ad")
    try:
        obs = _sanitize_obs_for_h5ad(dataset.obs)
        var = dataset.var.copy()
        if len(var) == dataset.n_vars:
            var.index = pd.Index(dataset.var_names, name="gene")
        else:
            var = pd.DataFrame(index=pd.Index(dataset.var_names, name="gene"))
        matrix = dataset.X.tocsr() if sp.issparse(dataset.X) else np.asarray(dataset.X, dtype=np.float32)
        adata = ad.AnnData(X=matrix, obs=obs, var=var)
        counts = dataset.counts.tocsr() if sp.issparse(dataset.counts) else np.asarray(dataset.counts, dtype=np.float32)
        adata.layers["counts"] = counts
        for key in ("X_pca", "X_emb"):
            if key in dataset.obsm:
                adata.obsm[key] = np.asarray(dataset.obsm[key], dtype=np.float32)
        adata.uns["dsh_singlecell"] = json.dumps(
            {
                "plugin_version": PLUGIN_VERSION,
                "dataset_id": dataset.id,
                "dataset_name": dataset.name,
                "source": dataset.source,
                "source_path": dataset.path,
                "steps": list(dataset.steps),
                "params": _jsonable(dataset.params),
                "engine": _engine_info(),
                "layout": dataset.layout,
                "exported_at": _now_iso(),
            },
            ensure_ascii=False,
        )
        adata.write_h5ad(str(path))
    except EngineError:
        raise
    except Exception as error:
        raise EngineError("bad_request", "导出 h5ad 失败：%s" % error)


def _markers_groups_for_export(dataset, top_n=10):
    if "clusters" not in dataset.obs.columns:
        raise EngineError("prerequisite", "导出 markers_csv 需要先完成 cluster 步骤")
    cached = dataset.markers
    if isinstance(cached, dict) and cached.get("groups"):
        return cached["groups"]
    return _ensure_cluster_markers(dataset, max(1, int(top_n)))


def _export_markers_csv(dataset, path):
    groups = _markers_groups_for_export(dataset)
    rows = []
    for group in groups:
        for record in group.get("genes", []):
            rows.append(
                {
                    "cluster": str(group.get("id")),
                    "gene": str(record.get("gene")),
                    "score": float(record.get("score", 0.0)),
                    "logfc": float(record.get("logfc", 0.0)),
                    "pval_adj": float(record.get("pval_adj", 1.0)),
                    "pct_in": float(record.get("pct_in", 0.0)),
                    "pct_out": float(record.get("pct_out", 0.0)),
                }
            )
    pd.DataFrame(rows).to_csv(path, index=False)


def _export_obs_csv(dataset, path):
    obs = dataset.obs.copy()
    obs.index = _unique_index(obs.index, "cell")
    obs.to_csv(path)


def _export_summary_json(dataset, path):
    cluster_summary = None
    categories, codes = dataset.cluster_codes()
    if categories is not None and codes is not None:
        sizes = np.bincount(codes, minlength=len(categories))
        cluster_summary = {
            "n_clusters": int(len(categories)),
            "sizes": {str(label): int(sizes[index]) for index, label in enumerate(categories)},
        }
    payload = {
        "dataset": {
            "id": dataset.id,
            "name": dataset.name,
            "path": dataset.path,
            "source": dataset.source,
            "n_obs": dataset.n_obs,
            "n_vars": dataset.n_vars,
            "steps": list(dataset.steps),
            "params": _jsonable(dataset.params),
            "warnings": list(dataset.warnings),
        },
        "filter": _jsonable(dataset.filter_info),
        "pca": {
            "n_pcs": dataset.uns.get("pca_n_pcs"),
            "variance_ratio": dataset.uns.get("pca_variance_ratio"),
            "solver": dataset.uns.get("pca_solver"),
        },
        "layout": dataset.layout,
        "clusters": cluster_summary,
        "cluster": _jsonable(dataset.uns.get("cluster")),
        "engine": _engine_info(),
        "generated_at": _now_iso(),
    }
    Path(path).write_text(json.dumps(_jsonable(payload), ensure_ascii=False, indent=2), encoding="utf-8")


def _build_protocol(dataset):
    manifest = build_repro_manifest(dataset)
    return {
        "protocol_version": 2,
        "plugin": {"name": "dsh-singlecell", "version": PLUGIN_VERSION},
        "generated_at": _now_iso(),
        "dataset": {
            "id": dataset.id,
            "name": dataset.name,
            "source": {
                "kind": dataset.source,
                "path": dataset.path,
                "demo_generator": "test/make_sample.py" if dataset.source == "demo" else None,
            },
            "version": int(dataset.version),
            "revision": int(dataset.revision),
            "n_obs": dataset.n_obs,
            "n_vars": dataset.n_vars,
            "obs_columns": [column["key"] for column in _obs_columns(dataset)],
            "steps": list(dataset.steps),
            "params": _jsonable(dataset.params),
            "matrix_kind": getattr(dataset, "matrix_kind", None),
            "counts_source": getattr(dataset, "counts_source", None),
            "counts_verified": bool(getattr(dataset, "counts_verified", True)),
            "warnings": list(dataset.warnings),
        },
        "engine": _engine_info(),
        "seed": SEED,
        "software": dict(PACKAGE_VERSIONS),
        "history": _jsonable(dataset.history),
        # EXPORT-01：Methods 与 protocol 共用同一份结构化执行记录
        "execution_record": _jsonable(execution_records(dataset)),
        # EXPORT-02：复现清单摘要（完整清单嵌入导出的 replay 脚本）
        "replay": {
            "manifest_version": REPLAY_MANIFEST_VERSION,
            "input": manifest.get("input"),
            "n_operations": len(manifest.get("operations") or []),
            "not_replayable": manifest.get("not_replayable"),
            "expectations": {
                "structural": (manifest.get("expectations") or {}).get("structural"),
                "hashes": (manifest.get("expectations") or {}).get("hashes"),
            },
            "reference_migration": manifest.get("reference_migration"),
        },
        "figures": _jsonable(
            [
                {
                    "kind": figure.get("kind"),
                    "file": figure.get("file"),
                    "width": figure.get("width"),
                    "height": figure.get("height"),
                    "bytes": figure.get("bytes"),
                    "params": figure.get("params"),
                }
                for figure in dataset.figures
            ]
        ),
    }


def _export_protocol(dataset, path):
    Path(path).write_text(json.dumps(_jsonable(_build_protocol(dataset)), ensure_ascii=False, indent=2), encoding="utf-8")


def _params_for_methods(dataset, step):
    params = dataset.params.get(step)
    if isinstance(params, dict) and params:
        return params
    for entry in reversed(dataset.history):
        if entry.get("step") == step and entry.get("ok"):
            return entry.get("params") or {}
    return {}


def _build_methods(dataset):
    """Methods 从结构化执行记录生成（EXPORT-01）：只描述真实成功、当前版本有效的步骤。"""
    records = execution_records(dataset)
    return "\n".join(_methods_lines(dataset, records))

def _fmt_version(value):
    return str(value) if value not in (None, "") else "未记录"


def _methods_software_text(software):
    software = software or {}
    return (
        "Python %s; NumPy %s; SciPy %s; scikit-learn %s; pandas %s; AnnData %s; matplotlib %s; "
        "scanpy %s; umap-learn %s; leidenalg %s"
        % (
            _fmt_version(software.get("python")),
            _fmt_version(software.get("numpy")),
            _fmt_version(software.get("scipy")),
            _fmt_version(software.get("sklearn")),
            _fmt_version(software.get("pandas")),
            _fmt_version(software.get("anndata")),
            _fmt_version(software.get("matplotlib")),
            _fmt_version(software.get("scanpy")),
            _fmt_version(software.get("umap")),
            _fmt_version(software.get("leidenalg")),
        )
    )


def _methods_params_text(record, skip=()):
    params = record.get("params_effective") or {}
    parts = []
    for key in sorted(params):
        if key in skip or key == "engine":
            continue
        value = params[key]
        if value in (None, "", {}, []):
            continue
        if isinstance(value, bool):
            text = "true" if value else "false"
        elif isinstance(value, (int, float)):
            text = _fmt_num(value)
        elif isinstance(value, (dict, list)):
            text = json.dumps(_jsonable(value), ensure_ascii=False)
        else:
            text = str(value)
        parts.append("%s = %s" % (key, text))
    return ", ".join(parts)


def _method_record_suffix(record):
    parts = []
    requested = record.get("requested_method")
    if requested not in (None, ""):
        parts.append("requested method: `%s`" % requested)
    parts.append("actual method: `%s`" % (record.get("actual_method") or "未记录"))
    if record.get("method_provenance"):
        parts.append("implementation: %s" % record.get("method_provenance"))
    fallback = record.get("fallback_reason")
    parts.append(
        "fallback reason: %s" % (fallback if fallback else "none (the requested method was executed)")
    )
    parts.append("seed = %s" % (record.get("seed") if record.get("seed") is not None else "未记录"))
    software = record.get("software") or {}
    parts.append(
        "software: dsh-singlecell v%s, Python %s, NumPy %s, scikit-learn %s, scipy %s"
        % (
            _fmt_version(software.get("dsh-singlecell") or PLUGIN_VERSION),
            _fmt_version(software.get("python")),
            _fmt_version(software.get("numpy")),
            _fmt_version(software.get("sklearn")),
            _fmt_version(software.get("scipy")),
        )
    )
    return "; ".join(parts)


def _methods_qc_section(dataset, step_records):
    qc_record = step_records.get("qc")
    filter_record = step_records.get("filter")
    if qc_record is None and filter_record is None:
        return []
    lines = ["### Quality control and filtering", ""]
    params = dict(((filter_record or {}).get("params_effective") or (qc_record or {}).get("params_effective") or {}))
    parts = []
    if params.get("min_genes") not in (None, ""):
        parts.append("at least %s detected genes per cell" % _fmt_num(params.get("min_genes")))
    if params.get("max_genes") not in (None, ""):
        parts.append("at most %s detected genes per cell" % _fmt_num(params.get("max_genes")))
    if params.get("min_counts") not in (None, ""):
        parts.append("at least %s total counts" % _fmt_num(params.get("min_counts")))
    if params.get("max_pct_mt") not in (None, ""):
        parts.append("mitochondrial counts below %s%% of total counts" % _fmt_num(params.get("max_pct_mt")))
    sentence = ""
    if qc_record is not None:
        sentence += (
            "QC metrics (detected genes, total counts and mitochondrial fraction per cell) were computed "
            "on the raw count matrix with the builtin implementation (actual method: `qc_metrics(builtin)`; "
            "requested method: `builtin`; fallback reason: none; seed = %s). "
            % (qc_record.get("seed") if qc_record.get("seed") is not None else "未记录")
        )
    sentence += "Cells were retained when they had " + ", ".join(parts) + "."
    if params.get("min_cells") not in (None, ""):
        sentence += " Genes detected in fewer than %s cells were discarded." % _fmt_num(params.get("min_cells"))
    if filter_record is not None:
        sentence += (
            " Thresholds were applied by the builtin filter (actual method: `threshold_filter(builtin)`; "
            "requested method: `threshold_filter`; fallback reason: none; seed = %s)."
            % (filter_record.get("seed") if filter_record.get("seed") is not None else "未记录")
        )
    filter_info = dataset.filter_info or {}
    if filter_info:
        sentence += " This retained %s of %s cells and %s of %s genes." % (
            _fmt_num(filter_info.get("n_obs_after")), _fmt_num(filter_info.get("n_obs_before")),
            _fmt_num(filter_info.get("n_vars_after")), _fmt_num(filter_info.get("n_vars_before")),
        )
    lines.append(sentence)
    return lines


def _methods_step_section(dataset, step, record):
    """单个流程步骤的 Methods 段落（内容全部来自该步骤的执行记录）。"""
    params = record.get("params_effective") or {}
    titles = {
        "normalize": "Normalisation",
        "hvg": "Highly variable genes",
        "pca": "Principal component analysis",
        "neighbors": "Nearest-neighbour graph",
        "cluster": "Clustering",
        "embed": "Two-dimensional embedding",
        "markers": "Marker genes",
        "annotate": "Cell-type annotation",
    }
    title = titles.get(step)
    if title is None:
        return []
    input_scale = record.get("input_scale") or "未记录"
    output_scale = record.get("output_scale") or "未记录"
    lines = ["### " + title, ""]
    if step == "normalize":
        sentence = (
            "Counts were scaled to a target library size of %s per cell and log1p-transformed "
            "(log(1 + x)); input matrix scale: %s; output matrix scale: %s. %s."
            % (_fmt_num(params.get("target_sum", 10000)), input_scale, output_scale,
               _method_record_suffix(record))
        )
    elif step == "hvg":
        selected = (dataset.uns.get("hvg") or {}).get("n_selected")
        sentence = (
            "Highly variable genes were selected with the Seurat-style mean-binned normalised "
            "dispersion method (20 equal-frequency mean bins; n_top_genes = %s%s; input matrix scale: %s). %s."
            % (_fmt_num(params.get("n_top_genes")),
               "; %s genes selected" % _fmt_num(selected) if selected is not None else "",
               input_scale, _method_record_suffix(record))
        )
    elif step == "pca":
        n_pcs = dataset.uns.get("pca_n_pcs") or params.get("n_pcs")
        pc1 = None
        ratios = dataset.uns.get("pca_variance_ratio") or []
        if ratios:
            pc1 = float(ratios[0]) * 100.0
        sentence = (
            "Principal component analysis was run on the scaled highly variable genes with %s components%s; "
            "input matrix scale: %s. %s."
            % (_fmt_num(n_pcs),
               "; PC1 explained %.1f%% of the variance" % pc1 if pc1 is not None else "",
               input_scale, _method_record_suffix(record))
        )
    elif step == "neighbors":
        neighbor_meta = dataset.uns.get("neighbors") or {}
        sentence = (
            "A k-nearest-neighbour graph (k = %s, distance metric = %s) was built in PCA space and "
            "symmetrised; input matrix scale: %s. %s."
            % (_fmt_num(neighbor_meta.get("n_neighbors", params.get("n_neighbors"))),
               neighbor_meta.get("metric", params.get("metric")), input_scale,
               _method_record_suffix(record))
        )
    elif step == "cluster":
        cluster_meta = dataset.uns.get("cluster") or {}
        n_clusters = cluster_meta.get("n_clusters")
        if n_clusters is None:
            categories, _codes = dataset.cluster_codes()
            n_clusters = len(categories) if categories is not None else None
        actual = str(record.get("actual_method") or "")
        if actual.lower().startswith("louvain"):
            headline = "Cells were clustered with **Louvain community detection (builtin pure-Python implementation)**"
        elif actual.lower().startswith("leiden"):
            headline = "Cells were clustered with **Leiden community detection (scanpy implementation)**"
        else:
            headline = "Cells were clustered with **%s**" % actual
        sentence = (
            "%s on the k-nearest-neighbour graph. %s; resolution = %s; random seed = %s; "
            "input matrix scale: %s. %s clusters were obtained."
            % (headline,
               _method_record_suffix(record),
               _fmt_num(cluster_meta.get("resolution", params.get("resolution"))),
               record.get("seed") if record.get("seed") is not None else "未记录",
               input_scale,
               _fmt_num(n_clusters) if n_clusters is not None else "未记录")
        )
    elif step == "embed":
        embed_meta = dataset.uns.get("embed") or {}
        layout = str(embed_meta.get("layout") or dataset.layout or params.get("layout") or "")
        if layout == "tsne":
            headline = "Two-dimensional visualisation was produced with t-SNE (perplexity = %s, PCA initialisation)" % _fmt_num(
                embed_meta.get("perplexity") or params.get("perplexity"))
        elif layout == "umap":
            headline = "Two-dimensional visualisation was produced with UMAP"
        elif layout.startswith("umap"):
            headline = "Two-dimensional visualisation was produced with %s" % layout
        else:
            headline = "Two-dimensional visualisation used the first two principal components"
        sentence = "%s. %s." % (headline, _method_record_suffix(record))
    elif step == "markers":
        sentence = (
            "Cluster marker genes were identified with one-vs-rest Welch t-tests on log-normalised "
            "expression, with p-values adjusted by the Benjamini-Hochberg procedure (top_n = %s, "
            "groupby = %s). %s."
            % (_fmt_num(params.get("top_n", 10)), params.get("groupby", "clusters"),
               _method_record_suffix(record))
        )
    else:  # annotate
        mapping = params.get("mapping") or {}
        pairs = ", ".join("%s -> %s" % (key, value) for key, value in list(mapping.items())[:20])
        sentence = "Cluster labels were manually curated as: %s. %s." % (
            pairs if pairs else "（空映射，cell_type 暂时保持簇号）",
            _method_record_suffix(record),
        )
    lines.append(sentence)
    return lines


def _methods_provenance(dataset, records):
    lines = ["### Data import and provenance", ""]
    import_record = next((entry for entry in records if entry.get("kind") == "open"), None)
    source_format = getattr(dataset, "source_format", None) or dataset.source or "unknown"
    full_obs = getattr(dataset, "full_obs", None)
    full_n_obs = int(full_obs.shape[0]) if full_obs is not None else dataset.n_obs
    full_n_vars = int(len(getattr(dataset, "full_var_names", None) or dataset.var_names))
    lines.append(
        "- Source: `%s` (format: %s; import record: %s)."
        % (dataset.path or "demo generator (test/make_sample.py)", source_format,
           "recorded" if import_record is not None else "未记录")
    )
    lines.append(
        "- Matrix scale of the current X: %s; counts layer source: %s (counts_verified=%s)."
        % (_matrix_scale(dataset), getattr(dataset, "counts_source", None) or "未记录",
           bool(getattr(dataset, "counts_verified", False)))
    )
    lines.append(
        "- Input size: %s cells × %s genes; current size: %s cells × %s genes."
        % (_fmt_num(full_n_obs), _fmt_num(full_n_vars), _fmt_num(dataset.n_obs), _fmt_num(dataset.n_vars))
    )
    if import_record is not None:
        args = import_record.get("args") or {}
        if args:
            preview = json.dumps(_jsonable(args), ensure_ascii=False)
            if len(preview) > 240:
                preview = preview[:240] + "…（完整参数见 protocol/复现清单）"
            lines.append("- Import arguments recorded: `%s`." % preview)
    return lines


def _methods_unrecorded_section(dataset, step_records):
    missing = [step for step in STEP_ORDER if step in dataset.steps and step not in step_records]
    if not missing:
        return []
    import_state = (getattr(dataset, "uns", {}) or {}).get("import_state") or {}
    details = import_state.get("details") or {}
    lines = ["### Steps without an execution record", ""]
    for step in missing:
        hint = "该步骤来自导入数据或对象自带结果，本插件没有成功执行记录，不能补写具体算法。"
        if details.get("source_class"):
            hint += " Imported object class: %s." % details.get("source_class")
        lines.append("- %s: 未记录 (not recorded). %s" % (step, hint))
    lines.append("")
    return lines


def _methods_manual_sections(dataset, records):
    lines = []
    manual = [entry for entry in records if entry.get("kind") in ("command", "subset")]
    if manual:
        lines.extend(["### Manual metadata edits", ""])
        for entry in manual:
            replayable = "yes" if entry.get("replayable", True) else "no (%s)" % (
                entry.get("not_replayable_reason") or "只能作为记录/示例")
            columns = ", ".join(entry.get("obs_columns") or []) or "none"
            lines.append(
                "- %s: `%s` — %s; affected obs columns: %s; replayable: %s."
                % (entry.get("at"), entry.get("command"),
                   entry.get("actual_method") or "未记录", columns, replayable)
            )
        lines.append("")
    external = [entry for entry in records if entry.get("kind") in ("analysis", "external")]
    if external:
        lines.extend(["### Analyses and external tools", ""])
        for entry in external:
            method = entry.get("actual_method") or "未记录"
            if entry.get("kind") == "external":
                lines.append(
                    "- %s: external tool `%s` — method: %s; 未能重放（外部 R/Python 依赖），"
                    "相关结果只能作为记录/示例." % (entry.get("at"), entry.get("command"), method)
                )
            else:
                lines.append(
                    "- %s: analysis `%s` — method: %s; the recorded call is re-invoked on replay, "
                    "but analysis-specific numeric outputs are not value-checked."
                    % (entry.get("at"), (entry.get("args") or {}).get("id"), method)
                )
        lines.append("")
    return lines


def _methods_lines(dataset, records):
    step_records = {}
    for entry in records:
        if entry.get("kind") == "step" and entry.get("step"):
            step_records[entry.get("step")] = entry
    step_software = [entry for entry in records if entry.get("kind") == "step"]
    software = ((step_software[-1].get("software") if step_software else None) or _record_software())
    lines = ["## Methods", ""]
    lines.append(
        "Single-cell RNA-seq data (%s cells × %s genes in the current matrix; dataset version v%s, revision %s) "
        "were analysed with **dsh-singlecell** v%s (%s; %s)."
        % (
            _fmt_num(dataset.n_obs), _fmt_num(dataset.n_vars), dataset.version, dataset.revision,
            PLUGIN_VERSION,
            "builtin NumPy/SciPy/scikit-learn engine" if ENGINE_MODE != "scanpy"
            else "scanpy %s backend" % _fmt_version(PACKAGE_VERSIONS.get("scanpy")),
            _methods_software_text(software),
        )
    )
    lines.append("")
    lines.extend(_methods_provenance(dataset, records))
    lines.append("")
    qc_lines = _methods_qc_section(dataset, step_records)
    if qc_lines:
        lines.extend(qc_lines)
        lines.append("")
    for step in STEP_ORDER:
        if step in ("qc", "filter"):
            continue
        record = step_records.get(step)
        if record is None:
            continue
        section = _methods_step_section(dataset, step, record)
        if section:
            lines.extend(section)
            lines.append("")
    lines.extend(_methods_manual_sections(dataset, records))
    lines.extend(_methods_unrecorded_section(dataset, step_records))
    lines.extend(["### Software and reproducibility", ""])
    lines.append(
        "All random seeds were fixed at %s. %s." % (SEED, _methods_software_text(software))
    )
    lines.append("")
    lines.append(
        "This section was generated from the structured execution record of dataset version v%s "
        "(revision %s); only steps that were actually executed successfully and still belong to the "
        "current data version are described. Values that could not be recovered from the record are "
        "marked \"未记录\" instead of being inferred from interface labels." % (dataset.version, dataset.revision)
    )
    return lines


def _replay_snapshot(dataset):
    """导出时的期望结果快照；重放后用同一函数计算并逐项比较。"""
    categories, codes = dataset.cluster_codes()
    x_fp = _matrix_fingerprint(dataset.X)
    counts_fp = _matrix_fingerprint(dataset.counts)
    pca = dataset.obsm.get("X_pca")
    pca_fp = _matrix_fingerprint(pca)
    embed = dataset.obsm.get("X_emb")
    structural = {
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
        "steps": list(dataset.steps),
        "matrix_scale": _matrix_scale(dataset),
        "counts_source": getattr(dataset, "counts_source", None),
        "counts_verified": bool(getattr(dataset, "counts_verified", False)),
        "n_clusters": int(len(categories)) if categories is not None else None,
        "layout": dataset.layout,
        "pca_n_pcs": dataset.uns.get("pca_n_pcs"),
        "cluster_method": (dataset.uns.get("cluster") or {}).get("method"),
        "filter": (
            {
                "n_obs_before": dataset.filter_info.get("n_obs_before"),
                "n_obs_after": dataset.filter_info.get("n_obs_after"),
                "n_vars_before": dataset.filter_info.get("n_vars_before"),
                "n_vars_after": dataset.filter_info.get("n_vars_after"),
            }
            if dataset.filter_info else None
        ),
    }
    hashes = {"X": x_fp["sha256"] if x_fp else None, "counts": counts_fp["sha256"] if counts_fp else None}
    if categories is not None:
        hashes["clusters"] = _obs_digest(dataset.obs["clusters"])
    for column in _tracked_obs_columns(dataset):
        hashes["obs:" + column] = _obs_digest(dataset.obs[column])
    cluster_sizes = []
    if codes is not None:
        cluster_sizes = sorted(int(value) for value in np.bincount(codes, minlength=len(categories)))
    markers = {}
    if isinstance(dataset.markers, dict):
        for group in dataset.markers.get("groups") or []:
            markers[str(group.get("id"))] = [
                [
                    str(record.get("gene")),
                    round(float(record.get("score") or 0.0), 10),
                    round(float(record.get("logfc") or 0.0), 10),
                ]
                for record in (group.get("genes") or [])[:5]
            ]
    numbers = {
        "cluster_sizes": cluster_sizes,
        "pca_variance_ratio": [round(float(value), 12) for value in (dataset.uns.get("pca_variance_ratio") or [])],
        "markers": markers,
        "embed": (
            {
                "shape": [int(value) for value in np.asarray(embed).shape],
                "sum": float(np.nansum(embed)),
                "min": float(np.nanmin(embed)),
                "max": float(np.nanmax(embed)),
            }
            if embed is not None and np.asarray(embed).ndim == 2 else None
        ),
    }
    return {
        "structural": structural,
        "hashes": hashes,
        "numbers": numbers,
        "fingerprints": {"X": x_fp, "counts": counts_fp, "X_pca": pca_fp},
    }


def _check_close(actual, expected, rel=1e-6, absolute=1e-9):
    try:
        actual = float(actual)
        expected = float(expected)
    except Exception:  # noqa: BLE001
        return False
    if not (math.isfinite(actual) and math.isfinite(expected)):
        return False
    return abs(actual - expected) <= absolute + rel * max(abs(actual), abs(expected))


def compare_replay_snapshots(expected, actual, *, exact=True, skipped_steps=()):
    """比较记录期望与重放结果；exact=False 时对浮点/矩阵使用预定义容差。"""
    checks = []

    def add(identifier, ok, expected_value=None, actual_value=None, mode="exact", detail=""):
        checks.append({
            "id": identifier,
            "ok": bool(ok),
            "mode": mode,
            "expected": _jsonable(expected_value),
            "actual": _jsonable(actual_value),
            "detail": detail,
            "required": True,
        })
        return bool(ok)

    expected_structural = expected.get("structural") or {}
    actual_structural = actual.get("structural") or {}
    skipped = set(str(step) for step in skipped_steps or ())
    for key in ("n_obs", "n_vars", "matrix_scale", "counts_source", "counts_verified",
                "layout", "pca_n_pcs", "cluster_method", "filter"):
        add("structural:" + key, expected_structural.get(key) == actual_structural.get(key),
            expected_structural.get(key), actual_structural.get(key))
    expected_steps = [str(step) for step in expected_structural.get("steps") or [] if str(step) not in skipped]
    add("structural:steps", expected_steps == [str(step) for step in actual_structural.get("steps") or []],
        expected_steps, actual_structural.get("steps"),
        detail=("跳过不可重放步骤：%s" % ", ".join(sorted(skipped))) if skipped else "")
    if skipped:
        add("structural:n_clusters", True, expected_structural.get("n_clusters"),
            actual_structural.get("n_clusters"), mode="record_only",
            detail="聚类步骤不可重放，簇数仅作记录")
    else:
        add("structural:n_clusters", expected_structural.get("n_clusters") == actual_structural.get("n_clusters"),
            expected_structural.get("n_clusters"), actual_structural.get("n_clusters"))

    expected_hashes = expected.get("hashes") or {}
    actual_hashes = actual.get("hashes") or {}
    expected_fp = expected.get("fingerprints") or {}
    actual_fp = actual.get("fingerprints") or {}
    expected_numbers = expected.get("numbers") or {}
    actual_numbers = actual.get("numbers") or {}

    for name in ("X", "counts"):
        if exact:
            add("matrix:%s" % name, expected_hashes.get(name) == actual_hashes.get(name),
                expected_hashes.get(name), actual_hashes.get(name), mode="exact",
                detail="逐值（float32 指纹）比较")
        else:
            left = expected_fp.get(name) or {}
            right = actual_fp.get(name) or {}
            ok = (
                left.get("shape") == right.get("shape")
                and left.get("nnz") == right.get("nnz")
                and _check_close(right.get("sum"), left.get("sum"), rel=1e-6)
                and _check_close(right.get("min"), left.get("min"), rel=1e-6)
                and _check_close(right.get("max"), left.get("max"), rel=1e-6)
            )
            add("matrix:%s" % name, ok, left.get("sha256"), right.get("sha256"), mode="tolerance",
                detail="依赖版本不同：按 shape/nnz/sum/min/max 容差比较（rel=1e-6）")

    for key, value in expected_hashes.items():
        if not key.startswith("obs:"):
            continue
        actual_value = actual_hashes.get(key)
        if actual_value is None and skipped:
            add(key, True, value, actual_value, mode="record_only",
                detail="相关步骤不可重放，该 obs 列只作为记录")
        else:
            add(key, value == actual_value, value, actual_value)

    if "clusters" in expected_hashes:
        if exact and not skipped:
            add("clusters:labels", expected_hashes.get("clusters") == actual_hashes.get("clusters"),
                expected_hashes.get("clusters"), actual_hashes.get("clusters"), mode="exact",
                detail="Louvain 为确定性实现，标签逐值比较")
        elif skipped:
            add("clusters:labels", True, "recorded", None, mode="record_only",
                detail="聚类步骤不可重放，标签只作为记录")
        else:
            left_sizes = expected_numbers.get("cluster_sizes") or []
            right_sizes = actual_numbers.get("cluster_sizes") or []
            add("clusters:labels", left_sizes == right_sizes,
                left_sizes, right_sizes, mode="tolerance",
                detail="依赖版本不同：只比较簇大小分布")

    expected_variance = expected_numbers.get("pca_variance_ratio") or []
    if expected_variance:
        actual_variance = actual_numbers.get("pca_variance_ratio") or []
        if exact:
            ok = len(actual_variance) == len(expected_variance) and all(
                _check_close(actual_variance[index], expected_variance[index], rel=0.0, absolute=0.0)
                for index in range(len(expected_variance))
            )
            mode = "exact"
        else:
            ok = len(actual_variance) == len(expected_variance) and all(
                _check_close(actual_variance[index], expected_variance[index], rel=1e-6)
                for index in range(len(expected_variance))
            )
            mode = "tolerance"
        add("pca:variance_ratio", ok, len(expected_variance), len(actual_variance), mode=mode)

    expected_markers = expected_numbers.get("markers") or {}
    actual_markers = actual_numbers.get("markers") or {}
    if expected_markers:
        expected_genes = {group: [item[0] for item in items] for group, items in expected_markers.items()}
        actual_genes = {group: [item[0] for item in items] for group, items in actual_markers.items()}
        add("markers:genes", expected_genes == actual_genes, expected_genes, actual_genes)
        score_ok = True
        for group, items in expected_markers.items():
            actual_items = actual_markers.get(group) or []
            if len(actual_items) != len(items):
                score_ok = False
                break
            for left, right in zip(items, actual_items):
                rel = 0.0 if exact else 1e-4
                if not (_check_close(right[1], left[1], rel=rel) and _check_close(right[2], left[2], rel=rel)):
                    score_ok = False
        add("markers:scores", score_ok, "gene scores", "gene scores",
            mode="exact" if exact else "tolerance")

    expected_embed = expected_numbers.get("embed")
    if expected_embed:
        actual_embed = actual_numbers.get("embed") or {}
        ok = (
            expected_embed.get("shape") == actual_embed.get("shape")
            and _check_close(actual_embed.get("sum"), expected_embed.get("sum"), rel=1e-5)
            and _check_close(actual_embed.get("min"), expected_embed.get("min"), rel=1e-5)
            and _check_close(actual_embed.get("max"), expected_embed.get("max"), rel=1e-5)
        )
        add("embed:coordinates", ok, expected_embed.get("shape"), actual_embed.get("shape"),
            mode="tolerance", detail="t-SNE/UMAP 预定义容差：shape 精确，统计量 rel=1e-5")
    return checks


def _operation_payload(dataset, record):
    kind = str(record.get("kind") or "")
    base = {
        "kind": kind,
        "at": record.get("at"),
        "actual_method": record.get("actual_method"),
        "replayable": bool(record.get("replayable", True)),
        "not_replayable_reason": record.get("not_replayable_reason"),
    }
    if kind == "open":
        base.update({
            "command": "open",
            "args": record.get("args") or {},
            "source_format": record.get("source_format"),
            "input_n_obs": record.get("input_n_obs"),
            "input_n_vars": record.get("input_n_vars"),
        })
        return base
    if kind == "step":
        replayable, reason = _step_replayable(record)
        base.update({
            "step": record.get("step"),
            "requested_method": record.get("requested_method"),
            "actual_method": record.get("actual_method"),
            "method_provenance": record.get("method_provenance"),
            "fallback_reason": record.get("fallback_reason"),
            "params_requested": record.get("params_requested") or {},
            "params_effective": record.get("params_effective") or {},
            "replay_params": record.get("replay_params") or {},
            "input_scale": record.get("input_scale"),
            "output_scale": record.get("output_scale"),
            "seed": record.get("seed"),
            "software": record.get("software"),
            "n_obs_after": record.get("n_obs_after"),
            "n_vars_after": record.get("n_vars_after"),
            "replayable": replayable,
            "not_replayable_reason": reason,
        })
        return base
    if kind in ("command", "analysis"):
        args = dict(record.get("args") or {})
        base.update({
            "command": "analysis" if kind == "analysis" else record.get("command"),
            "args": args,
            "obs_columns": record.get("obs_columns") or [],
            "replayable": bool(record.get("replayable", True)),
            "not_replayable_reason": record.get("not_replayable_reason"),
            "method_provenance": record.get("method_provenance"),
        })
        if kind == "analysis":
            base["analysis_id"] = args.get("id")
        return base
    base.update({
        "command": record.get("command"),
        "args": record.get("args") or {},
    })
    return base


def _step_replayable(record):
    actual = str(record.get("actual_method") or "").lower()
    provenance = str(record.get("method_provenance") or "").lower()
    # 只有真正的外部 Seurat 实现才不可重放；hvg 的 seurat_binned_dispersion 是内置实现
    if actual.startswith("seurat::") or "seurat 6.0" in provenance or "seurat6" in provenance:
        return False, "需要 Seurat 6.0 + R 桥（外部依赖）；该步骤只能作为记录/示例"
    if actual.startswith("leiden"):
        return False, "需要 scanpy + leidenalg（外部依赖），不能保证同算法结果"
    return True, None


def _input_manifest(dataset, records):
    import_record = next((entry for entry in records if entry.get("kind") == "open"), None)
    args = dict((import_record or {}).get("args") or {})
    if not args and getattr(dataset, "path", None):
        args = {"path": dataset.path}
    replayable = bool((import_record or {}).get("replayable", True))
    reason = (import_record or {}).get("not_replayable_reason")
    kind = "path"
    if args.get("url"):
        kind = "url"
        replayable, reason = False, "URL 输入依赖网络与远端文件；该产物只能作为记录/示例"
    elif args.get("demo"):
        kind = "demo"
    elif args.get("paths") or args.get("samples"):
        kind = "merged"
    if str(getattr(dataset, "source", "")) == "subset":
        kind = "subset"
        replayable, reason = False, "子集数据集需要先重放父数据集；该产物只能作为记录/示例"
    files = []
    if kind != "url":
        for path in _input_files(dataset, args):
            size = None
            try:
                size = int(Path(path).stat().st_size)
            except Exception:  # noqa: BLE001
                size = None
            files.append({"path": path, "sha256": _sha256_file(path), "bytes": size})
    if not files and kind == "path":
        replayable, reason = False, "没有记录到本地输入文件；该产物只能作为记录/示例"
    full_obs = getattr(dataset, "full_obs", dataset.obs)
    full_var_names = list(getattr(dataset, "full_var_names", None) or dataset.var_names)
    cells = [str(value) for value in full_obs.index]
    genes = [str(value) for value in full_var_names]
    return {
        "kind": kind,
        "open_args": _jsonable(args),
        "format": getattr(dataset, "source_format", None) or dataset.source,
        "path": dataset.path,
        "sha256": (files[0].get("sha256") if files else None),
        "demo_generator": "test/make_sample.py" if str(getattr(dataset, "source", "")) == "demo" else None,
        "files": files,
        "n_obs": len(cells),
        "n_vars": len(genes),
        "current_n_obs": int(dataset.n_obs),
        "current_n_vars": int(dataset.n_vars),
        "cells": {"n": len(cells), "sha256": _sha256_json(cells),
                  "head": cells[:5], "tail": cells[-5:]},
        "genes": {"n": len(genes), "sha256": _sha256_json(genes),
                  "head": genes[:5], "tail": genes[-5:]},
        "obs_columns": [str(column) for column in dataset.obs.columns],
        "samples": _jsonable((dataset.uns or {}).get("samples") or []),
        "groups": _summary_groups(dataset),
        "group_column": (dataset.uns or {}).get("group_column"),
        "matrix_kind": getattr(dataset, "matrix_kind", None),
        "counts_source": getattr(dataset, "counts_source", None),
        "counts_verified": bool(getattr(dataset, "counts_verified", True)),
        "replayable": replayable,
        "not_replayable_reason": reason,
    }


def _reference_migration_payload(dataset, records):
    step_records = {
        entry.get("step"): entry for entry in records
        if entry.get("kind") == "step" and entry.get("step")
    }
    differences = []
    cluster = step_records.get("cluster")
    if cluster is None:
        differences.append("- clustering: 没有执行记录（未记录），参考迁移脚本不包含聚类步骤。")
    else:
        actual = str(cluster.get("actual_method") or "")
        resolution = (cluster.get("params_effective") or {}).get("resolution")
        if actual.lower().startswith("louvain"):
            differences.append(
                "- clustering: 实际执行 = louvain(builtin)（内置纯 Python 实现，resolution=%s）；"
                "参考脚本调用 sc.tl.louvain（scanpy 的另一实现），簇编号/边界可能不同；"
                "sc.tl.leiden（Leiden 算法）与 Louvain 是不同算法，不能用来替代本次结果。" % _fmt_num(resolution)
            )
        elif actual.lower().startswith("leiden"):
            differences.append(
                "- clustering: 实际执行 = Leiden（scanpy %s，resolution=%s）；参考脚本同样调用 "
                "sc.tl.leiden，但依赖版本不同时结果可能不同。"
                % (_fmt_version((cluster.get("software") or {}).get("scanpy")), _fmt_num(resolution))
            )
        else:
            differences.append(
                "- clustering: 实际执行 = %s（外部实现），scanpy 无法等价复现。" % (actual or "未记录")
            )
    hvg = step_records.get("hvg")
    if hvg is not None:
        differences.append(
            "- hvg: 实际执行 = seurat_binned_dispersion(builtin)；参考脚本的 "
            "sc.pp.highly_variable_genes(flavor='seurat') 是近似等价实现，选出的基因可能不同。"
        )
    pca = step_records.get("pca")
    if pca is not None:
        differences.append(
            "- pca: 实际执行 = %s（固定 seed）；参考脚本的 sc.tl.pca 使用同一类随机化 SVD，"
            "但数值实现/依赖版本不同时解释方差可能略有差异。" % (pca.get("actual_method") or "未记录")
        )
    embed = step_records.get("embed")
    if embed is not None:
        differences.append(
            "- embedding: 实际执行 = %s；参考脚本的 sc.tl.tsne / sc.tl.umap 是不同实现，"
            "二维坐标不做逐值等价承诺。" % (embed.get("actual_method") or "未记录")
        )
    markers = step_records.get("markers")
    if markers is not None:
        differences.append(
            "- markers: 实际执行 = Welch t-test(builtin)+BH；参考脚本的 "
            "sc.tl.rank_genes_groups(method='t-test') 只做近似等价，P 值/排名可能不同。"
        )
    differences.append(
        "- filtering/counts: 参考脚本在 AnnData 上原地过滤，counts 层随切片同步（不回填过滤前矩阵）；"
        "数值与本次一致的前提是过滤阈值、矩阵方向与输入完全相同。"
    )
    return {
        "title": "参考迁移脚本（Reference Scanpy migration script）",
        "note": (
            "参考迁移脚本只是迁移参考，不是重放路径；标称可重放的是导出文件前半部分的 "
            "REPLAY_MANIFEST 驱动（走同一套 builtin 实现并经逐值/容差校验）。"
        ),
        "differences": differences,
    }


def build_repro_manifest(dataset):
    """导出/重放共用的复现清单：输入清单+校验值、操作序列、实际参数、期望结果。"""
    records = execution_records(dataset)
    operations = [_operation_payload(dataset, record) for record in records]
    manifest = {
        "manifest_version": REPLAY_MANIFEST_VERSION,
        "generated_at": _now_iso(),
        "plugin": {"name": "dsh-singlecell", "version": PLUGIN_VERSION},
        "seed": SEED,
        "software": {
            "engine": {
                "sc_engine_path": str(Path(__file__).resolve()),
                "engine_mode": ENGINE_MODE,
                "python": PY_VERSION,
                "plugin_version": PLUGIN_VERSION,
                "seed": SEED,
                "packages": dict(PACKAGE_VERSIONS),
            },
            "fingerprint": _software_fingerprint(),
        },
        "dataset": {
            "id": dataset.id,
            "name": dataset.name,
            "version": int(dataset.version),
            "revision": int(dataset.revision),
            "steps": list(dataset.steps),
            "n_obs": int(dataset.n_obs),
            "n_vars": int(dataset.n_vars),
        },
        "input": _input_manifest(dataset, records),
        "operations": operations,
        "expectations": _replay_snapshot(dataset),
        "reference_migration": _reference_migration_payload(dataset, records),
        "not_replayable": [
            {
                "kind": op.get("kind"),
                "step": op.get("step"),
                "command": op.get("command"),
                "reason": op.get("not_replayable_reason"),
            }
            for op in operations if op.get("replayable") is False
        ],
    }
    return manifest


def _operation_file_problem(op):
    """command/analysis 操作引用的本地文件缺失或被改动时，返回不可重放原因。"""
    args = op.get("args") or {}
    if not args.get("path"):
        return None
    expected = args.get("sha256")
    path = str(args.get("path"))
    if not Path(path).exists():
        return "操作引用的文件不存在：%s" % path
    if expected:
        actual = _sha256_file(path)
        if actual and actual != expected:
            return "操作引用的文件校验值已变化：%s" % path
    return None


def _replay_finish(out_dir, result):
    path = Path(out_dir) / "replay_manifest.json"
    result["out_files"] = dict(result.get("out_files") or {})
    result["out_files"]["replay_manifest.json"] = str(path)
    try:
        path.write_text(json.dumps(_jsonable(result), ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception as error:  # noqa: BLE001
        result["warnings"].append("无法写出 replay_manifest.json：%s" % error)
    return result


def run_replay(manifest, out_dir, *, reference_script=None, progress=None):
    """按清单重放：校验输入 → 用真实命令入口重跑 → 用同一快照函数逐项比较。

    返回 {"status", "exit_code", "checks", "skipped_operations", ...}。
    status: verified / verified_with_tolerance / partial（记录/示例）/ failed / input_error /
            record_only。
    """
    out = Path(out_dir)
    try:
        out.mkdir(parents=True, exist_ok=True)
    except Exception as error:  # noqa: BLE001
        raise EngineError("bad_request", "无法创建重放输出目录（%s）：%s" % (out, error))
    say = progress if callable(progress) else (lambda message: None)
    result = {
        "status": "failed",
        "exit_code": 1,
        "message": "",
        "checks": [],
        "skipped_operations": [],
        "failed_operations": [],
        "warnings": [],
        "out_files": {},
        "software": _software_fingerprint(),
        "engine": {"lib": str(Path(__file__).resolve()), "plugin_version": PLUGIN_VERSION},
        "generated_at": _now_iso(),
        "input": {},
    }

    def finish(status, exit_code, message):
        result["status"] = status
        result["exit_code"] = exit_code
        result["message"] = message
        say("replay status: %s（%s）" % (status, message))
        return _replay_finish(out, result)

    manifest = manifest if isinstance(manifest, dict) else {}
    version = manifest.get("manifest_version")
    if version != REPLAY_MANIFEST_VERSION:
        return finish("record_only", 3,
                      "复现清单版本不匹配（%s），该导出产物只能作为记录/示例" % version)
    fingerprint = manifest.get("software", {}).get("fingerprint") or {}
    exact = _software_fingerprint() == fingerprint
    if not exact:
        result["warnings"].append(
            "依赖/插件版本与原执行不同：矩阵与统计量改用预定义容差比较，方法名与规模仍要求一致"
        )

    input_info = manifest.get("input") or {}
    result["input"] = {
        "kind": input_info.get("kind"),
        "format": input_info.get("format"),
        "path": input_info.get("path"),
        "sha256": [(item.get("path"), item.get("sha256")) for item in (input_info.get("files") or [])],
        "replayable": input_info.get("replayable", True),
    }
    if not input_info.get("replayable", True):
        return finish("record_only", 3, "该导出产物只能作为记录/示例：%s"
                      % (input_info.get("not_replayable_reason") or "输入不可重放"))
    for file_info in input_info.get("files") or []:
        path = file_info.get("path")
        if not path:
            continue
        if not Path(path).exists():
            return finish("input_error", 2,
                          "输入文件不存在，拒绝重放；该导出产物只能作为记录/示例：%s" % path)
        expected_sha = file_info.get("sha256")
        actual_sha = _sha256_file(path)
        if expected_sha and actual_sha and expected_sha != actual_sha:
            return finish("input_error", 2,
                          "输入文件 SHA-256 与导出记录不一致，拒绝重放；该导出产物只能作为记录/示例：%s" % path)

    operations = list(manifest.get("operations") or [])
    open_op = next((op for op in operations if op.get("kind") == "open"), None)
    if open_op is None:
        return finish("record_only", 3, "复现清单没有导入记录，无法定位原始输入；该产物只能作为记录/示例")
    if open_op.get("replayable") is False:
        return finish("record_only", 3, "输入不可重放：%s；该导出产物只能作为记录/示例"
                      % (open_op.get("not_replayable_reason") or "未记录原因"))
    say("open  输入：%s（%s）" % (input_info.get("path"), input_info.get("format")))
    recorded_paths = {
        os.path.realpath(str(item.get("path")))
        for item in (input_info.get("files") or []) if item.get("path")
    }
    for key, existing in list(DATASETS.items()):
        existing_path = getattr(existing, "path", None)
        if existing_path and os.path.realpath(str(existing_path)) in recorded_paths:
            DATASETS.pop(key, None)
    try:
        opened = cmd_open(dict(open_op.get("args") or {}))
    except EngineError as error:
        status = "input_error" if error.code in ("unsupported", "bad_request", "prerequisite") else "failed"
        return finish(status, 2 if status == "input_error" else 1,
                      "重新导入输入失败（%s）：%s；该导出产物只能作为记录/示例" % (error.code, error.message))
    dataset_id = str((opened.get("dataset") or {}).get("id") or "")
    if dataset_id not in DATASETS:
        return finish("failed", 1, "重新导入输入后没有可用数据集，无法重放")
    replayed = DATASETS[dataset_id]
    initial_checks = []
    expected_input_n_obs = int(input_info.get("n_obs") or 0)
    expected_input_n_vars = int(input_info.get("n_vars") or 0)
    if expected_input_n_obs:
        initial_checks.append({
            "id": "open:n_obs", "ok": int(replayed.n_obs) == expected_input_n_obs,
            "mode": "exact", "expected": expected_input_n_obs, "actual": int(replayed.n_obs),
            "required": True, "detail": "从同一输入恢复的细胞数",
        })
    if expected_input_n_vars:
        initial_checks.append({
            "id": "open:n_vars", "ok": int(replayed.n_vars) == expected_input_n_vars,
            "mode": "exact", "expected": expected_input_n_vars, "actual": int(replayed.n_vars),
            "required": True, "detail": "从同一输入恢复的基因数",
        })

    skipped_steps = []
    skipped_step_set = set()
    analysis_replayed = False
    for index, op in enumerate(operations, start=1):
        if op is open_op:
            continue
        kind = op.get("kind")
        label = op.get("step") or op.get("command") or kind
        # 依赖步骤不可重放时，下游步骤只能跳过（否则会以 prerequisite 失败伪装成真实失败）
        blocked_by = []
        if kind == "step":
            blocked_by = [dep for dep in (STEP_DEPS.get(op.get("step")) or []) if dep in skipped_step_set]
        elif kind == "command" and op.get("command") == "markers":
            blocked_by = [dep for dep in ("normalize", "cluster") if dep in skipped_step_set]
        if blocked_by:
            reason = "依赖步骤不可重放：%s" % ", ".join(blocked_by)
            result["skipped_operations"].append({
                "index": index, "kind": kind, "step": op.get("step"),
                "command": op.get("command"), "reason": reason,
            })
            if op.get("step"):
                skipped_steps.append(op.get("step"))
                skipped_step_set.add(op.get("step"))
            say("skip  %s：%s" % (label, reason))
            continue
        if not op.get("replayable", True):
            reason = op.get("not_replayable_reason") or "标记为不可重放"
            result["skipped_operations"].append({
                "index": index, "kind": kind, "step": op.get("step"),
                "command": op.get("command"), "reason": reason,
            })
            if op.get("step"):
                skipped_steps.append(op.get("step"))
                skipped_step_set.add(op.get("step"))
            say("skip  %s：%s" % (label, reason))
            continue
        file_problem = _operation_file_problem(op)
        if file_problem:
            result["skipped_operations"].append({
                "index": index, "kind": kind, "step": op.get("step"),
                "command": op.get("command"), "reason": file_problem,
            })
            if op.get("step"):
                skipped_steps.append(op.get("step"))
                skipped_step_set.add(op.get("step"))
            say("skip  %s：%s" % (label, file_problem))
            continue
        try:
            if kind == "step":
                params = dict(op.get("replay_params") or op.get("params_effective") or {})
                cmd_run({"dataset": dataset_id, "step": op.get("step"), "params": params})
            elif kind == "analysis":
                args = dict(op.get("args") or {})
                args["dataset"] = dataset_id
                cmd_analysis(args)
                analysis_replayed = True
            elif kind == "command":
                command = str(op.get("command") or "")
                if command not in REPLAYABLE_COMMANDS:
                    raise EngineError("unsupported", "命令 %s 不在可重放白名单里" % command)
                args = dict(op.get("args") or {})
                args["dataset"] = dataset_id
                COMMANDS[command](args)
            else:
                raise EngineError("unsupported", "未知操作类型：%s" % kind)
        except EngineError as error:
            if error.code == "unsupported":
                result["skipped_operations"].append({
                    "index": index, "kind": kind, "step": op.get("step"),
                    "command": op.get("command"),
                    "reason": "重放环境缺少依赖：%s" % error.message,
                })
                if op.get("step"):
                    skipped_steps.append(op.get("step"))
                    skipped_step_set.add(op.get("step"))
                say("skip  %s：%s" % (label, error.message))
                continue
            if error.code == "prerequisite" and skipped_step_set:
                result["skipped_operations"].append({
                    "index": index, "kind": kind, "step": op.get("step"),
                    "command": op.get("command"),
                    "reason": "依赖步骤不可重放：%s" % error.message,
                })
                if op.get("step"):
                    skipped_steps.append(op.get("step"))
                    skipped_step_set.add(op.get("step"))
                say("skip  %s：%s" % (label, error.message))
                continue
            result["failed_operations"].append({
                "index": index, "kind": kind, "step": op.get("step"),
                "command": op.get("command"), "code": error.code,
                "message": error.message,
            })
            say("fail  %s：%s" % (label, error.message))
            continue
        replayed = DATASETS[dataset_id]
        say("ok    %s（%d cells × %d genes）" % (label, replayed.n_obs, replayed.n_vars))
        after_checks = []
        expected_n_obs = op.get("n_obs_after")
        expected_n_vars = op.get("n_vars_after")
        if expected_n_obs is not None:
            after_checks.append({
                "id": "op:%02d:%s:n_obs" % (index, label),
                "ok": int(replayed.n_obs) == int(expected_n_obs),
                "mode": "exact", "expected": int(expected_n_obs),
                "actual": int(replayed.n_obs), "required": True,
            })
        if expected_n_vars is not None:
            after_checks.append({
                "id": "op:%02d:%s:n_vars" % (index, label),
                "ok": int(replayed.n_vars) == int(expected_n_vars),
                "mode": "exact", "expected": int(expected_n_vars),
                "actual": int(replayed.n_vars), "required": True,
            })
        if kind == "step":
            step_records = [
                entry for entry in replayed.execution_record
                if entry.get("kind") == "step" and entry.get("step") == op.get("step")
            ]
            actual_method = step_records[-1].get("actual_method") if step_records else None
            after_checks.append({
                "id": "op:%02d:%s:actual_method" % (index, label),
                "ok": str(actual_method) == str(op.get("actual_method")),
                "mode": "exact", "expected": op.get("actual_method"),
                "actual": actual_method, "required": True,
                "detail": "重放必须走记录里的实际方法（含回退后的方法）",
            })
        result["checks"].extend(after_checks)

    snapshot = _replay_snapshot(DATASETS[dataset_id])
    comparison = compare_replay_snapshots(
        manifest.get("expectations") or {}, snapshot, exact=exact, skipped_steps=skipped_steps,
    )
    result["checks"] = initial_checks + result["checks"] + comparison
    result["actual_snapshot"] = {
        "structural": snapshot.get("structural"),
        "hashes": snapshot.get("hashes"),
        "numbers": snapshot.get("numbers"),
    }

    out_files = {}
    try:
        h5ad_path = out / "reproduced.h5ad"
        _export_h5ad(DATASETS[dataset_id], h5ad_path)
        out_files["h5ad"] = str(h5ad_path)
    except Exception as error:  # noqa: BLE001
        result["warnings"].append("无法写出 reproduced.h5ad：%s" % error)
    try:
        obs_path = out / "reproduced_obs.csv"
        _export_obs_csv(DATASETS[dataset_id], obs_path)
        out_files["obs_csv"] = str(obs_path)
    except Exception as error:  # noqa: BLE001
        result["warnings"].append("无法写出 reproduced_obs.csv：%s" % error)
    try:
        summary_path = out / "reproduced_summary.json"
        _export_summary_json(DATASETS[dataset_id], summary_path)
        out_files["summary_json"] = str(summary_path)
    except Exception as error:  # noqa: BLE001
        result["warnings"].append("无法写出 reproduced_summary.json：%s" % error)
    if reference_script:
        try:
            reference_path = out / "reference_scanpy_migration.py"
            reference_path.write_text(str(reference_script), encoding="utf-8")
            out_files["reference_scanpy"] = str(reference_path)
        except Exception as error:  # noqa: BLE001
            result["warnings"].append("无法写出参考迁移脚本：%s" % error)
    if analysis_replayed:
        result["warnings"].append(
            "native analysis 已按记录重新调用，但分析输出未逐值校验（只校验数据快照）"
        )
    result["out_files"] = out_files

    required = [check for check in result["checks"] if check.get("required", True)]
    bad = [check for check in required if not check.get("ok")]
    if result["failed_operations"]:
        return finish("failed", 1, "有可重放操作执行失败（%d 项），不能声称已重放"
                      % len(result["failed_operations"]))
    if bad:
        return finish("failed", 1, "重放结果与记录不一致（%d 项检查失败），不能声称已重放：%s"
                      % (len(bad), ", ".join(str(check.get("id")) for check in bad[:5])))
    if result["skipped_operations"]:
        return finish("partial", 3,
                      "部分操作不可重放（%d 项）；该导出产物只能作为记录/示例"
                      % len(result["skipped_operations"]))
    if exact:
        return finish("verified", 0, "已从同一输入按记录重放；确定性步骤逐值一致")
    return finish("verified_with_tolerance", 0,
                  "已从同一输入按记录重放；依赖版本不同，确定性步骤按预定义容差比较通过")


def _scanpy_load_block(input_info):
    """参考迁移脚本的输入读取段（按格式给出可执行代码或明确迁移说明）。"""
    kind = str(input_info.get("kind") or "path")
    fmt = str(input_info.get("format") or "")
    path = input_info.get("path")
    lines = []
    if kind == "demo":
        generator = input_info.get("demo_generator") or "test/make_sample.py"
        lines.append("# demo 输入：与引擎相同的生成器 %s" % generator)
        lines.append("import sys, pathlib")
        lines.append("sys.path.insert(0, str(pathlib.Path.cwd() / 'test'))")
        lines.append("from make_sample import demo_adata")
        lines.append("adata = demo_adata()")
    elif kind == "merged":
        lines.append("# 合并导入（%s 个文件）：参考脚本只演示第一个文件的读取；")
        lines.append("# 完整合并请按导出清单 input.files 逐个读取后取基因交集。")
        first = (input_info.get("files") or [{}])[0].get("path") or path
        if fmt == "h5ad" or str(first).endswith(".h5ad"):
            lines.append("adata = sc.read_h5ad(%r)" % first)
        else:
            lines.append("frame = pd.read_csv(%r, index_col=0)" % first)
            lines.append("adata = AnnData(X=frame.T.to_numpy(dtype=np.float32))")
    elif str(path or "").endswith(".h5ad") or fmt == "h5ad":
        lines.append("adata = sc.read_h5ad(%r)" % path)
        lines.append("if 'counts' in adata.layers:")
        lines.append("    adata.layers['counts'] = adata.layers['counts'].copy()")
        lines.append("else:")
        lines.append("    adata.layers['counts'] = adata.X.copy()")
    elif fmt in ("matrix_text", "csv", "tsv", "txt") or str(path or "").endswith((".csv", ".tsv", ".txt")):
        lines.append("# 文本矩阵：按导出记录的方向读取（原执行由导入器自动判断并可能转置）")
        lines.append("frame = pd.read_csv(%r, sep=',', index_col=0)" % path)
        lines.append("adata = AnnData(X=frame.T.to_numpy(dtype=np.float32))")
        lines.append("adata.layers['counts'] = adata.X.copy()")
    elif fmt == "rds" or fmt == "rdata":
        lines.append("# RDS/RData 输入：请先用 R（SeuratDisk / sceasy）转换为 h5ad，")
        lines.append("# 再把下面的读取语句换成 sc.read_h5ad('converted.h5ad')。")
        lines.append("raise SystemExit('参考迁移脚本需要先把 RDS/RData 转换为 h5ad')")
    else:
        lines.append("# 该格式（%s）请先用 scanpy 对应读取函数载入；输入路径：%r" % (fmt, path))
    return lines


def _reference_scanpy_script(manifest):
    """生成「参考迁移脚本」：明确标注不是重放路径，并列出与原执行的差异。"""
    input_info = manifest.get("input") or {}
    operations = manifest.get("operations") or []
    reference = manifest.get("reference_migration") or {}
    lines = []
    lines.append("#!/usr/bin/env python3")
    lines.append("# -*- coding: utf-8 -*-")
    lines.append("# =====================================================================")
    lines.append("# REFERENCE SCANPY MIGRATION SCRIPT (参考迁移脚本) — NOT the replay path")
    lines.append("# 该文件由 dsh-singlecell v%s 导出，仅作为迁移到 scanpy 的参考实现；" % PLUGIN_VERSION)
    lines.append("# 它**不是**本次分析的可重放路径（标称可重放的是导出文件前半部分的 REPLAY_MANIFEST 驱动）。")
    lines.append("# 与原执行方式的差异（逐条核对，不要据此声称结果一致）：")
    for difference in reference.get("differences") or []:
        lines.append("#   " + str(difference))
    lines.append("# =====================================================================")
    lines.append("")
    lines.append("import pandas as pd")
    lines.append("import scanpy as sc")
    lines.append("from anndata import AnnData")
    lines.append("")
    lines.append("SEED = %d" % int(manifest.get("seed") or 0))
    lines.append("")
    lines.append("# ---------------------------------------------------------------------")
    lines.append("# 1. 读取同一输入")
    lines.extend(_scanpy_load_block(input_info))
    lines.append("")
    step_index = 1
    for op in operations:
        if op.get("kind") != "step":
            continue
        step = str(op.get("step"))
        params = dict(op.get("params_effective") or {})
        replay_params = dict(op.get("replay_params") or {})
        lines.append("# ---------------------------------------------------------------------")
        lines.append("# %d. %s（实际方法：%s）" % (step_index, step, op.get("actual_method") or "未记录"))
        step_index += 1
        if step == "qc":
            lines.append("adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')")
            lines.append("sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)")
        elif step == "filter":
            keep = []
            if params.get("min_genes") not in (None, ""):
                keep.append("(adata.obs['n_genes_by_counts'] >= %s)" % _fmt_num(params.get("min_genes")))
            if params.get("max_genes") not in (None, ""):
                keep.append("(adata.obs['n_genes_by_counts'] <= %s)" % _fmt_num(params.get("max_genes")))
            if params.get("min_counts") not in (None, ""):
                keep.append("(adata.obs['total_counts'] >= %s)" % _fmt_num(params.get("min_counts")))
            if params.get("max_pct_mt") not in (None, ""):
                keep.append("(adata.obs['pct_counts_mt'] < %s)" % _fmt_num(params.get("max_pct_mt")))
            if keep:
                lines.append("keep = " + " & ".join(keep))
                lines.append("adata = adata[keep].copy()")
            if params.get("min_cells") not in (None, ""):
                lines.append("sc.pp.filter_genes(adata, min_cells=%s)" % _fmt_num(params.get("min_cells")))
            lines.append("# counts 层随 AnnData 同步切片，不回填过滤前矩阵")
        elif step == "normalize":
            lines.append("sc.pp.normalize_total(adata, target_sum=%s)" % _fmt_num(params.get("target_sum", 10000)))
            lines.append("sc.pp.log1p(adata)")
        elif step == "hvg":
            lines.append("sc.pp.highly_variable_genes(adata, n_top_genes=%s, flavor='seurat')"
                         % _fmt_num(params.get("n_top_genes", 2000)))
            lines.append("adata = adata[:, adata.var['highly_variable']].copy()")
        elif step == "pca":
            lines.append("sc.tl.pca(adata, n_comps=%s, svd_solver='randomized', random_state=SEED)"
                         % _fmt_num(params.get("n_pcs", 30)))
        elif step == "neighbors":
            lines.append("sc.pp.neighbors(adata, n_neighbors=%s, metric=%r, random_state=SEED)"
                         % (_fmt_num(params.get("n_neighbors", 15)), str(params.get("metric", "cosine"))))
        elif step == "cluster":
            actual = str(op.get("actual_method") or "")
            resolution = replay_params.get("resolution", params.get("resolution", 0.2))
            if actual.lower().startswith("leiden"):
                lines.append("sc.tl.leiden(adata, resolution=%s, random_state=SEED, key_added='clusters')"
                             % _fmt_num(resolution))
            elif actual.lower().startswith("louvain"):
                lines.append("# 差异：实际执行 = louvain(builtin)（纯 Python）；scanpy 用 sc.tl.louvain 是另一实现，")
                lines.append("# 簇编号可能不同；sc.tl.leiden（Leiden 算法）不能用来替代本次 Louvain 结果。")
                lines.append("sc.tl.louvain(adata, resolution=%s, random_state=SEED, key_added='clusters')"
                             % _fmt_num(resolution))
            else:
                lines.append("# 实际聚类方法为 %r，scanpy 无等价调用；请保留原工具结果。" % actual)
            if actual.lower().startswith(("leiden", "louvain")):
                lines.append("adata.obs['clusters'] = adata.obs['clusters'].astype(str)")
        elif step == "embed":
            layout = str((op.get("replay_params") or {}).get("layout") or "")
            if layout == "tsne":
                lines.append("sc.tl.tsne(adata, perplexity=%s, random_state=SEED)"
                             % _fmt_num(params.get("perplexity", 30)))
            elif layout == "umap":
                lines.append("sc.tl.umap(adata, random_state=SEED)")
            else:
                lines.append("# 实际回退到 PCA 前两维")
                lines.append("adata.obsm['X_emb'] = adata.obsm['X_pca'][:, :2]")
        elif step == "markers":
            lines.append("sc.tl.rank_genes_groups(adata, %r, method='t-test', n_genes=%s)"
                         % (str(params.get("groupby", "clusters")), _fmt_num(params.get("top_n", 10))))
        elif step == "annotate":
            mapping = params.get("mapping") or {}
            lines.append("annotation_mapping = %r" % (mapping,))
            lines.append("adata.obs['cell_type'] = adata.obs['clusters'].map(annotation_mapping).fillna(adata.obs['clusters'])")
        else:
            lines.append("# 未知步骤：%s" % step)
        lines.append("")
    lines.append("# ---------------------------------------------------------------------")
    lines.append("# 保存参考结果（不是重放结果）")
    lines.append("adata.write_h5ad('reference_scanpy_result.h5ad')")
    text = "\n".join(lines) + "\n"
    # 嵌入到驱动文件里的 r\"\"\" 字符串不能出现三引号
    return text.replace('"""', "'''")


def _build_script(dataset):
    """导出可重放驱动（前半部分）+ 参考迁移脚本（后半部分，明确标注）。"""
    manifest = build_repro_manifest(dataset)
    reference = _reference_scanpy_script(manifest)
    manifest_json = json.dumps(_jsonable(manifest), ensure_ascii=True, indent=1, sort_keys=False)
    if '"""' in manifest_json:  # 正常 JSON 不会出现；防御性兜底
        manifest_json = manifest_json.replace('"""', '\\"\\"\\"')
    lines = []
    lines.append("#!/usr/bin/env python3")
    lines.append("# -*- coding: utf-8 -*-")
    lines.append('"""Replayable workflow export — dsh-singlecell v%s' % PLUGIN_VERSION)
    lines.append("")
    lines.append("Dataset: %s (%s)" % (dataset.name, dataset.id))
    lines.append("Recorded: %s ; dataset version v%s, revision %s" % (manifest.get("generated_at"), dataset.version, dataset.revision))
    lines.append("")
    lines.append("This file has two clearly separated parts.")
    lines.append("")
    lines.append("[1] REPLAY DRIVER (the top part of this file). It replays the recorded workflow")
    lines.append("    through the SAME implementation that produced the original results")
    lines.append("    (dsh-singlecell lib/sc_engine.py):")
    lines.append("      * input files are verified by SHA-256 before use;")
    lines.append("      * every recorded operation is executed in order with its effective parameters;")
    lines.append("      * the recorded ACTUAL method is replayed (e.g. builtin Louvain, never Leiden);")
    lines.append("      * counts are sliced together with the filtered AnnData, never restored from a")
    lines.append("        pre-filter matrix;")
    lines.append("      * deterministic steps are compared value-by-value with the recorded expectations.")
    lines.append("")
    lines.append("[2] REFERENCE SCANPY MIGRATION SCRIPT (参考迁移脚本, at the bottom of this file as the")
    lines.append("    REFERENCE_SCANPY_SCRIPT constant). It is a migration reference, NOT the replay path,")
    lines.append("    and it lists item-by-item differences from the actual execution (see the manifest")
    lines.append("    section 'reference_migration'). It is also written to the replay output directory.")
    lines.append("")
    lines.append("Usage:")
    lines.append("    python3 <this file> --out replay_out")
    lines.append("    python3 <this file> --engine-lib /path/to/plugins/dsh-singlecell/lib")
    lines.append("    python3 <this file> --dump-manifest      # 打印记录清单 JSON")
    lines.append("    python3 <this file> --no-reference       # 不写参考迁移脚本")
    lines.append("")
    lines.append("Exit codes: 0 verified / verified_with_tolerance, 1 mismatch or failure,")
    lines.append("            2 input error, 3 record-only or partial (记录/示例), 4 engine unavailable.")
    lines.append('"""')
    lines.append("")
    lines.append("import argparse")
    lines.append("import importlib.util")
    lines.append("import json")
    lines.append("import os")
    lines.append("import sys")
    lines.append("from pathlib import Path")
    lines.append("")
    lines.append("REPLAY_MANIFEST = json.loads(r'''")
    lines.append(manifest_json)
    lines.append("''')")
    lines.append("")
    lines.append("")
    lines.append("def load_engine(explicit=None):")
    lines.append('    """按记录路径 / --engine-lib / DSH_SINGLECELL_LIB 加载同一套真实执行入口。"""')
    lines.append("    candidates = []")
    lines.append("    if explicit:")
    lines.append("        candidates.append(Path(explicit))")
    lines.append("    env = os.environ.get('DSH_SINGLECELL_LIB')")
    lines.append("    if env:")
    lines.append("        candidates.append(Path(env))")
    lines.append("    recorded = (REPLAY_MANIFEST.get('software') or {}).get('engine') or {}")
    lines.append("    if recorded.get('sc_engine_path'):")
    lines.append("        candidates.append(Path(recorded['sc_engine_path']))")
    lines.append("    for candidate in candidates:")
    lines.append("        path = candidate / 'sc_engine.py' if candidate.is_dir() else candidate")
    lines.append("        if path.exists():")
    lines.append("            spec = importlib.util.spec_from_file_location('dsh_singlecell_sc_engine', str(path))")
    lines.append("            module = importlib.util.module_from_spec(spec)")
    lines.append("            sys.modules['dsh_singlecell_sc_engine'] = module")
    lines.append("            spec.loader.exec_module(module)")
    lines.append("            return module")
    lines.append("    raise SystemExit('engine_unavailable: 找不到 dsh-singlecell lib/sc_engine.py；请用 --engine-lib 或环境变量 DSH_SINGLECELL_LIB 指定。没有执行入口时该导出物只能作为记录。')")
    lines.append("")
    lines.append("")
    lines.append("REFERENCE_SCANPY_SCRIPT = r'''")
    lines.append(reference.rstrip("\n"))
    lines.append("'''")
    lines.append("")
    lines.append("")
    lines.append("def main(argv=None):")
    lines.append("    parser = argparse.ArgumentParser(description='dsh-singlecell replay driver')")
    lines.append("    parser.add_argument('--out', default='replay_out')")
    lines.append("    parser.add_argument('--engine-lib', default=None)")
    lines.append("    parser.add_argument('--dump-manifest', action='store_true')")
    lines.append("    parser.add_argument('--no-reference', action='store_true')")
    lines.append("    args = parser.parse_args(argv)")
    lines.append("    if args.dump_manifest:")
    lines.append("        print(json.dumps(REPLAY_MANIFEST, ensure_ascii=False, indent=2))")
    lines.append("        return 0")
    lines.append("    os.environ.setdefault('MPLBACKEND', 'Agg')")
    lines.append("    engine_module = load_engine(args.engine_lib)")
    lines.append("    result = engine_module.run_replay(")
    lines.append("        REPLAY_MANIFEST, args.out,")
    lines.append("        reference_script=None if args.no_reference else REFERENCE_SCANPY_SCRIPT,")
    lines.append("        progress=lambda line: print(line, flush=True),")
    lines.append("    )")
    lines.append("    print('replay status: %s' % result.get('status'))")
    lines.append("    print('message: %s' % result.get('message'))")
    lines.append("    checks = result.get('checks') or []")
    lines.append("    failed = [item for item in checks if item.get('required', True) and not item.get('ok')]")
    lines.append("    print('checks: %d total, %d failed' % (len(checks), len(failed)))")
    lines.append("    for item in failed[:10]:")
    lines.append("        print('  FAIL %-42s expected=%r actual=%r' % (item.get('id'), item.get('expected'), item.get('actual')))")
    lines.append("    for item in result.get('skipped_operations') or []:")
    lines.append("        print('  SKIPPED %-30s %s' % (item.get('step') or item.get('command'), item.get('reason')))")
    lines.append("    print('outputs: %s' % json.dumps(result.get('out_files') or {}, ensure_ascii=False))")
    lines.append("    return int(result.get('exit_code', 1))")
    lines.append("")
    lines.append("")
    lines.append("if __name__ == '__main__':")
    lines.append("    raise SystemExit(main())")
    lines.append("")
    return "\n".join(lines)


def _export_methods(dataset, path):
    Path(path).write_text(_build_methods(dataset), encoding="utf-8")


def _build_report(dataset, report_path):
    parent = Path(report_path).parent
    lines = []
    lines.append("# 单细胞分析报告 · %s" % dataset.name)
    lines.append("")
    lines.append("- **数据集**：`%s`（%s）" % (dataset.id, dataset.name))
    if dataset.path:
        lines.append("- **来源**：`%s`（%s）" % (dataset.path, dataset.source))
    else:
        lines.append("- **来源**：内置 demo 生成器（test/make_sample.py）")
    lines.append("- **规模**：%d 个细胞 × %d 个基因" % (dataset.n_obs, dataset.n_vars))
    lines.append(
        "- **引擎**：%s（Python %s；dsh-singlecell v%s）"
        % ("scanpy" if ENGINE_MODE == "scanpy" else "builtin", PY_VERSION, PLUGIN_VERSION)
    )
    lines.append("- **生成时间**：%s" % _now_iso())
    lines.append("")

    lines.append("## 数据摘要")
    lines.append("")
    lines.append("| 项目 | 值 |")
    lines.append("| --- | --- |")
    lines.append("| 细胞数 | %d |" % dataset.n_obs)
    lines.append("| 基因数 | %d |" % dataset.n_vars)
    lines.append("| 已完成步骤 | %s |" % (" → ".join(dataset.steps) if dataset.steps else "无"))
    filter_info = dataset.filter_info or {}
    if filter_info:
        lines.append(
            "| 过滤保留 | %s / %s 细胞，%s / %s 基因 |"
            % (
                _fmt_num(filter_info.get("n_obs_after")), _fmt_num(filter_info.get("n_obs_before")),
                _fmt_num(filter_info.get("n_vars_after")), _fmt_num(filter_info.get("n_vars_before")),
            )
        )
    categories, codes = dataset.cluster_codes()
    if categories is not None and codes is not None:
        lines.append("| 聚类 | %d 个簇 |" % len(categories))
    lines.append("")

    methods = _build_methods(dataset)
    lines.append(methods)
    lines.append("")

    # 图表：生成并复制到报告同目录
    lines.append("## 图表")
    lines.append("")
    figure_specs = [
        ("embedding", {"color_by": {"kind": "cluster", "key": "clusters"}}, "聚类嵌入图"),
        ("qc_violin", {"metrics": list(QC_METRICS)}, "QC 指标分布"),
    ]
    figure_added = False
    for kind, params, caption in figure_specs:
        try:
            info = _make_plot(dataset, kind, params)
        except EngineError:
            continue
        source = Path(info["file"])
        try:
            shutil.copyfile(str(source), str(parent / source.name))
            relative = source.name
        except Exception:
            relative = str(source)
        lines.append("![%s](%s)" % (caption, relative))
        lines.append("")
        figure_added = True
    if not figure_added:
        lines.append("（尚未生成图：需要先完成 qc / embed 步骤）")
        lines.append("")

    # 标记基因表
    lines.append("## 各簇 top markers")
    lines.append("")
    try:
        groups = _markers_groups_for_export(dataset, 5)
    except EngineError:
        groups = []
    if groups:
        lines.append("| cluster | gene | logFC | adj.P | pct_in |")
        lines.append("| --- | --- | --- | --- | --- |")
        for group in groups:
            for record in group.get("genes", [])[:5]:
                lines.append(
                    "| %s | %s | %s | %s | %s |"
                    % (
                        group.get("id"),
                        record.get("gene"),
                        _fmt_num(record.get("logfc")),
                        _fmt_num(record.get("pval_adj")),
                        _fmt_num(record.get("pct_in")),
                    )
                )
    else:
        lines.append("（尚未计算标记基因）")
    lines.append("")

    # 自动注释建议
    lines.append("## 自动注释建议")
    lines.append("")
    try:
        auto = _autoannotate(dataset, "clusters")
        suggestions = auto.get("suggestions") or []
    except EngineError:
        suggestions = []
    if suggestions:
        lines.append("参考面板：`%s`（仅建议，未写入 obs）" % auto.get("panel", PANEL_ID))
        lines.append("")
        lines.append("| cluster | 建议类型 | score | margin | evidence |")
        lines.append("| --- | --- | --- | --- | --- |")
        for item in suggestions:
            lines.append(
                "| %s | %s | %s | %s | %s |"
                % (
                    item.get("cluster"),
                    item.get("label"),
                    _fmt_num(item.get("score")),
                    _fmt_num(item.get("margin")),
                    ", ".join(item.get("evidence") or []),
                )
            )
    else:
        lines.append("（尚未聚类或参考面板基因不足，无法自动注释）")
    lines.append("")

    # 自检
    audit = _audit(dataset)
    lines.append("## 自检 Audit")
    lines.append("")
    lines.append("**score = %d / 100（%s）**" % (audit["score"], audit["level"]))
    lines.append("")
    for check in audit["checks"]:
        marker = {"ok": "OK", "warn": "WARN", "fail": "FAIL"}.get(check["level"], check["level"])
        hint = ("；建议：" + check["hint"]) if check.get("hint") else ""
        lines.append("- [%s] **%s**：%s%s" % (marker, check["title"], check["detail"], hint))
    lines.append("")

    # 历史
    lines.append("## 分析历史")
    lines.append("")
    if dataset.history:
        lines.append("| 时间 | 步骤 | 参数 | 结果摘要 | 耗时 |")
        lines.append("| --- | --- | --- | --- | --- |")
        for entry in dataset.history:
            params_text = json.dumps(_jsonable(entry.get("params") or {}), ensure_ascii=False)
            if len(params_text) > 80:
                params_text = params_text[:77] + "..."
            lines.append(
                "| %s | %s | `%s` | %s | %s ms |"
                % (
                    entry.get("at", ""),
                    entry.get("step", ""),
                    params_text,
                    entry.get("summary", ""),
                    _fmt_num(entry.get("ms", 0)),
                )
            )
    else:
        lines.append("（还没有分析记录）")
    lines.append("")
    return "\n".join(lines)


def _export_report(dataset, path):
    content = _build_report(dataset, path)
    Path(path).write_text(content, encoding="utf-8")


def _export_script(dataset, path):
    Path(path).write_text(_build_script(dataset), encoding="utf-8")


EXPORT_WRITERS = {
    "h5ad": _export_h5ad,
    "markers_csv": _export_markers_csv,
    "obs_csv": _export_obs_csv,
    "summary_json": _export_summary_json,
    "protocol": _export_protocol,
    "methods": _export_methods,
    "report": _export_report,
    "script": _export_script,
}


def _export(dataset, what, user_path):
    if what not in EXPORT_WRITERS:
        raise EngineError(
            "bad_request", "未知导出格式：%s（可用：%s）" % (what, ", ".join(EXPORT_WHATS))
        )
    path = _resolve_export_path(dataset, what, user_path)
    EXPORT_WRITERS[what](dataset, path)
    if not path.exists():
        raise EngineError("internal", "导出失败：文件没有写成（%s）" % path)
    return {"file": str(path), "bytes": int(path.stat().st_size), "what": what}


# ===========================================================================
# 15. 协议命令实现
# ===========================================================================


def cmd_ping(args):
    packages = {
        "numpy": PACKAGE_VERSIONS.get("numpy"),
        "pandas": PACKAGE_VERSIONS.get("pandas"),
        "scipy": PACKAGE_VERSIONS.get("scipy"),
        "sklearn": PACKAGE_VERSIONS.get("sklearn"),
        "anndata": PACKAGE_VERSIONS.get("anndata"),
        "scanpy": PACKAGE_VERSIONS.get("scanpy"),
        "matplotlib": PACKAGE_VERSIONS.get("matplotlib"),
    }
    return {
        "python": PY_VERSION,
        "platform": sys.platform,
        "mode": ENGINE_MODE,
        "packages": packages,
        "layouts": list(LAYOUTS),
        "warnings": list(STARTUP_WARNINGS),
    }


def _open_payload(dataset):
    summary = _dataset_summary(dataset)
    payload = dict(summary)
    # Host（lib/impl.mjs）读取 result.dataset；同时顶层也保留 §2.3 的 dataset 字段
    payload["dataset"] = _dataset_summary(dataset)
    return payload


def _load_lib_module(name):
    """懒加载 lib/ 下的辅助模块（sc_io / analysis_registry）。

    文件缺失或语法错误时抛可读的 EngineError，而不是让 worker 启动即崩；
    等文件补齐后下一次调用会自动重试。
    """
    module_name = "_dsh_sc_lib_" + name
    cached = sys.modules.get(module_name)
    if cached is not None:
        return cached
    lib_dir = Path(__file__).resolve().parent
    path = lib_dir / (name + ".py")
    if not path.exists():
        raise EngineError("unsupported", "缺少模块文件：%s" % path)
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise EngineError("internal", "模块无法加载：%s" % path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as error:  # pragma: no cover - 开发期兜底
        sys.modules.pop(module_name, None)
        raise EngineError("internal", "模块加载失败（%s）：%s" % (path.name, error))
    return module


def _import_options(args):
    options = {}
    for key in ("format", "member", "object", "transpose", "max_bytes", "gse", "delimiter", "name"):
        if args.get(key) is not None:
            options[key] = args.get(key)
    return options


def _convert_import_error(error):
    code = str(getattr(error, "code", "internal") or "internal")
    if code not in EngineError.CODES:
        code = "internal"
    details = getattr(error, "details", None)
    details = dict(details) if isinstance(details, dict) else {}
    candidates = getattr(error, "candidates", None)
    if candidates:
        details["candidates"] = candidates
    return EngineError(code, str(error), details=details or None)


def _source_of(args):
    url = args.get("url")
    if url:
        return str(url), True
    path = args.get("path")
    if path:
        return str(Path(os.path.expanduser(str(path))).resolve()), False
    return None, False


def cmd_open(args):
    demo = bool(args.get("demo"))
    if demo:
        for dataset in DATASETS.values():
            if dataset.source == "demo":
                return _open_payload(dataset)
        dataset = _load_demo_dataset()
        _record_import(dataset, {"demo": True}, source_format="demo")
        DATASETS[dataset.id] = dataset
        return _open_payload(dataset)
    merge_paths = args.get("paths")
    merge_dir = bool(args.get("merge_dir"))
    merge_samples = isinstance(args.get("samples"), list) and bool(args.get("samples"))
    if merge_dir or merge_paths or merge_samples:
        return _open_merged(args)
    source, is_url = _source_of(args)
    if not source:
        raise EngineError("bad_request", "open 需要提供 path、url 或 demo=true")
    if not is_url:
        for dataset in DATASETS.values():
            if dataset.path and os.path.realpath(dataset.path) == os.path.realpath(source):
                return _open_payload(dataset)
    loader = _load_lib_module("sc_io")
    try:
        loaded = loader.load_source(source, workdir=str(WORKDIR), options=_import_options(args))
    except EngineError:
        raise
    except Exception as error:
        raise _convert_import_error(error)

    matrix = loaded.get("X")
    if matrix is None:
        raise EngineError("bad_request", "导入结果没有表达矩阵：%s" % source)
    counts = loaded.get("counts")
    var_names = list(loaded.get("var_names") or [])
    name = args.get("name") or loaded.get("name") or Path(str(loaded.get("source") or source)).name
    dataset = Dataset(
        _next_dataset_id(),
        name,
        str(loaded.get("source") or source),
        matrix,
        counts if counts is not None else matrix,
        var_names,
        loaded.get("obs"),
        var=loaded.get("var"),
        warnings=list(loaded.get("warnings") or []),
        source=str(loaded.get("format") or "auto"),
    )
    dataset.source_format = str(loaded.get("format") or "auto")
    dataset.notes = list(loaded.get("notes") or [])
    dataset.meta = dict(loaded.get("meta") or {})
    # SC-02：记录 counts 的来源与语义；未提供独立 counts 层时按数值判断是否可信
    meta = dataset.meta
    if meta.get("counts_verified") is not None:
        dataset.counts_verified = bool(meta.get("counts_verified"))
        dataset.counts_source = str(meta.get("counts_source") or "n/a")
    elif counts is not None:
        dataset.counts_verified = _integer_like(counts)
        dataset.counts_source = str(loaded.get("format") or "loader") + ":counts"
    else:
        dataset.counts_verified = _integer_like(matrix)
        dataset.counts_source = "assumed_from_X"
    dataset.matrix_kind = str(meta.get("matrix_kind") or ("counts" if _integer_like(matrix) else "lognorm"))
    if counts is None:
        if dataset.counts_verified:
            dataset.warnings.append(
                "该格式没有独立 counts 层，已把表达矩阵当作原始计数使用（数值像整数计数）"
            )
        else:
            dataset.warnings.append(
                "该格式没有独立 counts 层，且表达矩阵不像整数计数：counts 语义未验证，"
                "log-normalize 等依赖原始计数的步骤会被拒绝执行"
            )
    if isinstance(dataset.meta.get("warnings"), list):
        dataset.warnings.extend(str(item) for item in dataset.meta["warnings"])
    # 已跑过一半的数据：RDS/RData 总能检测；h5ad 需要带嵌入（obsm）才认作续跑
    _apply_import_state(dataset, loaded) if (
        str(loaded.get("format") or "") in ("rds", "rdata")
        or (loaded.get("import_state") and loaded.get("obsm"))
    ) else None
    _record_import(dataset, args, source_format=str(loaded.get("format") or "auto"))
    DATASETS[dataset.id] = dataset
    return _open_payload(dataset)


def _apply_import_state(dataset, loaded):
    """读取已跑过一半的数据：检测进度并写回 dataset（可续跑）。"""
    try:
        state_module = _load_lib_module("sc_state")
        state = state_module.apply(dataset, loaded or {})
    except Exception as error:  # noqa: BLE001 - 状态检测失败不影响导入
        dataset.warnings.append("已有进度检测失败：%s" % error)
        return None
    try:
        for line in state_module.summary_bullets(state):
            dataset.notes.append(line)
    except Exception:
        pass
    return state


def _open_selected(args):
    """按用户确认的样本表（含自定义样本名/大组）合并导入。"""
    loader = _load_lib_module("sc_io")
    meta = _load_lib_module("sc_metadata")
    entries = []
    for item in args.get("samples") or []:
        if not isinstance(item, dict) or not item.get("path"):
            continue
        options = {}
        fmt = str(item.get("format") or "")
        if item.get("trio_prefix") is not None:
            options["trio_prefix"] = item.get("trio_prefix")
        else:
            mapped = {
                "seurat_rds": "rds", "rdata": "rdata", "h5ad": "h5ad",
                "10x_h5": "10x_h5", "expr_matrix": "matrix_text",
                "geo_series": "geo_series", "loom": "loom",
            }.get(fmt)
            if mapped:
                options["format"] = mapped
        entries.append({"source": str(item["path"]), "name": item.get("name"),
                        "group": item.get("group"), "options": options})
    if not entries:
        raise EngineError("bad_request", "samples 不能为空：请在样本表里勾选参与分析的样本")
    try:
        merged = meta.merge_sources(loader, entries, str(WORKDIR), options=_import_options(args))
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))
    counts = merged.get("counts")
    name = args.get("name") or ("merged-%d-samples" % len(entries))
    dataset = Dataset(_next_dataset_id(), str(name), entries[0]["source"],
                      merged["X"], counts if counts is not None else merged["X"],
                      merged["var_names"], merged["obs"],
                      warnings=list(merged.get("warnings") or []), source="merged")
    dataset.source_format = "merged"
    dataset.notes = ["按样本表合并 %d 个样本（自定义样本名/大组已写入 obs）" % len(entries)]
    if isinstance(dataset.uns, dict):
        dataset.uns.update(merged.get("meta") or {})
    if counts is None:
        dataset.warnings.append("部分样本没有独立 counts 层，已把表达矩阵当作原始计数使用")
    _record_import(dataset, args, source_format="merged")
    DATASETS[dataset.id] = dataset
    return _open_payload(dataset)


def cmd_scan(args):
    """扫描 GEO 数据根目录，识别样本（对齐单细胞分析器输入逻辑）。"""
    scan = _load_lib_module("sc_scan")
    root = args.get("path") or args.get("root") or args.get("root_dir")
    if not root:
        raise EngineError("bad_request", "scan 需要 path（GEO 下载的数据根目录）")
    try:
        result = scan.scan_data_root(root, recursive=bool(args.get("recursive")))
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))
    if result.get("errors") and not result.get("samples"):
        raise EngineError("bad_request", "；".join(str(item) for item in result["errors"]))
    return result


def cmd_geo_lookup(args):
    scan = _load_lib_module("sc_scan")
    try:
        return scan.geo_lookup(args.get("gse") or args.get("id") or "")
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))


def _open_merged(args):
    """多文件 / 目录合并导入：每个文件 = 一个样本，基因取交集（CONTRACT §11）。"""
    if isinstance(args.get("samples"), list) and args.get("samples"):
        return _open_selected(args)
    loader = _load_lib_module("sc_io")
    meta = _load_lib_module("sc_metadata")
    raw = args.get("paths")
    if isinstance(raw, str):
        paths = [item.strip() for item in raw.split(os.pathsep) if item.strip()]
    elif isinstance(raw, list):
        paths = [str(item) for item in raw if item]
    else:
        paths = []
    if not paths:
        base = args.get("path")
        if not base:
            raise EngineError("bad_request", "合并导入需要 paths:[文件...] 或 path=目录 + merge_dir=true")
        base_path = Path(os.path.expanduser(str(base))).resolve()
        if base_path.is_dir():
            paths = meta.list_sources(base_path)
        else:
            paths = [str(base_path)]
    else:
        expanded = []
        for item in paths:
            candidate = Path(os.path.expanduser(str(item))).resolve()
            if candidate.is_dir():
                expanded.extend(meta.list_sources(candidate))
            else:
                expanded.append(str(candidate))
        paths = expanded
    if not paths:
        raise EngineError("bad_request", "没有可合并的数据文件")
    try:
        merged = meta.merge_sources(loader, paths, str(WORKDIR), options=_import_options(args))
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))
    name = args.get("name") or ("merged-%d-samples" % len(paths))
    counts = merged.get("counts")
    dataset = Dataset(
        _next_dataset_id(),
        str(name),
        str(paths[0]),
        merged["X"],
        counts if counts is not None else merged["X"],
        merged["var_names"],
        merged["obs"],
        warnings=list(merged.get("warnings") or []),
        source="merged",
    )
    dataset.source_format = "merged"
    dataset.notes = ["合并 %d 个样本文件（基因取交集，sample 列=原始文件名）" % len(paths)]
    if isinstance(dataset.uns, dict):
        dataset.uns.update(merged.get("meta") or {})
    if counts is None:
        dataset.warnings.append("部分样本没有独立 counts 层，已把表达矩阵当作原始计数使用")
    _record_import(dataset, args, source_format="merged")
    DATASETS[dataset.id] = dataset
    return _open_payload(dataset)


def cmd_inspect(args):
    source, _is_url = _source_of(args)
    if not source:
        raise EngineError("bad_request", "inspect 需要提供 path 或 url")
    loader = _load_lib_module("sc_io")
    try:
        result = loader.inspect_source(source, workdir=str(WORKDIR), options=_import_options(args))
    except EngineError:
        raise
    except Exception as error:
        raise _convert_import_error(error)
    payload = dict(result or {})
    payload["ok"] = True
    return payload


def cmd_formats(args):
    registry = _load_lib_module("analysis_registry")
    return {"formats": registry.formats()}


REPLACED_EXTERNAL = {
    "cellchat_exact": "ext_cellchat", "nichenet": "ext_nichenet",
    "pyscenic_exact": "ext_pyscenic", "scenic_exact": "ext_scenic_r",
    "monocle_exact": "ext_monocle3", "palantir_exact": "ext_palantir",
    "cellphonedb": "ext_cellphonedb", "rna_velocity": "ext_scvelo",
    "celloracle": "ext_celloracle", "regvelo": "ext_regvelo",
    "infercnv_exact": "ext_infercnv", "copykat": "ext_copykat",
    "soupx_decontx": "ext_soupx",
}


def cmd_analyses(args):
    registry = _load_lib_module("analysis_registry")
    external = _load_lib_module("sc_external")
    dataset = None
    try:
        dataset = _resolve_dataset(args, required=False)
    except EngineError:
        dataset = None
    ctx = registry.build_ctx(dataset, engine_mode=ENGINE_MODE) if dataset is not None else None
    payload = registry.catalog(ctx)
    try:
        tools = external.tools_payload(dataset)
    except Exception as error:  # noqa: BLE001 - 探测失败不影响内置分析目录
        tools = {"tools": [], "errors": [str(error)]}
    ext_items = []
    for item in tools.get("tools") or []:
        entry = dict(item)
        entry["native"] = True          # Host/UI 据此允许运行
        entry["external"] = True
        ext_items.append(entry)
    kept = []
    for item in payload.get("analyses") or []:
        if not item.get("native") and item.get("id") in REPLACED_EXTERNAL:
            continue
        kept.append(item)
    order = {name: index for index, name in enumerate(payload.get("groups") or [])}
    merged = kept + ext_items
    merged.sort(key=lambda item: (
        order.get(str(item.get("group")), 99),
        0 if item.get("native") else 1,
        str(item.get("id")),
    ))
    payload["analyses"] = merged
    payload["external"] = {
        "rscript": tools.get("rscript"),
        "python": tools.get("python"),
        "seurat_version": tools.get("seurat_version"),
        "seurat_lib": tools.get("seurat_lib"),
        "versions": tools.get("versions"),
        "errors": tools.get("errors") or [],
    }
    payload["dataset"] = dataset.id if dataset is not None else None
    return payload


def cmd_analysis(args):
    dataset = _resolve_dataset(args)
    analysis_id = str(args.get("id") or "").strip()
    if analysis_id.startswith("ext_"):
        run_args = dict(args)
        run_args["dataset"] = dataset.id
        run_args["id"] = analysis_id
        return cmd_external_run(run_args)
    if not analysis_id:
        raise EngineError("bad_request", "analysis 需要 id（先用 analyses 命令看可用分析）")
    params = args.get("params") or {}
    if not isinstance(params, dict):
        raise EngineError("bad_request", "params 必须是 JSON 对象")
    registry = _load_lib_module("analysis_registry")
    ctx = registry.build_ctx(dataset, engine_mode=ENGINE_MODE)
    figure_dir = WORKDIR / "analysis"
    try:
        figure_dir.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    try:
        payload = registry.run_analysis(analysis_id, ctx, params, str(figure_dir))
    except EngineError:
        raise
    except Exception as error:
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        message = getattr(error, "message", None) or str(error)
        raise EngineError(code, message, details=getattr(error, "details", None))
    payload = dict(payload)

    # 把分析产生的 obs 列写回数据集（供嵌入图着色 / 后续分析使用）
    obs_values = payload.pop("obs_values", None) or {}
    merged = []
    for column in payload.get("obs_columns") or []:
        if not isinstance(column, dict):
            continue
        key = str(column.get("key") or "").strip()
        if not key:
            continue
        values = column.pop("values", None)
        if values is None:
            values = obs_values.get(key)
        if values is None:
            continue
        try:
            if len(values) != dataset.n_obs:
                raise ValueError("长度 %d 与细胞数 %d 不一致" % (len(values), dataset.n_obs))
            series = pd.Series(list(values), index=dataset.obs.index)
            if column.get("kind") == "categorical":
                series = series.astype("object")
            else:
                series = pd.to_numeric(series, errors="coerce")
            dataset.obs[key] = series
            merged.append(key)
        except Exception as error:  # noqa: BLE001
            payload.setdefault("warnings", []).append("未能写回 obs 列 %s：%s" % (key, error))
    payload["obs_merged"] = merged
    payload["dataset"] = dataset.id
    payload["steps"] = list(dataset.steps)
    dataset.analyses_run.append({"id": analysis_id, "params": _jsonable(params), "at": _now_iso()})
    if len(dataset.analyses_run) > 50:
        dataset.analyses_run = dataset.analyses_run[-50:]
    _record_operation(
        dataset, "analysis", "analysis", {"id": analysis_id, "params": params},
        str(payload.get("method") or "未记录"),
        method_provenance=str(payload.get("method") or "") or None,
        replayable=True,
        obs_columns=[str(column) for column in merged],
        extra={"analysis_id": analysis_id, "analysis_title": payload.get("title")},
    )
    return payload


def cmd_list(args):
    return {"datasets": [_dataset_summary(dataset) for dataset in DATASETS.values()]}


def cmd_close(args):
    dataset_id = args.get("dataset")
    if dataset_id is None:
        raise EngineError("bad_request", "close 需要 dataset")
    key = str(dataset_id)
    if key not in DATASETS:
        raise EngineError("not_found", "数据集不存在：%s" % key)
    DATASETS.pop(key)
    return {"closed": True}


def cmd_preview(args):
    dataset = _resolve_dataset(args)
    genes = [str(gene) for gene in dataset.var_names[:50]]
    return {
        "genes": genes,
        "obs_columns": _obs_columns(dataset),
        "n_obs": dataset.n_obs,
        "n_vars": dataset.n_vars,
        "steps": list(dataset.steps),
        "current": dataset.steps[-1] if dataset.steps else None,
        "history": dataset.history[-20:],
        "mode": ENGINE_MODE,
        "source_format": getattr(dataset, "source_format", None) or dataset.source,
        "notes": list(getattr(dataset, "notes", []) or []),
        "analyses_run": list(getattr(dataset, "analyses_run", []) or []),
        "samples": list((getattr(dataset, "uns", {}) or {}).get("samples") or []),
        "groups": _summary_groups(dataset),
        "group_column": (getattr(dataset, "uns", {}) or {}).get("group_column"),
        "import_state": (getattr(dataset, "uns", {}) or {}).get("import_state"),
    }


def cmd_run(args):
    dataset = _resolve_dataset(args)
    step = str(args.get("step") or "")
    if step not in STEP_ORDER:
        raise EngineError("bad_request", "未知步骤：%s（可用：%s）" % (step, ", ".join(STEP_ORDER)))
    user_params = args.get("params") or {}
    started = time.perf_counter()
    try:
        params, summary = _execute_step(dataset, step, user_params)
    except EngineError as error:
        elapsed = (time.perf_counter() - started) * 1000.0
        dataset.push_history(_history_entry(dataset, step, user_params, False, elapsed, error=error.message))
        raise
    except Exception as error:  # pragma: no cover
        elapsed = (time.perf_counter() - started) * 1000.0
        dataset.push_history(_history_entry(dataset, step, user_params, False, elapsed, error=str(error)))
        raise
    elapsed = (time.perf_counter() - started) * 1000.0
    dataset.push_history(_history_entry(dataset, step, params, True, elapsed, summary=summary))
    return {
        "step": step,
        "summary": summary,
        "steps": list(dataset.steps),
        "current": step,
        "dataset": _dataset_summary(dataset),
        "history": dataset.history[-20:],
    }


def cmd_embedding(args):
    dataset = _resolve_dataset(args)
    color_by = args.get("color_by") or {"kind": "cluster", "key": "clusters"}
    max_points = args.get("max_points", 12000)
    try:
        max_points = int(max_points)
    except Exception:
        max_points = 12000
    return _embedding_payload(dataset, color_by, max_points=max_points)


def cmd_markers(args):
    dataset = _resolve_dataset(args)
    groupby = str(args.get("groupby") or "clusters")
    top_n = _as_int(args.get("top_n", 10), "top_n")
    top_n = max(1, min(200, top_n))
    method = str(args.get("method") or "t-test")
    if "clusters" not in dataset.obs.columns and groupby == "clusters":
        raise EngineError("prerequisite", "markers 需要先完成 cluster 步骤")
    if "normalize" not in dataset.steps:
        raise EngineError("prerequisite", "markers 需要先完成 normalize 步骤")
    groups = _compute_markers(dataset, groupby, top_n, method)
    dataset.uns["markers"] = {
        "groupby": groupby,
        "top_n": int(top_n),
        "requested_method": method,
        "actual_method": "Welch t-test(builtin, one-vs-rest)+Benjamini-Hochberg",
        "fallback_reason": None,
        "n_groups": len(groups),
    }
    _record_operation(
        dataset, "command", "markers", {"groupby": groupby, "top_n": top_n, "method": method},
        "Welch t-test(builtin, one-vs-rest)+Benjamini-Hochberg",
    )
    return {"groups": groups, "groupby": groupby, "top_n": top_n, "method": "t-test(builtin)"}


def cmd_samples_remove(args):
    """按样本移除细胞（真正从数据集里剔除，重置聚类/嵌入等下游结果）。"""
    dataset = _resolve_dataset(args)
    raw = args.get("samples")
    names = [str(item).strip() for item in (raw if isinstance(raw, list) else [raw]) if str(item).strip()]
    if not names:
        raise EngineError("bad_request", "samples 不能为空")
    if "sample" not in dataset.obs.columns:
        raise EngineError("bad_request", "该数据集没有 sample 列，无法按样本移除")
    labels = [str(value) for value in dataset.obs["sample"].tolist()]
    keep = np.array([label not in names for label in labels], dtype=bool)
    removed_cells = int((~keep).sum())
    if removed_cells == 0:
        raise EngineError("bad_request", "没有匹配到要移除的样本：%s" % ", ".join(names))
    if not keep.any():
        raise EngineError("bad_request", "不能移除全部样本")
    index = np.flatnonzero(keep)
    dataset.X = dataset.X[index]
    if getattr(dataset, "counts", None) is not None:
        dataset.counts = dataset.counts[index]
    dataset.obs = dataset.obs.iloc[index].copy()
    layers = getattr(dataset, "layers", None) or {}
    if isinstance(layers, dict):
        dataset.layers = {key: value[index] for key, value in layers.items() if hasattr(value, "shape") and value.shape[0] == len(keep)}
    dataset.obsm = {}
    for column in ("clusters", "cell_type", "_emb_x", "_emb_y", "leiden", "louvain", "seurat_clusters"):
        if column in dataset.obs.columns:
            dataset.obs = dataset.obs.drop(columns=[column])
    for key in ("cluster", "neighbors", "pca_engine", "pca_variance_ratio", "embedding", "import_state"):
        dataset.uns.pop(key, None)
    dataset.markers = None
    dataset.adj = None
    dataset.knn = None
    dataset.layout = None
    dataset.resolution_scan = None
    dataset.cell_type_default = True
    dataset.steps = [step for step in ("qc", "filter", "normalize", "hvg") if step in dataset.steps]
    samples = list(dataset.uns.get("samples") or [])
    if samples:
        dataset.uns["samples"] = [
            item for item in samples
            if str(item.get("id") or item.get("name")) not in names
        ]
    if "group" in dataset.obs.columns:
        try:
            groups = sorted({str(value) for value in dataset.obs["group"].dropna().tolist()})
            dataset.uns["groups"] = groups
        except Exception:
            pass
    dataset.warnings.append("已移除样本 %s（%d 个细胞）；聚类/嵌入等下游结果已重置，请重新运行降维聚类" % (
        ", ".join(names), removed_cells))
    dataset.touch("samples_remove")
    dataset.history.append("samples_remove: %s (%d cells)" % (", ".join(names), removed_cells))
    # 步骤被重置，旧执行记录随数据版本失效；本次移除本身要作为可重放操作记录下来
    _prune_execution_records(dataset)
    _record_operation(
        dataset, "command", "samples_remove", {"samples": names},
        "manual sample removal (%d cells removed)" % removed_cells,
        extra={"removed_cells": removed_cells, "removed_samples": names},
    )
    summary = _dataset_summary(dataset)
    return {"dataset": summary, "removed_cells": removed_cells, "removed_samples": names}


def cmd_subset(args):
    """提取小亚群：按簇（或任意 obs 列取值）子集化，生成新数据集用于二次降维/聚类。"""
    dataset = _resolve_dataset(args)
    raw = args.get("clusters")
    if raw is None:
        raw = args.get("labels")
    values = [str(item).strip() for item in (raw if isinstance(raw, list) else [raw]) if str(item).strip()]
    if not values:
        raise EngineError("bad_request", "clusters 不能为空（要提取的簇号列表，如 [\"3\",\"7\"]）")
    column = str(args.get("column") or args.get("cluster_column") or "clusters")
    if column not in dataset.obs.columns:
        available = [str(item) for item in dataset.obs.columns if dataset.obs[item].nunique() <= 200][:12]
        raise EngineError("bad_request", "obs 里没有 %s 列，无法提取小亚群（可用列：%s）" % (column, ", ".join(available)))
    labels = [str(item) for item in dataset.obs[column].tolist()]
    keep = np.array([label in values for label in labels], dtype=bool)
    subset_cells = int(keep.sum())
    if subset_cells == 0:
        raise EngineError("bad_request", "没有匹配到簇：%s（该列取值示例：%s）" % (
            ", ".join(values), ", ".join(sorted(set(labels))[:10])))
    if subset_cells == len(labels):
        raise EngineError("bad_request", "不能提取全部细胞：请只选部分簇")
    index = np.flatnonzero(keep)
    name = str(args.get("name") or "") or ("%s-亚群-%s" % (dataset.name, "-".join(values[:4])))
    subset = Dataset(
        _next_dataset_id(), name, dataset.path or dataset.source,
        dataset.X[index], dataset.counts[index], dataset.var_names,
        dataset.obs.iloc[index].copy(), var=dataset.var.copy(),
        warnings=list(dataset.warnings), source="subset",
    )
    subset.source_format = "subset"
    subset.meta = dict(getattr(dataset, "meta", {}) or {})
    subset.uns.update({
        "parent_dataset": dataset.id,
        "subset_clusters": values,
        "subset_column": column,
        "samples": list(dataset.uns.get("samples") or []),
        "group_column": dataset.uns.get("group_column"),
        "groups": list(dataset.uns.get("groups") or []),
    })
    subset.steps = [step for step in ("qc", "filter", "normalize", "hvg") if step in dataset.steps]
    subset.notes = ["从「%s」提取 %s 列簇 %s：%d 个细胞（父数据集 %s）" % (
        dataset.name, column, ", ".join(values), subset_cells, dataset.id)]
    _record_operation(
        subset, "subset", "subset",
        {"parent": dataset.id, "parent_path": dataset.path, "clusters": values,
         "column": column, "parent_version": int(dataset.version)},
        "subset extraction from parent dataset",
        replayable=False,
        not_replayable_reason="子集数据集需要先重放父数据集，当前导出物只能作为记录",
        extra={"parent": dataset.id, "subset_cells": subset_cells},
    )
    DATASETS[subset.id] = subset
    return {
        "dataset": _dataset_summary(subset),
        "parent": dataset.id,
        "parent_name": dataset.name,
        "subset_cells": subset_cells,
        "kept_clusters": values,
        "cluster_column": column,
    }


def cmd_plot(args):
    dataset = _resolve_dataset(args)
    kind = str(args.get("kind") or "")
    if kind not in PLOT_KINDS:
        raise EngineError("bad_request", "未知图类型：%s（可用：%s）" % (kind, ", ".join(PLOT_KINDS)))
    params = args.get("params") or {}
    info = _make_plot(dataset, kind, params)
    return {"file": info["file"], "width": info["width"], "height": info["height"], "bytes": info["bytes"]}


def cmd_export(args):
    dataset = _resolve_dataset(args)
    what = str(args.get("what") or "")
    if what not in EXPORT_EXTENSIONS:
        raise EngineError("bad_request", "未知导出格式：%s（可用：%s）" % (what, ", ".join(EXPORT_WHATS)))
    return _export(dataset, what, args.get("path"))


def cmd_suggest(args):
    dataset = None
    try:
        if args.get("dataset") is not None:
            dataset = _resolve_dataset(args)
        elif DATASETS:
            dataset = list(DATASETS.values())[-1]
    except EngineError:
        dataset = None
    return {"suggestions": _suggestions_for(dataset)}


def cmd_audit(args):
    dataset = _resolve_dataset(args)
    return _audit(dataset)


def cmd_autoannotate(args):
    dataset = _resolve_dataset(args)
    groupby = str(args.get("groupby") or "clusters")
    return _autoannotate(dataset, groupby)


def cmd_recipe(args):
    recipe_id = args.get("id")
    if recipe_id is None or str(recipe_id).strip() == "":
        return _recipes_payload()
    return _recipes_payload(str(recipe_id))


def cmd_resolution_scan(args):
    dataset = _resolve_dataset(args)
    resolutions = args.get("resolutions")
    return _resolution_scan(dataset, resolutions)


def _merge_obs_payload(dataset, payload):
    """把分析/外部工具写回 obs 的列合并进 dataset（按细胞名对齐）。"""
    obs_values = payload.pop("obs_values", None) or {}
    merged = []
    for column in payload.get("obs_columns") or []:
        if not isinstance(column, dict):
            continue
        key = str(column.get("key") or "").strip()
        if not key:
            continue
        values = column.pop("values", None)
        if values is None:
            values = obs_values.get(key)
        if values is None:
            continue
        try:
            if len(values) != dataset.n_obs:
                raise ValueError("长度 %d 与细胞数 %d 不一致" % (len(values), dataset.n_obs))
            series = pd.Series(list(values), index=dataset.obs.index)
            if column.get("kind") == "categorical":
                series = series.astype("object")
            else:
                series = pd.to_numeric(series, errors="coerce")
            dataset.obs[key] = series
            merged.append(key)
        except Exception as error:  # noqa: BLE001
            payload.setdefault("warnings", []).append("未能写回 obs 列 %s：%s" % (key, error))
    payload["obs_merged"] = merged
    payload["dataset"] = dataset.id
    payload["steps"] = list(dataset.steps)
    return payload


def cmd_external_tools(args):
    external = _load_lib_module("sc_external")
    dataset = None
    try:
        dataset = _resolve_dataset(args, required=False)
    except EngineError:
        dataset = None
    return external.tools_payload(dataset)


def cmd_external_run(args):
    dataset = _resolve_dataset(args)
    tool_id = str(args.get("id") or "").strip()
    if not tool_id:
        raise EngineError("bad_request", "external_run 需要 id（先用 external_tools 看可用工具）")
    params = args.get("params") or {}
    if not isinstance(params, dict):
        raise EngineError("bad_request", "params 必须是 JSON 对象")
    external = _load_lib_module("sc_external")
    try:
        payload = external.run_external(dataset, tool_id, params, str(WORKDIR))
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))
    payload = _merge_obs_payload(dataset, dict(payload))
    dataset.analyses_run.append({"id": tool_id, "params": _jsonable(params), "at": _now_iso(), "external": True})
    if len(dataset.analyses_run) > 50:
        dataset.analyses_run = dataset.analyses_run[-50:]
    _record_operation(
        dataset, "external", "external_run", {"id": tool_id, "params": params},
        str(payload.get("method") or "未记录"),
        method_provenance=str(payload.get("method") or ""),
        replayable=False,
        not_replayable_reason="外部工具依赖独立的 R / Python 环境，导出物只能作为记录/示例",
        obs_columns=[str(column) for column in (payload.get("obs_merged") or [])],
    )
    return payload


def _meta_module():
    return _load_lib_module("sc_metadata")


def _meta_call(func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        code = str(getattr(error, "code", "internal") or "internal")
        if code not in EngineError.CODES:
            code = "internal"
        raise EngineError(code, getattr(error, "message", None) or str(error))


def cmd_obs_list(args):
    dataset = _resolve_dataset(args)
    return _meta_call(_meta_module().obs_list, dataset)


def cmd_samples_update(args):
    dataset = _resolve_dataset(args)
    before_columns = set(dataset.obs.columns)
    result = _meta_call(
        _meta_module().samples_update, dataset,
        args.get("samples") or [], groups=args.get("groups"),
        sample_column=str(args.get("sample_column") or "sample"),
        group_column=str(args.get("group_column") or "group"),
    )
    dataset.touch("samples_update")
    _record_operation(
        dataset, "command", "samples_update",
        {"samples": args.get("samples") or [], "groups": args.get("groups"),
         "sample_column": args.get("sample_column"), "group_column": args.get("group_column")},
        "manual sample table edit (sample names / groups)",
        obs_columns=_new_obs_columns(dataset, before_columns),
    )
    result["dataset"] = _dataset_summary(dataset)
    return result


def cmd_samples_auto_group(args):
    dataset = _resolve_dataset(args)
    pattern = args.get("pattern")
    if not pattern:
        raise EngineError("bad_request", "auto_group 需要 pattern（正则，第 1 个捕获组=大组）")
    before_columns = set(dataset.obs.columns)
    result = _meta_call(_meta_module().samples_auto_group, dataset, pattern,
                        sample_column=str(args.get("sample_column") or "sample"))
    dataset.touch("samples_auto_group")
    _record_operation(
        dataset, "command", "samples_auto_group",
        {"pattern": pattern, "sample_column": args.get("sample_column")},
        "manual group assignment by regex (samples_auto_group)",
        obs_columns=_new_obs_columns(dataset, before_columns),
    )
    return result


def cmd_samples_import_map(args):
    dataset = _resolve_dataset(args)
    path = args.get("path")
    if not path:
        raise EngineError("bad_request", "import_map 需要 path（样本→大组 对照 CSV）")
    before_columns = set(dataset.obs.columns)
    result = _meta_call(_meta_module().samples_import_map, dataset, path,
                        sample_column=args.get("sample_column"), group_column=args.get("group_column"))
    if args.get("apply"):
        dataset.touch("samples_import_map")
        sha = _sha256_file(path)
        replayable = sha is not None
        _record_operation(
            dataset, "command", "samples_import_map",
            {"path": str(path), "sample_column": args.get("sample_column"),
             "group_column": args.get("group_column"), "apply": True, "sha256": sha},
            "manual group assignment from imported sample map",
            replayable=replayable,
            not_replayable_reason=None if replayable else "样本对照表不存在，只能作为记录",
            obs_columns=_new_obs_columns(dataset, before_columns),
        )
    result["dataset"] = _dataset_summary(dataset) if args.get("apply") else None
    return result


def cmd_obs_set_group(args):
    dataset = _resolve_dataset(args)
    result = _meta_call(_meta_module().obs_set_group, dataset, args.get("column") or "")
    dataset.touch("obs_set_group")
    _record_operation(
        dataset, "command", "obs_set_group", {"column": args.get("column") or ""},
        "manual group column selection (obs_set_group)",
        extra={"group_column": (dataset.uns or {}).get("group_column")},
    )
    result["dataset"] = _dataset_summary(dataset)
    return result


def cmd_obs_join(args):
    dataset = _resolve_dataset(args)
    path = args.get("path")
    if not path:
        raise EngineError("bad_request", "obs_join 需要 path（metadata CSV）")
    before_columns = set(dataset.obs.columns)
    result = _meta_call(_meta_module().obs_join, dataset, path, key=args.get("key"), cell_column=args.get("cell_column"))
    dataset.touch("obs_join")
    sha = _sha256_file(path)
    _record_operation(
        dataset, "command", "obs_join",
        {"path": str(path), "key": args.get("key"), "cell_column": args.get("cell_column"),
         "sha256": sha},
        "manual metadata join (obs_join)",
        replayable=sha is not None,
        not_replayable_reason=None if sha is not None else "metadata CSV 不存在，只能作为记录",
        obs_columns=_new_obs_columns(dataset, before_columns),
    )
    result["dataset"] = _dataset_summary(dataset)
    return result


def cmd_obs_derive(args):
    dataset = _resolve_dataset(args)
    before_columns = set(dataset.obs.columns)
    result = _meta_call(
        _meta_module().obs_derive, dataset,
        args.get("name") or "", source=str(args.get("source") or "cell_name"),
        pattern=args.get("pattern"), group=int(args.get("group") or 1), value=args.get("value"),
    )
    dataset.touch("obs_derive")
    _record_operation(
        dataset, "command", "obs_derive",
        {"name": args.get("name") or "", "source": str(args.get("source") or "cell_name"),
         "pattern": args.get("pattern"), "group": int(args.get("group") or 1),
         "value": args.get("value")},
        "manual metadata derivation (obs_derive)",
        obs_columns=_new_obs_columns(dataset, before_columns),
    )
    result["dataset"] = _dataset_summary(dataset)
    return result


def cmd_shutdown(args):
    return {"bye": True}


COMMANDS = {
    "ping": cmd_ping,
    "open": cmd_open,
    "inspect": cmd_inspect,
    "formats": cmd_formats,
    "analyses": cmd_analyses,
    "analysis": cmd_analysis,
    "external_tools": cmd_external_tools,
    "external_run": cmd_external_run,
    "obs_list": cmd_obs_list,
    "obs_join": cmd_obs_join,
    "obs_derive": cmd_obs_derive,
    "obs_set_group": cmd_obs_set_group,
    "samples_update": cmd_samples_update,
    "samples_auto_group": cmd_samples_auto_group,
    "samples_import_map": cmd_samples_import_map,
    "scan": cmd_scan,
    "geo_lookup": cmd_geo_lookup,
    "samples_remove": cmd_samples_remove,
    "subset": cmd_subset,
    "list": cmd_list,
    "close": cmd_close,
    "preview": cmd_preview,
    "run": cmd_run,
    "embedding": cmd_embedding,
    "markers": cmd_markers,
    "plot": cmd_plot,
    "export": cmd_export,
    "suggest": cmd_suggest,
    "audit": cmd_audit,
    "autoannotate": cmd_autoannotate,
    "recipe": cmd_recipe,
    "resolution_scan": cmd_resolution_scan,
    "shutdown": cmd_shutdown,
}


# ===========================================================================
# 16. 主循环（stdin 一行一个请求，stdout 一行一个响应）
# ===========================================================================


def _handle_request(request):
    if not isinstance(request, dict):
        raise EngineError("bad_request", "请求必须是 JSON 对象")
    command = request.get("cmd")
    if not command:
        raise EngineError("bad_request", "请求缺少 cmd 字段")
    args = request.get("args")
    if args is None:
        args = {}
    if not isinstance(args, dict):
        raise EngineError("bad_request", "args 必须是 JSON 对象")
    handler = COMMANDS.get(str(command))
    if handler is None:
        raise EngineError("bad_request", "未知命令：%s" % command)
    return handler(args)


def main(argv=None):
    _log(
        "worker ready: mode=%s python=%s workdir=%s layouts=%s"
        % (ENGINE_MODE, PY_VERSION, WORKDIR, ",".join(LAYOUTS))
    )
    for line in sys.stdin:
        text = line.strip()
        if not text:
            continue
        request_id = None
        shutdown = False
        try:
            try:
                request = json.loads(text)
            except json.JSONDecodeError:
                raise EngineError("bad_request", "请求不是合法 JSON")
            if isinstance(request, dict):
                request_id = request.get("id")
            command = str(request.get("cmd") or "") if isinstance(request, dict) else ""
            result = _handle_request(request)
            _emit({"id": request_id, "ok": True, "result": result})
            shutdown = command == "shutdown"
        except EngineError as error:
            payload = {"id": request_id, "ok": False, "error": {"code": error.code, "message": error.message}}
            if getattr(error, "details", None):
                payload["error"]["details"] = error.details
            _emit(payload)
        except Exception as error:  # pragma: no cover - 兜底：单条失败不影响进程
            trace = traceback.format_exc()
            _log("命令内部错误：%s\n%s" % (error, trace))
            _emit(
                {
                    "id": request_id,
                    "ok": False,
                    "error": {
                        "code": "internal",
                        "message": "引擎内部错误：%s" % error,
                        "trace": trace[-4000:],
                    },
                }
            )
        if shutdown:
            return 0
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:  # pragma: no cover
        raise SystemExit(0)
    except BrokenPipeError:  # pragma: no cover
        raise SystemExit(0)
