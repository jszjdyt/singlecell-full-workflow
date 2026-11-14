#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench 官方 Python 工具适配器公共库（CONTRACT §10.2 / §10.3）。

所有 ``lib/external/*.py`` 适配器共用本模块，负责：

1. CLI 约定：``<externalPython> <tool>.py --bundle DIR --out DIR --params JSON``；
2. bundle（§10.2）读取：``counts.mtx`` / ``data.mtx`` 是 **genes × cells** 的
   MatrixMarket，读出后转成 cells × genes；``features.tsv`` / ``barcodes.tsv`` /
   ``obs.csv`` / ``embed.csv`` / ``pca.csv`` / ``meta.json``；
3. 规模上限与确定性下采样（种子固定，warning 记录）；
4. ``result.json`` / ``obs_extra.csv`` 写出（键严格按 §10.3）；
5. 出图规范：白底、dpi≥140、图上文字英文；
6. 退出码：0 成功；**3 = 缺前置/数据不满足**（stderr 中文原因）；其它 = 失败。

设计原则：stdout 必须干净（所有日志走 stderr），任何"数据不够"的情况都要给
中文原因并以退出码 3 结束，绝不编造结果。
"""

from __future__ import annotations

import json
import os
import struct
import sys
import traceback
from pathlib import Path

# ---------------------------------------------------------------------------
# 退出码 / 错误
# ---------------------------------------------------------------------------

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_PREREQ = 3


class PrerequisiteError(Exception):
    """缺前置（数据/依赖/参考文件），对应退出码 3 + 中文 stderr。"""


class AdapterError(Exception):
    """适配器内部错误，对应退出码 1。"""


def log(message):
    """日志统一走 stderr，保证 stdout 干净。"""
    print(str(message), file=sys.stderr, flush=True)


def require_module(module_name, package_name=None):
    """导入依赖；失败抛 PrerequisiteError（退出码 3）而不是崩栈。"""
    try:
        module = __import__(module_name, fromlist=["__version__"])
    except Exception as error:  # pragma: no cover - 依赖探测
        package_name = package_name or module_name
        raise PrerequisiteError(
            "缺少 Python 包 %s（import %s 失败：%s）。请用 externalPython 安装："
            "pip install %s" % (package_name, module_name, error, package_name)
        )
    return module


def version_of(module, fallback="?"):
    try:
        return str(getattr(module, "__version__"))
    except Exception:
        try:
            from importlib.metadata import version as _v

            return str(_v(getattr(module, "__name__", "")))
        except Exception:
            return fallback


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description="dsh-singlecell-workbench external adapter (CONTRACT §10.3)")
    parser.add_argument("--bundle", required=True, help="bundle 目录（§10.2）")
    parser.add_argument("--out", required=True, help="输出目录（result.json/图/obs_extra.csv）")
    parser.add_argument("--params", default="{}", help="JSON 参数")
    parsed = parser.parse_args(argv)
    try:
        params = json.loads(parsed.params) if parsed.params else {}
    except Exception as error:
        raise PrerequisiteError("--params 不是合法 JSON：%s" % error)
    if not isinstance(params, dict):
        raise PrerequisiteError("--params 必须是 JSON 对象")
    bundle = Path(os.path.abspath(os.path.expanduser(parsed.bundle)))
    out = Path(os.path.abspath(os.path.expanduser(parsed.out)))
    if not bundle.is_dir():
        raise PrerequisiteError("bundle 目录不存在：%s" % bundle)
    out.mkdir(parents=True, exist_ok=True)
    return bundle, out, params


def run_adapter(adapter):
    """统一的 CLI 包装：解析参数 → 运行 → 写 result.json；异常映射到退出码。"""
    import time

    started = time.time()
    try:
        bundle, out, params = parse_args()
        log("adapter start: bundle=%s out=%s" % (bundle, out))
        # 第三方库（pyscenic/df2regulons 等）会用 print() 往 stdout 写日志；
        # 为了让 stdout 保持"只有一行结果摘要"的干净状态，这里把 stdout 重定向到 stderr。
        import contextlib

        with contextlib.redirect_stdout(sys.stderr):
            payload = adapter(bundle, out, params)
        payload.setdefault("warnings", [])
        payload.setdefault("notes", [])
        payload = write_result(out, payload)
        payload["elapsed_ms"] = int(round((time.time() - started) * 1000))
        log("adapter done in %.1fs" % (time.time() - started))
        # stdout 只输出一行机器可读摘要，便于上游/测试确认
        print(json.dumps({"ok": True, "id": payload.get("id"), "elapsed_ms": payload["elapsed_ms"]}, ensure_ascii=False))
        return 0
    except PrerequisiteError as error:
        sys.stderr.write("\n[prerequisite] %s\n" % str(error).strip())
        sys.stderr.flush()
        return EXIT_PREREQ
    except SystemExit:
        raise
    except Exception as error:  # pragma: no cover - 兜底
        log("[error] 适配器运行失败：%s" % error)
        traceback.print_exc(file=sys.stderr)
        return EXIT_ERROR


# ---------------------------------------------------------------------------
# bundle 读取（§10.2）
# ---------------------------------------------------------------------------


def _read_lines(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return [line.rstrip("\r") for line in text.splitlines() if line.strip() != ""]


def _mmread_matrix(path):
    """读 MatrixMarket（genes × cells）→ cells × genes 的 csr float32。"""
    import numpy as np
    from scipy import sparse
    from scipy.io import mmread

    try:
        matrix = mmread(str(path))
    except Exception as error:
        raise PrerequisiteError("无法读取矩阵 %s：%s" % (path, error))
    if sparse.issparse(matrix):
        matrix = matrix.tocsr().transpose().tocsr()
    else:
        matrix = sparse.csr_matrix(np.asarray(matrix).T)
    matrix = matrix.astype("float32")
    return matrix


def read_mtx_matrix(path):
    """读 MatrixMarket（genes × cells）→ cells × genes 的 csr float32（公开给 scVelo 等层文件用）。"""
    return _mmread_matrix(Path(path))


class Bundle(object):
    """§10.2 bundle 的只读封装。"""

    def __init__(self, root):
        self.root = Path(root)
        self.features = _read_lines(self.root / "features.tsv") if (self.root / "features.tsv").exists() else []
        self.barcodes = _read_lines(self.root / "barcodes.tsv") if (self.root / "barcodes.tsv").exists() else []
        self.meta = {}
        meta_path = self.root / "meta.json"
        if meta_path.exists():
            try:
                self.meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:
                self.meta = {}
        self.obs = self._read_obs()
        if not self.features or not self.barcodes:
            raise PrerequisiteError("bundle 缺少 features.tsv 或 barcodes.tsv（§10.2 约定）")
        self._cache = {}

    def _read_obs(self):
        import pandas as pd

        path = self.root / "obs.csv"
        if not path.exists():
            return pd.DataFrame(index=[str(item) for item in self.barcodes]).rename_axis("cell")
        try:
            frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
        except Exception as error:
            raise PrerequisiteError("无法读取 obs.csv：%s" % error)
        if frame.shape[1] == 0:
            return pd.DataFrame(index=[str(item) for item in self.barcodes]).rename_axis("cell")
        first = frame.columns[0]
        frame[first] = frame[first].astype(str)
        frame = frame.set_index(first)
        frame.index = [str(item) for item in frame.index]
        # 与 barcodes 对齐（barcode 是权威顺序）
        frame = frame.reindex([str(item) for item in self.barcodes])
        frame = frame.fillna("")
        return frame

    @property
    def n_obs(self):
        return len(self.barcodes)

    @property
    def n_vars(self):
        return len(self.features)

    @property
    def x_kind(self):
        return str(self.meta.get("X_kind") or "counts")

    def matrix(self, kind="data"):
        """kind: data（当前矩阵，可能 lognorm）| counts（原始 counts）。"""
        key = "matrix:%s" % kind
        if key in self._cache:
            return self._cache[key]
        filename = "data.mtx" if kind == "data" else "counts.mtx"
        path = self.root / filename
        if not path.exists():
            if kind == "data":
                return self.matrix("counts")
            raise PrerequisiteError("bundle 缺少 %s（§10.2 约定）" % filename)
        matrix = _mmread_matrix(path)
        if matrix.shape != (self.n_obs, self.n_vars):
            raise PrerequisiteError(
                "矩阵 %s 形状 %s 与 barcodes/features（%d × %d）不一致"
                % (filename, tuple(matrix.shape), self.n_obs, self.n_vars)
            )
        self._cache[key] = matrix
        return matrix

    def embed(self):
        import pandas as pd

        path = self.root / "embed.csv"
        if not path.exists():
            return None
        try:
            frame = pd.read_csv(path, index_col=0)
        except Exception:
            return None
        frame = frame.reindex([str(item) for item in self.barcodes])
        if frame.shape[1] < 2:
            return None
        values = frame.iloc[:, :2].to_numpy(dtype="float64")
        return values

    def pca(self):
        import pandas as pd

        path = self.root / "pca.csv"
        if not path.exists():
            return None
        try:
            frame = pd.read_csv(path, index_col=0)
        except Exception:
            return None
        frame = frame.reindex([str(item) for item in self.barcodes])
        values = frame.to_numpy(dtype="float64")
        if values.ndim != 2 or values.shape[1] < 2:
            return None
        return values

    def obs_column(self, column, default=None):
        """返回与 barcodes 等长的字符串列（不存在时返回 default 或 None）。"""
        if column and column in self.obs.columns:
            return [str(item) for item in self.obs[column].tolist()]
        return default

    def obs_columns_like(self, candidates):
        """在 obs 中找第一个存在的候选列名。"""
        for column in candidates:
            if column and column in self.obs.columns:
                return column
        return None

    def to_adata(self, kind="data", subset=None, genes=None):
        """转成 anndata.AnnData（X=cells×genes，稠密 float32）。"""
        require_module("anndata", "anndata")
        import anndata as ad
        import numpy as np
        import pandas as pd

        matrix = self.matrix(kind)
        if subset is not None:
            matrix = matrix[subset]
        if genes is not None:
            matrix = matrix[:, genes]
        if hasattr(matrix, "toarray"):
            values = matrix.toarray()
        else:
            values = np.asarray(matrix)
        values = np.ascontiguousarray(values.astype("float32"))
        obs = self.obs.iloc[subset] if subset is not None else self.obs
        features = self.features
        if genes is not None:
            features = [features[index] for index in genes]
        barcodes = [str(self.barcodes[index]) for index in (range(self.n_obs) if subset is None else subset)]
        adata = ad.AnnData(
            X=values,
            obs=pd.DataFrame(index=pd.Index(barcodes, name="cell")),
            var=pd.DataFrame(index=pd.Index([str(item) for item in features], name="gene")),
        )
        for column in obs.columns:
            try:
                adata.obs[str(column)] = obs[column].reindex(adata.obs.index).astype(str).to_numpy()
            except Exception:
                pass
        return adata


# ---------------------------------------------------------------------------
# 规模控制（确定性下采样）
# ---------------------------------------------------------------------------


def cap_cells(bundle, max_cells, seed=0, warnings=None, label=""):
    """细胞数超限时确定性下采样；未超限返回 None，否则返回升序细胞下标数组。"""
    import numpy as np

    if not max_cells or bundle.n_obs <= int(max_cells):
        return None
    chooser = np.random.default_rng(int(seed))
    cells = np.sort(chooser.choice(bundle.n_obs, size=int(max_cells), replace=False))
    message = "为控制规模对细胞确定性下采样 %d → %d（%s）" % (bundle.n_obs, int(max_cells), label or "官方工具")
    if warnings is not None:
        warnings.append(message)
    log(message)
    return cells


def select_genes(matrix, features, max_genes, keep_genes=None, seed=0, warnings=None):
    """按方差选基因（可选强制保留 keep_genes），返回排序后的基因下标。"""
    import numpy as np

    n_vars = len(features)
    if not max_genes or n_vars <= int(max_genes):
        return None
    keep_idx = []
    known = {}
    if keep_genes:
        for index, name in enumerate(features):
            known.setdefault(str(name), index)
        for name in keep_genes:
            if str(name) in known:
                keep_idx.append(known[str(name)])
    keep_idx = list(dict.fromkeys(keep_idx))
    if hasattr(matrix, "toarray"):
        values = np.asarray(matrix.toarray(), dtype="float64")
    else:
        values = np.asarray(matrix, dtype="float64")
    if values.shape[1] != n_vars:
        raise AdapterError("select_genes 的矩阵列数与 features 不一致")
    if values.shape[0] > 2000:
        values = values[:2000]
    variance = values.var(axis=0)
    order = np.argsort(variance)[::-1]
    picked = list(dict.fromkeys(keep_idx + [int(index) for index in order.tolist()]))[: int(max_genes)]
    picked = np.sort(np.asarray(picked, dtype="int64"))
    message = "为控制规模取高方差基因 %d → %d（保留 TF 等必需基因）" % (n_vars, len(picked))
    if warnings is not None:
        warnings.append(message)
    log(message)
    return picked


def subset_arrays(bundle, cell_idx=None, gene_idx=None):
    """按细胞/基因子集返回常用数据的 dict。"""
    cell_idx = None if cell_idx is None else list(cell_idx)
    gene_idx = None if gene_idx is None else list(gene_idx)
    barcodes = [str(item) for item in bundle.barcodes]
    features = [str(item) for item in bundle.features]
    if cell_idx is not None:
        barcodes = [barcodes[index] for index in cell_idx]
    if gene_idx is not None:
        features = [features[index] for index in gene_idx]
    obs = bundle.obs.iloc[cell_idx] if cell_idx is not None else bundle.obs
    embed = bundle.embed()
    pca = bundle.pca()
    if embed is not None and cell_idx is not None:
        embed = embed[cell_idx]
    if pca is not None and cell_idx is not None:
        pca = pca[cell_idx]
    return {"cells": cell_idx, "genes": gene_idx, "barcodes": barcodes, "features": features,
            "obs": obs, "embed": embed, "pca": pca}


# ---------------------------------------------------------------------------
# result.json / obs_extra.csv（§10.3）
# ---------------------------------------------------------------------------


def png_size(path):
    try:
        with open(path, "rb") as handle:
            head = handle.read(24)
        if len(head) >= 24 and head[:8] == b"\x89PNG\r\n\x1a\n":
            return struct.unpack(">II", head[16:24])
    except Exception:
        pass
    return (0, 0)


def write_result(out, payload):
    """写 <out>/result.json，并补齐 §10.3 固定键。"""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    result = dict(payload or {})
    summary = result.get("summary")
    if not isinstance(summary, dict):
        summary = {}
    summary.setdefault("bullets", [])
    summary.setdefault("metrics", {})
    result["summary"] = summary
    figures = []
    for figure in result.get("figures") or []:
        if not isinstance(figure, dict) or not figure.get("file"):
            continue
        file = Path(str(figure["file"]))
        if not file.is_absolute():
            file = (out / file).resolve()
        if not file.exists():
            log("warning: 图不存在，已跳过：%s" % file)
            continue
        entry = dict(figure)
        entry["file"] = str(file)
        width, height = png_size(str(file))
        entry.setdefault("width", int(width))
        entry.setdefault("height", int(height))
        figures.append(entry)
    result["figures"] = figures
    tables = []
    for table in result.get("tables") or []:
        if not isinstance(table, dict):
            continue
        entry = dict(table)
        entry.setdefault("id", "table")
        entry.setdefault("title", entry["id"])
        entry.setdefault("columns", [])
        entry.setdefault("rows", [])
        tables.append(entry)
    result["tables"] = tables
    result.setdefault("warnings", [])
    result.setdefault("notes", [])
    if result.get("obs_file"):
        result["obs_file"] = str(result["obs_file"])
    path = out / "result.json"
    path.write_text(json.dumps(result, ensure_ascii=False, default=_json_default), encoding="utf-8")
    return result


def _json_default(value):
    try:
        import numpy as np

        if isinstance(value, (np.integer,)):
            return int(value)
        if isinstance(value, (np.floating,)):
            return float(value) if np.isfinite(value) else None
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, (np.bool_,)):
            return bool(value)
    except Exception:
        pass
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        import math

        return value if math.isfinite(value) else None
    return str(value)


def write_obs_extra(out, barcodes, columns, filename="obs_extra.csv"):
    """写 obs_extra.csv：第一列细胞名 + 新增列（§10.3）。返回绝对路径。"""
    import pandas as pd

    out = Path(out)
    frame = pd.DataFrame({"cell": [str(item) for item in barcodes]})
    for key, values in columns.items():
        frame[str(key)] = list(values)
    path = out / filename
    frame.to_csv(path, index=False, encoding="utf-8")
    return str(path.resolve())


def table(columns, rows, table_id=None, title=None, note=None):
    """构造 §10.3 表格：columns=[{"key","label","type"}]，rows=[{...}]。"""
    entry = {
        "id": table_id or "table",
        "title": title or (table_id or "Table"),
        "columns": columns,
        "rows": rows,
    }
    if note:
        entry["note"] = note
    return entry


def col(key, label, ptype="string"):
    return {"key": key, "label": label, "type": ptype}


def figure(out, fig, slug, figure_id=None, title=None, dpi=140):
    """保存 matplotlib figure 为白底 PNG（英文图注），返回 §10.3 figure 条目。"""
    import matplotlib.pyplot as plt

    out = Path(out)
    path = out / ("%s.png" % slug)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(path, dpi=int(dpi), facecolor="white", edgecolor="white", bbox_inches="tight")
    plt.close(fig)
    width, height = png_size(str(path))
    return {"id": figure_id or slug, "title": title or slug, "file": str(path.resolve()),
            "width": int(width), "height": int(height)}


def setup_matplotlib():
    """统一出图规范：Agg 后端、白底、dpi≥140、英文文字。"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "figure.dpi": 140,
        "savefig.dpi": 140,
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "DejaVu Sans",
    })
    return plt


# ---------------------------------------------------------------------------
# 参考文件探测（SCENIC motif 库 / CellOracle base GRN 等）
# ---------------------------------------------------------------------------

HOME = Path(os.path.expanduser("~"))


def first_existing(paths):
    for item in paths:
        if not item:
            continue
        path = Path(os.path.expanduser(str(item)))
        if path.exists():
            return path
    return None


def cache_data_dir():
    path = HOME / ".dsh" / "sc-plugin-pydata"
    try:
        path.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    return path


SCENIC_KNOWN_DIRS = [
    str(HOME / "Desktop" / "单细胞工作站" / "SCINCE尝试" / "human" / "python用"),
    str(HOME / "Desktop" / "单细胞工作站" / "SCINCE尝试" / "human"),
    str(HOME / "cisTarget_databases"),
    str(HOME / ".dsh" / "sc-plugin-pydata" / "cistarget"),
]


def find_scenic_db(params=None):
    """探测 SCENIC motif 库：返回 (db_files, annotation_file, db_dir)。

    优先级：params.db_dir / params.dbs / params.scenic_db → DS_SC_SCENIC_DB →
    已知用户目录 → 插件缓存目录。找不到抛 PrerequisiteError（退出码 3）。
    """
    params = params or {}
    candidates = []
    explicit = params.get("db_dir") or params.get("scenic_db") or params.get("db")
    if explicit:
        candidates.append(str(explicit))
    dbs = params.get("dbs")
    if isinstance(dbs, list):
        files = [Path(os.path.expanduser(str(item))) for item in dbs if item]
        existing = [item for item in files if item.is_file()]
        if existing:
            annotation = _find_annotation(params, [str(item.parent) for item in existing] + SCENIC_KNOWN_DIRS)
            if annotation is None:
                raise PrerequisiteError(
                    "找到 motif 库 feather 但缺少 motif annotation 表（*motifs*.tbl）；"
                    "请通过 params.annotations 指定，或放到 %s" % (cache_data_dir() / "cistarget")
                )
            return [str(item) for item in existing], str(annotation), str(existing[0].parent)
    env_db = os.environ.get("DS_SC_SCENIC_DB")
    if env_db:
        candidates.append(env_db)
    candidates.extend(SCENIC_KNOWN_DIRS)
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(os.path.expanduser(str(candidate)))
        if path.is_file() and path.suffix == ".feather":
            files = [path]
            search_dirs = [str(path.parent)]
        elif path.is_dir():
            files = sorted(path.glob("*.genes_vs_motifs.rankings.feather"))
            if not files:
                files = sorted(path.glob("*rankings*.feather"))
            search_dirs = [str(path)]
        else:
            continue
        if not files:
            continue
        annotation = _find_annotation(params, search_dirs + SCENIC_KNOWN_DIRS)
        if annotation is None:
            raise PrerequisiteError(
                "找到 SCENIC motif 库目录 %s，但缺少 motif annotation 表（motifs-v10*.tbl）。"
                "请从 https://resources.aertslab.org/cistarget/motif2tf/ 下载后放到该目录或 %s，"
                "或用 params.annotations 指定。" % (path, cache_data_dir() / "cistarget")
            )
        return [str(item) for item in files], str(annotation), str(path if path.is_dir() else path.parent)
    raise PrerequisiteError(
        "未找到 SCENIC motif 库（*.genes_vs_motifs.rankings.feather）。"
        "请用 params.db_dir 指定目录（如 ~/scenic_databases），"
        "或设置环境变量 DS_SC_SCENIC_DB。"
    )


def _find_annotation(params, dirs):
    explicit = params.get("annotations") or params.get("annotation") or params.get("motif_annotations")
    if explicit:
        path = Path(os.path.expanduser(str(explicit)))
        if path.is_file():
            return path
        raise PrerequisiteError("params.annotations 指定的文件不存在：%s" % path)
    env_annotation = os.environ.get("DS_SC_SCENIC_ANNOTATIONS")
    if env_annotation and Path(os.path.expanduser(env_annotation)).is_file():
        return Path(os.path.expanduser(env_annotation))
    patterns = ["*motifs*.tbl", "motifs-v*.tbl", "*.tbl"]
    for directory in dirs:
        if not directory:
            continue
        path = Path(os.path.expanduser(str(directory)))
        if not path.is_dir():
            continue
        for pattern in patterns:
            matches = sorted(path.glob(pattern))
            if matches:
                return matches[0]
    return None


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------


def json_safe(value):
    return json.loads(json.dumps(value, default=_json_default))


def to_builtin(value):
    """numpy/pandas 标量 → Python 原生类型（result.json 用）。"""
    import numpy as np

    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def top_values(series, n=30):
    """按绝对值取 top n（Series 索引=名称），返回 [(name, value)]。"""
    import numpy as np

    values = np.asarray(series, dtype="float64")
    order = np.argsort(np.abs(values))[::-1][: int(n)]
    return [(str(series.index[index]), float(values[index])) for index in order]
