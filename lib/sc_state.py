#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""已跑过一半的数据：判断跑到哪一步，并把结果接回工作台状态（可续跑）。

支持来源：
  - Seurat RDS / RData（R 桥已导出 obs + pca.csv + embed.csv + variable_features.txt + import_state）
  - h5ad（obsm['X_pca']/['X_umap']、var['highly_variable']、uns['rank_genes_groups']）

输出统一映射到本插件 10 步流程：
  qc → filter → normalize → hvg → pca → neighbors → cluster → embed → markers → annotate
已完成的步骤会被写进 dataset.steps，UI 流程卡直接显示"已完成"，
后续步骤（如 markers/annotate/embedding 着色）可以在此基础上继续跑。
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

STEP_ORDER = ["qc", "filter", "normalize", "hvg", "pca", "neighbors", "cluster", "embed", "markers", "annotate"]

CLUSTER_NAME_HINTS = (
    "clusters", "cluster", "seurat_clusters", "leiden", "louvain", "snn_res",
    "clusters_res", "kmeans", "subcluster", "cluster_r",
)
CELLTYPE_NAME_HINTS = (
    "celltype", "cell_type", "cell types", "annotation", "cell_annotation",
    "celltype_major", "celltype_minor", "celltype_fine", "cell_type_major", "type", "label",
)
# 不是细胞类型的常见元数据列（避免把样本/批次当注释）
METADATA_EXCLUDE = (
    "orig.ident", "orig_ident", "sample", "batch", "group", "condition", "treatment",
    "donor", "patient", "subject", "sex", "age", "stage", "timepoint", "library",
    "lane", "replicate", "stim", "genotype", "disease", "cluster", "seurat", "res",
)


def _as_str_list(values):
    return [str(item) for item in (values or []) if str(item).strip()]


def _categorical_columns(obs):
    out = []
    for name in obs.columns:
        series = obs[name]
        if pd.api.types.is_bool_dtype(series):
            continue
        if pd.api.types.is_numeric_dtype(series):
            continue
        try:
            unique = int(series.nunique(dropna=True))
        except Exception:
            continue
        if 1 < unique <= 200:
            out.append((str(name), unique))
    return out


def _column_unique(series):
    try:
        return int(series.nunique(dropna=True))
    except Exception:
        return 0


def _pick_cluster_column(obs):
    # 1) 名字命中聚类关键词（数值型也认：Seurat 的 seurat_clusters/RNA_snn_res.* 常是数值）
    for name in obs.columns:
        lower = str(name).lower()
        if not any(hint in lower for hint in CLUSTER_NAME_HINTS):
            continue
        series = obs[name]
        if pd.api.types.is_bool_dtype(series):
            continue
        unique = _column_unique(series)
        if 1 < unique <= 200:
            return str(name)
    categorical = _categorical_columns(obs)
    names = [name for name, _unique in categorical]
    lowered = {name.lower(): name for name in names}
    for hint in CLUSTER_NAME_HINTS:
        for lower, original in lowered.items():
            if hint in lower:
                return original
    # 兜底：取值像数字簇号（0/1/2…）且数量 2–100
    for name, unique in categorical:
        if unique < 2 or unique > 100:
            continue
        values = [str(v) for v in obs[name].dropna().unique().tolist()[:20]]
        if values and all(re.fullmatch(r"\d+", value) for value in values):
            return name
    return None


def _pick_celltype_column(obs, cluster_column=None):
    categorical = _categorical_columns(obs)
    names = [name for name, _unique in categorical]
    lowered = {name.lower(): name for name in names}
    for hint in CELLTYPE_NAME_HINTS:
        for lower, original in lowered.items():
            if hint in lower and original != cluster_column:
                if any(token in lower for token in METADATA_EXCLUDE):
                    continue
                return original
    # 兜底：字符串型（非纯数字）分类列，排除样本/批次等元数据列
    for name, unique in categorical:
        if name == cluster_column:
            continue
        lower = str(name).lower()
        if any(token in lower for token in METADATA_EXCLUDE):
            continue
        values = [str(v) for v in obs[name].dropna().unique().tolist()[:20]]
        if len(values) >= 2 and any(not re.fullmatch(r"\d+", value) for value in values):
            return name
    return None


def _x_is_lognormalized(dataset):
    matrix = getattr(dataset, "X", None)
    if matrix is None:
        return False
    try:
        if hasattr(matrix, "data"):
            sample = np.asarray(matrix.data[: min(2000, matrix.nnz)], dtype=np.float64)
        else:
            array = np.asarray(matrix)
            sample = array.reshape(-1)[:2000].astype(np.float64)
    except Exception:
        return False
    sample = sample[np.isfinite(sample)]
    if sample.size == 0:
        return False
    if np.all(np.abs(sample - np.round(sample)) < 1e-6):
        return False
    return float(np.max(sample)) < 60.0


def _has_counts_layer(dataset):
    layers = getattr(dataset, "layers", None) or {}
    return "counts" in layers


def _obs_qc_columns(obs):
    patterns = ("n_genes", "n_counts", "pct_counts", "percent.mt", "ncount", "nfeature", "qc")
    found = []
    for name in obs.columns:
        lower = str(name).lower()
        if any(pattern in lower for pattern in patterns):
            found.append(str(name))
    return found


def detect(dataset, loaded=None):
    """返回 {"steps": [...], "current": str|None, "details": {...}, "mapped": {...}}。"""
    loaded = loaded or {}
    obs = dataset.obs
    details = {}
    detected = set()

    # QC / filter：导入的数据默认已质控过滤（RDS/RData 基本都是分析中途存档）
    qc_columns = _obs_qc_columns(obs)
    details["qc_columns"] = qc_columns
    detected.update({"qc", "filter"})

    # normalize
    lognorm = _x_is_lognormalized(dataset)
    details["x_lognormalized"] = bool(lognorm)
    if lognorm:
        detected.add("normalize")

    # hvg
    variable_features = _as_str_list(loaded.get("var_features") or [])
    var = getattr(dataset, "var", None)
    if not variable_features and isinstance(var, pd.DataFrame) and "highly_variable" in var.columns:
        try:
            mask = np.asarray(var["highly_variable"].to_numpy(), dtype=bool)
            variable_features = [str(gene) for gene, keep in zip(dataset.var_names, mask) if keep]
        except Exception:
            variable_features = []
    details["n_variable_features"] = len(variable_features)
    if variable_features:
        detected.add("hvg")
        if isinstance(var, pd.DataFrame) and "highly_variable" not in var.columns:
            mask = np.isin(np.asarray(dataset.var_names, dtype=object), np.asarray(variable_features, dtype=object))
            dataset.var["highly_variable"] = pd.Series(mask, index=dataset.var.index)

    import_state = loaded.get("import_state") or {}
    reductions = _as_str_list(import_state.get("reductions") or [])
    graphs = _as_str_list(import_state.get("graphs") or [])
    details["reductions"] = reductions
    details["graphs"] = graphs
    details["source_class"] = import_state.get("class")

    # embeddings（R 桥 / h5ad 提前导出）
    obsm = loaded.get("obsm") or {}
    pca = obsm.get("X_pca")
    embed = obsm.get("X_emb")
    embed_kind = loaded.get("embed_kind") or _embed_name(reductions)
    if pca is not None:
        try:
            matrix = np.asarray(pca, dtype=np.float32)
            if matrix.ndim == 2 and matrix.shape[0] == dataset.n_obs:
                dataset.obsm["X_pca"] = matrix
                detected.add("pca")
                details["n_pcs"] = int(matrix.shape[1])
        except Exception:
            pass
    if "pca" in [name.lower() for name in reductions]:
        detected.add("pca")
        details.setdefault("n_pcs", None)
        if "pca" not in details.get("reductions", []):
            details["pca_from_reduction"] = True
    if embed is not None:
        try:
            matrix = np.asarray(embed, dtype=np.float32)
            if matrix.ndim == 2 and matrix.shape[0] == dataset.n_obs and matrix.shape[1] >= 2:
                dataset.obsm["X_emb"] = matrix[:, :2]
                dataset.obs["_emb_x"] = matrix[:, 0]
                dataset.obs["_emb_y"] = matrix[:, 1]
                dataset.layout = (embed_kind or "umap") + "(imported)"
                detected.add("embed")
                details["embed_kind"] = embed_kind or "unknown"
        except Exception:
            pass
    if embed_kind and embed is None:
        details["embed_reduction"] = embed_kind

    # cluster
    cluster_column = _pick_cluster_column(obs)
    cluster_count = None
    if cluster_column:
        labels = [str(v) for v in obs[cluster_column].tolist()]
        cluster_count = len(set(labels))
        details["cluster_column"] = cluster_column
        details["n_clusters"] = cluster_count
        detected.update({"pca", "neighbors", "cluster"})
        if "clusters" not in obs.columns or cluster_column != "clusters":
            if "clusters" not in obs.columns:
                dataset.obs["clusters"] = pd.Series(labels, index=obs.index, dtype=object)
        if graphs:
            detected.add("neighbors")

    # markers（部分对象把结果存 misc/uns；R 桥给出 rank_genes 标志）
    if import_state.get("rank_genes"):
        detected.add("markers")
        details["rank_genes"] = True
    uns = getattr(dataset, "uns", None) or {}
    if isinstance(uns, dict):
        for key in ("rank_genes_groups", "markers", "deg"):
            if key in uns:
                detected.add("markers")
                details["markers_key"] = key

    # cell type / annotate
    celltype_column = _pick_celltype_column(obs, cluster_column)
    if celltype_column:
        labels = [str(v) for v in obs[celltype_column].tolist()]
        if "cell_type" not in obs.columns or celltype_column != "cell_type":
            if "cell_type" not in obs.columns:
                dataset.obs["cell_type"] = pd.Series(labels, index=obs.index, dtype=object)
        dataset.cell_type_default = False
        detected.add("annotate")
        details["celltype_column"] = celltype_column
        details["n_celltypes"] = len(set(labels))

    # 只有出现"过程产物"才算续跑；仅有标准化（很多 h5ad 已 log 化）不算跑过一半
    evidence = detected & {"hvg", "pca", "neighbors", "cluster", "embed", "markers", "annotate"}
    if not evidence:
        return {
            "steps": [],
            "current": None,
            "details": dict(details, note="未检测到已有分析步骤（按原始数据从头开始）"),
            "mapped": {"cluster_column": None, "celltype_column": None, "n_clusters": None,
                       "n_variable_features": len(variable_features), "has_pca": False, "has_embed": False},
        }
    # 步骤取到"检测到的最高步"的连续前缀（保证依赖完整，可继续往下跑）
    highest = -1
    for index, step in enumerate(STEP_ORDER):
        if step in detected:
            highest = index
    steps = [step for step in STEP_ORDER[: highest + 1]] if highest >= 0 else []
    return {
        "steps": steps,
        "current": steps[-1] if steps else None,
        "details": details,
        "mapped": {
            "cluster_column": cluster_column,
            "celltype_column": celltype_column,
            "n_clusters": cluster_count,
            "n_variable_features": len(variable_features),
            "has_pca": "pca" in detected,
            "has_embed": "embed" in detected,
        },
    }


def _embed_name(reductions):
    lowered = {str(name).lower(): str(name) for name in reductions or []}
    for candidate in ("umap", "tsne", "tsne_2d", "harmony_umap"):
        if candidate in lowered:
            return lowered[candidate]
    for lower, original in lowered.items():
        if "umap" in lower or "tsne" in lower:
            return original
    return None


def apply(dataset, loaded=None):
    """把检测到的进度写回 dataset，返回状态摘要（用于 open payload / UI 提示）。"""
    state = detect(dataset, loaded)
    if isinstance(getattr(dataset, "uns", None), dict):
        dataset.uns["import_state"] = state
    dataset.steps = list(state["steps"])
    if state["current"]:
        dataset.params = dict(getattr(dataset, "params", {}) or {})
        dataset.params["resumed_from"] = state["current"]
    return state


def summary_bullets(state):
    details = state.get("details") or {}
    bullets = []
    if state.get("steps"):
        bullets.append("检测到已有进度：已完成 %s，可从「%s」之后继续" % (
            " → ".join(state["steps"]), state["current"]))
    else:
        bullets.append("未检测到已有分析步骤，将从头开始")
    mapped = state.get("mapped") or {}
    if mapped.get("cluster_column"):
        bullets.append("已有聚类列：%s（%s 个簇）" % (mapped["cluster_column"], mapped.get("n_clusters")))
    if mapped.get("celltype_column"):
        bullets.append("已有细胞类型列：%s" % mapped["celltype_column"])
    if mapped.get("has_embed"):
        bullets.append("已有二维嵌入（%s），嵌入图可直接显示" % (details.get("embed_kind") or "umap"))
    if mapped.get("n_variable_features"):
        bullets.append("已有高变基因：%d 个" % mapped["n_variable_features"])
    return bullets
