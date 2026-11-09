#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""官方 Palantir 适配器（CONTRACT §10.3）。

Palantir 官方流程：
  diffusion maps（palantir.utils.run_diffusion_maps）
  → multiscale space（determine_multiscale_space）
  → 伪时间/终端状态/命运概率（palantir.core.run_palantir）。

CLI：
  <externalPython> palantir.py --bundle DIR --out DIR --params JSON

params（可选）：
  root_cell      根细胞名（barcodes.tsv 中的名字）；留空自动选"表达基因数最多"的细胞
                 （CytoTRACE 式启发，写在 notes 里）
  n_components   扩散分量数（默认 10）
  knn            kNN 数（默认 30）
  num_waypoints  路标点数量（默认 500，越小越快）
  n_jobs         并行数（默认 4）
  n_pcs          用多少主成分（默认 30；bundle 有 pca.csv 时直接用）
  max_cells      规模上限（默认 5000，超限确定性下采样）

输出：<out>/result.json + 伪时间着色图 + 簇级伪时间分布图 + obs_extra.csv
（palantir_pseudotime / palantir_entropy / palantir_fate_<终态>）。
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

_EXT_DIR = str(Path(__file__).resolve().parent)
sys.path.insert(0, _EXT_DIR)
import pyext_common as C  # noqa: E402

# 适配器文件名与第三方包同名（pyscenic.py / cellphonedb.py / palantir.py …），
# 运行脚本时 Python 会把本目录放进 sys.path[0]，若不移除，import palantir 会导入
# 适配器自身而不是官方包（自遮蔽）。pyext_common 已导入，随后移除本目录。
while _EXT_DIR in sys.path:
    sys.path.remove(_EXT_DIR)


def _int_param(params, key, default):
    value = params.get(key)
    if value is None or value == "":
        return int(default)
    return int(float(value))


def _build_lognorm(bundle):
    """返回 (cells×genes 的 log-normalized 稠密矩阵, 是否内部做过标准化)。"""
    import numpy as np

    matrix = bundle.matrix("data")
    if bundle.x_kind == "lognorm":
        return matrix, False
    counts = bundle.matrix("counts")
    row_sum = np.asarray(counts.sum(axis=1)).ravel()
    row_sum[row_sum <= 0] = 1.0
    scaled = counts.multiply((1e4 / row_sum)[:, None]).tocsr()
    scaled.data = np.log1p(scaled.data)
    return scaled.astype("float32"), True


def run(bundle_dir, out, params):
    started = time.time()
    out = Path(out)
    warnings = []
    notes = []

    palantir = C.require_module("palantir", "palantir")
    # palantir 的 __init__ 不自动导入子模块，必须显式 import utils/core
    import palantir.core
    import palantir.utils
    C.require_module("scanpy", "scanpy")
    C.require_module("anndata", "anndata")
    import anndata as ad
    import numpy as np
    import pandas as pd

    n_components = _int_param(params, "n_components", 10)
    knn = _int_param(params, "knn", 30)
    num_waypoints = _int_param(params, "num_waypoints", 500)
    n_jobs = _int_param(params, "n_jobs", 4)
    n_pcs = _int_param(params, "n_pcs", 30)
    max_cells = _int_param(params, "max_cells", 5000)
    seed = _int_param(params, "seed", 42)

    bundle = C.Bundle(bundle_dir)
    if bundle.n_obs < 50 or bundle.n_vars < 50:
        raise C.PrerequisiteError("Palantir 至少需要 50 个细胞 × 50 个基因（当前 %d × %d）"
                                  % (bundle.n_obs, bundle.n_vars))
    cell_idx = C.cap_cells(bundle, max_cells, seed=seed, warnings=warnings, label="Palantir")
    matrix, internal_norm = _build_lognorm(bundle)
    if internal_norm:
        warnings.append("bundle 的 X 不是 log-normalized，适配器内部按 log1p(CP10K) 处理")
    if cell_idx is not None:
        matrix = matrix[cell_idx]
    genes = [str(item) for item in bundle.features]
    barcodes = [str(item) for item in bundle.barcodes]
    if cell_idx is not None:
        barcodes = [barcodes[index] for index in cell_idx]
    obs = bundle.obs.iloc[cell_idx] if cell_idx is not None else bundle.obs

    adata = ad.AnnData(
        X=np.asarray(matrix.toarray() if hasattr(matrix, "toarray") else matrix, dtype="float32"),
        obs=pd.DataFrame(index=pd.Index(barcodes, name="cell")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    for column in obs.columns:
        try:
            adata.obs[str(column)] = obs[column].reindex(adata.obs.index).astype(str).to_numpy()
        except Exception:
            pass

    # ---- PCA：优先用 bundle 的 pca.csv，否则用 scanpy 现场算 -----------------
    pca = bundle.pca()
    if pca is not None:
        pca = pca[cell_idx] if cell_idx is not None else pca
        if pca.shape[0] == adata.n_obs and pca.shape[1] >= 2:
            adata.obsm["X_pca"] = np.asarray(pca[:, : min(n_pcs, pca.shape[1])], dtype="float32")
            notes.append("使用 bundle 的 pca.csv（%d 个主成分）" % adata.obsm["X_pca"].shape[1])
    if "X_pca" not in adata.obsm:
        import scanpy as sc

        n_comps = int(min(n_pcs, adata.n_vars - 1, adata.n_obs - 1, 50))
        if n_comps < 2:
            raise C.PrerequisiteError("基因/细胞太少，无法计算 PCA（n_comps=%d）" % n_comps)
        sc.pp.pca(adata, n_comps=n_comps, random_state=seed)
        notes.append("bundle 无 pca.csv，用 scanpy PCA 现算 %d 个主成分" % n_comps)

    # ---- 根细胞 -----------------------------------------------------------
    root_cell = params.get("root_cell")
    if root_cell:
        root_cell = str(root_cell)
        if root_cell not in set(adata.obs_names):
            raise C.PrerequisiteError(
                "params.root_cell=%s 不在 bundle 的 barcodes 中（示例：%s）"
                % (root_cell, ", ".join(list(adata.obs_names)[:3]))
            )
        root_note = "根细胞由 params.root_cell 指定：%s" % root_cell
    else:
        n_genes_per_cell = np.asarray((adata.X > 0).sum(axis=1)).ravel()
        root_cell = str(adata.obs_names[int(np.argmax(n_genes_per_cell))])
        root_note = ("自动选根细胞 %s（表达基因数最多，CytoTRACE 式启发；"
                     "可用 params.root_cell 覆盖）" % root_cell)
    notes.append(root_note)

    # ---- 官方 Palantir 流程 ------------------------------------------------
    t0 = time.time()
    try:
        C.log("Palantir: 计算 diffusion maps…")
        dm_res = palantir.utils.run_diffusion_maps(
            adata, n_components=int(min(n_components, adata.n_obs - 1)), knn=knn, seed=seed
        )
        C.log("Palantir: multiscale space…")
        palantir.utils.determine_multiscale_space(adata)
        C.log("Palantir: 伪时间/终端状态…")
        palantir.core.run_palantir(
            adata, early_cell=root_cell, num_waypoints=num_waypoints,
            n_jobs=n_jobs, seed=seed,
        )
    except Exception as error:
        raise C.PrerequisiteError("Palantir 运行失败（数据可能不满足扩散图要求）：%s" % error)
    palantir_seconds = time.time() - t0
    notes.append("Palantir %s：diffusion maps + run_palantir 耗时 %.1fs"
                 % (C.version_of(palantir), palantir_seconds))

    if "palantir_pseudotime" not in adata.obs.columns:
        raise C.PrerequisiteError("Palantir 没有返回伪时间（palantir_pseudotime）")
    pseudotime = pd.to_numeric(adata.obs["palantir_pseudotime"], errors="coerce").to_numpy(dtype="float64")
    entropy = pd.to_numeric(adata.obs.get("palantir_entropy", pd.Series(index=adata.obs_names)), errors="coerce").to_numpy(dtype="float64")
    fate = adata.obsm.get("palantir_fate_probabilities")
    terminal_cells = []
    if fate is not None and getattr(fate, "shape", (0, 0))[1] > 0:
        terminal_cells = [str(item) for item in getattr(fate, "columns", [])]
    if not terminal_cells:
        warnings.append("Palantir 未检测到终端状态（数据更像离散细胞类型而非连续发育轨迹）；"
                        "仍输出了伪时间与熵，但没有命运概率列。")

    # ---- 图 1：伪时间着色嵌入图 --------------------------------------------
    plt = C.setup_matplotlib()
    embed = bundle.embed()
    if embed is not None:
        embed = embed[cell_idx] if cell_idx is not None else embed
    if embed is None or embed.shape[0] != adata.n_obs:
        import scanpy as sc

        C.log("Palantir: bundle 无 embed.csv，重算 UMAP 用于出图…")
        sc.pp.neighbors(adata, n_neighbors=int(min(30, adata.n_obs - 1)), use_rep="X_pca", random_state=seed)
        sc.tl.umap(adata, random_state=seed)
        embed = np.asarray(adata.obsm["X_umap"])
        notes.append("bundle 无 embed.csv，用 scanpy UMAP 现算嵌入用于出图")
    fig, axis = plt.subplots(figsize=(6.4, 5.2), dpi=140)
    scatter = axis.scatter(embed[:, 0], embed[:, 1], c=pseudotime, s=6, cmap="viridis", linewidths=0)
    colorbar = fig.colorbar(scatter, ax=axis, shrink=0.8)
    colorbar.set_label("Palantir pseudotime", fontsize=8)
    if terminal_cells:
        positions = {str(name): index for index, name in enumerate(adata.obs_names)}
        terminal_xy = [embed[positions[cell]] for cell in terminal_cells if cell in positions]
        if terminal_xy:
            terminal_xy = np.asarray(terminal_xy)
            axis.scatter(terminal_xy[:, 0], terminal_xy[:, 1], s=90, marker="*",
                         facecolor="white", edgecolor="black", linewidths=1.0, label="Terminal state")
            axis.legend(loc="best", fontsize=7, frameon=False)
    root_index = list(adata.obs_names).index(str(root_cell))
    axis.scatter(embed[root_index, 0], embed[root_index, 1], s=110, marker="X",
                 facecolor="red", edgecolor="white", linewidths=0.8, label="Root cell")
    axis.set_xlabel("UMAP1", fontsize=9)
    axis.set_ylabel("UMAP2", fontsize=9)
    axis.set_title("Palantir pseudotime (root: %s)" % root_cell, fontsize=10)
    figures = [C.figure(out, fig, "palantir_pseudotime", "pseudotime_embedding",
                        "Palantir pseudotime on embedding")]

    # ---- 图 2：簇级伪时间分布 ----------------------------------------------
    cluster_col = bundle.obs_columns_like(["clusters", "cluster", "cell_type"])
    if cluster_col:
        labels = obs[cluster_col].astype(str).to_numpy()
        clusters = []
        for label in labels:
            if label not in clusters:
                clusters.append(label)
        clusters.sort()
        grouped = []
        for cluster in clusters:
            mask = labels == cluster
            grouped.append(pseudotime[mask])
        fig, axis = plt.subplots(figsize=(max(5.0, 0.45 * len(clusters) + 2.0), 4.0), dpi=140)
        axis.boxplot(grouped, tick_labels=clusters, showfliers=False, widths=0.6)
        axis.set_ylabel("Palantir pseudotime", fontsize=9)
        axis.set_xlabel("Cluster (%s)" % cluster_col, fontsize=9)
        axis.set_title("Pseudotime distribution per cluster", fontsize=10)
        axis.tick_params(axis="x", labelrotation=45)
        figures.append(C.figure(out, fig, "palantir_pseudotime_by_cluster", "pseudotime_by_cluster",
                                "Palantir pseudotime per cluster"))
    else:
        fig, axis = plt.subplots(figsize=(5.2, 3.8), dpi=140)
        axis.hist(pseudotime, bins=40, color="#4C72B0", edgecolor="white", linewidth=0.3)
        axis.set_xlabel("Palantir pseudotime", fontsize=9)
        axis.set_ylabel("Number of cells", fontsize=9)
        axis.set_title("Palantir pseudotime distribution", fontsize=10)
        figures.append(C.figure(out, fig, "palantir_pseudotime_hist", "pseudotime_hist",
                                "Palantir pseudotime distribution"))

    # ---- 伪时间相关基因表（Pearson）----------------------------------------
    values = np.asarray(adata.X, dtype="float64")
    time_centered = pseudotime - pseudotime.mean()
    time_norm = np.sqrt((time_centered ** 2).sum()) or 1.0
    expression = values - values.mean(axis=0, keepdims=True)
    expression_norm = np.sqrt((expression ** 2).sum(axis=0))
    expression_norm[expression_norm == 0] = 1.0
    correlations = (time_centered @ expression) / (time_norm * expression_norm)
    order = np.argsort(correlations)[::-1]
    gene_rows = []
    for index in list(order[:25]) + list(order[-25:][::-1]):
        gene_rows.append({"gene": genes[int(index)],
                          "correlation": round(float(correlations[int(index)]), 4)})
    gene_rows.sort(key=lambda item: -abs(item["correlation"]))

    # ---- obs_extra ---------------------------------------------------------
    obs_columns = {
        "palantir_pseudotime": [None if not np.isfinite(value) else round(float(value), 6) for value in pseudotime],
        "palantir_entropy": [None if not np.isfinite(value) else round(float(value), 6) for value in entropy],
    }
    for column in terminal_cells[:8]:
        sanitized = "palantir_fate_" + str(column).replace("-", "_").replace(" ", "_")
        try:
            probabilities = np.asarray(fate[column], dtype="float64")
            obs_columns[sanitized] = [None if not np.isfinite(value) else round(float(value), 6)
                                      for value in probabilities]
        except Exception:
            continue
    obs_file = C.write_obs_extra(out, list(adata.obs_names), obs_columns)

    # ---- 终端状态表 --------------------------------------------------------
    tables = [C.table(
        [C.col("gene", "Gene"), C.col("correlation", "Pearson r with pseudotime", "number")],
        gene_rows, "pseudotime_genes", "与伪时间最相关的基因（Pearson）",
        note="基于 log-normalized 表达（bundle data.mtx 或内部 log1p CP10K）",
    )]
    if terminal_cells:
        terminal_rows = []
        for column in terminal_cells[:20]:
            try:
                probabilities = np.asarray(fate[column], dtype="float64")
                terminal_rows.append({
                    "terminal_state": str(column),
                    "n_cells_prob_gt_0.5": int((probabilities > 0.5).sum()),
                    "max_probability": round(float(np.nanmax(probabilities)), 4),
                })
            except Exception:
                continue
        tables.append(C.table(
            [C.col("terminal_state", "Terminal state"), C.col("n_cells_prob_gt_0.5", "Cells (p>0.5)", "number"),
             C.col("max_probability", "Max probability", "number")],
            terminal_rows, "terminal_states", "Palantir 终端状态",
        ))

    total_seconds = time.time() - started
    bullets = [
        "Palantir %s 官方流程：%d 细胞 × %d 基因，扩散分量 %d"
        % (C.version_of(palantir), adata.n_obs, adata.n_vars, min(n_components, adata.n_obs - 1)),
        "根细胞：%s；伪时间范围 %.3f–%.3f" % (root_cell, float(np.nanmin(pseudotime)), float(np.nanmax(pseudotime))),
        "终端状态 %d 个；obs_extra 写入 %d 个细胞级指标" % (len(terminal_cells), len(obs_columns)),
        "总耗时 %.1f 秒" % total_seconds,
    ]
    metrics = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_components": int(min(n_components, adata.n_obs - 1)),
        "n_terminal_states": int(len(terminal_cells)),
        "pseudotime_min": round(float(np.nanmin(pseudotime)), 5),
        "pseudotime_max": round(float(np.nanmax(pseudotime)), 5),
        "root_cell": root_cell,
    }
    return {
        "id": "ext_palantir",
        "title": "Palantir 官方轨迹",
        "method": C.version_of(palantir),
        "summary": {"bullets": bullets, "metrics": metrics},
        "tables": tables,
        "figures": figures,
        "obs_file": obs_file,
        "warnings": warnings,
        "notes": notes,
    }


def main():
    return C.run_adapter(run)


if __name__ == "__main__":
    sys.exit(main())
