#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""官方 scVelo 适配器（CONTRACT §10.3）。

scVelo 做 RNA velocity **必须**要有 spliced / unspliced 两层计数。当前工作台的
bundle（§10.2）只导出 counts/data，没有层数据时本适配器不会编造结果，而是
**退出码 3 + 中文说明**（引擎转成 prerequisite）。

层文件约定（扩展 §10.2，兼容 Lead 的导出器）：
  <bundle>/spliced.mtx   + <bundle>/unspliced.mtx      （genes × cells，MatrixMarket）
  或 <bundle>/layers/spliced.mtx + <bundle>/layers/unspliced.mtx

CLI：
  <externalPython> scvelo.py --bundle DIR --out DIR --params JSON

params（可选）：
  mode            stochastic（默认）| deterministic | dynamical
  n_top_genes     用于 velocity 的高变基因数（默认 2000）
  n_pcs           默认 30
  n_neighbors     默认 30
  min_shared_counts  默认 20
  n_jobs          默认 4
  max_cells       规模上限（默认 5000，超限确定性下采样）

输出：velocity 图 + obs_extra.csv（velocity_pseudotime / velocity_confidence /
latent_time 等）+ result.json。
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


def _layer_path(bundle_dir, name):
    root = Path(bundle_dir)
    for candidate in (root / ("%s.mtx" % name), root / "layers" / ("%s.mtx" % name)):
        if candidate.is_file():
            return candidate
    return None


def run(bundle_dir, out, params):
    started = time.time()
    out = Path(out)
    warnings = []
    notes = []

    scv = C.require_module("scvelo", "scvelo")
    C.require_module("anndata", "anndata")
    C.require_module("scanpy", "scanpy")

    mode = str(params.get("mode") or "stochastic").lower()
    if mode not in ("stochastic", "deterministic", "dynamical"):
        raise C.PrerequisiteError("params.mode 只支持 stochastic/deterministic/dynamical（收到 %s）" % mode)
    n_top_genes = _int_param(params, "n_top_genes", 2000)
    n_pcs = _int_param(params, "n_pcs", 30)
    n_neighbors = _int_param(params, "n_neighbors", 30)
    min_shared = _int_param(params, "min_shared_counts", 20)
    n_jobs = _int_param(params, "n_jobs", 4)
    max_cells = _int_param(params, "max_cells", 5000)

    # ---- 层数据检查（关键：没有就退出码 3，绝不假造）----------------------
    spliced_path = _layer_path(bundle_dir, "spliced")
    unspliced_path = _layer_path(bundle_dir, "unspliced")
    if spliced_path is None or unspliced_path is None:
        missing = []
        if spliced_path is None:
            missing.append("spliced")
        if unspliced_path is None:
            missing.append("unspliced")
        raise C.PrerequisiteError(
            "bundle 缺少 %s 层（§10.2 只导出 counts/data，没有剪接/未剪接计数）。"
            "RNA velocity 需要 spliced/unspliced 两种计数：请在数据导入时提供 velocyto / "
            "STARsolo Velocyto / kb-python 的 spliced+unspliced 矩阵，并由 Host 以 "
            "<bundle>/spliced.mtx 与 <bundle>/unspliced.mtx（genes×cells）导出后再运行。"
            "在没有该层数据时 scVelo 无法给出可信速度；可改用本插件的 trajectory_dpt/"
            "trajectory_paga 做发育轨迹（不依赖剪接层）。" % "、".join(missing)
        )

    bundle = C.Bundle(bundle_dir)
    cell_idx = C.cap_cells(bundle, max_cells, seed=_int_param(params, "seed", 42),
                           warnings=warnings, label="scVelo")
    spliced = C.read_mtx_matrix(spliced_path)
    unspliced = C.read_mtx_matrix(unspliced_path)
    if spliced.shape != (bundle.n_obs, bundle.n_vars):
        raise C.PrerequisiteError(
            "spliced.mtx 形状 %s 与 bundle（%d × %d 细胞×基因）不一致"
            % (tuple(spliced.shape), bundle.n_obs, bundle.n_vars)
        )
    if unspliced.shape != spliced.shape:
        raise C.PrerequisiteError("unspliced.mtx 形状 %s 与 spliced.mtx %s 不一致"
                                  % (tuple(unspliced.shape), tuple(spliced.shape)))
    if cell_idx is not None:
        spliced = spliced[cell_idx]
        unspliced = unspliced[cell_idx]
    genes = [str(item) for item in bundle.features]
    barcodes = [str(item) for item in bundle.barcodes]
    if cell_idx is not None:
        barcodes = [barcodes[index] for index in cell_idx]
    obs = bundle.obs.iloc[cell_idx] if cell_idx is not None else bundle.obs

    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.AnnData(
        X=spliced.copy(),
        obs=pd.DataFrame(index=pd.Index(barcodes, name="cell")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.layers["spliced"] = spliced.copy()
    adata.layers["unspliced"] = unspliced.copy()
    for column in obs.columns:
        try:
            adata.obs[str(column)] = obs[column].reindex(adata.obs.index).astype(str).to_numpy()
        except Exception:
            pass
    embed = bundle.embed()
    if embed is not None:
        embed = embed[cell_idx] if cell_idx is not None else embed
        if embed.shape[0] == adata.n_obs:
            adata.obsm["X_umap"] = np.asarray(embed, dtype="float32")
            notes.append("使用 bundle 的 embed.csv 作为 2D 嵌入（basis=umap）")
    notes.append("scVelo %s（mode=%s）" % (C.version_of(scv), mode))
    notes.append("层数据：%s / %s（%d 细胞 × %d 基因）" %
                 (spliced_path.name, unspliced_path.name, adata.n_obs, adata.n_vars))

    try:
        # 官方预处理流程（scVelo 0.3.x：filter_and_normalize → 可选 HVG → moments）
        scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared)
        if adata.n_vars < 50:
            raise C.PrerequisiteError(
                "spliced/unspliced 过滤后只剩 %d 个基因（<50），数据过于稀疏，scVelo 无法估计速度"
                % adata.n_vars
            )
        if n_top_genes and adata.n_vars > int(n_top_genes):
            import scanpy as sc

            # scVelo 0.3.x 的 filter_and_normalize 没有 n_top_genes，HVG 用 scanpy 官方实现
            sc.pp.highly_variable_genes(adata, n_top_genes=int(n_top_genes), flavor="seurat")
            adata = adata[:, adata.var["highly_variable"]].copy()
            notes.append("HVG：保留 %d 个高变基因（scanpy seurat flavor）" % adata.n_vars)
        if adata.n_vars < 50:
            raise C.PrerequisiteError("HVG 过滤后只剩 %d 个基因（<50），scVelo 无法估计速度" % adata.n_vars)
        n_pcs_used = int(min(n_pcs, max(2, adata.n_vars - 1), max(2, adata.n_obs - 1)))
        scv.pp.moments(adata, n_pcs=n_pcs_used, n_neighbors=int(min(n_neighbors, max(2, adata.n_obs - 1))))
        if mode == "dynamical":
            scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
            scv.tl.velocity(adata, mode="dynamical")
            scv.tl.velocity_graph(adata)
            scv.tl.latent_time(adata)
        else:
            scv.tl.velocity(adata, mode=mode)
            scv.tl.velocity_graph(adata)
            scv.tl.velocity_pseudotime(adata)
        scv.tl.velocity_confidence(adata)
        if "X_umap" not in adata.obsm:
            scv.tl.umap(adata)
            notes.append("bundle 无 embed.csv，scVelo 重算了 UMAP")
    except C.PrerequisiteError:
        raise
    except Exception as error:
        raise C.PrerequisiteError("scVelo 计算失败（数据可能不满足 velocity 要求）：%s" % error)

    # ---- 图：velocity stream ----
    plt = C.setup_matplotlib()
    color_key = None
    for candidate in ("clusters", "cluster", "cell_type"):
        if candidate in adata.obs.columns:
            color_key = candidate
            break
    figures = []
    try:
        fig, axis = plt.subplots(figsize=(6.8, 5.4), dpi=140)
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color=color_key, ax=axis, show=False,
            legend_loc="right margin" if color_key else "none",
            title="scVelo RNA velocity (%s)" % mode,
        )
        figures.append(C.figure(out, fig, "scvelo_velocity_stream", "velocity_stream",
                                "scVelo RNA velocity stream"))
    except Exception as error:
        warnings.append("velocity stream 图渲染失败（%s），改用散点图" % error)
        plt.close("all")
        fig, axis = plt.subplots(figsize=(6.2, 5.0), dpi=140)
        coords = adata.obsm["X_umap"]
        axis.scatter(coords[:, 0], coords[:, 1], s=4, c="0.6", linewidths=0)
        axis.set_xlabel("UMAP1")
        axis.set_ylabel("UMAP2")
        axis.set_title("scVelo velocity (%s)" % mode, fontsize=10)
        figures.append(C.figure(out, fig, "scvelo_velocity_stream", "velocity_stream",
                                "scVelo RNA velocity"))

    # 第二张图：velocity pseudotime / latent time
    time_key = "latent_time" if "latent_time" in adata.obs.columns else (
        "velocity_pseudotime" if "velocity_pseudotime" in adata.obs.columns else None)
    if time_key:
        try:
            fig, axis = plt.subplots(figsize=(6.0, 4.8), dpi=140)
            coords = adata.obsm["X_umap"]
            values = adata.obs[time_key].to_numpy(dtype="float64")
            scatter = axis.scatter(coords[:, 0], coords[:, 1], c=values, s=5, cmap="viridis", linewidths=0)
            colorbar = fig.colorbar(scatter, ax=axis, shrink=0.8)
            colorbar.set_label(time_key.replace("_", " "), fontsize=8)
            axis.set_xlabel("UMAP1")
            axis.set_ylabel("UMAP2")
            axis.set_title("scVelo %s" % time_key.replace("_", " "), fontsize=10)
            figures.append(C.figure(out, fig, "scvelo_%s" % time_key, time_key,
                                    "scVelo %s" % time_key.replace("_", " ")))
        except Exception as error:
            warnings.append("时间轴图渲染失败：%s" % error)

    # ---- 驱动基因表 ----
    var = adata.var.copy()
    driver_rows = []
    score_column = None
    for candidate in ("fit_likelihood", "velocity_gamma", "fit_alpha", "velocity_genes"):
        if candidate in var.columns:
            score_column = candidate
            break
    if score_column is not None:
        if var[score_column].dtype == bool:
            ranking = var[var[score_column]].index.tolist()
        else:
            ranking = var[score_column].abs().sort_values(ascending=False).index.tolist()
        for gene in ranking[:40]:
            record = var.loc[gene]
            row = {"gene": str(gene), "velocity_gene": bool(record.get("velocity_genes", score_column == "velocity_genes"))}
            for key in ("fit_likelihood", "velocity_gamma", "fit_alpha", "fit_beta", "fit_gamma"):
                if key in var.columns:
                    try:
                        value = float(record[key])
                        row[key] = round(value, 4) if np.isfinite(value) else None
                    except Exception:
                        row[key] = None
            driver_rows.append(row)

    # ---- obs_extra ----
    obs_columns = {}
    for key in ("velocity_pseudotime", "latent_time", "velocity_confidence",
                "velocity_length", "velocity_confidence_transition"):
        if key in adata.obs.columns:
            values = adata.obs[key].tolist()
            cleaned = []
            for value in values:
                try:
                    number = float(value)
                    cleaned.append(round(number, 6) if np.isfinite(number) else None)
                except Exception:
                    cleaned.append(None)
            obs_columns[key] = cleaned
    obs_file = C.write_obs_extra(out, list(adata.obs_names), obs_columns) if obs_columns else None

    total_seconds = time.time() - started
    n_velocity_genes = int(adata.var["velocity_genes"].sum()) if "velocity_genes" in adata.var.columns else 0
    bullets = [
        "scVelo %s 官方 velocity（mode=%s，%d 细胞 × %d 基因）"
        % (C.version_of(scv), mode, adata.n_obs, adata.n_vars),
        "velocity 基因 %d 个；输出 %d 个细胞级指标到 obs_extra" % (n_velocity_genes, len(obs_columns)),
        "总耗时 %.1f 秒" % total_seconds,
    ]
    metrics = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_velocity_genes": n_velocity_genes,
        "mode": mode,
    }
    for key in ("velocity_confidence", "velocity_pseudotime"):
        if key in adata.obs.columns:
            try:
                values = pd.to_numeric(adata.obs[key], errors="coerce")
                metrics["mean_%s" % key] = round(float(values.mean()), 5)
            except Exception:
                pass
    return {
        "id": "ext_scvelo",
        "title": "scVelo 官方 RNA velocity",
        "method": C.version_of(scv),
        "summary": {"bullets": bullets, "metrics": metrics},
        "tables": [C.table(
            [C.col("gene", "Gene"), C.col("velocity_gene", "Velocity gene"),
             C.col("fit_likelihood", "Fit likelihood", "number"),
             C.col("velocity_gamma", "Velocity gamma", "number"),
             C.col("fit_alpha", "Fit alpha", "number")],
            driver_rows, "velocity_genes", "Top velocity 驱动基因",
        )] if driver_rows else [],
        "figures": figures,
        "obs_file": obs_file,
        "warnings": warnings,
        "notes": notes,
    }


def main():
    return C.run_adapter(run)


if __name__ == "__main__":
    sys.exit(main())
