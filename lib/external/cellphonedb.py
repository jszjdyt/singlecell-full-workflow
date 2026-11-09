#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""官方 CellPhoneDB 适配器（CONTRACT §10.3）。

使用 CellPhoneDB 5.0.x 官方 Python API 的 statistical_analysis
（permutation test：打乱细胞类型标签估计 ligand-receptor 表达的显著性）。
官方数据库 cellphonedb.zip（v5.0.0）由安装阶段下载到插件缓存目录。

CLI：
  <externalPython> cellphonedb.py --bundle DIR --out DIR --params JSON

params（可选）：
  cell_type_col  细胞类型列（默认 clusters，其次 cluster/cell_type）
  n_perms        置换次数（对应官方 iterations，默认 100）
  pvalue         显著性阈值（默认 0.05）
  threshold      表达百分比阈值（官方默认 0.1）
  min_cells      每个细胞类型最少细胞数（默认 10，少于则从分析中剔除）
  n_jobs         threads（默认 4）
  db_path        官方 cellphonedb.zip 路径（默认自动探测插件缓存/用户目录）
  max_cells      规模上限（默认 5000，超限确定性下采样）
  seed           debug_seed（默认 42，确定性）

输出：<out>/result.json + 互作热图 + 互作图/表格；无 obs_extra（CellPhoneDB 是簇级分析）。
"""

from __future__ import annotations

import inspect
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


def _float_param(params, key, default):
    value = params.get(key)
    if value is None or value == "":
        return float(default)
    return float(value)


CPDB_DB_CANDIDATES = [
    C.HOME / ".dsh" / "sc-plugin-pydata" / "cellphonedb" / "cellphonedb.zip",
    C.HOME / "cellphonedb_data" / "cellphonedb.zip",
    C.HOME / "cpdb" / "cellphonedb.zip",
    C.HOME / ".cpdb" / "cellphonedb.zip",
    C.HOME / ".cellphonedb" / "cellphonedb.zip",
]


def _find_database(params):
    explicit = (params.get("db_path") or params.get("cpdb_file_path")
                or os.environ.get("DS_SC_CELLPHONEDB_DB")
                or os.environ.get("CELLPHONEDB_DATABASE_PATH"))
    candidates = []
    if explicit:
        explicit = os.path.expanduser(str(explicit))
        candidates.append(explicit)
        # 也允许直接给数据库目录
        candidates.append(os.path.join(explicit, "cellphonedb.zip"))
    candidates.extend(str(item) for item in CPDB_DB_CANDIDATES)
    for candidate in candidates:
        path = Path(os.path.expanduser(candidate))
        if path.is_file() and path.suffix == ".zip":
            return path
    try:
        import cellphonedb

        package_root = Path(cellphonedb.__file__).resolve().parent
        for found in sorted(package_root.rglob("cellphonedb*.zip")):
            return found
    except Exception:
        pass
    raise C.PrerequisiteError(
        "未找到 CellPhoneDB 官方数据库（cellphonedb.zip）。安装命令：\n"
        "  /opt/miniconda3/bin/python3 -c \"from cellphonedb.utils import db_utils; "
        "db_utils.download_database(os.path.expanduser('~/.dsh/sc-plugin-pydata/cellphonedb'), 'v5.0.0')\"\n"
        "或用 params.db_path 指定 zip 路径 / 设置环境变量 DS_SC_CELLPHONEDB_DB。"
    )


def _pick_cluster_column(bundle, params):
    column = bundle.obs_columns_like([
        str(params.get("cell_type_col") or ""), "clusters", "cluster", "cell_type",
    ])
    if column is None:
        raise C.PrerequisiteError(
            "bundle 的 obs.csv 里找不到细胞类型列（尝试了 cell_type_col/clusters/cluster/cell_type），"
            "CellPhoneDB 需要先完成聚类。"
        )
    return column


CPDB_META_COLUMNS = {
    "id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b",
    "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin",
    "directionality", "classification", "rank", "is_complex", "directionality",
}


def _to_matrix(frame):
    """CellPhoneDB 5.x 的返回表是"元数据列 + 细胞类型对列"的宽表（v4 是 interacting_pair 作索引）。"""
    import numpy as np
    import pandas as pd

    if not isinstance(frame, pd.DataFrame):
        return None, []
    pair_columns = [str(column) for column in frame.columns if str(column) not in CPDB_META_COLUMNS]
    if "interacting_pair" in frame.columns:
        matrix = frame.set_index("interacting_pair")[pair_columns]
    else:
        matrix = frame[pair_columns]
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix.columns = pair_columns
    return matrix, pair_columns


def _format_tables(results, pvalue_threshold):
    """从 CellPhoneDB 返回 dict 里整理 means / pvalues / significant_means。"""
    import numpy as np
    import pandas as pd

    def _get(key):
        value = results.get(key)
        return value if isinstance(value, pd.DataFrame) else None

    pvalues = _get("pvalues")
    means = _get("means")
    if pvalues is None:
        raise C.AdapterError("CellPhoneDB 返回结果缺少 pvalues 表")
    pvalues, pair_columns = _to_matrix(pvalues)
    means, _ = _to_matrix(means) if means is not None else (None, [])
    if pvalues is None or pvalues.shape[0] == 0:
        raise C.AdapterError("CellPhoneDB pvalues 表为空")
    if means is None or means.shape[0] != pvalues.shape[0]:
        means = pvalues.copy() * np.nan
    means = means.reindex(index=pvalues.index, columns=pair_columns)
    means = means.apply(pd.to_numeric, errors="coerce")
    interaction_names = [str(item) for item in pvalues.index]

    min_pvalues = pvalues.min(axis=1)
    n_significant_pairs = (pvalues < float(pvalue_threshold)).sum(axis=1)
    max_means = means.max(axis=1)
    rows = []
    order = np.lexsort((-max_means.to_numpy(dtype="float64"),
                        n_significant_pairs.to_numpy(), min_pvalues.to_numpy()))
    for position in order:
        if len(rows) >= 60:
            break
        rows.append({
            "interacting_pair": interaction_names[int(position)],
            "min_pvalue": round(float(min_pvalues.iloc[int(position)]), 4),
            "n_significant_pairs": int(n_significant_pairs.iloc[int(position)]),
            "max_mean": round(float(max_means.iloc[int(position)]), 4),
        })
    cell_pair_rows = []
    for column in pair_columns:
        cell_pair_rows.append({
            "cell_pair": str(column),
            "n_significant_interactions": int((pvalues[column] < float(pvalue_threshold)).sum()),
        })
    cell_pair_rows.sort(key=lambda item: -item["n_significant_interactions"])
    significant_mask = (pvalues < float(pvalue_threshold)).to_numpy()
    metrics = {
        "n_interactions_tested": int(pvalues.shape[0]),
        "n_cell_type_pairs": int(pvalues.shape[1]),
        "n_significant_entries": int(significant_mask.sum()),
        "n_significant_interactions": int((significant_mask.sum(axis=1) > 0).sum()),
        "n_significant_cell_pairs": len(pair_columns),
    }
    return {
        "pvalues": pvalues, "means": means, "interaction_rows": rows,
        "cell_pair_rows": cell_pair_rows, "metrics": metrics,
        "interaction_names": interaction_names, "pair_columns": pair_columns,
    }


def run(bundle_dir, out, params):
    started = time.time()
    out = Path(out)
    warnings = []
    notes = []

    try:
        from cellphonedb.src.core.methods import cpdb_statistical_analysis_method as cpdb_method
    except Exception as error:
        raise C.PrerequisiteError("import CellPhoneDB 失败：%s" % error)
    import cellphonedb
    import numpy as np
    import pandas as pd

    n_perms = _int_param(params, "n_perms", 100)
    pvalue_threshold = _float_param(params, "pvalue", 0.05)
    threshold = _float_param(params, "threshold", 0.1)
    min_cells = _int_param(params, "min_cells", 10)
    n_jobs = _int_param(params, "n_jobs", 4)
    max_cells = _int_param(params, "max_cells", 5000)
    seed = _int_param(params, "seed", 42)
    version = C.version_of(cellphonedb)

    database = _find_database(params)
    notes.append("CellPhoneDB %s（官方 statistical_analysis permutation test）" % version)
    notes.append("官方数据库：%s" % database)
    notes.append("参数：iterations=%d, pvalue=%.3f, threshold=%.2f, threads=%d, seed=%d"
                 % (n_perms, pvalue_threshold, threshold, n_jobs, seed))

    bundle = C.Bundle(bundle_dir)
    if bundle.n_obs < 20:
        raise C.PrerequisiteError("CellPhoneDB 至少需要 20 个细胞（当前 %d）" % bundle.n_obs)
    cluster_column = _pick_cluster_column(bundle, params)
    cell_idx = C.cap_cells(bundle, max_cells, seed=seed, warnings=warnings, label="CellPhoneDB")
    counts = bundle.matrix("counts")
    if cell_idx is not None:
        counts = counts[cell_idx]
    genes = [str(item) for item in bundle.features]
    barcodes = [str(item) for item in bundle.barcodes]
    if cell_idx is not None:
        barcodes = [barcodes[index] for index in cell_idx]
    obs = bundle.obs.iloc[cell_idx] if cell_idx is not None else bundle.obs
    clusters = np.asarray([str(item) for item in obs[cluster_column].tolist()], dtype=object)

    # 每个细胞类型至少 min_cells 个细胞（官方要求 / 稳定性）
    counts_by_cluster = pd.Series(clusters).value_counts()
    kept_clusters = [str(item) for item, number in counts_by_cluster.items() if int(number) >= min_cells]
    if not kept_clusters:
        raise C.PrerequisiteError(
            "没有细胞类型满足 min_cells=%d（各类型细胞数：%s）"
            % (min_cells, ", ".join("%s=%d" % (item, int(number)) for item, number in counts_by_cluster.items()))
        )
    dropped = [str(item) for item in counts_by_cluster.index if str(item) not in set(kept_clusters)]
    if dropped:
        warnings.append("细胞数 < %d 的细胞类型已从 CellPhoneDB 分析中剔除：%s"
                        % (min_cells, ", ".join(dropped)))
    keep_mask = np.asarray([label in set(kept_clusters) for label in clusters])
    counts = counts[keep_mask]
    barcodes = [barcodes[index] for index, keep in enumerate(keep_mask.tolist()) if keep]
    clusters = clusters[keep_mask]
    if len(barcodes) < 20:
        raise C.PrerequisiteError("过滤后只剩 %d 个细胞，CellPhoneDB 无法分析" % len(barcodes))
    notes.append("分析规模：%d 细胞 × %d 基因 × %d 个细胞类型（列 %s）"
                 % (len(barcodes), len(genes), len(kept_clusters), cluster_column))
    # 只保留 CellPhoneDB 官方数据库里出现的基因：既保证结果不变（DB 外的基因会被官方
    # 实现丢弃），又能把真实数据（~2 万基因）的输入文件从几百 MB 降到几 MB。
    try:
        from cellphonedb.utils import db_utils as cpdb_db_utils

        _interactions, db_genes, _cc, _ce, _syn, _r2t = cpdb_db_utils.get_interactions_genes_complex(
            str(database)
        )
        allowed = set()
        for column in ("hgnc_symbol", "gene_name"):
            if column in db_genes.columns:
                allowed |= set(str(value) for value in db_genes[column].dropna().tolist())
        if 50 <= len(allowed):
            keep_positions = [index for index, gene in enumerate(genes) if gene in allowed]
            if 50 <= len(keep_positions) < len(genes):
                counts = counts[:, keep_positions]
                genes = [genes[index] for index in keep_positions]
                notes.append("按 CellPhoneDB 官方数据库基因表过滤 counts：%d 个候选基因 → 保留 %d 个"
                             % (len(allowed), len(genes)))
    except Exception as error:
        warnings.append("按 CellPhoneDB 数据库过滤基因失败（%s），使用全部基因（输入文件会更大）" % error)

    # ---- 官方要求的两张输入表 ---------------------------------------------
    meta_frame = pd.DataFrame({"Cell": barcodes, "cell_type": clusters})
    meta_path = out / "cellphonedb_meta.txt"
    meta_frame.to_csv(meta_path, sep="\t", index=False, encoding="utf-8")
    counts_dense = np.asarray(counts.toarray() if hasattr(counts, "toarray") else counts)
    counts_frame = pd.DataFrame(counts_dense.T, index=genes, columns=barcodes)
    counts_frame.index.name = "Gene"
    counts_path = out / "cellphonedb_counts.txt"
    counts_frame.to_csv(counts_path, sep="\t", encoding="utf-8")
    notes.append("输入文件：%s（%d 细胞）、%s（%d 基因 × %d 细胞）"
                 % (meta_path.name, len(barcodes), counts_path.name, len(genes), len(barcodes)))

    # ---- 调用官方 Python API（按签名过滤参数，兼容 5.0.x）-----------------
    output_dir = out / "cellphonedb_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    call_kwargs = {
        "cpdb_file_path": str(database),
        "meta_file_path": str(meta_path),
        "counts_file_path": str(counts_path),
        "counts_data": "hgnc_symbol",
        "output_path": str(output_dir),
        "iterations": int(n_perms),
        "threshold": float(threshold),
        "threads": int(n_jobs),
        "debug_seed": int(seed),
        "result_precision": 3,
        "pvalue": float(pvalue_threshold),
        "separator": "|",
        "debug": False,
    }
    try:
        signature = inspect.signature(cpdb_method.call)
        accepted = {name for name in signature.parameters}
        filtered = {key: value for key, value in call_kwargs.items() if key in accepted}
        dropped_kwargs = sorted(set(call_kwargs) - set(filtered))
        if dropped_kwargs:
            warnings.append("当前 CellPhoneDB 版本不接受参数：%s（已自动忽略）" % ", ".join(dropped_kwargs))
        C.log("CellPhoneDB: 运行官方 statistical_analysis（iterations=%d, threads=%d）…"
              % (n_perms, n_jobs))
        results = cpdb_method.call(**filtered)
    except C.PrerequisiteError:
        raise
    except Exception as error:
        message = str(error)
        if "human" in message.lower() or "AllCountsFiltered" in type(error).__name__:
            raise C.PrerequisiteError("CellPhoneDB 过滤后没有可用基因（需人类基因 symbol 数据）：%s" % message)
        raise C.PrerequisiteError("CellPhoneDB statistical_analysis 运行失败：%s" % message)
    if not isinstance(results, dict):
        raise C.AdapterError("CellPhoneDB 返回结果不是 dict（%s）" % type(results))
    elapsed_cpdb = time.time() - started

    parsed = _format_tables(results, pvalue_threshold)
    metrics = parsed["metrics"]
    metrics.update({
        "n_cells": int(len(barcodes)),
        "n_genes": int(len(genes)),
        "n_cell_types": int(len(kept_clusters)),
        "iterations": int(n_perms),
    })
    # 落一份完整表格方便复用
    try:
        parsed["pvalues"].to_csv(out / "cellphonedb_pvalues.tsv", sep="\t", encoding="utf-8")
        parsed["means"].to_csv(out / "cellphonedb_means.tsv", sep="\t", encoding="utf-8")
    except Exception:
        pass

    # ---- 热图：top 显著互作 × 细胞类型对 -----------------------------------
    plt = C.setup_matplotlib()
    figures = []
    pvalues = parsed["pvalues"]
    means = parsed["means"]
    pair_columns = parsed["pair_columns"]
    significant_order = []
    pvalues_values = pvalues.to_numpy(dtype="float64")
    means_values = means.to_numpy(dtype="float64")
    for position in np.lexsort((-np.nanmax(means_values, axis=1),
                                np.nansum(pvalues_values < float(pvalue_threshold), axis=1),
                                np.nanmin(pvalues_values, axis=1))):
        significant_order.append(int(position))
    top_positions = [position for position in significant_order
                     if np.nanmin(pvalues_values[position]) < float(pvalue_threshold)][:25]
    if not top_positions:
        top_positions = significant_order[:25]
        warnings.append("没有 p < %.3f 的显著互作（置换次数 %d），热图改用 p 值最小的 %d 个互作"
                        % (pvalue_threshold, n_perms, len(top_positions)))
    matrix = -np.log10(np.clip(pvalues_values[top_positions], 1e-4, 1.0))
    row_labels = [str(pvalues.index[position]) for position in top_positions]
    fig, axis = plt.subplots(figsize=(max(7.0, 0.34 * len(pair_columns) + 3.0),
                                      max(4.5, 0.24 * len(row_labels) + 1.6)), dpi=140)
    image = axis.imshow(matrix, aspect="auto", cmap="magma_r", vmin=0.0,
                        vmax=float(max(2.0, np.nanmax(matrix) if np.isfinite(np.nanmax(matrix)) else 2.0)))
    axis.set_yticks(range(len(row_labels)))
    axis.set_yticklabels(row_labels, fontsize=6)
    axis.set_xticks(range(len(pair_columns)))
    axis.set_xticklabels(pair_columns, rotation=90, fontsize=6)
    axis.set_xlabel("Cell type pair", fontsize=9)
    axis.set_ylabel("Interacting pair", fontsize=9)
    axis.set_title("CellPhoneDB significant interactions (-log10 p-value)", fontsize=10)
    colorbar = fig.colorbar(image, ax=axis, shrink=0.7)
    colorbar.set_label("-log10(p-value)", fontsize=8)
    figures.append(C.figure(out, fig, "cellphonedb_interactions", "interaction_heatmap",
                            "CellPhoneDB interactions heatmap"))

    # 第二张图：各细胞类型对的显著互作数量
    pair_names = [item["cell_pair"] for item in parsed["cell_pair_rows"]]
    pair_counts = [item["n_significant_interactions"] for item in parsed["cell_pair_rows"]]
    fig, axis = plt.subplots(figsize=(max(6.0, 0.35 * len(pair_names) + 2.4), 3.8), dpi=140)
    axis.bar(range(len(pair_names)), pair_counts, color="#4C72B0", edgecolor="white", linewidth=0.3)
    axis.set_xticks(range(len(pair_names)))
    axis.set_xticklabels(pair_names, rotation=90, fontsize=7)
    axis.set_ylabel("Significant interactions", fontsize=9)
    axis.set_title("Significant interactions per cell type pair", fontsize=10)
    figures.append(C.figure(out, fig, "cellphonedb_pair_counts", "pair_counts",
                            "Significant interactions per cell type pair"))

    total_seconds = time.time() - started
    bullets = [
        "CellPhoneDB %s 官方统计（permutation test）：%d 细胞 × %d 细胞类型 × %d 个互作"
        % (version, len(barcodes), len(kept_clusters), metrics["n_interactions_tested"]),
        "显著（p < %.3f）互作 %d 个 / 显著细胞类型对组合 %d 个"
        % (pvalue_threshold, metrics["n_significant_interactions"], metrics["n_significant_entries"]),
        "置换 %d 次；输出互作表 + 热图（%s）" % (n_perms, database.name),
        "总耗时 %.1f 秒（CellPhoneDB %.1f 秒）" % (total_seconds, elapsed_cpdb),
    ]
    return {
        "id": "ext_cellphonedb",
        "title": "CellPhoneDB 官方通讯",
        "method": version,
        "summary": {"bullets": bullets, "metrics": metrics},
        "tables": [
            C.table(
                [C.col("interacting_pair", "Interacting pair"),
                 C.col("min_pvalue", "Min p-value", "number"),
                 C.col("n_significant_pairs", "Significant cell pairs", "number"),
                 C.col("max_mean", "Max mean expression", "number")],
                parsed["interaction_rows"], "interactions",
                "显著互作（按最小 p 值排序，最多 60 条）",
                note="p 值来自 CellPhoneDB 官方 permutation test；pair 用 | 分隔细胞类型",
            ),
            C.table(
                [C.col("cell_pair", "Cell type pair"),
                 C.col("n_significant_interactions", "Significant interactions", "number")],
                parsed["cell_pair_rows"], "cell_pairs", "各细胞类型对的显著互作数量",
            ),
        ],
        "figures": figures,
        "obs_file": None,
        "warnings": warnings,
        "notes": notes,
    }


def main():
    return C.run_adapter(run)


if __name__ == "__main__":
    sys.exit(main())
