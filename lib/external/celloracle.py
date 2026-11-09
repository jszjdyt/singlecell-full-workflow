#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""官方 CellOracle 适配器 —— 虚拟敲除 / 过表达（CONTRACT §10.3）。

官方流程（CellOracle 0.20）：
  1. import_anndata_as_raw_count：载入原始 counts + 聚类列 + 2D 嵌入；
  2. import_TF_data：载入 human promoter base GRN（hg38 gimmemotifsv5）；
  3. perform_PCA → knn_imputation（kNN 平滑）；
  4. fit_GRN_for_simulation（按簇拟合 GRN）；
  5. simulate_shift（perturb_condition={target: 0.0} 敲除 / 高表达值过表达）；
  6. estimate_transition_prob → calculate_embedding_shift → calculate_grid_arrows。

CLI：
  <externalPython> celloracle.py --bundle DIR --out DIR --params JSON

params：
  target        目标基因（必填；缺省自动选 base GRN 覆盖的最高表达基因）
  mode          knockout（默认）| overexpress
  n_top_genes   参与分析的高变基因数（默认 2000；bundle 基因更少则全用）
  base_grn      base GRN parquet 路径（默认自动探测 ~/celloracle_data/... 或用官方下载器）
  alpha         岭回归正则强度（官方推荐 10，默认 10）
  n_propagation 信号传播步数（默认 3）
  k_imputation  kNN 平滑的 k（留空 = 2.5% 细胞数，官方推荐）
  n_jobs        并行数（默认 4）
  max_cells     规模上限（默认 5000，超限确定性下采样）

输出：<out>/result.json + 扰动箭头图 + 受影响基因条形图 + obs_extra.csv
（celloracle_shift_magnitude / celloracle_shift_x / celloracle_shift_y /
celloracle_delta_<target>）。
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


def _float_param(params, key, default):
    value = params.get(key)
    if value is None or value == "":
        return float(default)
    return float(value)


CELLORACLE_CACHE = C.HOME / "celloracle_data"
HUMAN_BASE_GRN_DEFAULT = (
    CELLORACLE_CACHE / "promoter_base_GRN"
    / "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"
)


def _load_base_grn(params, notes, warnings):
    """返回 (base_grn DataFrame, 来源说明)。探测顺序：params → 缓存 → 官方下载器 → toy GRN。"""
    import pandas as pd

    explicit = params.get("base_grn") or os.environ.get("DS_SC_CELLORACLE_BASE_GRN")
    candidates = []
    if explicit:
        candidates.append(str(explicit))
    candidates.append(str(HUMAN_BASE_GRN_DEFAULT))
    candidates.append(str(CELLORACLE_CACHE / "promoter_base_GRN"
                          / "hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"))
    for candidate in candidates:
        path = Path(os.path.expanduser(candidate))
        if path.is_file():
            try:
                frame = pd.read_parquet(path)
            except Exception as error:
                warnings.append("base GRN parquet 读取失败（%s）：%s" % (path, error))
                continue
            if "peak_id" in frame.columns and "gene_short_name" in frame.columns:
                notes.append("base GRN：%s（%d 行 × %d 列）" % (path, frame.shape[0], frame.shape[1]))
                return frame, "human promoter base GRN（%s）" % path.name
            warnings.append("base GRN 文件缺 peak_id/gene_short_name 列：%s" % path)
    # 官方下载器（会联网）
    try:
        import celloracle as co

        C.log("CellOracle: 本地无 base GRN，尝试官方 data.load_human_promoter_base_GRN() 下载…")
        frame = co.data.load_human_promoter_base_GRN(version="hg38_gimmemotifsv5_fpr2")
        if frame is not None and "peak_id" in frame.columns:
            notes.append("base GRN：官方下载器下载的 human hg38 gimmemotifsv5 fpr2")
            return frame, "human promoter base GRN（官方下载）"
    except Exception as error:
        warnings.append("官方下载 human base GRN 失败：%s" % error)
    # 最后兜底：用本插件 TF 列表 + 共表达生成 toy GRN（仅供功能演示，warnings 明确标注）
    warnings.append(
        "未能获得官方 human base GRN，改用插件内置 toy base GRN（TF 列表 × 高变基因共表达代理，"
        "仅用于功能演示；结果的生物学解释需谨慎，建议下载官方 base GRN 后重跑）。"
    )
    return None, "toy base GRN（插件内置代理，非官方）"


def _toy_base_grn(genes, tf_list, adata_for_corr):
    """构造最小 base GRN（DataFrame: peak_id/gene_short_name/<TF> 0/1）。

    仅当官方 base GRN 不可用时使用：把与每个基因相关性最高的少量 TF 视为潜在调控。
    """
    import numpy as np
    import pandas as pd

    gene_set = set(genes)
    tfs = [tf for tf in (tf_list or []) if tf in gene_set]
    if not tfs:
        tfs = list(genes)[:100]
    values = np.asarray(adata_for_corr, dtype="float64")
    gene_index = {gene: index for index, gene in enumerate(genes)}
    tf_indices = [gene_index[tf] for tf in tfs]
    tf_matrix = values[:, tf_indices]
    # 相关（基因 × TF）
    tf_centered = tf_matrix - tf_matrix.mean(axis=0, keepdims=True)
    tf_norm = np.sqrt((tf_centered ** 2).sum(axis=0))
    tf_norm[tf_norm == 0] = 1.0
    centered = values - values.mean(axis=0, keepdims=True)
    norm = np.sqrt((centered ** 2).sum(axis=0))
    norm[norm == 0] = 1.0
    correlation = (centered / norm[:, None]).T @ (tf_centered / tf_norm[None, :])
    frame = pd.DataFrame(0, index=range(len(genes)), columns=["peak_id", "gene_short_name"] + tfs)
    frame["peak_id"] = ["toy_peak_%d" % index for index in range(len(genes))]
    frame["gene_short_name"] = genes
    top_n = max(3, min(20, len(tfs)))
    for gene_index_pos in range(len(genes)):
        order = np.argsort(np.abs(correlation[gene_index_pos]))[::-1][:top_n]
        for position in order:
            frame.iloc[gene_index_pos, 2 + int(position)] = 1
    return frame


def run(bundle_dir, out, params):
    started = time.time()
    out = Path(out)
    warnings = []
    notes = []

    # ---- 依赖（genomepy/gimmemotifs 等 import 级依赖缺失 → 退出码 3）--------
    try:
        import celloracle as co
    except Exception as error:
        raise C.PrerequisiteError(
            "import celloracle 失败：%s。CellOracle 的 import 级依赖（genomepy/gimmemotifs/"
            "velocyto/louvain 等）可能未装齐，请检查 /opt/miniconda3 环境。" % error
        )
    C.require_module("anndata", "anndata")
    C.require_module("scanpy", "scanpy")
    import anndata as ad
    import numpy as np
    import pandas as pd

    mode = str(params.get("mode") or "knockout").lower()
    if mode not in ("knockout", "overexpress"):
        raise C.PrerequisiteError("params.mode 只支持 knockout / overexpress（收到 %s）" % mode)
    n_top_genes = _int_param(params, "n_top_genes", 2000)
    alpha = _float_param(params, "alpha", 10.0)
    n_propagation = _int_param(params, "n_propagation", 3)
    n_jobs = _int_param(params, "n_jobs", 4)
    max_cells = _int_param(params, "max_cells", 5000)
    seed = _int_param(params, "seed", 42)
    np.random.seed(seed)
    notes.append("CellOracle %s（%s，官方 in silico perturbation）" % (C.version_of(co), mode))

    # ---- 读 bundle --------------------------------------------------------
    bundle = C.Bundle(bundle_dir)
    if bundle.n_obs < 100:
        raise C.PrerequisiteError("CellOracle 至少需要 100 个细胞（当前 %d）" % bundle.n_obs)
    cluster_col = bundle.obs_columns_like(["clusters", "cluster", "cell_type"])
    if cluster_col is None:
        raise C.PrerequisiteError(
            "bundle 的 obs.csv 里没有聚类列（clusters/cluster/cell_type），"
            "CellOracle 需要先完成聚类（clusters）"
        )
    cell_idx = C.cap_cells(bundle, max_cells, seed=seed, warnings=warnings, label="CellOracle")
    counts = bundle.matrix("counts")
    if cell_idx is not None:
        counts = counts[cell_idx]
    genes = [str(item) for item in bundle.features]
    barcodes = [str(item) for item in bundle.barcodes]
    if cell_idx is not None:
        barcodes = [barcodes[index] for index in cell_idx]
    obs = bundle.obs.iloc[cell_idx] if cell_idx is not None else bundle.obs

    # 基因选择：CellOracle 推荐 1000–3000 个高变基因；base GRN 覆盖的 TF 尽量保留
    gene_idx = C.select_genes(counts, genes, n_top_genes, seed=seed, warnings=warnings)
    if gene_idx is not None:
        counts = counts[:, gene_idx]
        genes = [genes[index] for index in gene_idx]
    if len(genes) < 200:
        raise C.PrerequisiteError("可用于 CellOracle 的基因只有 %d 个（至少 200 个）" % len(genes))
    counts_dense = np.asarray(counts.toarray() if hasattr(counts, "toarray") else counts)
    if counts_dense.min() < 0:
        raise C.PrerequisiteError("counts.mtx 含负值，不是原始 counts，无法做 CellOracle GRN 拟合")
    counts_dense = np.rint(counts_dense).astype("float32")
    notes.append("分析规模：%d 细胞 × %d 基因（原始 counts），聚类列 %s（%d 类）"
                 % (len(barcodes), len(genes), cluster_col, int(obs[cluster_col].nunique())))

    # ---- AnnData 准备 ------------------------------------------------------
    embedding = bundle.embed()
    if embedding is not None:
        embedding = embedding[cell_idx] if cell_idx is not None else embedding
    if embedding is None or embedding.shape[0] != len(barcodes):
        C.log("CellOracle: bundle 无 embed.csv，用 scanpy 计算 2D UMAP 作为嵌入…")
        import scanpy as sc

        temp = ad.AnnData(X=counts_dense, obs=pd.DataFrame(index=pd.Index(barcodes, name="cell")),
                          var=pd.DataFrame(index=pd.Index(genes, name="gene")))
        sc.pp.normalize_total(temp, target_sum=1e4)
        sc.pp.log1p(temp)
        sc.pp.pca(temp, n_comps=int(min(30, temp.n_vars - 1, temp.n_obs - 1)), random_state=seed)
        sc.pp.neighbors(temp, n_neighbors=30, random_state=seed)
        sc.tl.umap(temp, random_state=seed)
        embedding = np.asarray(temp.obsm["X_umap"], dtype="float32")
        warnings.append("bundle 无 embed.csv，CellOracle 内部用 scanpy UMAP 计算 2D 嵌入")
        notes.append("2D 嵌入由适配器内部计算（scanpy UMAP）")
    else:
        notes.append("使用 bundle 的 embed.csv 作为 2D 嵌入（X_umap）")

    adata = ad.AnnData(
        X=counts_dense,
        obs=pd.DataFrame(index=pd.Index(barcodes, name="cell")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.obs[cluster_col] = pd.Categorical(obs[cluster_col].astype(str).to_numpy())
    adata.obsm["X_umap"] = np.asarray(embedding, dtype="float32")
    # CellOracle 的 base GRN 使用 hgnc symbol；这里 var_names 就是 symbol
    oracle = co.Oracle()
    try:
        oracle.import_anndata_as_raw_count(
            adata=adata, cluster_column_name=cluster_col, embedding_name="X_umap"
        )
    except Exception as error:
        raise C.PrerequisiteError("CellOracle 导入数据失败：%s" % error)

    # ---- base GRN ----------------------------------------------------------
    base_grn, base_source = _load_base_grn(params, notes, warnings)
    if base_grn is None:
        lib_dir = Path(__file__).resolve().parents[1]
        try:
            if str(lib_dir) not in sys.path:
                sys.path.insert(0, str(lib_dir))
            from biopriors import TF_LIST  # type: ignore

            tf_list = list(TF_LIST)
        except Exception:
            tf_list = None
        base_grn = _toy_base_grn(genes, tf_list, adata.X)
        notes.append("toy base GRN：%d 基因 × %d 候选 TF（%s）"
                     % (base_grn.shape[0], base_grn.shape[1] - 2, base_source))
    try:
        oracle.import_TF_data(TF_info_matrix=base_grn)
    except Exception as error:
        raise C.PrerequisiteError("CellOracle 载入 base GRN 失败：%s" % error)
    n_tf_pairs = sum(len(value) for value in getattr(oracle, "TFdict", {}).values())
    notes.append("base GRN 覆盖：%d 个基因、%d 个 TF-基因潜在调控关系（来源：%s）"
                 % (len(getattr(oracle, "TFdict", {})), n_tf_pairs, base_source))
    if n_tf_pairs == 0:
        raise C.PrerequisiteError(
            "base GRN 与表达矩阵基因名的交集为空（base GRN 用 hgnc symbol，bundle 基因名可能不是 symbol），"
            "CellOracle 无法构建 GRN。"
        )

    # ---- 目标基因（CellOracle 官方限定：目标必须是 base GRN 中有 motif 信息的 TF）----
    gene_index = {gene: index for index, gene in enumerate(genes)}
    target = str(params.get("target") or "").strip()
    if target and target not in gene_index:
        raise C.PrerequisiteError(
            "params.target=%s 不在 bundle 的基因列表里（CellOracle 只能扰动表达矩阵里的基因；"
            "虚拟敲除官方实现要求目标是 base GRN 中有 motif 信息的调控因子/TF）" % target
        )

    # ---- 官方 CellOracle 流程 ---------------------------------------------
    t0 = time.time()
    C.log("CellOracle: PCA + kNN 平滑…")
    try:
        oracle.perform_PCA()
        k = params.get("k_imputation")
        if k:
            k = int(k)
        else:
            k = int(max(10, min(200, 0.025 * oracle.adata.n_obs)))
        oracle.knn_imputation(
            n_pca_dims=int(min(50, oracle.pcs.shape[1])), k=int(k),
            balanced=True, b_sight=int(min(k * 8, oracle.adata.n_obs - 1)),
            b_maxl=int(min(k * 4, oracle.adata.n_obs - 1)), n_jobs=n_jobs,
        )
    except Exception as error:
        raise C.PrerequisiteError("CellOracle PCA/kNN 平滑失败：%s" % error)
    notes.append("kNN 平滑：k=%d, n_pca_dims=%d" % (k, int(min(50, oracle.pcs.shape[1]))))

    C.log("CellOracle: 拟合 GRN（按簇，alpha=%s）…" % alpha)
    try:
        oracle.fit_GRN_for_simulation(GRN_unit="cluster", alpha=float(alpha), verbose_level=0)
    except Exception as error:
        raise C.PrerequisiteError("CellOracle GRN 拟合失败：%s" % error)

    # 目标基因：显式指定，或自动选 GRN 中最活跃的调控因子
    all_regulators = set(getattr(oracle, "all_regulatory_genes_in_TFdict", []))
    active_regulators = set(getattr(oracle, "active_regulatory_genes", []))
    if not target:
        candidates = [gene for gene in genes if gene in active_regulators]
        if not candidates:
            raise C.PrerequisiteError(
                "当前数据拟合出的 GRN 中没有可扰动的调控因子（TF）。CellOracle 只能对 base GRN 中"
                "有 motif 信息的 TF 做 in silico 敲除；请检查基因名是否为 human symbol，或换一份"
                "覆盖更多真实基因的数据。"
            )
        mean_expression = np.asarray(counts_dense.mean(axis=0)).ravel()
        candidates.sort(key=lambda gene: -mean_expression[gene_index[gene]])
        target = candidates[0]
        warnings.append("未指定 params.target，自动选 GRN 中最活跃的调控因子 %s" % target)
        notes.append("自动目标基因：%s（可用 params.target 指定）" % target)
    if target not in all_regulators:
        examples = [gene for gene in genes if gene in all_regulators][:8]
        raise C.PrerequisiteError(
            "params.target=%s 不是 CellOracle base GRN 中的调控因子（TF，或没有 motif 信息），"
            "官方 simulate_shift 拒绝执行。可扰动的候选基因示例：%s"
            % (target, "、".join(examples) or "（本数据无）")
        )
    if target not in active_regulators:
        examples = [gene for gene in genes if gene in active_regulators][:8]
        raise C.PrerequisiteError(
            "params.target=%s 在 base GRN 中是 TF，但在当前数据拟合出的 GRN 中调控连接不足，"
            "官方 simulate_shift 拒绝执行。可扰动的候选基因示例：%s"
            % (target, "、".join(examples) or "（本数据无）")
        )

    target_index = gene_index[target]
    target_mean = float(counts_dense[:, target_index].mean())
    notes.append("目标基因：%s（原始 counts 均值 %.2f）" % (target, target_mean))

    if mode == "knockout":
        perturb_value = 0.0
        condition_text = "敲除（表达置 0）"
    else:
        target_imputed = np.asarray(oracle.adata.layers["imputed_count"])[:, target_index]
        perturb_value = float(max(1.0, np.percentile(target_imputed, 99)))
        condition_text = "过表达（置为 imputed 表达 99 分位 %.3f）" % perturb_value
    notes.append("扰动条件：%s = %s（%s）" % (target, perturb_value, condition_text))

    C.log("CellOracle: simulate_shift（n_propagation=%d）…" % n_propagation)
    try:
        oracle.simulate_shift(
            perturb_condition={target: perturb_value},
            GRN_unit="cluster", n_propagation=n_propagation, clip_delta_X=True,
        )
    except Exception as error:
        raise C.PrerequisiteError("CellOracle simulate_shift 失败：%s" % error)
    if "delta_X" not in oracle.adata.layers:
        raise C.PrerequisiteError("CellOracle 没有生成 delta_X（模拟失败）")
    delta_x = np.asarray(oracle.adata.layers["delta_X"], dtype="float64")
    mean_delta = delta_x.mean(axis=0)
    mean_abs_delta = np.abs(delta_x).mean(axis=0)

    C.log("CellOracle: 计算嵌入位移（transition / embedding shift / grid arrows）…")
    try:
        oracle.estimate_transition_prob(
            n_neighbors=int(min(200, max(10, len(barcodes) // 10))),
            knn_random=True, sampled_fraction=0.3, n_jobs=n_jobs,
        )
        oracle.calculate_embedding_shift(sigma_corr=0.05)
        oracle.calculate_grid_arrows(
            smooth=0.8, steps=(40, 40),
            n_neighbors=int(min(200, max(10, len(barcodes) // 10))), n_jobs=n_jobs,
        )
    except Exception as error:
        raise C.PrerequisiteError("CellOracle 嵌入位移计算失败：%s" % error)
    flow_seconds = time.time() - t0
    delta_embedding = np.asarray(getattr(oracle, "delta_embedding", np.zeros((len(barcodes), 2))))
    if delta_embedding.shape != (len(barcodes), 2):
        delta_embedding = np.zeros((len(barcodes), 2))
    shift_magnitude = np.linalg.norm(delta_embedding, axis=1)
    notes.append("GRN 拟合 + 模拟 + 嵌入位移：%.1fs" % flow_seconds)

    # ---- 图 1：扰动箭头图（官方 grid flow：flow_grid + flow 向量场）--------
    plt = C.setup_matplotlib()
    figures = []
    flow = getattr(oracle, "flow", None)
    grid = getattr(oracle, "flow_grid", None)
    try:
        if flow is None or grid is None or len(grid) == 0:
            raise RuntimeError("CellOracle 没有输出 grid flow（flow_grid/flow）")
        flow = np.asarray(flow, dtype="float64")
        grid = np.asarray(grid, dtype="float64")
        fig, axis = plt.subplots(figsize=(6.6, 5.6), dpi=140)
        colors = np.asarray(getattr(oracle, "colorandum", np.array(["#888888"] * len(barcodes))))
        axis.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=4, alpha=0.5, linewidths=0)
        # 只在中位网格间距 2.5 倍内有细胞的区域画箭头，避免空白区伪箭头
        from sklearn.neighbors import NearestNeighbors

        grid_x = np.unique(grid[:, 0])
        grid_y = np.unique(grid[:, 1])
        spacing = float(min(np.mean(np.diff(grid_x)), np.mean(np.diff(grid_y))))
        neighbour = NearestNeighbors(n_neighbors=int(min(10, len(barcodes)))).fit(embedding)
        distance, _ = neighbour.kneighbors(grid)
        mask = distance[:, -1] <= 2.5 * spacing
        if mask.sum() < 10:
            mask = np.ones(len(grid), dtype=bool)
        magnitude = np.linalg.norm(flow, axis=1)
        reference = float(np.nanpercentile(magnitude[mask], 99)) if magnitude[mask].any() else 1.0
        scale_factor = 0.6 * spacing / max(reference, 1e-12)
        colors_arrows = magnitude[mask]
        quiver = axis.quiver(
            grid[mask, 0], grid[mask, 1],
            flow[mask, 0] * scale_factor, flow[mask, 1] * scale_factor,
            colors_arrows, angles="xy", scale_units="xy", scale=1.0,
            cmap="viridis", width=0.003, edgecolor="k", linewidth=0.1,
        )
        colorbar = fig.colorbar(quiver, ax=axis, shrink=0.75)
        colorbar.set_label("Simulated shift magnitude", fontsize=8)
        axis.set_title("CellOracle in silico %s: %s (%d grid arrows)"
                       % ("knockout" if mode == "knockout" else "overexpression", target, int(mask.sum())),
                       fontsize=10)
        axis.set_xlabel("UMAP1", fontsize=9)
        axis.set_ylabel("UMAP2", fontsize=9)
        figures.append(C.figure(out, fig, "celloracle_quiver", "perturbation_quiver",
                                "CellOracle in silico perturbation flow"))
    except Exception as error:
        warnings.append("CellOracle grid-flow 图失败（%s），改用细胞级位移箭头" % error)
        plt.close("all")
        fig, axis = plt.subplots(figsize=(6.6, 5.6), dpi=140)
        axis.scatter(embedding[:, 0], embedding[:, 1], c="lightgray", s=5, linewidths=0)
        magnitude = np.linalg.norm(delta_embedding, axis=1)
        reference = float(np.nanpercentile(magnitude, 99)) if len(magnitude) else 1.0
        scale_factor = 1.2 / max(reference, 1e-12)
        axis.quiver(embedding[:, 0], embedding[:, 1],
                    delta_embedding[:, 0] * scale_factor, delta_embedding[:, 1] * scale_factor,
                    magnitude, angles="xy", scale_units="xy", scale=1.0, cmap="viridis",
                    width=0.0025, edgecolor="k", linewidth=0.1)
        axis.set_title("CellOracle in silico %s: %s" % (mode, target), fontsize=10)
        axis.set_xlabel("UMAP1", fontsize=9)
        axis.set_ylabel("UMAP2", fontsize=9)
        figures.append(C.figure(out, fig, "celloracle_quiver", "perturbation_quiver",
                                "CellOracle in silico perturbation flow"))

    # ---- 图 2：受影响基因条形图 --------------------------------------------
    order = np.argsort(mean_abs_delta)[::-1]
    top_indices = [int(index) for index in order[:25]]
    fig, axis = plt.subplots(figsize=(6.0, max(4.0, 0.24 * len(top_indices) + 1.4)), dpi=140)
    selected = top_indices[::-1]
    values = [float(mean_delta[index]) for index in selected]
    colors = ["#C44E52" if value < 0 else "#4C72B0" for value in values]
    axis.barh([genes[index] for index in selected], values, color=colors, edgecolor="white", linewidth=0.3)
    axis.axvline(0, color="black", linewidth=0.6)
    axis.set_xlabel("Mean expression shift (simulated - wild type)", fontsize=9)
    axis.set_title("Top genes affected by %s %s" %
                   ("knockout" if mode == "knockout" else "overexpression", target), fontsize=10)
    axis.tick_params(axis="y", labelsize=7)
    figures.append(C.figure(out, fig, "celloracle_affected_genes", "affected_genes",
                            "Top affected genes (CellOracle)"))

    # ---- 表：受影响基因 ----------------------------------------------------
    affected_rows = []
    for index in order[:60]:
        index = int(index)
        affected_rows.append({
            "gene": genes[index],
            "mean_delta": round(float(mean_delta[index]), 5),
            "mean_abs_delta": round(float(mean_abs_delta[index]), 5),
            "direction": "down" if mean_delta[index] < 0 else "up",
        })
    # 簇级位移汇总
    cluster_rows = []
    labels = obs[cluster_col].astype(str).to_numpy()
    for cluster in sorted(set(labels)):
        mask = labels == cluster
        cluster_rows.append({
            "cluster": str(cluster),
            "n_cells": int(mask.sum()),
            "mean_shift_magnitude": round(float(shift_magnitude[mask].mean()), 5),
        })

    # ---- obs_extra ---------------------------------------------------------
    delta_target = delta_x[:, target_index].tolist()
    obs_columns = {
        "celloracle_shift_magnitude": [round(float(value), 6) for value in shift_magnitude],
        "celloracle_shift_x": [round(float(value), 6) for value in delta_embedding[:, 0]],
        "celloracle_shift_y": [round(float(value), 6) for value in delta_embedding[:, 1]],
        "celloracle_delta_%s" % target: [round(float(value), 6) for value in delta_target],
    }
    obs_file = C.write_obs_extra(out, barcodes, obs_columns)

    total_seconds = time.time() - started
    n_affected = int((mean_abs_delta > 0.05).sum())
    bullets = [
        "CellOracle %s 官方 in silico %s：目标基因 %s（%s）"
        % (C.version_of(co), "敲除" if mode == "knockout" else "过表达", target, condition_text),
        "GRN：%d 细胞 × %d 基因，base GRN 覆盖 %d 个基因 / %d 条潜在调控（%s）"
        % (len(barcodes), len(genes), len(getattr(oracle, "TFdict", {})), n_tf_pairs, base_source),
        "模拟传播 %d 步；|平均表达位移| > 0.05 的基因 %d 个；obs_extra 写入 %d 个细胞级指标"
        % (n_propagation, n_affected, len(obs_columns)),
        "总耗时 %.1f 秒（其中 GRN+模拟+嵌入位移 %.1f 秒）"
        % (total_seconds, flow_seconds),
    ]
    metrics = {
        "n_cells": int(len(barcodes)),
        "n_genes": int(len(genes)),
        "n_tfs_in_grn": int(len(getattr(oracle, "TFdict", {}))),
        "n_tf_gene_links": int(n_tf_pairs),
        "target": target,
        "mode": mode,
        "n_affected_genes": n_affected,
        "mean_shift_magnitude": round(float(shift_magnitude.mean()), 5),
    }
    return {
        "id": "ext_celloracle",
        "title": "CellOracle 官方虚拟敲除",
        "method": C.version_of(co),
        "summary": {"bullets": bullets, "metrics": metrics},
        "tables": [
            C.table(
                [C.col("gene", "Gene"), C.col("mean_delta", "Mean delta", "number"),
                 C.col("mean_abs_delta", "Mean |delta|", "number"), C.col("direction", "Direction")],
                affected_rows, "affected_genes", "受影响基因（按 |平均表达位移| 排序）",
                note="delta_X = simulated - wildtype（CellOracle 官方 simulate_shift 输出，log 表达空间）",
            ),
            C.table(
                [C.col("cluster", "Cluster"), C.col("n_cells", "Cells", "number"),
                 C.col("mean_shift_magnitude", "Mean shift magnitude", "number")],
                cluster_rows, "cluster_shift", "各簇平均嵌入位移强度",
            ),
        ],
        "figures": figures,
        "obs_file": obs_file,
        "warnings": warnings,
        "notes": notes,
    }


def main():
    return C.run_adapter(run)


if __name__ == "__main__":
    sys.exit(main())
