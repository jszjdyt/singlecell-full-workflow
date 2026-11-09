#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""官方 pySCENIC 适配器（CONTRACT §10.3）。

流程 = pySCENIC 官方三件套（与 ``pyscenic grn → ctx → aucell`` 相同的实现）：
  1. GRNBoost2（arboreto）推断 TF → 靶基因共表达网络；
  2. cisTarget motif 富集（ctxcore.rnkdb + pyscenic.prune.prune2df）得到 regulon；
  3. AUCell 计算每个细胞每个 regulon 的活性。

CLI：
  <externalPython> pyscenic.py --bundle DIR --out DIR --params JSON

params（均可选）：
  n_workers       并行进程数（默认 4）
  n_top_genes     参与分析的基因数上限（按方差选，默认 2000；bundle 更小则全用）
  db_dir          SCENIC motif 库目录 / 单个 .feather 文件（默认自动探测，含 §10.4 配置）
  dbs             显式 feather 文件列表
  annotations     motif annotation 表（*motifs-v10*.tbl）；默认在 db 目录/插件缓存中探测
  tf_list         TF 列表文本文件（每行一个符号）；默认读 lib/biopriors.py 的 TF_LIST
  n_pcs / seed    预留；seed 默认 42（确定性）
  rank_threshold / auc_threshold / nes_threshold   ctx 官方默认 5000 / 0.05 / 3.0
  max_cells       规模上限（默认 3000，超限确定性下采样）
  min_regulons    少于该数量视为数据不满足（默认 1）

输出：<out>/result.json + 图 + obs_extra.csv（regulon_<TF> 活性列）。
缺 spliced 等无关，本适配器只需要 counts/data 矩阵；缺 motif 库/annotation → 退出码 3。
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


# ---------------------------------------------------------------------------
# 参数与参考文件
# ---------------------------------------------------------------------------


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


def _load_tf_list(params, warnings, notes):
    """TF 列表：params.tf_list 文件 → lib/biopriors.py 的 TF_LIST → 全基因（warning）。"""
    explicit = params.get("tf_list")
    if explicit:
        path = Path(os.path.expanduser(str(explicit)))
        if not path.is_file():
            raise C.PrerequisiteError("params.tf_list 指定的 TF 列表文件不存在：%s" % path)
        names = [line.strip() for line in path.read_text(encoding="utf-8", errors="replace").splitlines()]
        names = [name for name in names if name and not name.startswith("#")]
        notes.append("TF 列表来自文件：%s（%d 个）" % (path, len(names)))
        return names
    lib_dir = Path(__file__).resolve().parents[1]
    try:
        if str(lib_dir) not in sys.path:
            sys.path.insert(0, str(lib_dir))
        from biopriors import TF_LIST  # type: ignore

        notes.append("TF 列表来自 lib/biopriors.py TF_LIST（%d 个）" % len(TF_LIST))
        return list(TF_LIST)
    except Exception as error:
        warnings.append("未能读取 lib/biopriors.py 的 TF_LIST（%s），改为把全部基因当作候选 TF" % error)
        return None


# ---------------------------------------------------------------------------
# GRN 推断：照抄 pySCENIC 自带的官方 multiprocessing 流程
# （pyscenic/cli/arboreto_with_multiprocessing.py，绕开 arboreto 0.1.6 与
#  dask 2025.x 的不兼容：新版 dask 不允许对空 delayed 列表调用 from_delayed）
# ---------------------------------------------------------------------------


def _grn_worker(target_gene_index, gene_names, ex_matrix, tf_matrix,
                tf_matrix_gene_names, seed):
    from arboreto.core import EARLY_STOP_WINDOW_LENGTH, SGBM_KWARGS, infer_partial_network

    return infer_partial_network(
        regressor_type="GBM",
        regressor_kwargs=SGBM_KWARGS,
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=gene_names[target_gene_index],
        target_gene_expression=ex_matrix[:, target_gene_index],
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=seed,
    )


def _run_grn_multiprocessing(expression, tf_names, n_workers, seed):
    """返回 TF-target 网络 DataFrame（列：TF/target/importance，与 pyscenic grn 输出一致）。"""
    from functools import partial
    from multiprocessing import Pool

    import pandas as pd
    from arboreto.algo import _prepare_input
    from arboreto.core import target_gene_indices, to_tf_matrix

    ex_matrix, gene_names, tf_used = _prepare_input(
        expression, list(expression.columns), list(tf_names)
    )
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(ex_matrix, gene_names, tf_used)
    worker = partial(
        _grn_worker, gene_names=gene_names, ex_matrix=ex_matrix,
        tf_matrix=tf_matrix, tf_matrix_gene_names=tf_matrix_gene_names, seed=seed,
    )
    indices = list(target_gene_indices(gene_names, target_genes="all"))
    frames = []
    if int(n_workers) > 1:
        try:
            with Pool(int(n_workers)) as pool:
                for chunk in pool.imap(worker, indices, chunksize=1):
                    if chunk is not None and len(chunk) > 0:
                        frames.append(chunk)
        except Exception as error:
            C.log("GRN 多进程失败（%s），回退串行" % error)
            frames = []
    if not frames:
        for index in indices:
            chunk = worker(index)
            if chunk is not None and len(chunk) > 0:
                frames.append(chunk)
    if not frames:
        raise C.PrerequisiteError("GRNBoost2 没有推断出任何 TF-靶基因边（数据可能过于稀疏或方差不足）")
    network = pd.concat(frames, axis=0, ignore_index=True)
    return network


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------


def _install_numpy_alias_compat():
    """pySCENIC 0.12.1 内部仍使用 numpy<1.24 的别名（np.object / np.float 等）。

    numpy 2.x 已移除这些别名（会抛 AttributeError）。这里只在适配器进程内把别名补回
    numpy 模块（不改 site-packages、不影响其它工具），使官方 pySCENIC 0.12.1 能在
    numpy 2.x 上运行。
    """
    import numpy as np

    aliases = {
        "object": object, "float": float, "int": int, "bool": bool,
        "str": str, "unicode": str, "complex": complex,
    }
    for name, builtin_type in aliases.items():
        try:
            if not hasattr(np, name):
                setattr(np, name, builtin_type)
        except Exception:
            pass


def run(bundle_dir, out, params):
    started = time.time()
    out = Path(out)
    warnings = []
    notes = []

    # ---- 依赖（缺依赖 → 退出码 3）----------------------------------------
    _install_numpy_alias_compat()
    pyscenic = C.require_module("pyscenic", "pyscenic")
    C.require_module("ctxcore", "ctxcore")
    C.require_module("arboreto", "arboreto")
    import numpy as np
    import pandas as pd

    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.aucell import aucell
    from pyscenic.prune import prune2df
    from pyscenic.transform import df2regulons
    from pyscenic.utils import modules_from_adjacencies

    seed = _int_param(params, "seed", 42)
    n_workers = max(1, _int_param(params, "n_workers", 4))
    max_cells = _int_param(params, "max_cells", 3000)
    n_top_genes = _int_param(params, "n_top_genes", 2000)
    rank_threshold = _int_param(params, "rank_threshold", 5000)
    auc_threshold = _float_param(params, "auc_threshold", 0.05)
    nes_threshold = _float_param(params, "nes_threshold", 3.0)
    min_regulons = _int_param(params, "min_regulons", 1)

    notes.append("pySCENIC %s（GRNBoost2 + cisTarget ctx + AUCell，官方 Python 实现）" %
                 C.version_of(pyscenic))
    notes.append("参数：n_workers=%d, n_top_genes=%d, rank_threshold=%d, auc_threshold=%.2f, "
                 "nes_threshold=%.1f, seed=%d" % (n_workers, n_top_genes, rank_threshold,
                                                  auc_threshold, nes_threshold, seed))

    # ---- motif 库 / annotation（缺 → 退出码 3）---------------------------
    db_files, annotation_file, db_dir = C.find_scenic_db(params)
    notes.append("motif 库目录：%s（%d 个 .feather）" % (db_dir, len(db_files)))
    notes.append("motif annotation：%s" % annotation_file)
    for db in db_files:
        notes.append("  - %s" % Path(db).name)

    # ---- TF 列表 ----------------------------------------------------------
    tf_candidates = _load_tf_list(params, warnings, notes)

    # ---- 读 bundle --------------------------------------------------------
    bundle = C.Bundle(bundle_dir)
    if bundle.n_vars < 50:
        raise C.PrerequisiteError("bundle 只有 %d 个基因，无法做 pySCENIC（至少需要 50 个）" % bundle.n_vars)
    cell_idx = C.cap_cells(bundle, max_cells, seed=seed, warnings=warnings, label="pySCENIC")
    counts = bundle.matrix("counts")
    data = bundle.matrix("data")
    if cell_idx is not None:
        counts = counts[cell_idx]
        data = data[cell_idx]
    genes = [str(item) for item in bundle.features]
    barcodes = [str(item) for item in bundle.barcodes]
    if cell_idx is not None:
        barcodes = [barcodes[index] for index in cell_idx]
    obs = bundle.obs.iloc[cell_idx] if cell_idx is not None else bundle.obs
    if data is not None and bundle.x_kind != "lognorm":
        # 没有 log 标准化时，按官方推荐做 log1p(CPM/1e4) 以便相关系数稳健
        warnings.append("bundle 的 X 不是 log-normalized，适配器内部按 log1p(CP10K) 处理")
        row_sum = np.asarray(counts.sum(axis=1)).ravel()
        row_sum[row_sum <= 0] = 1.0
        scaled = counts.multiply((1e4 / row_sum)[:, None]).tocsr()
        scaled.data = np.log1p(scaled.data)
        data = scaled.astype("float32")

    # ---- 选基因（保留 TF）-------------------------------------------------
    tf_names = None
    if tf_candidates:
        tf_set = set(str(item) for item in tf_candidates)
        tf_present = [gene for gene in genes if gene in tf_set]
        if len(tf_present) < 10:
            warnings.append("表达矩阵只覆盖 TF 列表的 %d 个基因，改为把全部基因当作候选 TF" % len(tf_present))
            tf_names = None
        else:
            tf_names = tf_present
            notes.append("TF 列表命中表达矩阵：%d 个（列表共 %d 个）" % (len(tf_names), len(tf_set)))
    gene_idx = C.select_genes(data, genes, n_top_genes, keep_genes=tf_candidates, seed=seed, warnings=warnings)
    if gene_idx is not None:
        counts = counts[:, gene_idx]
        data = data[:, gene_idx]
        genes = [genes[index] for index in gene_idx]
        if tf_names is not None:
            gene_set = set(genes)
            tf_names = [gene for gene in tf_names if gene in gene_set]
    if len(genes) < 50:
        raise C.PrerequisiteError("可用于 pySCENIC 的基因只有 %d 个（至少 50 个）" % len(genes))
    if not tf_names:
        tf_names = list(genes)  # 官方允许 'all'，这里显式传全部基因
        warnings.append("未提供有效 TF 列表，把全部 %d 个基因当作候选 TF（运行会变慢）" % len(genes))
    notes.append("分析规模：%d 细胞 × %d 基因 × %d 候选 TF" % (len(barcodes), len(genes), len(tf_names)))

    expression = pd.DataFrame(
        np.asarray(counts.toarray(), dtype="float32") if hasattr(counts, "toarray") else np.asarray(counts),
        index=pd.Index(barcodes, name="cell"),
        columns=pd.Index(genes, name="gene"),
    )
    log_expression = pd.DataFrame(
        np.asarray(data.toarray(), dtype="float32") if hasattr(data, "toarray") else np.asarray(data),
        index=pd.Index(barcodes, name="cell"),
        columns=pd.Index(genes, name="gene"),
    )

    # ---- 1) GRNBoost2（官方 multiprocessing 流程）--------------------------
    t0 = time.time()
    C.log("pySCENIC: GRNBoost2（%d cells × %d genes, %d TFs, %d workers）…" %
          (expression.shape[0], expression.shape[1], len(tf_names), n_workers))
    network = _run_grn_multiprocessing(expression, tf_names, n_workers, seed)
    network = network.dropna().reset_index(drop=True)
    if network.empty:
        raise C.PrerequisiteError("GRNBoost2 没有推断出任何 TF-靶基因边，数据可能过于稀疏或方差不足")
    grn_seconds = time.time() - t0
    notes.append("GRNBoost2：%d 条 TF→靶基因边，耗时 %.1fs" % (len(network), grn_seconds))

    # ---- 2) cisTarget motif 富集（ctx）------------------------------------
    t0 = time.time()
    C.log("pySCENIC: 生成共表达模块并做 cisTarget motif 富集…")
    dbs = []
    try:
        for db_file in db_files:
            dbs.append(RankingDatabase(fname=db_file, name=os.path.splitext(os.path.basename(db_file))[0]))
    except Exception as error:
        raise C.PrerequisiteError("无法打开 SCENIC motif 数据库（%s）：%s" % (db_files, error))

    # motif 库只对"有 promoter ranking 的基因"有统计功效。合成数据/靶向面板里常有大量
    # 无注释的假基因（如 demo 的 BG####），若直接建模块会让 cisTarget 富集完全没有功效。
    # 这里统计 motif 库覆盖的基因集合，必要时把模块构建限制到已注释基因（真实全转录组
    # 数据几乎全部有注释，此步骤不改变行为）。
    db_gene_sets = []
    for database in dbs:
        try:
            genes_in_db = database.genes
            db_gene_sets.append(set(str(item) for item in genes_in_db))
        except Exception as error:
            warnings.append("读取 motif 库基因名失败（%s），跳过模块基因限制" % error)
    if db_gene_sets:
        annotated_genes = set()
        for gene_set in db_gene_sets:
            annotated_genes |= gene_set
        kept_genes = [gene for gene in genes if gene in annotated_genes]
        notes.append("motif 库覆盖表达矩阵基因：%d / %d" % (len(kept_genes), len(genes)))
        if len(kept_genes) < 100:
            raise C.PrerequisiteError(
                "表达矩阵只有 %d 个基因在 SCENIC motif 库中有注释（<100），cisTarget 富集没有"
                "统计功效；请确认基因名是 human symbol，或换用覆盖度更高的数据。"
                % len(kept_genes)
            )
        if len(kept_genes) < 0.9 * len(genes):
            warnings.append(
                "表达矩阵 %d 个基因中只有 %d 个在 SCENIC motif 库有注释（可能是合成/面板数据）；"
                "为让 cisTarget 富集有统计功效，共表达模块限制到这些已注释基因（全转录组真实数据不受影响）。"
                % (len(genes), len(kept_genes))
            )
            network_for_modules = network[network["target"].isin(set(kept_genes))]
            expression_for_modules = log_expression[kept_genes]
            notes.append("模块构建：网络边限制到已注释基因（%d → %d 条），表达矩阵取 %d 个已注释基因"
                         % (len(network), len(network_for_modules), len(kept_genes)))
        else:
            network_for_modules = network
            expression_for_modules = log_expression
    else:
        network_for_modules = network
        expression_for_modules = log_expression
    modules = list(modules_from_adjacencies(
        network_for_modules, expression_for_modules,
        thresholds=(0.75, 0.90),
        top_n_targets=(50,),
        top_n_regulators=(5, 10, 50),
        min_genes=20,
        keep_only_activating=True,
    ))
    if not modules:
        raise C.PrerequisiteError(
            "GRNBoost2 得到 %d 条边，但没能生成任何共表达模块（<20 基因）；"
            "可能是基因数太少或细胞类型结构太弱，pySCENIC 无法继续。" % len(network_for_modules)
        )
    notes.append("共表达模块：%d 个（官方默认阈值 0.75/0.90, top50 targets）" % len(modules))
    try:
        motif_df = prune2df(
            dbs, modules, str(annotation_file),
            rank_threshold=rank_threshold,
            auc_threshold=auc_threshold,
            nes_threshold=nes_threshold,
            motif_similarity_fdr=0.001,
            orthologuous_identity_threshold=0.0,
            client_or_address="custom_multiprocessing",
            num_workers=n_workers,
            module_chunksize=100,
        )
    except Exception as error:
        raise C.PrerequisiteError("cisTarget motif 富集失败：%s" % error)
    ctx_seconds = time.time() - t0
    if motif_df is None or len(motif_df) == 0:
        raise C.PrerequisiteError(
            "cisTarget 未找到任何显著富集 motif（%d 个模块）；pySCENIC 无法生成 regulon。"
            "可尝试换用 10kbp 数据库或放宽 nes_threshold。" % len(modules)
        )
    regulons = list(df2regulons(motif_df))
    notes.append("cisTarget ctx：%d 个显著 motif 富集，得到 %d 个 regulon，耗时 %.1fs" %
                 (len(motif_df), len(regulons), ctx_seconds))
    if len(regulons) < min_regulons:
        raise C.PrerequisiteError("只得到 %d 个 regulon（< %d），数据不满足 pySCENIC 的最低要求"
                                  % (len(regulons), min_regulons))

    # ---- 3) AUCell ---------------------------------------------------------
    t0 = time.time()
    C.log("pySCENIC: AUCell 计算 regulon 活性…")
    regulon_signatures = [reg for reg in regulons if len(reg.genes) > 0]
    if not regulon_signatures:
        raise C.PrerequisiteError("所有 regulon 的靶基因都是空的，无法做 AUCell")
    auc_mtx = aucell(
        log_expression, regulon_signatures,
        auc_threshold=auc_threshold,
        noweights=False,
        seed=seed,
        num_workers=max(1, min(n_workers, 4)),
    )
    auc_seconds = time.time() - t0
    notes.append("AUCell：%d 细胞 × %d regulon 活性，耗时 %.1fs" %
                 (auc_mtx.shape[0], auc_mtx.shape[1], auc_seconds))

    # ---- 结果整理 ----------------------------------------------------------
    def _regulon_tf(name):
        return str(name).replace("(+)", "").replace("(-)", "")

    regulon_mean = auc_mtx.mean(axis=0).sort_values(ascending=False)
    rows = []
    name_by_regulon = {}
    for regulon in regulon_signatures:
        name = str(regulon.name)
        name_by_regulon[name] = regulon
        tf_name = _regulon_tf(name)
        subset = motif_df[motif_df["TF"] == tf_name] if "TF" in motif_df.columns else motif_df.iloc[0:0]
        if len(subset) > 0:
            best_index = subset["NES"].idxmax() if "NES" in subset.columns else subset.index[0]
            record = subset.loc[best_index]
            motif = str(record.get("motif", ""))
            nes = float(record.get("NES", float("nan")))
            auc = float(record.get("AUC", float("nan")))
        else:
            motif, nes, auc = "", float("nan"), float("nan")
        rows.append({
            "regulon": name,
            "tf": tf_name,
            "n_targets": int(len(regulon.genes)),
            "best_motif": motif,
            "nes": round(nes, 3) if np.isfinite(nes) else None,
            "auc": round(auc, 4) if np.isfinite(auc) else None,
            "mean_activity": round(float(regulon_mean.get(name, 0.0)), 5),
        })
    rows.sort(key=lambda item: (-item["n_targets"], item["regulon"]))
    regulon_rows = rows[:100]

    edge_rows = []
    for regulon in regulon_signatures:
        tf_name = _regulon_tf(regulon.name)
        genes_of = [str(gene) for gene in regulon.genes]
        weights = []
        gene2weight = getattr(regulon, "gene2weight", None)
        if isinstance(gene2weight, dict) and gene2weight:
            weights = [float(gene2weight.get(gene, 0.0)) for gene in genes_of]
        elif hasattr(regulon, "weights") and regulon.weights is not None:
            try:
                weights = [float(weight) for weight in regulon.weights]
                if len(weights) != len(genes_of):
                    weights = [0.0] * len(genes_of)
            except Exception:
                weights = [0.0] * len(genes_of)
        else:
            weights = [0.0] * len(genes_of)
        order = np.argsort(np.abs(weights))[::-1]
        for index in order[:5]:
            edge_rows.append({"tf": tf_name, "target": genes_of[index],
                              "weight": round(float(weights[index]), 5)})
    edge_rows.sort(key=lambda item: -abs(item["weight"]))
    edge_rows = edge_rows[:150]

    # 图：每簇 regulon 平均活性热图（有 clusters 列时），否则 top regulon 细胞级热图
    plt = C.setup_matplotlib()
    cluster_col = bundle.obs_columns_like(["clusters", "cluster", "cell_type", str(params.get("cluster_col") or "")])
    cluster_labels = obs[cluster_col].astype(str).tolist() if cluster_col else None
    top_names = list(regulon_mean.index[: min(25, len(regulon_mean))])
    fig, axis = plt.subplots(figsize=(max(7.0, 0.32 * len(top_names) + 2.5), 4.6), dpi=140)
    if cluster_labels:
        clusters = sorted(set(cluster_labels))
        matrix = np.zeros((len(clusters), len(top_names)), dtype="float64")
        for row_index, cluster in enumerate(clusters):
            mask = np.asarray([label == cluster for label in cluster_labels])
            if mask.any():
                matrix[row_index] = auc_mtx.loc[mask, top_names].to_numpy(dtype="float64").mean(axis=0)
        std = matrix.std(axis=0)
        std[std == 0] = 1.0
        zscore = (matrix - matrix.mean(axis=0)) / std
        image = axis.imshow(zscore, aspect="auto", cmap="RdBu_r", vmin=-2.0, vmax=2.0)
        axis.set_yticks(range(len(clusters)))
        axis.set_yticklabels(clusters, fontsize=8)
        axis.set_title("pySCENIC regulon activity (mean AUC per cluster, z-score)", fontsize=10)
    else:
        matrix = auc_mtx[top_names].to_numpy(dtype="float64").T
        std = matrix.reshape(matrix.shape[0], -1).std(axis=1)
        std[std == 0] = 1.0
        zscore = (matrix - matrix.mean(axis=1, keepdims=True)) / std[:, None]
        image = axis.imshow(zscore, aspect="auto", cmap="RdBu_r", vmin=-2.0, vmax=2.0)
        axis.set_yticks(range(len(top_names)))
        axis.set_yticklabels(top_names, fontsize=7)
        axis.set_title("pySCENIC regulon activity (AUC per cell)", fontsize=10)
    axis.set_xticks(range(len(top_names)))
    axis.set_xticklabels(top_names, rotation=90, fontsize=7)
    axis.set_xlabel("Regulon", fontsize=9)
    axis.set_ylabel("Cluster" if cluster_labels else "Regulon", fontsize=9)
    colorbar = fig.colorbar(image, ax=axis, shrink=0.75)
    colorbar.set_label("z-score", fontsize=8)
    figures = [C.figure(out, fig, "pyscenic_regulon_heatmap", "regulon_heatmap",
                        "pySCENIC regulon activity heatmap")]

    # obs_extra：按方差取 top regulon 活性列
    activity_columns = {}
    variance = auc_mtx.var(axis=0).sort_values(ascending=False)
    for name in list(variance.index[:60]):
        sanitized = "regulon_" + _regulon_tf(name).replace("-", "_").replace(" ", "_")
        if sanitized in activity_columns:
            continue
        activity_columns[sanitized] = [round(float(value), 6) for value in auc_mtx[name].tolist()]
    obs_file = C.write_obs_extra(out, barcodes, activity_columns) if activity_columns else None

    total_seconds = time.time() - started
    bullets = [
        "pySCENIC %s 官方流程：GRNBoost2 推断 %d 条 TF→靶基因边（%d 个候选 TF）"
        % (C.version_of(pyscenic), len(network), len(tf_names)),
        "cisTarget motif 富集（%d 个 motif 库）得到 %d 个 regulon，平均靶基因 %d 个"
        % (len(db_files), len(regulons), int(round(np.mean([len(reg.genes) for reg in regulon_signatures]))) if regulon_signatures else 0),
        "AUCell 输出 %d 细胞 × %d regulon 活性矩阵；obs_extra 写入 %d 个 regulon 活性列"
        % (auc_mtx.shape[0], auc_mtx.shape[1], len(activity_columns)),
        "总耗时 %.1f 秒（GRN %.1fs + ctx %.1fs + AUCell %.1fs）"
        % (total_seconds, grn_seconds, ctx_seconds, auc_seconds),
    ]
    metrics = {
        "n_cells": int(auc_mtx.shape[0]),
        "n_genes": int(len(genes)),
        "n_tfs": int(len(tf_names)),
        "n_edges": int(len(network)),
        "n_modules": int(len(modules)),
        "n_regulons": int(len(regulons)),
        "n_motif_hits": int(len(motif_df)),
        "n_obs_columns": int(len(activity_columns)),
    }
    return {
        "id": "ext_pyscenic",
        "title": "pySCENIC 官方转录调控",
        "method": C.version_of(pyscenic),
        "summary": {"bullets": bullets, "metrics": metrics},
        "tables": [
            C.table(
                [C.col("regulon", "Regulon"), C.col("tf", "TF"), C.col("n_targets", "Targets", "number"),
                 C.col("best_motif", "Best motif"), C.col("nes", "NES", "number"),
                 C.col("auc", "AUC", "number"), C.col("mean_activity", "Mean AUC", "number")],
                regulon_rows, "regulons", "Regulons（按靶基因数排序）",
                note="regulon 的靶基因由 cisTarget motif 富集剪枝得到；NES/AUC 取该 TF 最优 motif 的值",
            ),
            C.table(
                [C.col("tf", "TF"), C.col("target", "Target"), C.col("weight", "Weight", "number")],
                edge_rows, "edges", "Top TF→target 边（含 regulon 权重）",
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
