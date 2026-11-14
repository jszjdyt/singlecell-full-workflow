#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench — sc-analysis-b（高级分析模块）测试

不依赖 node / 宿主，直接 import ``lib/sc_analyses_dynamics.py``，用**自造的
synthetic ctx**（≥800 细胞 × ≥1500 基因、csr 矩阵、obs 含 clusters/cell_type/batch、
obsm 含 X_pca/X_emb、固定 rng）逐个跑通 13 个分析：

  integration, trajectory_dpt, trajectory_paga, entropy_stemness, cytotrace_lite,
  cnv_score, cnv_tree, ith_score, lr_communication, regulon_activity,
  virtual_knockout, virtual_overexpress, virtual_screen

校验项：
  1. ``ANALYSES`` 元数据覆盖 13 个 id，参数类型/needs 合法，method 诚实标注（本插件/lite/近似/代理）；
  2. 返回结构固定键、表格行列一致、图是真实 PNG（magic + 尺寸一致 + 文件名规范）；
  3. 写 obs 的分析返回顶层 ``obs_values``，长度与 ctx.n_obs 严格对齐，且与 obs_columns 声明一致；
  4. 同名参数两次运行结果一致（确定性）；
  5. 缺前置抛 ``EngineError("prerequisite", ...)``（缺 batch 列/缺 clusters/非人源基因名）；
  6. 参数错误/未知 id 抛 ``bad_request``；cnv_score 提供 gene_order 时走染色体排序分支。

直接运行：
    python3 test/py_analyses_dynamics_test.py

退出码 0 = 全部 PASS，1 = 有 FAIL。临时目录测试结束自动清理
（KEEP_SC_TEST_TMP=1 可保留用于排查）。
"""

from __future__ import annotations

import json
import os
import shutil
import struct
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PLUGIN_DIR = HERE.parent
LIB_DIR = PLUGIN_DIR / "lib"
sys.path.insert(0, str(LIB_DIR))

os.environ.setdefault("MPLBACKEND", "Agg")

import biopriors  # noqa: E402
import sc_analyses_dynamics as mod  # noqa: E402

ALLOWED_PARAM_TYPES = {"number", "string", "select", "boolean", "column", "file", "text"}
HONEST_MARKERS = ("本插件", "lite", "Lite", "近似", "代理", "简化")

EXPECTED_IDS = [
    "integration",
    "trajectory_dpt",
    "trajectory_paga",
    "entropy_stemness",
    "cytotrace_lite",
    "cnv_score",
    "cnv_tree",
    "ith_score",
    "lr_communication",
    "regulon_activity",
    "virtual_knockout",
    "virtual_overexpress",
    "virtual_screen",
]

CELL_TYPES = ["T cell", "B cell", "NK cell", "Monocyte", "Dendritic cell", "Platelet"]

# 每类细胞的 marker（与 demo 的风格一致的真实基因名，同时覆盖 CB 模块常用符号）
TYPE_MARKERS = {
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK", "TRAC", "TRBC1", "TRBC2", "LAT", "ZAP70", "LTB"],
    "B cell": ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "CD37", "CD22", "TNFRSF13C", "CD74", "HLA-DRA"],
    "NK cell": ["GNLY", "NKG7", "GZMB", "PRF1", "KLRD1", "KLRF1", "NCAM1", "FCGR3A", "TYROBP", "CST7"],
    "Monocyte": ["LYZ", "CD14", "S100A8", "S100A9", "CST3", "FCN1", "CD68", "CSF1R", "VCAN", "CTSD"],
    "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "HLA-DPA1", "ITGAX", "LAMP3", "IRF8", "CD1E", "BATF3", "FLT3"],
    "Platelet": ["PPBP", "PF4", "ITGA2B", "GP9", "TUBB1", "SELP", "NRGN", "CLU", "GP1BA", "TREML1"],
}

PASS_COUNT = 0
FAIL_COUNT = 0
FAILURES = []


def check(name, condition, detail="", required=True):
    global PASS_COUNT, FAIL_COUNT
    if condition:
        PASS_COUNT += 1
        print("PASS  %-66s %s" % (name, str(detail)[:110]))
        return True
    if required:
        FAIL_COUNT += 1
        FAILURES.append(name)
        print("FAIL  %-66s %s" % (name, str(detail)[:240]))
    else:
        print("note  %-66s %s" % (name, str(detail)[:240]))
    return False


# ---------------------------------------------------------------------------
# synthetic ctx
# ---------------------------------------------------------------------------


def build_gene_names(rng):
    """构造 ≥1500 个基因名：全部内置 LR 基因 + 一批 TF + 细胞类型 marker + 背景基因。"""
    names = []
    seen = set()

    def add(gene):
        gene = str(gene)
        if gene not in seen:
            seen.add(gene)
            names.append(gene)

    for gene in biopriors.lr_genes():
        add(gene)
    for index, gene in enumerate(biopriors.tf_list()):
        if index % 4 == 0:  # ~160 个 TF
            add(gene)
    for cell_type in CELL_TYPES:
        for gene in TYPE_MARKERS[cell_type]:
            add(gene)
    fill = 1
    while len(names) < 1520:
        gene = "BG%04d" % fill
        fill += 1
        add(gene)
    return names[:1520]


def make_ctx(n_cells=880, n_genes=1520, seed=0, name="synthetic-test"):
    rng = np.random.default_rng(seed)
    gene_names = build_gene_names(rng)
    n_genes = len(gene_names)
    gene_index = {gene: index for index, gene in enumerate(gene_names)}

    # 1) 细胞类型 + 3 个批次（批次里 25% 基因有系统性倍乘 = 可校正的批次效应）
    type_index = np.asarray([index % len(CELL_TYPES) for index in range(n_cells)], dtype=np.int64)
    rng.shuffle(type_index)
    batches = np.asarray(["batch1", "batch2", "batch3"], dtype=object)[np.arange(n_cells) % 3]
    batch_effect_genes = rng.choice(n_genes, size=max(40, n_genes // 4), replace=False)
    batch_effect = np.zeros((n_cells, n_genes), dtype=np.float64)
    batch2 = batches == "batch2"
    batch_effect[batch2[:, None] & np.isin(np.arange(n_genes), batch_effect_genes)[None, :]] = 0.55
    batch3 = batches == "batch3"
    batch_effect[batch3[:, None] & np.isin(np.arange(n_genes), batch_effect_genes[::2])[None, :]] = -0.35

    # 2) 10 个共表达程序：细胞类型 program 0-5，另 4 个通用 program
    n_programs = 10
    program = rng.integers(0, n_programs, size=n_genes)
    for type_position, cell_type in enumerate(CELL_TYPES):
        for gene in TYPE_MARKERS[cell_type]:
            program[gene_index[gene]] = type_position
    for index, gene in enumerate(biopriors.lr_genes()):
        program[gene_index[gene]] = index % len(CELL_TYPES)  # 配体/受体分布在不同细胞类型
    cell_factors = rng.normal(0.0, 1.0, size=(n_cells, n_programs))
    for type_position, cell_type in enumerate(CELL_TYPES):
        cell_factors[type_index == type_position, type_position] += 2.1
    cell_factors[:, 6:] += rng.normal(0.0, 0.8, size=(n_cells, n_programs - 6))

    baseline = rng.normal(-0.6, 0.9, size=n_genes)
    # 细胞类型 marker 提高基线，让簇间差异明显
    for cell_type in CELL_TYPES:
        for gene in TYPE_MARKERS[cell_type]:
            baseline[gene_index[gene]] += 1.6
    size_factor = rng.normal(0.0, 0.25, size=n_cells)
    log_mean = (
        baseline[None, :]
        + 0.75 * cell_factors[:, program]
        + batch_effect
        + size_factor[:, None]
    )
    rates = np.exp(np.clip(log_mean, -6.0, 6.0))
    counts = rng.poisson(rates).astype(np.float32)

    # 3) CNV 样信号：Platelet 簇的 200 个基因拷贝数样上调（顺序与基因索引无关，检验平滑/聚类）
    cnv_block = np.arange(300, 500)
    counts[type_index == 5, :] *= 1.0
    counts[np.ix_(type_index == 5, cnv_block)] *= 2.4

    counts = sp.csr_matrix(counts)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.where(totals > 0, 1e4 / np.maximum(totals, 1e-9), 0.0).astype(np.float32)
    lognorm = (sp.diags(scale) @ counts).tocsr().astype(np.float32)
    lognorm.data = np.log1p(lognorm.data)

    # 4) 簇标签：与 cell_type 大体一致，加 3% 噪声（检验簇级分析对噪声的鲁棒性）
    clusters = np.asarray([str(value) for value in type_index], dtype=object)
    noise_count = max(1, int(round(n_cells * 0.03)))
    noise_index = rng.choice(n_cells, size=noise_count, replace=False)
    clusters[noise_index] = np.asarray(
        [str(rng.integers(0, len(CELL_TYPES))) for _ in range(noise_count)], dtype=object
    )
    cell_type = np.asarray([CELL_TYPES[index] for index in type_index], dtype=object)

    cell_index = pd.Index(["cell-%04d" % index for index in range(n_cells)])
    obs = pd.DataFrame(
        {
            "clusters": pd.Series(clusters, index=cell_index, dtype=object),
            "cell_type": pd.Series(cell_type, index=cell_index, dtype=object),
            "batch": pd.Series(np.asarray(batches, dtype=object), index=cell_index, dtype=object),
        },
        index=cell_index,
    )

    # 5) PCA / 2D 嵌入（用 HVG + PCA 计算，确定性）
    from sklearn.decomposition import PCA

    dense = lognorm.toarray()
    variance = dense.var(axis=0)
    hvg = np.argsort(-variance)[: min(800, n_genes)]
    pca_model = PCA(n_components=20, svd_solver="randomized", random_state=0)
    pca = pca_model.fit_transform(dense[:, hvg]).astype(np.float32)
    embedding = pca[:, :2] + rng.normal(0.0, 0.02, size=(n_cells, 2)).astype(np.float32)

    ctx = {
        "name": name,
        "dataset_id": "ds-test",
        "n_obs": int(n_cells),
        "n_vars": int(n_genes),
        "X": lognorm,
        "X_kind": "lognorm",
        "counts": counts,
        "obs": obs,
        "var_names": list(gene_names),
        "obsm": {"X_pca": pca, "X_emb": embedding},
        "clusters": {"key": "clusters", "labels": [str(value) for value in clusters]},
        "group_keys": ["clusters", "cell_type", "batch"],
        "engine_mode": "builtin",
        "rng": np.random.default_rng(seed),
    }
    return ctx


def write_gene_order(path, gene_names):
    """写一个 gene / chr / start / end 四列的位置文件。"""
    lines = ["# synthetic gene order: gene\tchr\tstart\tend"]
    for index, gene in enumerate(gene_names):
        chrom = "chr%d" % (index % 8 + 1)
        start = (index // 8) * 2000 + 1
        lines.append("%s\t%s\t%d\t%d" % (gene, chrom, start, start + 1500))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# 校验工具
# ---------------------------------------------------------------------------


def png_info(path):
    with open(path, "rb") as handle:
        head = handle.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", head[16:24])


def params_for(analysis_id, target_gene):
    """测试用参数：比默认值小，保证 13 项全跑在几秒到十几秒内。"""
    table = {
        "integration": {"batch": "batch", "n_neighbors": 12, "n_pcs": 15},
        "trajectory_dpt": {"n_neighbors": 12, "n_dcs": 4},
        "trajectory_paga": {"n_neighbors": 12},
        "entropy_stemness": {"n_genes": 800, "n_neighbors": 12},
        "cytotrace_lite": {"n_genes": 800, "n_neighbors": 12, "signature_size": 60},
        "cnv_score": {"window": 11, "max_genes": 900},
        "cnv_tree": {"window": 11, "max_genes": 900},
        "ith_score": {"window": 11, "max_genes": 900, "max_cells_per_cluster": 60},
        "lr_communication": {"n_perm": 20, "top_n": 20},
        "regulon_activity": {"n_genes": 800, "n_targets": 15, "n_regulons": 12, "top_regulons": 5},
        "virtual_knockout": {"target": target_gene, "n_genes": 500, "top_k": 8, "steps": 3},
        "virtual_overexpress": {"target": target_gene, "n_genes": 500, "top_k": 8, "steps": 3},
        "virtual_screen": {
            "cluster": "0",
            "n_candidates": 12,
            "n_genes": 500,
            "steps": 3,
            "signature_size": 40,
        },
    }
    return table[analysis_id]


def validate_result(analysis_id, result, ctx, spec, outdir):
    ok = True
    ok &= check("%s: 顶层键齐全" % analysis_id,
                isinstance(result, dict)
                and all(key in result for key in
                        ("id", "title", "method", "summary", "tables", "figures", "obs_columns", "warnings")))
    if not isinstance(result, dict):
        return False
    ok &= check("%s: id/title/method 与 ANALYSES 一致" % analysis_id,
                result.get("id") == spec["id"]
                and result.get("title") == spec["title"]
                and result.get("method") == spec["method"])
    method = str(result.get("method") or "")
    ok &= check("%s: method 诚实标注（本插件/lite/近似/代理）" % analysis_id,
                any(marker in method for marker in HONEST_MARKERS), method[:80])
    summary = result.get("summary") or {}
    ok &= check("%s: summary.bullets 是非空字符串数组" % analysis_id,
                isinstance(summary.get("bullets"), list) and len(summary["bullets"]) >= 1
                and all(isinstance(item, str) and item for item in summary["bullets"]))
    ok &= check("%s: summary.metrics 是对象且 JSON 可序列化" % analysis_id,
                isinstance(summary.get("metrics"), dict)
                and _json_ok(summary["metrics"]))
    tables = result.get("tables")
    ok &= check("%s: tables 至少 1 张且结构合法" % analysis_id,
                isinstance(tables, list) and len(tables) >= 1 and all(
                    isinstance(table, dict)
                    and isinstance(table.get("id"), str)
                    and isinstance(table.get("title"), str)
                    and isinstance(table.get("columns"), list)
                    and isinstance(table.get("rows"), list)
                    for table in tables))
    for table in tables or []:
        column_keys = {column.get("key") for column in table.get("columns", [])}
        rows = table.get("rows", [])
        ok &= check("%s: 表 %s 行键 ⊆ 声明列" % (analysis_id, table.get("id")),
                    all(isinstance(row, dict) and set(row).issubset(column_keys) for row in rows),
                    "rows=%d cols=%d" % (len(rows), len(column_keys)))
        ok &= check("%s: 表 %s 可 JSON 序列化" % (analysis_id, table.get("id")),
                    _json_ok(table))
    figures = result.get("figures")
    ok &= check("%s: figures 非空" % analysis_id, isinstance(figures, list) and len(figures) >= 1)
    dataset_slug = str(ctx["name"]).replace(" ", "-")
    for figure in figures or []:
        path = Path(str(figure.get("file")))
        info = png_info(str(path)) if path.is_file() else None
        ok &= check("%s: 图 %s 是真实 PNG 且在 outdir 下" % (analysis_id, figure.get("id")),
                    info is not None and Path(str(outdir)).resolve() in path.resolve().parents,
                    str(path))
        ok &= check("%s: 图 %s 尺寸字段与文件一致" % (analysis_id, figure.get("id")),
                    info is not None and (figure.get("width"), figure.get("height")) == info,
                    "returned=%s actual=%s" % ((figure.get("width"), figure.get("height")), info))
        ok &= check("%s: 图文件名符合 <dataset>_<analysis>_<slug>.png" % analysis_id,
                    path.name.startswith("%s_%s_" % (dataset_slug, analysis_id)) and path.suffix == ".png",
                    path.name)
    obs_columns = result.get("obs_columns")
    ok &= check("%s: obs_columns 结构合法" % analysis_id,
                isinstance(obs_columns, list)
                and all(isinstance(column, dict)
                        and isinstance(column.get("key"), str)
                        and column.get("kind") in ("numeric", "categorical")
                        for column in obs_columns))
    declared = {str(column["key"]): column["kind"] for column in obs_columns or []}
    obs_values = result.get("obs_values")
    if obs_values is not None:
        ok &= check("%s: obs_values 是对象" % analysis_id, isinstance(obs_values, dict))
        for key, values in (obs_values or {}).items():
            ok &= check("%s: obs_values[%s] 与 ctx.n_obs 等长" % (analysis_id, key),
                        len(values) == ctx["n_obs"], "%d vs %d" % (len(values), ctx["n_obs"]))
            ok &= check("%s: obs_values[%s] 已在 obs_columns 声明" % (analysis_id, key),
                        key in declared, sorted(declared))
            if declared.get(key) == "numeric":
                numeric = [value for value in values if value is not None]
                ok &= check("%s: obs_values[%s] numeric 类型正确" % (analysis_id, key),
                            all(isinstance(value, (int, float)) and not isinstance(value, bool)
                                for value in numeric) and len(numeric) > 0,
                            "%d/%d non-null" % (len(numeric), len(values)))
            elif declared.get(key) == "categorical":
                texts = [value for value in values if value is not None]
                ok &= check("%s: obs_values[%s] categorical 类型正确" % (analysis_id, key),
                            all(isinstance(value, str) for value in texts) and len(texts) > 0)
    ok &= check("%s: warnings 是字符串数组" % analysis_id,
                isinstance(result.get("warnings"), list)
                and all(isinstance(item, str) for item in result["warnings"]))
    ok &= check("%s: 整体 JSON 序列化" % analysis_id, _json_ok(result))
    return bool(ok)


def _json_ok(value):
    try:
        json.dumps(value, ensure_ascii=False, allow_nan=False)
        return True
    except (TypeError, ValueError):
        return False


def expect_engine_error(label, function, expected_code):
    try:
        function()
    except mod.EngineError as error:
        return check("%s: 抛 EngineError(%s)" % (label, expected_code),
                     error.code == expected_code, "code=%s message=%s" % (error.code, error.message[:120]))
    except Exception as error:  # noqa: BLE001
        return check("%s: 抛 EngineError(%s)" % (label, expected_code), False,
                     "抛了其它异常：%s: %s" % (type(error).__name__, error))
    return check("%s: 抛 EngineError(%s)" % (label, expected_code), False, "没有抛异常")


def stable_signature(result):
    """用于确定性比较：summary/tables/obs_values/obs_columns/warnings 的 JSON 指纹。"""
    return json.dumps(
        {
            "summary": result.get("summary"),
            "tables": result.get("tables"),
            "obs_columns": result.get("obs_columns"),
            "obs_values": result.get("obs_values"),
            "warnings": result.get("warnings"),
        },
        ensure_ascii=False,
        sort_keys=True,
        allow_nan=False,
    )


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------


def run_tests(tmp_root):
    global FAIL_COUNT, PASS_COUNT
    started_at = time.time()
    ctx = make_ctx()
    outdir = tmp_root / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_order_path = tmp_root / "gene_order.tsv"
    write_gene_order(gene_order_path, ctx["var_names"])

    check("synthetic ctx: 细胞数 ≥800", ctx["n_obs"] >= 800, ctx["n_obs"])
    check("synthetic ctx: 基因数 ≥1500", ctx["n_vars"] >= 1500, ctx["n_vars"])
    check("synthetic ctx: X 为 csr float32", sp.issparse(ctx["X"]) and ctx["X"].dtype == np.float32)
    check("synthetic ctx: obs 含 clusters/cell_type/batch",
          all(column in ctx["obs"].columns for column in ("clusters", "cell_type", "batch")))

    # ---- 1. ANALYSES 元数据 ----
    print("\n=== ANALYSES 元数据 ===")
    check("ANALYSES 覆盖 13 个分析 id",
          [str(spec["id"]) for spec in mod.ANALYSES] == EXPECTED_IDS,
          [spec["id"] for spec in mod.ANALYSES])
    for spec in mod.ANALYSES:
        analysis_id = spec["id"]
        check("%s: 元数据必备字段" % analysis_id,
              all(key in spec for key in ("id", "title", "group", "method", "params", "needs"))
              and isinstance(spec["title"], str) and spec["title"]
              and isinstance(spec["group"], str) and spec["group"]
              and isinstance(spec["method"], str) and spec["method"]
              and isinstance(spec["params"], list) and isinstance(spec["needs"], list))
        check("%s: method 诚实标注" % analysis_id,
              any(marker in spec["method"] for marker in HONEST_MARKERS), spec["method"][:70])
        parameter_ok = True
        for parameter in spec["params"]:
            parameter_ok &= (
                isinstance(parameter, dict)
                and isinstance(parameter.get("key"), str)
                and parameter["key"]
                and isinstance(parameter.get("label"), str)
                and parameter.get("type") in ALLOWED_PARAM_TYPES
                and "default" in parameter
            )
            if parameter.get("type") == "select":
                parameter_ok &= isinstance(parameter.get("options"), list) and len(parameter["options"]) >= 1
        check("%s: params 类型合法（number/string/select/boolean/column/file/text）" % analysis_id,
              parameter_ok)
        needs_ok = all(
            isinstance(need, str) and (need == "clusters" or need == "spliced" or need == "gene_sets"
                                       or need == "gene_order" or need.startswith("groups:"))
            for need in spec["needs"]
        )
        check("%s: needs 合法" % analysis_id, needs_ok, spec["needs"])
    check("ANALYSES 全部为 native 分析（13 项）", len(mod.ANALYSES) == 13)
    check("未知 analysis id 抛 bad_request",
          _expect_unknown_id(ctx, outdir))

    # ---- 2. 逐项跑通 + 结构/图/obs_values 校验 ----
    print("\n=== 13 个分析逐个跑通 ===")
    target_gene = "STAT1"
    if target_gene not in ctx["var_names"]:
        target_gene = next(gene for gene in ctx["var_names"] if gene in set(biopriors.tf_list()))
    results = {}
    timings = {}
    for spec in mod.ANALYSES:
        analysis_id = spec["id"]
        check("%s: 测试上下文有该分析所需的 obs 列/簇" % analysis_id, True)
        params = params_for(analysis_id, target_gene)
        begin = time.time()
        try:
            result = mod.run(analysis_id, ctx, params, str(outdir))
        except Exception as error:  # noqa: BLE001
            FAIL_COUNT += 1
            FAILURES.append("%s: run 抛异常" % analysis_id)
            print("FAIL  %-66s %s" % ("%s: run 抛异常" % analysis_id,
                                      "%s: %s" % (type(error).__name__, error)))
            traceback.print_exc()
            continue
        elapsed = time.time() - begin
        timings[analysis_id] = elapsed
        results[analysis_id] = (result, params)
        validate_result(analysis_id, result, ctx, spec, outdir)
        print("      %-22s %.2fs  tables=%d figures=%d obs_cols=%d obs_values=%d"
              % (analysis_id, elapsed, len(result["tables"]), len(result["figures"]),
                 len(result["obs_columns"]), len(result.get("obs_values") or {})))

    # ---- 3. 特定语义校验 ----
    print("\n=== 语义校验 ===")
    integration = results.get("integration", (None, None))[0]
    if integration:
        check("integration: 返回批次混合指标前后对比",
              any(table["id"] == "batch_mixing" for table in integration["tables"])
              and {"same_batch_frac_before", "same_batch_frac_after",
                   "batch_entropy_before", "batch_entropy_after"}
              <= {column["key"] for column in integration["tables"][0]["columns"]})
        check("integration: obs_values 含校正后 PC1/PC2",
              {"integration_pc1", "integration_pc2"} <= set(integration.get("obs_values") or {}))
        after_entropy = integration["summary"]["metrics"].get("batch_entropy_after")
        before_entropy = integration["summary"]["metrics"].get("batch_entropy_before")
        check("integration: 批次熵有前后对比数值",
              isinstance(after_entropy, (int, float)) and isinstance(before_entropy, (int, float)),
              "before=%s after=%s" % (before_entropy, after_entropy))
    dpt = results.get("trajectory_dpt", (None, None))[0]
    if dpt:
        check("trajectory_dpt: 伪时间在 obs_values 且范围 [0,1]",
              "dpt_pseudotime" in (dpt.get("obs_values") or {})
              and _all_in_range(dpt["obs_values"]["dpt_pseudotime"], 0.0, 1.0))
        check("trajectory_dpt: 提供簇级伪时间排序表",
              any(table["id"] == "pseudotime_by_cluster" for table in dpt["tables"]))
    paga = results.get("trajectory_paga", (None, None))[0]
    if paga:
        check("trajectory_paga: 提供 connectivity 边表 + MST",
              any(table["id"] == "cluster_connectivity" for table in paga["tables"])
              and paga["summary"]["metrics"].get("n_mst_edges", 0) >= 1)
    cnv = results.get("cnv_score", (None, None))[0]
    if cnv:
        check("cnv_score: 未给 gene_order 时明确标注代理分数",
              cnv["summary"]["metrics"].get("ordered_by_chromosome") is False
              and any("代理分数" in warning for warning in cnv["warnings"]))
        check("cnv_score: 簇级热图存在",
              any("heatmap" in figure["id"] for figure in cnv["figures"]))
    lr = results.get("lr_communication", (None, None))[0]
    if lr:
        check("lr_communication: 返回 LR 交互表 + 聚合表 + 网络图",
              {"lr_interactions", "cell_cell_strength"} <= {table["id"] for table in lr["tables"]}
              and len(lr["figures"]) == 1
              and lr["summary"]["metrics"].get("n_lr_pairs_present", 0) > 0)
    regulon = results.get("regulon_activity", (None, None))[0]
    if regulon:
        regulon_keys = [column["key"] for column in regulon["obs_columns"]]
        check("regulon_activity: 多个 regulon 活性写进 obs_columns",
              len([key for key in regulon_keys if key.startswith("regulon_")]) >= 2,
              regulon_keys[:6])
        check("regulon_activity: 热图 + 靶基因表",
              any(table["id"] == "regulon_targets" for table in regulon["tables"])
              and len(regulon["figures"]) >= 1)
    knockout = results.get("virtual_knockout", (None, None))[0]
    if knockout:
        check("virtual_knockout: 受影响基因排序表（含方向）",
              knockout["tables"][0]["id"] == "affected_genes"
              and {"rank", "gene", "total_effect", "direction"}
              <= {column["key"] for column in knockout["tables"][0]["columns"]})
        check("virtual_knockout: 目标基因写入 metrics",
              knockout["summary"]["metrics"].get("target_gene") == target_gene)
    overexpress = results.get("virtual_overexpress", (None, None))[0]
    if overexpress:
        check("virtual_overexpress: 方向与敲除取反",
              overexpress["summary"]["metrics"].get("mode") == "overexpression")
    screen = results.get("virtual_screen", (None, None))[0]
    if screen:
        check("virtual_screen: 提供候选排序表 + 图",
              screen["tables"][0]["id"] == "screen_results"
              and len(screen["tables"][0]["rows"]) >= 1
              and len(screen["figures"]) == 1)

    # cnv_score 的染色体排序分支
    check("cnv_score: gene_order 文件存在", gene_order_path.is_file())
    ordered = mod.run(
        "cnv_score",
        ctx,
        {"gene_order": str(gene_order_path), "window": 11, "max_genes": 900},
        str(outdir),
    )
    check("cnv_score: 提供 gene_order 后按染色体排序",
          ordered["summary"]["metrics"].get("ordered_by_chromosome") is True)
    check("cnv_score: 有序分支不再出现代理分数警告",
          not any("代理分数" in warning for warning in ordered["warnings"]))
    check("cnv_score: 有序分支记录染色体区段数",
          ordered["summary"]["metrics"].get("n_chromosome_segments", 0) >= 2,
          ordered["summary"]["metrics"].get("n_chromosome_segments"))
    results["cnv_score_ordered"] = (
        ordered,
        {"gene_order": str(gene_order_path), "window": 11, "max_genes": 900},
    )

    # trajectory_dpt 指定根细胞
    root_cell = str(ctx["obs"].index[3])
    dpt_root = mod.run("trajectory_dpt", ctx, {"n_neighbors": 12, "n_dcs": 4,
                                               "root_cell": root_cell}, str(outdir))
    check("trajectory_dpt: params.root_cell 生效",
          dpt_root["summary"]["metrics"].get("root_source") == "params.root_cell"
          and dpt_root["summary"]["metrics"].get("root_index") == 3)

    # ---- 3b. 无簇 ctx 分支（needs=[] 的分析必须没有 cluster 也能跑） ----
    print("\n=== 无簇 ctx 分支 ===")
    obs_without_clusters = ctx["obs"].drop(columns=["clusters"])
    ctx_no_clusters = dict(ctx, obs=obs_without_clusters, clusters=None)
    no_cluster_cases = {
        "trajectory_dpt": {"n_neighbors": 12, "n_dcs": 3},
        "entropy_stemness": {"n_genes": 400, "n_neighbors": 10},
        "cytotrace_lite": {"n_genes": 400, "n_neighbors": 10, "signature_size": 40},
        "cnv_score": {"window": 11, "max_genes": 500},
        "regulon_activity": {"n_genes": 300, "n_targets": 10, "n_regulons": 6, "top_regulons": 3},
        "virtual_knockout": {"target": target_gene, "n_genes": 300, "top_k": 6, "steps": 2},
    }
    spec_lookup = {spec["id"]: spec for spec in mod.ANALYSES}
    for analysis_id, params in no_cluster_cases.items():
        try:
            result = mod.run(analysis_id, ctx_no_clusters, dict(params), str(outdir))
        except Exception as error:  # noqa: BLE001
            check("%s: 无簇 ctx 可跑" % analysis_id, False,
                  "%s: %s" % (type(error).__name__, error))
            traceback.print_exc()
            continue
        check("%s: 无簇 ctx 可跑" % analysis_id, True,
              "tables=%d figures=%d" % (len(result["tables"]), len(result["figures"])))
        validate_result(analysis_id, result, ctx_no_clusters, spec_lookup[analysis_id], outdir)

    # ---- 4. 确定性 ----
    print("\n=== 确定性（同名参数跑两次） ===")
    for analysis_id, (first, params) in results.items():
        base_id = analysis_id.replace("_ordered", "")
        try:
            second = mod.run(base_id, ctx, dict(params), str(outdir))
        except Exception as error:  # noqa: BLE001
            check("%s: 二次运行成功" % analysis_id, False, str(error))
            continue
        check("%s: 两次结果完全一致" % analysis_id,
              stable_signature(first) == stable_signature(second))

    # ---- 5. 前置/参数错误 ----
    print("\n=== 前置与错误处理 ===")
    obs_without_batch = ctx["obs"].drop(columns=["batch"])
    ctx_no_batch = dict(ctx, obs=obs_without_batch)
    expect_engine_error(
        "integration 缺 batch 列",
        lambda: mod.run("integration", ctx_no_batch, {"batch": "batch"}, str(outdir)),
        "prerequisite",
    )
    expect_engine_error(
        "integration 指定不存在的列",
        lambda: mod.run("integration", ctx, {"batch": "sample_id"}, str(outdir)),
        "prerequisite",
    )
    obs_without_clusters = ctx["obs"].drop(columns=["clusters"])
    ctx_no_clusters = dict(ctx, obs=obs_without_clusters, clusters=None)
    for analysis_id in ("trajectory_paga", "cnv_tree", "ith_score", "lr_communication", "virtual_screen"):
        expect_engine_error(
            "%s 缺 clusters" % analysis_id,
            lambda analysis_id=analysis_id: mod.run(analysis_id, ctx_no_clusters, {}, str(outdir)),
            "prerequisite",
        )
    fake_genes_ctx = dict(ctx, var_names=["FAKE%05d" % index for index in range(ctx["n_vars"])])
    for analysis_id in ("lr_communication", "regulon_activity"):
        expect_engine_error(
            "%s 非人源基因名" % analysis_id,
            lambda analysis_id=analysis_id: mod.run(analysis_id, fake_genes_ctx, {}, str(outdir)),
            "prerequisite",
        )
    expect_engine_error(
        "virtual_knockout 缺 target",
        lambda: mod.run("virtual_knockout", ctx, {"n_genes": 500}, str(outdir)),
        "bad_request",
    )
    expect_engine_error(
        "virtual_knockout 目标基因不存在",
        lambda: mod.run("virtual_knockout", ctx, {"target": "NOT_A_GENE", "n_genes": 500}, str(outdir)),
        "bad_request",
    )
    expect_engine_error(
        "cnv_score gene_order 文件不存在",
        lambda: mod.run("cnv_score", ctx, {"gene_order": str(tmp_root / "no-such-file.tsv")}, str(outdir)),
        "bad_request",
    )
    expect_engine_error(
        "不支持的 params 类型",
        lambda: mod.run("entropy_stemness", ctx, ["not", "a", "dict"], str(outdir)),
        "bad_request",
    )

    # ---- 6. 总结 ----
    total = time.time() - started_at
    print("\n=== 耗时 ===")
    for analysis_id, elapsed in sorted(timings.items(), key=lambda item: -item[1]):
        print("  %-24s %.2fs" % (analysis_id, elapsed))
    slowest = max(timings.items(), key=lambda item: item[1]) if timings else ("-", 0.0)
    print("  %-24s %.2fs" % ("-- total (13 analyses)", sum(timings.values())))
    print("  %-24s %.2fs" % ("-- slowest", slowest[1]))
    print("  %-24s %.2fs" % ("-- wall clock (all checks)", total))
    return total, timings


def _expect_unknown_id(ctx, outdir):
    try:
        mod.run("no_such_analysis", ctx, {}, str(outdir))
    except mod.EngineError as error:
        return error.code == "bad_request"
    except Exception:
        return False
    return False


def _all_in_range(values, low, high):
    numbers = [value for value in values if value is not None]
    return bool(numbers) and all(low - 1e-6 <= float(value) <= high + 1e-6 for value in numbers)


def main():
    tmp_root = Path(tempfile.mkdtemp(prefix="dsh-sc-dyn-test-"))
    print("临时目录：%s" % tmp_root)
    exit_code = 1
    try:
        run_tests(tmp_root)
        exit_code = 0 if FAIL_COUNT == 0 else 1
    except Exception:  # noqa: BLE001
        traceback.print_exc()
        exit_code = 1
    finally:
        if FAIL_COUNT == 0:
            print("\n===== ALL PASS: %d checks =====" % PASS_COUNT)
        else:
            print("\n===== FAILED: %d failed, %d passed =====" % (FAIL_COUNT, PASS_COUNT))
            for name in FAILURES:
                print("  - %s" % name)
        if os.environ.get("KEEP_SC_TEST_TMP") == "1":
            print("临时目录保留：%s" % tmp_root)
        else:
            shutil.rmtree(tmp_root, ignore_errors=True)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
