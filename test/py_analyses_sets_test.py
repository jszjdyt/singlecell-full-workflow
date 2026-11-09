#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sc-analysis-a 模块测试：自造 ctx，逐个跑通 §9.3 的 14 个分析。

运行：``python3 test/py_analyses_sets_test.py``（不依赖 node / 插件安装）
覆盖：
  1. ANALYSES 覆盖 §9.3 sc-analysis-a 的 14 个 id，元数据结构合法；
  2. 14 个分析逐个跑通，结果结构（§9.2 固定键）与 table/figure/obs 结构合法；
  3. 图存在、是 PNG magic、宽高 > 200；
  4. obs_values 与 ctx 细胞数等长、kind 匹配（数值列允许 null）；
  5. 同参数两次结果（tables/summary）完全一致（确定性）；
  6. 缺 group 列时 proportion_test/de_group/volcano/proportion_sankey 抛 prerequisite；
  7. 用户 GMT 文件通路（ORA）可跑通；未知 id / 非法参数给可读错误。
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import shutil
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
LIB = HERE.parent / "lib"
MPL_DIR = tempfile.mkdtemp(prefix="sc-analyses-sets-mpl-")
os.environ["MPLBACKEND"] = "Agg"
os.environ["MPLCONFIGDIR"] = MPL_DIR
sys.dont_write_bytecode = True


def _load(module_name, filename):
    spec = importlib.util.spec_from_file_location(module_name, str(LIB / filename))
    if spec is None or spec.loader is None:
        raise RuntimeError("无法加载 %s" % filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GENE_SETS = _load("_test_gene_sets_builtin", "gene_sets_builtin.py")
MODULE = _load("_test_sc_analyses_sets", "sc_analyses_sets.py")

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
CELL_TYPES = ["T cell", "B cell", "NK cell", "Monocyte"]
MARKER_GENES = {
    "T cell": ["CD3D", "CD3E", "CD3G", "IL7R", "LCK", "CD28"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD19", "IGHM", "BANK1"],
    "NK cell": ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "KLRF1"],
    "Monocyte": ["LYZ", "CD14", "FCGR3A", "S100A8", "S100A9", "CTSS"],
}
EFFECT_GENES = ["IFIT1", "ISG15", "STAT1", "MX1", "IFI44L", "OAS1"]


def _gene_names(n_vars):
    """真实人类 symbol 轮转取前 1100 个（保证内置通路集有交集），其余合成基因补足。"""
    groups = [list(genes) for genes in GENE_SETS.BUILTIN_GMT.values()]
    groups.append(list(GENE_SETS.CELL_CYCLE_S))
    groups.append(list(GENE_SETS.CELL_CYCLE_G2M))
    ordered = []
    seen = set()
    for index in range(max(len(group) for group in groups)):
        for group in groups:
            if index < len(group):
                gene = group[index]
                if gene not in seen:
                    seen.add(gene)
                    ordered.append(gene)
    names = ordered[:1100]
    synthetic_index = 0
    while len(names) < n_vars:
        names.append("GENE_%04d" % synthetic_index)
        synthetic_index += 1
    return names[:n_vars]


def _log_normalize(counts):
    row_sum = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.zeros_like(row_sum)
    np.divide(10000.0, row_sum, out=scale, where=row_sum > 0)
    matrix = (sp.diags(scale) @ counts).tocsr()
    matrix.data = np.log1p(matrix.data)
    matrix.eliminate_zeros()
    return matrix.astype(np.float32)


def build_ctx(n_cells=900, n_vars=1600, seed=0):
    """自造带真实基因结构 + 4 簇/4 细胞类型/2 组/数值元数据的分析上下文。"""
    rng = np.random.default_rng(seed)
    genes = _gene_names(n_vars)
    index = {gene: position for position, gene in enumerate(genes)}
    cluster_sizes = [0] * 4
    remaining = n_cells
    for cluster in range(3):
        size = max(40, int(round(n_cells * 0.25)))
        cluster_sizes[cluster] = size
        remaining -= size
    cluster_sizes[3] = remaining
    cluster_labels = np.repeat(np.arange(4), cluster_sizes)
    cell_type = np.array([CELL_TYPES[cluster] for cluster in cluster_labels], dtype=object)
    healthy_probability = [0.62, 0.55, 0.50, 0.32]
    random_draw = rng.random(n_cells)
    group = np.array(["healthy"] * n_cells, dtype=object)
    for cluster_index in range(4):
        mask = cluster_labels == cluster_index
        group[mask] = np.where(
            random_draw[mask] < healthy_probability[cluster_index], "healthy", "disease")

    counts = rng.poisson(1.2, size=(n_cells, n_vars)).astype(np.float64)
    for cell_type_name, markers in MARKER_GENES.items():
        mask = cell_type == cell_type_name
        for gene in markers:
            position = index.get(gene)
            if position is not None:
                counts[mask, position] += rng.poisson(9.0, size=int(mask.sum()))
    effector_mask = (cell_type == "T cell") & (group == "disease")
    for gene in EFFECT_GENES:
        position = index.get(gene)
        if position is not None:
            counts[effector_mask, position] += rng.poisson(7.0, size=int(effector_mask.sum()))
    cell_cycle_genes = list(GENE_SETS.CELL_CYCLE_S[:20]) + list(GENE_SETS.CELL_CYCLE_G2M[:20])
    for gene in cell_cycle_genes:
        position = index.get(gene)
        if position is not None:
            counts[:, position] += rng.poisson(1.0, size=n_cells)
    counts += rng.poisson(0.6, size=(n_cells, n_vars))
    counts = sp.csr_matrix(np.maximum(counts, 0.0).astype(np.float32))
    matrix = _log_normalize(counts)

    obs = pd.DataFrame(
        {
            "clusters": [str(value) for value in cluster_labels.tolist()],
            "cell_type": [str(value) for value in cell_type.tolist()],
            "group": [str(value) for value in group.tolist()],
            "age": np.round(rng.normal(50.0, 15.0, size=n_cells), 2),
            "batch": np.where(rng.random(n_cells) < 0.5, "batch1", "batch2"),
        },
        index=["cell_%04d" % (index + 1) for index in range(n_cells)],
    )
    from sklearn.decomposition import PCA

    dense = np.asarray(matrix.toarray(), dtype=np.float64)
    pca = PCA(n_components=20, svd_solver="randomized", random_state=0)
    scores = pca.fit_transform(dense)
    embedding = scores[:, :2].copy()
    return {
        "dataset_id": "ds-test",
        "name": "synthetic-test",
        "n_obs": n_cells,
        "n_vars": n_vars,
        "X": matrix,
        "X_kind": "lognorm",
        "counts": counts,
        "obs": obs,
        "var_names": genes,
        "obsm": {"X_pca": scores, "X_emb": embedding},
        "clusters": {"key": "clusters", "labels": obs["clusters"].tolist()},
        "group_keys": ["clusters", "cell_type", "group", "batch"],
        "engine_mode": "builtin",
        "rng": np.random.default_rng(0),
    }


# ---------------------------------------------------------------------------
# 校验工具
# ---------------------------------------------------------------------------

FAILURES = []
PASS_COUNT = 0


def check(condition, message, failures=None):
    global PASS_COUNT
    if condition:
        PASS_COUNT += 1
        return True
    (FAILURES if failures is None else failures).append(message)
    return False


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _validate_figure(figure, failures, label):
    check(isinstance(figure, dict), "%s: figure 不是 dict" % label, failures)
    for key in ("id", "title", "file", "width", "height"):
        check(key in figure, "%s: figure 缺键 %s" % (label, key), failures)
    path = Path(str(figure.get("file") or ""))
    check(path.is_absolute(), "%s: figure 路径不是绝对路径" % label, failures)
    check(path.exists() and path.is_file(), "%s: figure 文件不存在：%s" % (label, path), failures)
    if path.exists():
        with open(path, "rb") as handle:
            header = handle.read(8)
        check(header == PNG_MAGIC, "%s: figure 不是 PNG magic：%s" % (label, path), failures)
        check(int(figure.get("width") or 0) > 200 and int(figure.get("height") or 0) > 200,
              "%s: figure 宽高过小：%sx%s" % (label, figure.get("width"), figure.get("height")),
              failures)


def _validate_table(table, failures, label):
    check(isinstance(table, dict), "%s: table 不是 dict" % label, failures)
    for key in ("id", "title", "columns", "rows"):
        check(key in table, "%s: table 缺键 %s" % (label, key), failures)
    columns = table.get("columns")
    check(isinstance(columns, list) and len(columns) > 0, "%s: table.columns 为空" % label, failures)
    keys = set()
    kinds = {}
    for column in columns or []:
        check(isinstance(column, dict), "%s: column 不是 dict" % label, failures)
        for key in ("key", "label", "type"):
            check(key in column, "%s: column 缺键 %s" % (label, key), failures)
        check(column.get("type") in ("string", "number"),
              "%s: column type 非法：%r" % (label, column.get("type")), failures)
        keys.add(str(column.get("key")))
        kinds[str(column.get("key"))] = column.get("type")
    rows = table.get("rows")
    check(isinstance(rows, list), "%s: table.rows 不是 list" % label, failures)
    for row in rows or []:
        check(isinstance(row, dict), "%s: row 不是 dict" % label, failures)
        for key in keys:
            check(key in row, "%s: row 缺列 %s" % (label, key), failures)
        for key, kind in kinds.items():
            if key not in row:
                continue
            value = row[key]
            if kind == "number":
                check(value is None or _is_number(value),
                      "%s: 数值列 %s 值非法：%r" % (label, key, value), failures)
            else:
                check(value is None or isinstance(value, str),
                      "%s: 字符串列 %s 值非法：%r" % (label, key, value), failures)
    try:
        json.dumps(table, ensure_ascii=False, allow_nan=False)
    except Exception as error:
        failures.append("%s: table 不能 JSON 序列化：%s" % (label, error))


def _validate_result(result, analysis_id, ctx, failures):
    label = "analysis:%s" % analysis_id
    if not isinstance(result, dict):
        failures.append("%s: 返回值不是 dict：%r" % (label, type(result)))
        return
    for key in ("id", "title", "group", "method", "summary", "tables", "figures",
                "obs_columns", "warnings"):
        check(key in result, "%s: 缺少固定键 %s" % (label, key), failures)
    for key in ("title", "group", "method"):
        value = result.get(key)
        check(isinstance(value, str) and value.strip() != "", "%s: %s 为空" % (label, key), failures)
    check(result.get("id") == analysis_id, "%s: id 不一致：%r" % (label, result.get("id")), failures)
    summary = result.get("summary")
    check(isinstance(summary, dict), "%s: summary 不是 dict" % label, failures)
    if isinstance(summary, dict):
        bullets = summary.get("bullets")
        check(isinstance(bullets, list) and len(bullets) > 0,
              "%s: summary.bullets 为空" % label, failures)
        check(all(isinstance(item, str) for item in bullets or []),
              "%s: bullets 含非字符串" % label, failures)
        check(isinstance(summary.get("metrics"), dict), "%s: metrics 不是 dict" % label, failures)
    tables = result.get("tables")
    check(isinstance(tables, list) and len(tables) > 0, "%s: tables 为空" % label, failures)
    for table in tables or []:
        _validate_table(table, failures, label)
    figures = result.get("figures")
    check(isinstance(figures, list) and len(figures) > 0, "%s: figures 为空" % label, failures)
    for figure in figures or []:
        _validate_figure(figure, failures, label)
    obs_columns = result.get("obs_columns")
    obs_values = result.get("obs_values")
    check(isinstance(obs_columns, list), "%s: obs_columns 不是 list" % label, failures)
    check(isinstance(obs_values, dict), "%s: obs_values 不是 dict" % label, failures)
    for column in obs_columns or []:
        check(isinstance(column, dict), "%s: obs_columns 条目不是 dict" % label, failures)
        key = str(column.get("key") or "")
        kind = column.get("kind")
        check(key != "", "%s: obs 列名为空" % label, failures)
        check(kind in ("categorical", "numeric"), "%s: obs kind 非法：%r" % (label, kind), failures)
        check(key in (obs_values or {}), "%s: obs_values 缺少列 %s" % (label, key), failures)
        values = (obs_values or {}).get(key)
        if values is None:
            continue
        check(isinstance(values, list), "%s: obs_values[%s] 不是 list" % (label, key), failures)
        check(len(values) == int(ctx["n_obs"]),
              "%s: obs_values[%s] 长度 %d != 细胞数 %d" % (
                  label, key, len(values), int(ctx["n_obs"])), failures)
        if kind == "categorical":
            check(all(isinstance(value, str) for value in values),
                  "%s: 分类 obs 值含非字符串" % label, failures)
        else:
            check(all(value is None or _is_number(value) for value in values),
                  "%s: 数值 obs 值含非数值" % label, failures)
            check(all(value is None or not (isinstance(value, float) and not math.isfinite(value))
                      for value in values),
                  "%s: 数值 obs 值含 NaN/Inf" % label, failures)
    warnings = result.get("warnings")
    check(isinstance(warnings, list) and all(isinstance(item, str) for item in warnings),
          "%s: warnings 不是字符串 list" % label, failures)
    try:
        json.dumps(result, ensure_ascii=False, allow_nan=False)
    except Exception as error:
        failures.append("%s: 结果不能严格 JSON 序列化：%s" % (label, error))


EXPECTED_IDS = sorted([
    "cell_cycle", "subcluster", "proportion_test", "de_group", "volcano",
    "proportion_sankey", "featureplot", "pathway_ora", "gsea_prerank", "gsva_like",
    "auc_score", "cnmf_programs", "wgcna_lite", "geneset_clinical",
])

PARAM_CASES = {
    "cell_cycle": {},
    "subcluster": {"cluster": "0", "resolution": 0.8, "n_top_genes": 300,
                   "n_neighbors": 10, "n_pcs": 10},
    "proportion_test": {},
    "de_group": {"top_n": 20},
    "volcano": {},
    "proportion_sankey": {},
    "featureplot": {},
    "pathway_ora": {},
    "gsea_prerank": {"n_perm": 100},
    "gsva_like": {"n_sets": 6},
    "auc_score": {"n_sets": 6},
    "cnmf_programs": {"k": 4, "n_top_genes": 300, "max_iter": 100},
    "wgcna_lite": {"n_top_genes": 200, "min_module_size": 5},
    "geneset_clinical": {},
}

OBS_EXPECTATIONS = {
    "cell_cycle": {"phase", "S_score", "G2M_score"},
    "subcluster": {"subcluster"},
}
OBS_PREFIXES = {
    "gsva_like": "gsva_",
    "auc_score": "auc_",
    "cnmf_programs": "nmf_",
    "wgcna_lite": "wgcna_",
}


def main():
    print("=" * 78)
    print("sc-analysis-a 模块测试（§9.3 的 14 个分析）")
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. 元数据校验
    # ------------------------------------------------------------------
    print("\n[1] ANALYSES 元数据")
    metadata_failures = []
    analyses = list(getattr(MODULE, "ANALYSES", []) or [])
    ids = [str(item.get("id")) for item in analyses if isinstance(item, dict)]
    check(len(analyses) == 14, "ANALYSES 数量不是 14：%d" % len(analyses), metadata_failures)
    check(sorted(ids) == EXPECTED_IDS,
          "ANALYSES ids 与 §9.3 不一致：%r" % (sorted(ids),), metadata_failures)
    allowed_param_types = {"number", "string", "select", "boolean", "column", "file", "text"}
    allowed_needs = {"clusters", "gene_sets", "gene_order", "spliced"}
    for item in analyses:
        entry_id = str(item.get("id"))
        for key in ("id", "title", "group", "method", "params", "needs"):
            check(key in item, "%s: 元数据缺少 %s" % (entry_id, key), metadata_failures)
        check(isinstance(item.get("title"), str) and item["title"].strip() != "",
              "%s: title 为空" % entry_id, metadata_failures)
        check(isinstance(item.get("group"), str) and item["group"].strip() != "",
              "%s: group 为空" % entry_id, metadata_failures)
        check(isinstance(item.get("method"), str) and item["method"].strip() != "",
              "%s: method 为空" % entry_id, metadata_failures)
        params = item.get("params")
        check(isinstance(params, list), "%s: params 不是 list" % entry_id, metadata_failures)
        for param in params or []:
            check(isinstance(param, dict), "%s: param 不是 dict" % entry_id, metadata_failures)
            check(str(param.get("key") or "") != "", "%s: param.key 为空" % entry_id, metadata_failures)
            check(str(param.get("label") or "") != "", "%s: param.label 为空" % entry_id, metadata_failures)
            check(param.get("type") in allowed_param_types,
                  "%s: param %s type 非法：%r" % (entry_id, param.get("key"), param.get("type")),
                  metadata_failures)
            if param.get("type") == "select":
                options = param.get("options")
                check(isinstance(options, list) and len(options) > 0,
                      "%s: select 参数 %s 缺 options" % (entry_id, param.get("key")),
                      metadata_failures)
                check(param.get("default") in options,
                      "%s: select 参数 %s default 不在 options" % (entry_id, param.get("key")),
                      metadata_failures)
        for need in item.get("needs") or []:
            text = str(need)
            check(text.startswith("groups:") or text in allowed_needs,
                  "%s: needs 项非法：%r" % (entry_id, text), metadata_failures)
    check(hasattr(MODULE, "run") and callable(MODULE.run), "模块缺少 run()", metadata_failures)
    check(hasattr(MODULE, "EngineError"), "模块缺少 EngineError", metadata_failures)
    if metadata_failures:
        FAILURES.extend("[metadata] " + item for item in metadata_failures)
        print("    FAIL (%d)" % len(metadata_failures))
        for message in metadata_failures:
            print("      -", message)
    else:
        print("    PASS  (14 个 id / params / needs 结构合法)")

    # ------------------------------------------------------------------
    # 2. 自造 ctx 并逐个跑分析（每个跑两遍验证确定性）
    # ------------------------------------------------------------------
    print("\n[2] 自造 ctx（900 细胞 × 1600 基因，log 标准化 + 4 簇/4 类型/2 组）")
    outdir = Path(tempfile.mkdtemp(prefix="sc-analyses-sets-out-"))
    ctx = build_ctx()
    print("    ctx: n_obs=%d n_vars=%d clusters=%d cell_types=%d groups=%d nnz=%d" % (
        ctx["n_obs"], ctx["n_vars"], ctx["obs"]["clusters"].nunique(),
        ctx["obs"]["cell_type"].nunique(), ctx["obs"]["group"].nunique(), ctx["X"].nnz))

    timings = {}
    results = {}
    analysis_failures = []
    print("\n[3] 逐个分析（跑两遍校验确定性）")
    for analysis_id in EXPECTED_IDS:
        params = PARAM_CASES.get(analysis_id, {})
        failures = []
        started = time.perf_counter()
        try:
            first = MODULE.run(analysis_id, ctx, params, outdir)
        except Exception as error:
            FAILURES.append("analysis:%s 第一次运行异常：%s\n%s" % (
                analysis_id, error, traceback.format_exc(limit=6)))
            print("    %-20s FAIL  运行异常：%s" % (analysis_id, error))
            continue
        elapsed = time.perf_counter() - started
        try:
            second = MODULE.run(analysis_id, ctx, params, outdir)
        except Exception as error:
            FAILURES.append("analysis:%s 第二次运行异常：%s\n%s" % (
                analysis_id, error, traceback.format_exc(limit=6)))
            print("    %-20s FAIL  第二次运行异常：%s" % (analysis_id, error))
            continue
        timings[analysis_id] = elapsed
        results[analysis_id] = first
        _validate_result(first, analysis_id, ctx, failures)
        # 确定性：tables/summary/obs_columns 两次完全一致
        try:
            deterministic = (
                json.dumps(first.get("tables"), ensure_ascii=False, sort_keys=True, allow_nan=False)
                == json.dumps(second.get("tables"), ensure_ascii=False, sort_keys=True, allow_nan=False)
                and json.dumps(first.get("summary"), ensure_ascii=False, sort_keys=True, allow_nan=False)
                == json.dumps(second.get("summary"), ensure_ascii=False, sort_keys=True, allow_nan=False)
                and json.dumps(first.get("obs_columns"), ensure_ascii=False, sort_keys=True)
                == json.dumps(second.get("obs_columns"), ensure_ascii=False, sort_keys=True)
            )
        except Exception:
            deterministic = False
        check(deterministic, "%s: 同参数两次结果不一致（非确定性）" % analysis_id, failures)
        for figure in second.get("figures") or []:
            check(Path(str(figure.get("file"))).exists(),
                  "%s: 第二次运行图文件缺失" % analysis_id, failures)
        # obs 写回期望
        obs_keys = {str(column.get("key")) for column in first.get("obs_columns") or []}
        expected = OBS_EXPECTATIONS.get(analysis_id)
        if expected is not None:
            check(obs_keys == expected,
                  "%s: obs 列 %r != 期望 %r" % (analysis_id, sorted(obs_keys), sorted(expected)),
                  failures)
        prefix = OBS_PREFIXES.get(analysis_id)
        if prefix is not None:
            check(len(obs_keys) > 0 and all(key.startswith(prefix) for key in obs_keys),
                  "%s: obs 列应以前缀 %s 开头：%r" % (analysis_id, prefix, sorted(obs_keys)),
                  failures)
        if analysis_id not in OBS_EXPECTATIONS and analysis_id not in OBS_PREFIXES:
            check(len(obs_keys) == 0,
                  "%s: 不应写 obs 列却写了 %r" % (analysis_id, sorted(obs_keys)), failures)
        if failures:
            analysis_failures.extend(failures)
            print("    %-20s FAIL  %.2fs" % (analysis_id, elapsed))
            for message in failures:
                print("      -", message)
        else:
            n_tables = len(first.get("tables") or [])
            n_figures = len(first.get("figures") or [])
            print("    %-20s PASS  %.2fs  tables=%d figures=%d warnings=%d" % (
                analysis_id, elapsed, n_tables, n_figures, len(first.get("warnings") or [])))
    FAILURES.extend(analysis_failures)

    # ------------------------------------------------------------------
    # 4. 缺 group 列 → prerequisite
    # ------------------------------------------------------------------
    print("\n[4] 缺前置数据抛 prerequisite")
    prerequisite_failures = []
    ctx_no_group = dict(ctx)
    ctx_no_group["obs"] = ctx["obs"].drop(columns=["group"])
    ctx_no_group["group_keys"] = ["clusters", "cell_type", "batch"]
    for analysis_id in ("proportion_test", "de_group", "volcano", "proportion_sankey"):
        try:
            MODULE.run(analysis_id, ctx_no_group, {"group": "group"}, outdir)
            prerequisite_failures.append("%s: 缺 group 列没有抛 prerequisite" % analysis_id)
        except MODULE.EngineError as error:
            check(getattr(error, "code", None) == "prerequisite",
                  "%s: 错误码不是 prerequisite：%r" % (analysis_id, getattr(error, "code", None)),
                  prerequisite_failures)
        except Exception as error:  # noqa: BLE001
            prerequisite_failures.append("%s: 抛了非 EngineError：%r" % (analysis_id, error))
    try:
        MODULE.run("subcluster", ctx, {"cluster": "999"}, outdir)
        prerequisite_failures.append("subcluster: 不存在的簇没有抛 prerequisite")
    except MODULE.EngineError as error:
        check(getattr(error, "code", None) == "prerequisite",
              "subcluster: 错误码不是 prerequisite：%r" % getattr(error, "code", None),
              prerequisite_failures)
    if prerequisite_failures:
        FAILURES.extend("[prerequisite] " + item for item in prerequisite_failures)
        print("    FAIL (%d)" % len(prerequisite_failures))
        for message in prerequisite_failures:
            print("      -", message)
    else:
        print("    PASS  (5 个缺前置场景均抛 EngineError(prerequisite))")

    # ------------------------------------------------------------------
    # 5. 非法参数 / 未知 id
    # ------------------------------------------------------------------
    print("\n[5] 非法输入的错误码")
    error_failures = []
    try:
        MODULE.run("no_such_analysis", ctx, {}, outdir)
        error_failures.append("未知 id 没有抛 not_found")
    except MODULE.EngineError as error:
        check(getattr(error, "code", None) == "not_found",
              "未知 id 错误码不是 not_found：%r" % getattr(error, "code", None), error_failures)
    try:
        MODULE.run("featureplot", ctx, {"kind": "bogus"}, outdir)
        error_failures.append("非法 select 值没有抛 bad_request")
    except MODULE.EngineError as error:
        check(getattr(error, "code", None) == "bad_request",
              "非法 select 错误码不是 bad_request：%r" % getattr(error, "code", None), error_failures)
    try:
        MODULE.run("cell_cycle", ctx, {"scoring": 123}, outdir)
        error_failures.append("非法数值参数没有抛 bad_request")
    except MODULE.EngineError as error:
        check(getattr(error, "code", None) == "bad_request",
              "非法数值错误码不是 bad_request：%r" % getattr(error, "code", None), error_failures)
    if error_failures:
        FAILURES.extend("[error-code] " + item for item in error_failures)
        print("    FAIL (%d)" % len(error_failures))
        for message in error_failures:
            print("      -", message)
    else:
        print("    PASS  (not_found / bad_request 均可读)")

    # ------------------------------------------------------------------
    # 6. 用户 GMT 文件（ORA）
    # ------------------------------------------------------------------
    print("\n[6] 用户 GMT 文件通路分析")
    gmt_failures = []
    gmt_path = Path(tempfile.mkdtemp(prefix="sc-analyses-sets-gmt-")) / "custom.gmt"
    with open(gmt_path, "w", encoding="utf-8") as handle:
        handle.write("T_cell_markers\tsynthetic\t%s\n" % "\t".join(MARKER_GENES["T cell"]))
        handle.write("B_cell_markers\tsynthetic\t%s\n" % "\t".join(MARKER_GENES["B cell"]))
        handle.write("Interferon\t%s\n" % "\t".join(EFFECT_GENES))
    try:
        gmt_result = MODULE.run(
            "pathway_ora", ctx, {"source": "gmt", "gmt_file": str(gmt_path)}, outdir)
        _validate_result(gmt_result, "pathway_ora", ctx, gmt_failures)
        terms = {row["term"] for row in (gmt_result.get("tables") or [{}])[0].get("rows", [])}
        check(terms, "GMT ORA 没有返回任何富集条目", gmt_failures)
    except Exception as error:
        gmt_failures.append("GMT ORA 运行异常：%s\n%s" % (error, traceback.format_exc(limit=6)))
    try:
        gsva_gmt = MODULE.run(
            "gsva_like", ctx, {"gmt_file": str(gmt_path), "n_sets": 3}, outdir)
        _validate_result(gsva_gmt, "gsva_like", ctx, gmt_failures)
    except Exception as error:
        gmt_failures.append("GMT GSVA 运行异常：%s\n%s" % (error, traceback.format_exc(limit=6)))
    if gmt_failures:
        FAILURES.extend("[gmt] " + item for item in gmt_failures)
        print("    FAIL (%d)" % len(gmt_failures))
        for message in gmt_failures:
            print("      -", message)
    else:
        print("    PASS  (用户 .gmt 可用于 ORA / GSVA)")

    # ------------------------------------------------------------------
    # 7. 汇总
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    if timings:
        slowest = max(timings.items(), key=lambda item: item[1])
        print("性能：最慢 %s %.2fs；合计 %.2fs（14 项 × 2 遍 + 图）" % (
            slowest[0], slowest[1], sum(timings.values())))
        for analysis_id, seconds in sorted(timings.items(), key=lambda item: -item[1])[:5]:
            print("      %-20s %.2fs" % (analysis_id, seconds))
    if FAILURES:
        print("\n结果：FAIL（%d 条失败）" % len(FAILURES))
        for message in FAILURES:
            print("  -", message)
        print("PASS 断言数：%d" % PASS_COUNT)
        return 1
    print("\n结果：全部 PASS（%d 条断言，14 个分析）" % PASS_COUNT)
    return 0


if __name__ == "__main__":
    exit_code = 1
    try:
        exit_code = main()
    finally:
        shutil.rmtree(MPL_DIR, ignore_errors=True)
        sys.exit(exit_code)
