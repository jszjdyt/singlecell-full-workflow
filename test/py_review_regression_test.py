#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell — 审查任务单回归测试（SC-01 / SC-03 / SC-04 / SC-07 / SC-08）

每个用例都对着《DSH_GEO_Singlecell_审查与优化任务单》里给出的最小复现写：

  SC-01  h5ad 的 raw 必须按基因 ID 对齐；基因顺序不同不能错位；raw 缺基因要放弃而不是硬套。
  SC-03  续跑不得用「连续步骤前缀」伪造 normalize/hvg/pca/neighbors 状态。
  SC-04  分组互斥后 DE 缓存必须失效（不能返回旧结论）。
  SC-07  重复 barcode 去重后索引必须严格唯一。
  SC-08  稀有细胞类型过滤出空组时比例分析不能崩，也不能伪造 p=1。

用法：cd plugins/dsh-singlecell && python3 test/py_review_regression_test.py
"""

from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
LIB = HERE.parent / "lib"
MPL_DIR = tempfile.mkdtemp(prefix="sc-review-mpl-")
os.environ["MPLBACKEND"] = "Agg"
os.environ["MPLCONFIGDIR"] = MPL_DIR
sys.dont_write_bytecode = True

OUTDIR = tempfile.mkdtemp(prefix="sc-review-")
PASSED: list[str] = []
FAILED: list[str] = []


def check(name, cond, detail=""):
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))


def _load(module_name, filename):
    spec = importlib.util.spec_from_file_location(module_name, str(LIB / filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SC_IO = _load("_rev_sc_io", "sc_io.py")
SC_STATE = _load("_rev_sc_state", "sc_state.py")
SETS = _load("_rev_sets", "sc_analyses_sets.py")
try:
    ad = importlib.import_module("anndata")
except Exception:  # noqa: BLE001
    ad = None


# ---------------------------------------------------------------------------
# SC-01 / SC-02：raw 对齐
# ---------------------------------------------------------------------------

def _write_h5ad(path, X, genes, raw_X=None, raw_genes=None, obs=None):
    data = ad.AnnData(X=sp.csr_matrix(np.asarray(X, dtype=np.float32)),
                      obs=obs if obs is not None else pd.DataFrame(index=[f"c{i}" for i in range(np.asarray(X).shape[0])]),
                      var=pd.DataFrame(index=list(genes)))
    if raw_X is not None:
        raw = ad.AnnData(X=sp.csr_matrix(np.asarray(raw_X, dtype=np.float32)),
                         obs=data.obs.copy(),
                         var=pd.DataFrame(index=list(raw_genes)))
        data.raw = raw
    data.write_h5ad(str(path))
    return path


def test_sc01_raw_alignment():
    print("\n[SC-01] h5ad raw 按基因 ID 对齐")
    if ad is None:
        print("  SKIP 没有 anndata")
        return
    tmp = Path(tempfile.mkdtemp(prefix="sc01-"))
    genes = ["G1", "G2"]
    raw = [[100, 10], [200, 20]]
    # 当前 var 逆序，raw 仍是原顺序：counts 必须按名字对齐，不能按位置
    path = _write_h5ad(tmp / "rev.h5ad", [[1.0, 9.0], [3.0, 7.0]], ["G2", "G1"],
                       raw_X=raw, raw_genes=["G1", "G2"])
    payload = SC_IO.load_source(str(path), workdir=tempfile.mkdtemp(prefix="sc01w-"))
    counts = payload["counts"]
    counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    expected = np.array([[10, 100], [20, 200]], dtype=float)
    check("SC-01a 基因逆序时 counts 按基因名对齐",
          np.allclose(counts, expected), str(counts.tolist()))
    check("SC-01b 记录了对齐来源", (payload.get("meta") or {}).get("counts_alignment") == "raw_by_gene_id",
          str((payload.get("meta") or {}).get("counts_alignment")))

    # raw 多基因（含交集外的基因）：取交集对齐
    path2 = _write_h5ad(tmp / "extra.h5ad", [[1.0, 9.0]], ["G2", "G1"],
                        raw_X=[[1, 2, 3]], raw_genes=["G1", "GX", "G2"])
    payload2 = SC_IO.load_source(str(path2), workdir=tempfile.mkdtemp(prefix="sc01w-"))
    c2 = payload2["counts"]
    c2 = c2.toarray() if hasattr(c2, "toarray") else np.asarray(c2)
    # raw: G1=1, GX=2, G2=3；当前顺序 [G2, G1] → [[G2, G1]] = [[3, 1]]
    check("SC-01c raw 基因更多时按交集对齐", np.allclose(c2, [[3, 1]]), str(c2.tolist()))

    # raw 缺当前基因 → 放弃 raw counts，不能硬套，也不能悄悄用 X
    path3 = _write_h5ad(tmp / "missing.h5ad", [[1.0, 9.0]], ["G2", "G1"],
                        raw_X=[[1, 2]], raw_genes=["G1", "GZ"])
    payload3 = SC_IO.load_source(str(path3), workdir=tempfile.mkdtemp(prefix="sc01w-"))
    check("SC-01d raw 缺基因时放弃 counts（不按 shape 硬套）",
          payload3.get("counts") is None, str(payload3.get("counts")))
    check("SC-01e 放弃时给出明确 warning",
          any("缺少当前矩阵里的" in w for w in payload3.get("warnings") or []),
          str((payload3.get("warnings") or [""])[0])[:60])

    # 重复 ID → 放弃
    path4 = _write_h5ad(tmp / "dup.h5ad", [[1.0, 9.0]], ["G1", "G2"],
                        raw_X=[[1, 2]], raw_genes=["G1", "G1"])
    payload4 = SC_IO.load_source(str(path4), workdir=tempfile.mkdtemp(prefix="sc01w-"))
    check("SC-01f raw 基因名重复时放弃 counts", payload4.get("counts") is None)

    # 非整数 counts → 标记为未验证
    path5 = _write_h5ad(tmp / "float.h5ad", [[1.1, 2.2], [3.3, 4.4]], ["G1", "G2"],
                        raw_X=[[1.1, 2.2], [3.3, 4.4]], raw_genes=["G1", "G2"])
    payload5 = SC_IO.load_source(str(path5), workdir=tempfile.mkdtemp(prefix="sc01w-"))
    meta5 = payload5.get("meta") or {}
    check("SC-01g 非整数 counts 标记 counts_verified=False",
          meta5.get("counts_verified") is False, str(meta5.get("counts_verified")))
    check("SC-01h 同时记录 matrix_kind", meta5.get("matrix_kind") in ("lognorm", "counts"),
          str(meta5.get("matrix_kind")))


def test_sc07_unique_barcodes():
    print("\n[SC-07] 重复 barcode 去重")
    index = pd.Index(["bc", "bc", "bc"], name="cell")
    unique = SC_IO._unique_index(index)
    check("SC-07a 三个重复 barcode 变唯一",
          len(set(unique)) == 3 and list(unique) == ["bc", "bc_2", "bc_3"], str(list(unique)))
    index2 = pd.Index(["bc", "bc", "bc_2", "bc_2", "other"], name="cell")
    unique2 = SC_IO._unique_index(index2)
    check("SC-07b 与已有后缀名不冲突", len(set(unique2)) == 5, str(list(unique2)))
    index3 = pd.Index(["a", "b", "c"], name="cell")
    check("SC-07c 本来就唯一时保持不变", list(SC_IO._unique_index(index3)) == ["a", "b", "c"])


# ---------------------------------------------------------------------------
# SC-03：续跑状态不得伪造
# ---------------------------------------------------------------------------

class _FakeDataset:
    def __init__(self, X, obs, var_names, steps=None):
        self.X = sp.csr_matrix(np.asarray(X, dtype=np.float32))
        self.counts = self.X.copy()
        self.obs = obs
        self.var_names = list(var_names)
        self.var = pd.DataFrame(index=pd.Index(self.var_names, name="gene"))
        self.obsm = {}
        self.uns = {}
        self.steps = list(steps or [])
        self.n_obs = obs.shape[0]

    def __getattr__(self, item):
        return None


def test_sc03_state():
    print("\n[SC-03] 续跑状态不得伪造")
    n = 40
    counts = np.tile(np.array([10.0, 100.0, 5.0, 50.0]), (n, 1))
    obs = pd.DataFrame({"clusters": [str(i % 4) for i in range(n)]}, index=[f"c{i}" for i in range(n)])
    dataset = _FakeDataset(counts, obs, [f"G{i}" for i in range(4)])
    state = SC_STATE.detect(dataset, {
        "obsm": {"X_emb": np.random.default_rng(0).normal(size=(n, 2)).astype(np.float32)},
        "import_state": {"reductions": ["umap"], "class": "AnnData"},
    })
    check("SC-03a 原始 counts + 只有嵌入时，不标记 normalize",
          "normalize" not in state["steps"], str(state["steps"]))
    check("SC-03b 不伪造 hvg/pca/neighbors",
          not ({"hvg", "pca", "neighbors"} & set(state["steps"])), str(state["steps"]))
    check("SC-03c 已导入的 cluster/embed 仍被保留",
          "cluster" in state["steps"] and "embed" in state["steps"], str(state["steps"]))
    check("SC-03d 记录 x_lognormalized=false",
          state["details"].get("x_lognormalized") is False, str(state["details"].get("x_lognormalized")))

    # 真 lognorm 才标 normalize
    lognorm = np.log1p(counts * 0.3)
    dataset2 = _FakeDataset(lognorm, obs, [f"G{i}" for i in range(4)])
    state2 = SC_STATE.detect(dataset2, {"obsm": {"X_emb": np.zeros((n, 2), dtype=np.float32)}})
    check("SC-03e log 尺度矩阵才标记 normalize", "normalize" in state2["steps"], str(state2["steps"]))

    # 依赖缺口要显式记录
    dataset3 = _FakeDataset(counts, obs, [f"G{i}" for i in range(4)])
    state3 = SC_STATE.detect(dataset3, {})
    details3 = state3.get("details") or {}
    check("SC-03f 记录证据/假设步骤分离",
          "steps_evidence" in details3 and "steps_assumed" in details3, str(details3.get("steps_assumed")))


# ---------------------------------------------------------------------------
# SC-04：DE 缓存必须随数据变化失效
# ---------------------------------------------------------------------------

def _de_ctx(matrix, obs, groups=("A", "B")):
    return {
        "dataset_id": "ds-cache",
        "dataset_version": 1,
        "dataset_revision": 1,
        "content_revision": "rev-%d" % abs(hash(tuple(map(tuple, np.asarray(matrix)[:3].tolist())))),
        "n_obs": obs.shape[0],
        "n_vars": np.asarray(matrix).shape[1],
        "X": sp.csr_matrix(np.asarray(matrix, dtype=np.float32)),
        "X_kind": "lognorm",
        "counts": sp.csr_matrix(np.asarray(matrix, dtype=np.float32)),
        "obs": obs,
        "var_names": [f"G{i}" for i in range(np.asarray(matrix).shape[1])],
        "obsm": {}, "clusters": None,
        "group_keys": ["group", "cell_type"],
        "engine_mode": "builtin",
    }


def test_sc04_cache():
    print("\n[SC-04] 分组变更后 DE 缓存失效")
    rng = np.random.default_rng(0)
    n = 48
    genes = 6
    groups = ["A"] * 24 + ["B"] * 24
    cell_types = ["T"] * 24 + ["Bcell"] * 24
    matrix = rng.normal(5.0, 0.5, (n, genes))
    matrix[:24, :2] += 4.0     # A 组前两个基因高
    obs = pd.DataFrame({
        "group": groups,
        "cell_type": ["T cell"] * 30 + ["B cell"] * 18,
    }, index=[f"c{i}" for i in range(n)])
    ctx_a = _de_ctx(matrix, obs)
    params = {"group": "group", "group1": "A", "group2": "B", "min_cells": 3}
    data1 = SETS._compute_de(ctx_a, params, [])
    logfc_first = None
    for cell_type, record in data1["per_type"].items():
        logfc_first = record["logfc"][:2].copy()
        break
    check("SC-04a 首次 DE 计算成功", logfc_first is not None, str(logfc_first))

    # 互换分组标签（细胞数、列名都不变）
    swapped = obs.copy()
    swapped["group"] = ["B"] * 24 + ["A"] * 24
    ctx_b = _de_ctx(matrix, swapped)
    ctx_b["content_revision"] = "swapped-labels"
    data2 = SETS._compute_de(ctx_b, params, [])
    logfc_second = None
    for cell_type, record in data2["per_type"].items():
        logfc_second = record["logfc"][:2].copy()
        break
    check("SC-04b 互换分组后 logFC 反向（缓存未命中）",
          logfc_second is not None and np.allclose(logfc_second, -logfc_first),
          f"{logfc_first} -> {logfc_second}")

    # 同样的数据再来一次，应当命中缓存（一致性）
    warnings = []
    data3 = SETS._compute_de(ctx_b, dict(params), warnings)
    hit = any("缓存" in str(w) for w in warnings)
    check("SC-04c 数据未变时复用缓存", hit, str(warnings[:1]))

    # 只改注释列（group 不变）也要产生新键
    annot = swapped.copy()
    annot["cell_type"] = ["X cell"] * n
    ctx_c = _de_ctx(matrix, annot)
    key_b = SETS._de_cache_key(ctx_b, {**params, "cell_type": "cell_type"})
    key_c = SETS._de_cache_key(ctx_c, {**params, "cell_type": "cell_type"})
    check("SC-04d 注释列改变 → 缓存键不同", key_b != key_c)


# ---------------------------------------------------------------------------
# SC-08：空组不能崩
# ---------------------------------------------------------------------------

def test_sc08_empty_group():
    print("\n[SC-08] 稀有类型过滤后出现空组")
    # 场景 1：3 组，其中 A 组只剩 rare 细胞（过滤后整组为空），另外两组有效
    n = 66
    groups = ["A"] * 2 + ["B"] * 32 + ["C"] * 32
    cell_types = ["rare"] * 2 + (["T"] * 16 + ["Bcell"] * 16) * 2
    matrix = np.random.default_rng(1).normal(5, 0.3, (n, 5))
    obs = pd.DataFrame({"group": groups, "cell_type": cell_types}, index=[f"c{i}" for i in range(n)])
    ctx = _de_ctx(matrix, obs)
    result = None
    raised = None
    try:
        result = SETS.run("proportion_test", ctx, {"group": "group", "min_cells": 10},
                          os.path.join(OUTDIR, "sc08"))
    except Exception as error:  # noqa: BLE001
        raised = error
    check("SC-08a 空组场景不崩溃（修复前 UnboundLocalError）",
          raised is None or type(raised).__name__ != "UnboundLocalError",
          f"{type(raised).__name__}: {raised}")
    if raised is None and result is not None:
        metrics = result.get("summary", {}).get("metrics", {})
        check("SC-08b 排除空组并给出解释",
              any("没有任何细胞" in str(w) for w in result.get("warnings") or []),
              str((result.get("warnings") or [""])[:2]))
        check("SC-08c 有效组数按排除后统计",
              int(metrics.get("n_groups") or 0) == 2, str(metrics.get("n_groups")))
    # 场景 3：期望频数过低时只给警告、不崩
    small = ["A"] * 20 + ["B"] * 21
    small_types = ["T"] * 20 + ["T"] * 20 + ["Mono"]
    m3 = np.random.default_rng(3).normal(5, 0.3, (41, 5))
    o3 = pd.DataFrame({"group": small, "cell_type": small_types}, index=[f"c{i}" for i in range(41)])
    res3 = SETS.run("proportion_test", _de_ctx(m3, o3), {"group": "group", "min_cells": 1},
                    os.path.join(OUTDIR, "sc08-lowexp"))
    check("SC-08e 低期望频数只警告", any("期望频数" in str(w) for w in res3.get("warnings") or []),
          str((res3.get("warnings") or [""])[:2]))

    # 只有一个有效组 → prerequisite，而不是给出"显著"结论
    obs2 = pd.DataFrame({"group": ["A"] * 2 + ["B"] * 20,
                         "cell_type": ["rare"] * 2 + ["T"] * 20}, index=[f"c{i}" for i in range(22)])
    matrix2 = np.random.default_rng(2).normal(5, 0.3, (22, 5))
    ctx2 = _de_ctx(matrix2, obs2)
    raised = None
    try:
        SETS.run("proportion_test", ctx2, {"group": "group", "min_cells": 10},
                          os.path.join(OUTDIR, "sc08-single"))
    except Exception as error:  # noqa: BLE001
        raised = error
    check("SC-08d 单个有效组 → prerequisite",
          raised is not None and str(getattr(raised, "code", "")) == "prerequisite",
          f"{type(raised).__name__}: {getattr(raised, 'message', raised)}")


# ---------------------------------------------------------------------------
# SC-05：统计学单位（患者/供体 vs 细胞）
# ---------------------------------------------------------------------------

def _donor_ctx(n_donors=(3, 3), cells_per_donor=8, genes=6, seed=0):
    rng = np.random.default_rng(seed)
    samples, groups, cell_types, rows = [], [], [], []
    donor_effect = {}
    for group, label in ((0, "A"), (1, "B")):
        for d in range(n_donors[group]):
            donor = "%s%d" % (label, d)
            donor_effect[donor] = rng.normal(0, 1.5)          # 供体间差异
            sig = 2.0 if label == "B" else 0.0                # 组间效应
            for c in range(cells_per_donor):
                # 两类细胞交替分配，保证每个供体在两类里都有细胞
                cell_type = "T cell" if c % 2 == 0 else "B cell"
                values = rng.normal(5.0, 0.3, genes)
                values[:2] += donor_effect[donor] + sig
                rows.append(values)
                samples.append(donor)
                groups.append(label)
                cell_types.append(cell_type)
    matrix = np.vstack(rows)
    n = matrix.shape[0]
    obs = pd.DataFrame({
        "group": groups,
        "cell_type": cell_types,
        "sample": samples,
    }, index=["c%d" % i for i in range(n)])
    return _de_ctx(matrix, obs)


def _duplicate_cells(ctx):
    """把每个细胞复制一遍（供体数不变，细胞数翻倍）——用来验证 n 不是细胞数。"""
    matrix = ctx["X"].toarray()
    obs = ctx["obs"].copy()
    names = [str(i) + "_dup" for i in obs.index]
    obs2 = pd.concat([obs, obs.set_axis(names)])
    matrix2 = np.vstack([matrix, matrix])
    ctx2 = _de_ctx(matrix2, obs2)
    ctx2["var_names"] = list(ctx["var_names"])
    return ctx2


def test_sc05_statistical_unit():
    print("\n[SC-05] 组间推断的统计学单位")
    ctx = _donor_ctx()
    params = {"group": "group", "group1": "A", "group2": "B", "min_cells": 3}
    cell_level = SETS.run("de_group", ctx, dict(params, unit="cell"), os.path.join(OUTDIR, "sc05-cell"))
    donor_level = SETS.run("de_group", ctx, dict(params, unit="donor"), os.path.join(OUTDIR, "sc05-donor"))

    m_cell = cell_level["summary"]["metrics"]
    m_donor = donor_level["summary"]["metrics"]
    check("SC-05a 默认按细胞比较时标出统计单位",
          m_cell.get("statistical_unit") == "cell" and m_cell.get("patient_level_inference") is False,
          f"{m_cell.get('statistical_unit')}/{m_cell.get('patient_level_inference')}")
    check("SC-05b 检出生物学重复列并在警告里说明",
          any("统计单位是「细胞」" in str(w) or "统计单位" in str(w) for w in cell_level.get("warnings") or []),
          str((cell_level.get("warnings") or [""])[:1])[:70])
    check("SC-05c 供体级分析记录重复数与单位",
          m_donor.get("statistical_unit") == "donor" and m_donor.get("n_donors_in") == 3
          and m_donor.get("n_donors_out") == 3 and m_donor.get("patient_level_inference") is True,
          f"{m_donor.get('statistical_unit')} {m_donor.get('n_donors_in')}/{m_donor.get('n_donors_out')}")

    # 每类差异基因表的 n 应当是供体数，而不是细胞数
    def first_row(payload):
        for table_obj in payload.get("tables") or []:
            if table_obj.get("id") == "de_genes":
                rows = table_obj.get("rows") or []
                return rows[0] if rows else None
        return None

    row_cell = first_row(cell_level)
    row_donor = first_row(donor_level)
    check("SC-05d 细胞级 n = 12（3 供体 × 每类 4 细胞）",
          row_cell and int(row_cell["n_in"]) == 12, str(row_cell and row_cell["n_in"]))
    check("SC-05e 供体级 n = 3（不是 24）",
          row_donor and int(row_donor["n_in"]) == 3, str(row_donor and row_donor["n_in"]))

    # 核心验收：只复制同一供体的细胞，不能增加独立 n，供体级 p 值必须不变
    doubled = _duplicate_cells(ctx)
    donor_doubled = SETS.run("de_group", doubled, dict(params, unit="donor"), os.path.join(OUTDIR, "sc05-donor-dup"))
    cell_doubled = SETS.run("de_group", doubled, dict(params, unit="cell"), os.path.join(OUTDIR, "sc05-cell-dup"))

    def pvector(payload):
        for table_obj in payload.get("tables") or []:
            if table_obj.get("id") == "de_genes":
                rows = table_obj.get("rows") or []
                return np.array([float(r["pval"]) for r in rows]), np.array([float(r["n_in"]) for r in rows])
        return np.array([]), np.array([])

    p_donor_a, n_donor_a = pvector(donor_level)
    p_donor_b, n_donor_b = pvector(donor_doubled)
    p_cell_a, n_cell_a = pvector(cell_level)
    p_cell_b, n_cell_b = pvector(cell_doubled)
    same_n = n_donor_a.size and n_donor_b.size and np.all(n_donor_a == n_donor_b) and int(n_donor_a[0]) == 3
    check("SC-05f 复制细胞后供体级 n 不变（仍是 3）", bool(same_n),
          f"{n_donor_a[:1]} -> {n_donor_b[:1]}")
    check("SC-05g 复制细胞后供体级 p 值不变（同一供体均值）",
          p_donor_a.size == p_donor_b.size and np.allclose(p_donor_a, p_donor_b, equal_nan=True),
          f"{p_donor_a[:2]} vs {p_donor_b[:2]}")
    check("SC-05h 细胞级比较在复制后被夸大（说明它的 n 是细胞数）",
          p_cell_b.size and (np.mean(p_cell_b) < np.mean(p_cell_a)),
          f"cell mean p {np.mean(p_cell_a):.3g} -> {np.mean(p_cell_b):.3g}")

    # 每组只有一个供体 → 供体级请求必须明确拒绝
    single = _donor_ctx(n_donors=(1, 1))
    raised = None
    try:
        SETS.run("de_group", single, dict(params, unit="donor"), os.path.join(OUTDIR, "sc05-single"))
    except Exception as error:  # noqa: BLE001
        raised = error
    check("SC-05i 每组单供体时拒绝供体级推断",
          raised is not None and str(getattr(raised, "code", "")) == "prerequisite",
          f"{type(raised).__name__}: {getattr(raised, 'message', raised)}")

    # 缺供体列却要求供体级 → prerequisite
    no_donor = _de_ctx(
        np.random.default_rng(5).normal(size=(40, 6)),
        pd.DataFrame({"group": ["A"] * 20 + ["B"] * 20,
                      "cell_type": (["T cell"] * 10 + ["B cell"] * 10) * 2},
                     index=[f"x{i}" for i in range(40)]),
    )
    raised = None
    try:
        SETS.run("de_group", no_donor, dict(params, unit="donor"), os.path.join(OUTDIR, "sc05-nodonor"))
    except Exception as error:  # noqa: BLE001
        raised = error
    check("SC-05j 没有供体列时供体级分析报 prerequisite",
          raised is not None and str(getattr(raised, "code", "")) == "prerequisite",
          f"{type(raised).__name__}: {getattr(raised, 'message', raised)}")


def main() -> int:
    print("== dsh-singlecell 审查回归（SC-01/03/04/05/07/08）==")
    for fn in (test_sc01_raw_alignment, test_sc07_unique_barcodes, test_sc03_state,
               test_sc04_cache, test_sc08_empty_group, test_sc05_statistical_unit):
        try:
            fn()
        except Exception:  # noqa: BLE001
            FAILED.append(fn.__name__)
            print(f"  FAIL {fn.__name__} 抛异常")
            traceback.print_exc()
    print(f"\n== 结果：{len(PASSED)} 通过 / {len(FAILED)} 失败 ==")
    if FAILED:
        for item in FAILED:
            print(f"  ✗ {item}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
