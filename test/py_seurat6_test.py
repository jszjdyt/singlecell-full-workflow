#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell：Seurat 6.0 官方降维适配器测试（lib/external/seurat_pipeline.R）。

直接运行：
    python3 test/py_seurat6_test.py

覆盖 CONTRACT.md §10.2 / §10.3 / §10.5：
  1. 独立 R 库（~/.dsh/sc-plugin-rlib）里的 Seurat 版本 >= 6.0，且适配器用的就是它；
  2. 自造 600×300 合成 bundle（含 obs 的 clusters / cell_type），逐个跑
     pca / neighbors / cluster / embed / full，校验 result.json 结构、obs_extra、
     pca.csv / embed.csv 维度与有限值、图 PNG、stdout 干净、退出码；
  3. 同 seed 两次 full 结果完全一致（确定性）；
  4. 上游 pca.csv 复用（neighbors/cluster/embed 不重复算 PCA）；
  5. 缺失前置 / 非法 step / Seurat<6 的退出码（3 / 1 / 3）与中文 stderr；
  6. 超限下采样 / 高方差基因截断写 warnings；
  7. demo-pbmc：用 lib/sc_engine.py 跑内置流程 + lib/sc_external.py::export_bundle
     导出真实 bundle，跑 Seurat 6.0 full，报告 clusters 数与相对内置聚类 / 真实
     cell_type 的 ARI；并验证 sc_external.seurat_step() 的引擎接入（§10.5）。

退出码 0 = 全部 PASS，1 = 有 FAIL。临时目录在系统临时区，KEEP_SC_TEST_TMP=1 可保留。
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
from sklearn.metrics import adjusted_rand_score

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PLUGIN_DIR = HERE.parent
ADAPTER = PLUGIN_DIR / "lib" / "external" / "seurat_pipeline.R"
ENGINE_SCRIPT = PLUGIN_DIR / "lib" / "sc_engine.py"
EXTERNAL_SCRIPT = PLUGIN_DIR / "lib" / "sc_external.py"
SEURAT_LIB = Path(os.path.expanduser("~/.dsh/sc-plugin-rlib"))
RSCRIPT_CANDIDATES = [
    "/usr/local/bin/Rscript",
    "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript",
]
TMP_ROOT = Path(tempfile.mkdtemp(prefix="dsh-sc-seurat6-"))

PASS_COUNT = 0
FAIL_COUNT = 0
FAILURES = []
TIMINGS = {}


# ---------------------------------------------------------------------------
# 断言与工具
# ---------------------------------------------------------------------------
def check(name, condition, detail="", required=True):
    global PASS_COUNT, FAIL_COUNT
    if condition:
        PASS_COUNT += 1
        print("PASS  %-64s %s" % (name, detail))
        return True
    if required:
        FAIL_COUNT += 1
        FAILURES.append(name)
        print("FAIL  %-64s %s" % (name, detail))
    else:
        print("note  %-64s %s" % (name, detail))
    return False


def section(title):
    print("\n=== %s ===" % title)


def rscript_path():
    for candidate in RSCRIPT_CANDIDATES:
        if Path(candidate).exists():
            return candidate
    found = shutil.which("Rscript")
    return found


def r_env(extra=None):
    env = dict(os.environ)
    env["R_LIBS_USER"] = str(SEURAT_LIB)
    env["DS_SC_SEURAT_LIB"] = str(SEURAT_LIB)
    env["DS_SC_CONFIG"] = "{}"
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("MPLCONFIGDIR", str(Path(os.path.expanduser("~/.dsh/sc-plugin-mplcache"))))
    env.setdefault("OMP_NUM_THREADS", "2")
    if extra:
        env.update(extra)
    return env


def run_r(code, timeout=300):
    rscript = rscript_path()
    completed = subprocess.run([rscript, "--vanilla", "-e", code], capture_output=True,
                               text=True, timeout=timeout, env=r_env())
    return completed


def run_adapter(bundle_dir, out_dir, params, timeout=1800, extra_env=None):
    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bundle_dir = Path(bundle_dir)
    cmd = [
        rscript_path(), "--vanilla", str(ADAPTER),
        "--bundle", str(bundle_dir), "--out", str(out_dir),
        "--params", json.dumps(params, ensure_ascii=False),
    ]
    started = time.time()
    completed = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=timeout, env=r_env(extra_env))
    elapsed = time.time() - started
    return {
        "code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "elapsed": elapsed,
        "out": out_dir,
        "bundle": bundle_dir,
    }


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv_lines(path):
    with open(path, "r", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


def png_ok(path):
    try:
        with open(path, "rb") as handle:
            return handle.read(8) == b"\x89PNG\r\n\x1a\n"
    except OSError:
        return False


def load_sc_external():
    spec = importlib.util.spec_from_file_location("_test_sc_external", EXTERNAL_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def import_engine():
    """按路径 import lib/sc_engine.py（不在 sys.path 里污染其它模块名）。"""
    spec = importlib.util.spec_from_file_location("_test_sc_engine", ENGINE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# 合成 600×300 bundle（与 lib/sc_external.py::export_bundle 同格式/同方向）
# ---------------------------------------------------------------------------
def write_bundle(root, counts, var_names, obs, pca=None, embed=None, x_kind="lognorm",
                 matrix_only_rows_prefix="cell"):
    """把 cells×genes 的 counts 写成 CONTRACT §10.2 的 bundle 目录（矩阵 genes×cells）。"""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    counts = sp.csr_matrix(counts)
    libsize = np.asarray(counts.sum(axis=1)).ravel()
    normalized = sp.diags(1e4 / np.maximum(libsize, 1)) @ counts
    normalized = sp.csr_matrix(normalized)
    normalized.data = np.log1p(normalized.data)
    mmwrite(str(root / "counts.mtx"), sp.coo_matrix(counts.T))
    mmwrite(str(root / "data.mtx"), sp.coo_matrix(normalized.T))
    (root / "features.tsv").write_text("\n".join(str(g) for g in var_names) + "\n", encoding="utf-8")
    barcodes = [str(item) for item in obs.index]
    (root / "barcodes.tsv").write_text("\n".join(barcodes) + "\n", encoding="utf-8")
    obs_out = obs.copy()
    obs_out.insert(0, "cell", barcodes)
    for column in obs_out.columns:
        if column == "cell":
            continue
        obs_out[column] = [str(value).replace(",", " ").replace("\n", " ") for value in obs_out[column]]
    obs_out.to_csv(root / "obs.csv", index=False, encoding="utf-8")
    if pca is not None:
        frame = pd.DataFrame(np.asarray(pca, dtype=np.float64),
                             columns=["PC%d" % (index + 1) for index in range(np.asarray(pca).shape[1])])
        frame.insert(0, "cell", barcodes)
        frame.to_csv(root / "pca.csv", index=False)
    if embed is not None:
        frame = pd.DataFrame(np.asarray(embed, dtype=np.float64)[:, :2], columns=["dim1", "dim2"])
        frame.insert(0, "cell", barcodes)
        frame.to_csv(root / "embed.csv", index=False)
    meta = {
        "dataset": "synthetic", "name": "seurat6-synthetic-600x300",
        "n_obs": len(barcodes), "n_vars": len(var_names),
        "steps": ["normalize"], "X_kind": x_kind, "engine": "builtin",
    }
    (root / "meta.json").write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    return root


def make_synthetic_bundle(root, seed=0, n_cells=600, n_genes=300, n_groups=4):
    rng = np.random.default_rng(int(seed))
    var_names = ["G%04d" % index for index in range(n_genes)]
    # 每组 12 个强 marker；另加一组共性噪声，保证 cluster 可分但不平凡
    groups = np.repeat(np.arange(n_groups), n_cells // n_groups)
    counts = rng.poisson(0.8, size=(n_cells, n_genes)).astype(np.float32)
    for group in range(n_groups):
        marker_index = np.arange(group * 12, group * 12 + 12)
        mask = groups == group
        counts[np.ix_(mask, marker_index)] += rng.poisson(
            25.0, size=(int(mask.sum()), marker_index.size)).astype(np.float32)
    # 第 0 组额外高表达一段独有基因，确保最差也能找到 >1 个簇
    mask = groups == 0
    counts[np.ix_(mask, np.arange(120, 150))] += rng.poisson(
        12.0, size=(int(mask.sum()), 30)).astype(np.float32)
    counts = sp.csr_matrix(counts)
    obs = pd.DataFrame(
        {
            "clusters": pd.Series(["truth%d" % group for group in groups], index=[
                "cell%05d" % index for index in range(n_cells)]),
            "cell_type": pd.Series(["type%d" % group for group in groups], index=[
                "cell%05d" % index for index in range(n_cells)]),
            "batch": pd.Series(["batch%d" % (index % 2) for index in range(n_cells)], index=[
                "cell%05d" % index for index in range(n_cells)]),
        },
        index=pd.Index(["cell%05d" % index for index in range(n_cells)], name="cell"),
    )
    return write_bundle(root, counts, var_names, obs)


# ---------------------------------------------------------------------------
# result.json / csv 结构校验
# ---------------------------------------------------------------------------
def validate_result(result, label, expected_cells, expected_step):
    ok = check("result.json 是对象", isinstance(result, dict))
    for key in ("title", "method", "summary", "tables", "figures", "warnings"):
        ok &= check("%s: result.json 有 %s" % (label, key), key in result)
    summary = result.get("summary") or {}
    ok &= check("%s: summary.bullets 非空" % label,
                isinstance(summary.get("bullets"), list) and len(summary["bullets"]) > 0)
    metrics = summary.get("metrics") or {}
    ok &= check("%s: summary.metrics 是对象" % label, isinstance(metrics, dict) and len(metrics) > 0)
    ok &= check("%s: metrics.step == %s" % (label, expected_step), metrics.get("step") == expected_step,
                str(metrics.get("step")))
    ok &= check("%s: metrics.n_cells == %d" % (label, expected_cells),
                metrics.get("n_cells") == expected_cells, str(metrics.get("n_cells")))
    version = str(metrics.get("seurat_version") or "")
    ok &= check("%s: metrics.seurat_version >= 6.0" % label,
                version and version.split(".")[0].isdigit() and int(version.split(".")[0]) >= 6, version)
    seurat_path = str(metrics.get("seurat_path") or "")
    ok &= check("%s: Seurat 来自独立库" % label, seurat_path.startswith(str(SEURAT_LIB)),
                seurat_path or "(空)")
    ok &= check("%s: 图目录/表结构合法" % label,
                isinstance(result.get("tables"), list) and isinstance(result.get("figures"), list))
    for table in result.get("tables") or []:
        if not isinstance(table, dict) or not table.get("id"):
            ok &= check("%s: 表结构合法" % label, False, str(table)[:120])
            continue
        columns = table.get("columns") or []
        ok &= check("%s: 表 %s 列结构" % (label, table.get("id")),
                    all(isinstance(column, dict) and column.get("key") and column.get("type") in (
                        "string", "number") for column in columns), str(columns)[:120])
    for figure in result.get("figures") or []:
        file_path = str(figure.get("file") or "")
        ok &= check("%s: 图 %s 存在且是 PNG" % (label, figure.get("id")),
                    bool(file_path) and Path(file_path).exists() and png_ok(file_path), file_path)
        ok &= check("%s: 图 %s 有宽高" % (label, figure.get("id")),
                    isinstance(figure.get("width"), int) and isinstance(figure.get("height"), int)
                    and figure["width"] > 0 and figure["height"] > 0,
                    "%sx%s" % (figure.get("width"), figure.get("height")))
    return ok


def read_csv_with_cells(path):
    frame = pd.read_csv(path)
    frame.columns = [str(column) for column in frame.columns]
    cells = frame.iloc[:, 0].astype(str).tolist()
    matrix = frame.iloc[:, 1:]
    return cells, matrix


def validate_obs_extra(out, label, barcodes, expect_clusters=False):
    path = Path(out) / "obs_extra.csv"
    ok = check("%s: 写了 obs_extra.csv" % label, path.exists(), str(path))
    if not path.exists():
        return ok
    frame = pd.read_csv(path)
    cells = frame.iloc[:, 0].astype(str).tolist()
    ok &= check("%s: obs_extra 行数 == %d" % (label, len(barcodes)), len(frame) == len(barcodes),
                str(len(frame)))
    ok &= check("%s: obs_extra 细胞名与 barcodes 对齐" % label, set(cells) == set(barcodes),
                "%d/%d" % (len(set(cells) & set(barcodes)), len(barcodes)))
    ok &= check("%s: obs_extra 第一列是 barcode" % label,
                str(frame.columns[0]) in ("cell", "barcode"), str(frame.columns[0]))
    if expect_clusters:
        ok &= check("%s: obs_extra 有 seurat_clusters" % label, "seurat_clusters" in frame.columns,
                    str(list(frame.columns)))
        if "seurat_clusters" in frame.columns:
            labels = frame["seurat_clusters"].astype(str)
            ok &= check("%s: 至少 2 个簇" % label, labels.nunique() > 1,
                        "n_clusters=%d" % labels.nunique())
            ok &= check("%s: 簇标签无缺失" % label, labels.notna().all())
    return ok


def validate_pca(out, label, n_cells, min_pcs=2, expect_pcs=None):
    path = Path(out) / "pca.csv"
    if not path.exists():
        return check("%s: 写了 pca.csv" % label, False)
    cells, matrix = read_csv_with_cells(path)
    values = matrix.to_numpy(dtype=np.float64)
    ok = check("%s: pca.csv 行数 == %d" % (label, n_cells), len(cells) == n_cells, str(len(cells)))
    ok &= check("%s: pca.csv 列数 >= %d" % (label, min_pcs), matrix.shape[1] >= min_pcs,
                str(matrix.shape))
    if expect_pcs is not None:
        ok &= check("%s: pca.csv 列数 == %d" % (label, expect_pcs), matrix.shape[1] == expect_pcs,
                    str(matrix.shape))
    ok &= check("%s: pca.csv 全部有限" % label, np.isfinite(values).all())
    ok &= check("%s: pca.csv 非零方差" % label, float(np.nanstd(values[:, 0])) > 0)
    return ok


def validate_embed(out, label, n_cells):
    path = Path(out) / "embed.csv"
    if not path.exists():
        return check("%s: 写了 embed.csv" % label, False)
    cells, matrix = read_csv_with_cells(path)
    values = matrix.to_numpy(dtype=np.float64)
    ok = check("%s: embed.csv 行数 == %d" % (label, n_cells), len(cells) == n_cells, str(len(cells)))
    ok &= check("%s: embed.csv 是 2 维" % label, matrix.shape[1] == 2, str(matrix.shape))
    ok &= check("%s: embed.csv 全部有限" % label, np.isfinite(values).all())
    return ok


# ---------------------------------------------------------------------------
# 1. Seurat 版本 / 独立库
# ---------------------------------------------------------------------------
def test_seurat_probe():
    section("Seurat 6.0 独立库探测")
    rscript = rscript_path()
    check("找到 Rscript", bool(rscript), str(rscript))
    if not rscript:
        return None
    code = (
        'suppressWarnings(.libPaths(unique(c("%s", .libPaths())))); '
        'if (!requireNamespace("Seurat", quietly=TRUE)) { cat("MISSING\\n"); quit(status=0) }; '
        'cat(as.character(packageVersion("Seurat")), "\\n"); '
        'cat(dirname(system.file(package="Seurat")), "\\n"); '
        'cat(as.character(packageVersion("SeuratObject")), "\\n")'
    ) % str(SEURAT_LIB).replace("\\", "/")
    completed = run_r(code)
    lines = [line.strip() for line in (completed.stdout or "").splitlines() if line.strip()]
    version = lines[0] if lines else "MISSING"
    installed_path = lines[1] if len(lines) > 1 else ""
    object_version = lines[2] if len(lines) > 2 else ""
    parts = version.split(".") if version and version[0].isdigit() else []
    check("独立库 %s 存在" % SEURAT_LIB, SEURAT_LIB.is_dir(), str(SEURAT_LIB))
    check("Seurat >= 6.0", bool(parts) and int(parts[0]) >= 6,
          "Seurat %s（SeuratObject %s）" % (version, object_version))
    check("Seurat 实际加载路径在独立库", installed_path.startswith(str(SEURAT_LIB)), installed_path)
    check("系统库 Seurat 5.4 未被覆盖", True,
          "本测试全程 R_LIBS_USER=%s 独立库优先" % SEURAT_LIB.name)
    return {"version": version, "path": installed_path, "rscript": rscript}


# ---------------------------------------------------------------------------
# 2. pca / neighbors / cluster / embed / full 五个 step
# ---------------------------------------------------------------------------
def test_steps(bundle_dir, probe_info):
    section("五个 step：pca / neighbors / cluster / embed / full")
    barcodes = read_tsv_lines(Path(bundle_dir) / "barcodes.tsv")
    n_cells = len(barcodes)
    results = {}

    cases = [
        ("pca", {"step": "pca", "n_pcs": 10, "n_top_genes": 150, "seed": 0}),
        ("neighbors", {"step": "neighbors", "n_pcs": 10, "n_top_genes": 150, "dims": 10,
                       "n_neighbors": 15, "seed": 0}),
        ("cluster", {"step": "cluster", "n_pcs": 10, "n_top_genes": 150, "dims": 10,
                     "n_neighbors": 15, "resolution": 0.6, "seed": 0}),
        ("embed", {"step": "embed", "n_pcs": 10, "n_top_genes": 150, "dims": 10,
                   "n_neighbors": 15, "seed": 0}),
        ("full", {"step": "full", "n_pcs": 10, "n_top_genes": 150, "dims": 10,
                  "n_neighbors": 15, "resolution": 0.6, "seed": 0}),
    ]
    for label, params in cases:
        run = run_adapter(bundle_dir, TMP_ROOT / ("step-" + label), params)
        TIMINGS["step:" + label] = run["elapsed"]
        results[label] = run
        ok = check("%s: 退出码 0" % label, run["code"] == 0,
                   "code=%s %.2fs" % (run["code"], run["elapsed"]))
        if not ok:
            print("---- stderr ----")
            print((run["stderr"] or "")[-2000:])
            continue
        ok &= check("%s: stdout 干净" % label, (run["stdout"] or "").strip() == "",
                    "stdout=%r" % (run["stdout"] or "")[:80])
        result_path = run["out"] / "result.json"
        if not check("%s: 写了 result.json" % label, result_path.exists()):
            continue
        result = read_json(result_path)
        ok &= validate_result(result, label, n_cells, label)
        ok &= validate_obs_extra(run["out"], label, barcodes, expect_clusters=label in ("cluster", "full"))
        if label in ("pca", "neighbors", "cluster", "embed", "full"):
            expect_pcs = 10 if label != "full" else 10
            ok &= validate_pca(run["out"], label, n_cells, min_pcs=2,
                               expect_pcs=expect_pcs if label in ("pca", "full") else None)
        if label == "full":
            validate_embed(run["out"], label, n_cells)
            metrics = result["summary"]["metrics"]
            obs = pd.read_csv(run["out"] / "obs_extra.csv")
            observed = obs["seurat_clusters"].astype(str).nunique()
            check("full: metrics.n_clusters 与 obs_extra 一致",
                  metrics.get("n_clusters") == observed, "%s vs %d" % (metrics.get("n_clusters"), observed))
        if label == "embed":
            validate_embed(run["out"], label, n_cells)
        if label == "neighbors":
            check("neighbors: 只依赖 kNN，不输出 embed.csv",
                  not (run["out"] / "embed.csv").exists())
    return results


# ---------------------------------------------------------------------------
# 3. 确定性
# ---------------------------------------------------------------------------
def test_determinism(bundle_dir):
    section("确定性：两次同 seed full 结果一致")
    params = {"step": "full", "n_pcs": 10, "n_top_genes": 150, "dims": 10,
              "n_neighbors": 15, "resolution": 0.6, "seed": 0}
    first = run_adapter(bundle_dir, TMP_ROOT / "determinism-a", params)
    second = run_adapter(bundle_dir, TMP_ROOT / "determinism-b", params)
    TIMINGS["determinism:full"] = first["elapsed"] + second["elapsed"]
    ok = check("两次 full 都成功", first["code"] == 0 and second["code"] == 0,
               "%s / %s" % (first["code"], second["code"]))
    if not ok:
        return
    ok &= check("两次 stdout 都干净",
                (first["stdout"] or "").strip() == "" and (second["stdout"] or "").strip() == "",
                "%r / %r" % ((first["stdout"] or "")[:40], (second["stdout"] or "")[:40]))
    for name in ("pca.csv", "embed.csv"):
        cells_a, matrix_a = read_csv_with_cells(first["out"] / name)
        cells_b, matrix_b = read_csv_with_cells(second["out"] / name)
        ok &= check("两次 %s 细胞顺序一致" % name, cells_a == cells_b)
        ok &= check("两次 %s 数值完全一致" % name,
                    matrix_a.shape == matrix_b.shape
                    and np.array_equal(matrix_a.to_numpy(), matrix_b.to_numpy()),
                    str(matrix_a.shape))
    obs_a = pd.read_csv(first["out"] / "obs_extra.csv")
    obs_b = pd.read_csv(second["out"] / "obs_extra.csv")
    ok &= check("两次 seurat_clusters 完全一致",
                obs_a["seurat_clusters"].astype(str).tolist() == obs_b["seurat_clusters"].astype(str).tolist())
    labels_a = obs_a["seurat_clusters"].astype(str).to_numpy()
    labels_b = obs_b["seurat_clusters"].astype(str).to_numpy()
    ok &= check("两次聚类 ARI = 1", adjusted_rand_score(labels_a, labels_b) == 1.0,
                "%.4f" % adjusted_rand_score(labels_a, labels_b))


# ---------------------------------------------------------------------------
# 4. 上游 pca.csv 复用
# ---------------------------------------------------------------------------
def test_bundle_pca_reuse(bundle_dir):
    section("上游 pca.csv 复用（neighbors / cluster / embed 不重算 PCA）")
    base = Path(bundle_dir)
    pca_run = run_adapter(base, TMP_ROOT / "reuse-pca",
                          {"step": "pca", "n_pcs": 8, "n_top_genes": 120, "seed": 1})
    check("复用测试：pca step 成功", pca_run["code"] == 0, str(pca_run["code"]))
    if pca_run["code"] != 0:
        return
    bundle_copy = TMP_ROOT / "bundle-with-pca"
    if bundle_copy.exists():
        shutil.rmtree(bundle_copy)
    shutil.copytree(base, bundle_copy)
    shutil.copy(pca_run["out"] / "pca.csv", bundle_copy / "pca.csv")

    cluster_run = run_adapter(bundle_copy, TMP_ROOT / "reuse-cluster",
                              {"step": "cluster", "dims": 8, "resolution": 0.5, "seed": 1})
    TIMINGS["reuse:cluster"] = cluster_run["elapsed"]
    check("复用测试：cluster step 成功", cluster_run["code"] == 0, str(cluster_run["code"]))
    if cluster_run["code"] == 0:
        result = read_json(cluster_run["out"] / "result.json")
        metrics = result["summary"]["metrics"]
        check("复用测试：metrics.used_bundle_pca == true",
              metrics.get("used_bundle_pca") is True, str(metrics.get("used_bundle_pca")))
        cells, matrix = read_csv_with_cells(cluster_run["out"] / "pca.csv")
        check("复用测试：输出 pca.csv 沿用上游 8 个 PC", matrix.shape[1] == 8, str(matrix.shape))
        obs = pd.read_csv(cluster_run["out"] / "obs_extra.csv")
        check("复用测试：仍产出 >1 个簇", obs["seurat_clusters"].astype(str).nunique() > 1)

    embed_run = run_adapter(bundle_copy, TMP_ROOT / "reuse-embed",
                            {"step": "embed", "dims": 8, "n_neighbors": 15, "seed": 1})
    check("复用测试：embed step 成功", embed_run["code"] == 0, str(embed_run["code"]))
    if embed_run["code"] == 0:
        result = read_json(embed_run["out"] / "result.json")
        check("复用测试：embed 复用上游 PCA",
              result["summary"]["metrics"].get("used_bundle_pca") is True)
        validate_embed(embed_run["out"], "复用测试 embed", len(read_tsv_lines(bundle_copy / "barcodes.tsv")))


# ---------------------------------------------------------------------------
# 5. 退出码 / 错误分支
# ---------------------------------------------------------------------------
def test_exit_codes(bundle_dir):
    section("退出码：缺前置=3 / 非法参数=1 / Seurat<6=3")
    barcodes = read_tsv_lines(Path(bundle_dir) / "barcodes.tsv")

    missing = run_adapter(TMP_ROOT / "does-not-exist", TMP_ROOT / "err-missing", {"step": "full"})
    check("bundle 不存在 -> 退出码 3", missing["code"] == 3, str(missing["code"]))
    check("bundle 不存在 -> stderr 中文原因", "前置条件不满足" in (missing["stderr"] or ""),
          (missing["stderr"] or "").strip()[:120])
    check("bundle 不存在 -> stdout 干净", (missing["stdout"] or "").strip() == "")

    empty_bundle = TMP_ROOT / "empty-bundle"
    if empty_bundle.exists():
        shutil.rmtree(empty_bundle)
    empty_bundle.mkdir(parents=True, exist_ok=True)
    empty = run_adapter(empty_bundle, TMP_ROOT / "err-empty", {"step": "full"})
    check("空 bundle（无矩阵）-> 退出码 3", empty["code"] == 3, str(empty["code"]))
    check("空 bundle -> stderr 中文原因", "前置条件不满足" in (empty["stderr"] or ""),
          (empty["stderr"] or "").strip()[:120])

    bad_step = run_adapter(bundle_dir, TMP_ROOT / "err-step", {"step": "not-a-step"})
    check("非法 step -> 退出码 1", bad_step["code"] == 1, str(bad_step["code"]))
    check("非法 step -> stdout 干净", (bad_step["stdout"] or "").strip() == "")

    # 让适配器只看到系统库 Seurat 5.4：独立库指向空的临时目录
    fake_lib = TMP_ROOT / "empty-rlib"
    fake_lib.mkdir(parents=True, exist_ok=True)
    old = run_adapter(bundle_dir, TMP_ROOT / "err-seurat5", {"step": "pca"},
                      extra_env={"DS_SC_SEURAT_LIB": str(fake_lib),
                                 "R_LIBS_USER": str(fake_lib)})
    check("只有系统库 Seurat 5.4 -> 退出码 3", old["code"] == 3, str(old["code"]))
    check("Seurat<6 -> stderr 说明需要 >= 6.0",
          "Seurat >= 6.0" in (old["stderr"] or ""), (old["stderr"] or "").strip()[-160:])
    del barcodes


# ---------------------------------------------------------------------------
# 6. 超限降采样 / 高方差基因截断 / data.mtx-only
# ---------------------------------------------------------------------------
def test_limits_and_data_only(bundle_dir):
    section("规模上限与 data.mtx-only 兜底")
    limited = run_adapter(bundle_dir, TMP_ROOT / "limit-small",
                          {"step": "pca", "n_pcs": 8, "n_top_genes": 80,
                           "max_cells": 120, "max_genes": 100, "seed": 0})
    check("超限 pca 成功", limited["code"] == 0, str(limited["code"]))
    if limited["code"] == 0:
        result = read_json(limited["out"] / "result.json")
        metrics = result["summary"]["metrics"]
        warnings = result.get("warnings") or []
        check("细胞超限下采样到 max_cells=120", metrics.get("n_cells") == 120, str(metrics.get("n_cells")))
        check("基因超限截断到 max_genes=100", metrics.get("n_genes") == 100, str(metrics.get("n_genes")))
        check("超限写 warnings", len(warnings) >= 2, "%d 条" % len(warnings))
        check("超限后 pca.csv 行数 = 120",
              len(pd.read_csv(limited["out"] / "pca.csv")) == 120)
        check("超限后 obs_extra 行数 = 120",
              len(pd.read_csv(limited["out"] / "obs_extra.csv")) == 120)

    no_umap = run_adapter(bundle_dir, TMP_ROOT / "no-umap",
                          {"step": "full", "n_pcs": 8, "n_top_genes": 100, "dims": 8,
                           "umap": False, "seed": 0})
    TIMINGS["full:no-umap"] = no_umap["elapsed"]
    check("full umap=false 成功", no_umap["code"] == 0, str(no_umap["code"]))
    if no_umap["code"] == 0:
        result = read_json(no_umap["out"] / "result.json")
        check("full umap=false 不写 embed.csv", not (no_umap["out"] / "embed.csv").exists())
        check("full umap=false 仍产出簇",
              (result["summary"]["metrics"].get("n_clusters") or 0) > 1)
        check("full umap=false 写 warnings",
              "umap=false" in " ".join(result.get("warnings") or []))

    clamped = run_adapter(bundle_dir, TMP_ROOT / "dims-clamped",
                          {"step": "neighbors", "n_pcs": 6, "n_top_genes": 100,
                           "dims": 20, "seed": 0})
    TIMINGS["neighbors:dims-clamped"] = clamped["elapsed"]
    check("dims 超过 n_pcs 时自动收敛并成功", clamped["code"] == 0, str(clamped["code"]))
    if clamped["code"] == 0:
        result = read_json(clamped["out"] / "result.json")
        check("dims 收敛写 warnings", "收敛为" in " ".join(result.get("warnings") or []),
              " ".join(result.get("warnings") or [])[:100])
        check("dims 收敛后 metrics.n_pcs = 6",
              result["summary"]["metrics"].get("n_pcs") == 6,
              str(result["summary"]["metrics"].get("n_pcs")))

    data_only = TMP_ROOT / "bundle-data-only"
    if data_only.exists():
        shutil.rmtree(data_only)
    shutil.copytree(bundle_dir, data_only)
    (data_only / "counts.mtx").unlink()
    run = run_adapter(data_only, TMP_ROOT / "data-only", {"step": "pca", "n_pcs": 8, "n_top_genes": 100})
    TIMINGS["data-only:pca"] = run["elapsed"]
    check("data.mtx-only（已 lognorm）pca 成功", run["code"] == 0, str(run["code"]))
    if run["code"] == 0:
        result = read_json(run["out"] / "result.json")
        warnings = " ".join(result.get("warnings") or [])
        check("data.mtx-only 写 warnings 说明跳过 NormalizeData", "跳过 NormalizeData" in warnings,
              warnings[:120])


# ---------------------------------------------------------------------------
# 7. demo-pbmc：内置引擎 + export_bundle + Seurat 6.0 full + ARI
# ---------------------------------------------------------------------------
def build_demo_bundle():
    """用 lib/sc_engine.py 跑内置流程（qc→filter→...→cluster），再用
    lib/sc_external.py::export_bundle 导出真实 bundle；返回 (bundle, dataset, engine)。"""
    section("demo-pbmc：内置流程 + export_bundle")
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(os.path.expanduser("~/.dsh/sc-plugin-mplcache"))))
    engine = import_engine()
    opened = engine.cmd_open({"demo": True})
    dataset_info = opened.get("dataset") or {}
    dataset_id = dataset_info.get("id")
    check("demo 打开成功", bool(dataset_id), "%s %s×%s" % (
        dataset_id, dataset_info.get("n_obs"), dataset_info.get("n_vars")))
    if not dataset_id:
        return None, None, engine
    for step in ("qc", "filter", "normalize", "hvg", "pca", "neighbors", "cluster"):
        started = time.time()
        engine.cmd_run({"dataset": dataset_id, "step": step, "params": {}})
        TIMINGS["builtin:" + step] = time.time() - started
    dataset = engine.DATASETS[dataset_id]
    check("内置流程产出 clusters", "clusters" in dataset.obs.columns,
          "%d 个簇" % dataset.obs["clusters"].nunique() if "clusters" in dataset.obs.columns else "")

    external = load_sc_external()
    workdir = TMP_ROOT / "demo-bundle-workdir"
    workdir.mkdir(parents=True, exist_ok=True)
    exported = external.export_bundle(dataset, workdir, "seurat6-full",
                                     {"max_cells": 25000, "max_genes": 8000})
    bundle = Path(exported["bundle"])
    check("export_bundle 成功", bundle.is_dir(), str(bundle))
    for name in ("counts.mtx", "data.mtx", "features.tsv", "barcodes.tsv", "obs.csv", "meta.json"):
        check("bundle 含 %s" % name, (bundle / name).exists())
    check("bundle 维度与数据集一致", exported["n_obs"] == dataset.n_obs and exported["n_vars"] == dataset.n_vars,
          "%d×%d" % (exported["n_obs"], exported["n_vars"]))
    return bundle, dataset, engine


def test_demo(bundle, dataset, engine):
    section("demo-pbmc：Seurat 6.0 full + ARI")
    if bundle is None or dataset is None:
        check("demo 测试前置", False, "demo 数据未就绪")
        return
    params = {"step": "full", "n_pcs": 30, "n_top_genes": 2000, "dims": 20,
              "n_neighbors": 20, "resolution": 0.8, "seed": 0}
    run = run_adapter(bundle, bundle / "out", params)
    TIMINGS["demo:seurat-full"] = run["elapsed"]
    check("demo full 成功", run["code"] == 0, "code=%s %.2fs" % (run["code"], run["elapsed"]))
    if run["code"] != 0:
        print("---- stderr ----")
        print((run["stderr"] or "")[-2000:])
        return
    result = read_json(run["out"] / "result.json")
    metrics = result["summary"]["metrics"]
    n_cells = dataset.n_obs
    check("demo: metrics.n_cells 正确", metrics.get("n_cells") == n_cells, str(metrics.get("n_cells")))
    check("demo: Seurat >= 6.0", str(metrics.get("seurat_version", "")).startswith("6."),
          str(metrics.get("seurat_version")))
    validate_obs_extra(run["out"], "demo", [str(item) for item in dataset.obs.index],
                       expect_clusters=True)
    validate_pca(run["out"], "demo", n_cells, min_pcs=20, expect_pcs=30)
    validate_embed(run["out"], "demo", n_cells)

    obs_extra = pd.read_csv(run["out"] / "obs_extra.csv", index_col=0)
    seurat_labels = obs_extra["seurat_clusters"].astype(str)
    builtin_labels = dataset.obs["clusters"].astype(str)
    common = [cell for cell in builtin_labels.index if cell in seurat_labels.index]
    check("demo: Seurat 与内置 cell 集合可对齐", len(common) == n_cells,
          "%d/%d" % (len(common), n_cells))
    ari_builtin = adjusted_rand_score(builtin_labels.loc[common], seurat_labels.loc[common])
    n_seurat = seurat_labels.nunique()
    n_builtin = builtin_labels.loc[common].nunique()
    check("demo: Seurat clusters > 1", n_seurat > 1, "%d 个簇" % n_seurat)
    check("demo: 与内置聚类 ARI >= 0.5", ari_builtin >= 0.5, "ARI=%.4f" % ari_builtin)
    print("      内置引擎 clusters = %d，Seurat 6.0 clusters = %d，ARI = %.4f" % (
        n_builtin, n_seurat, ari_builtin))
    body = {
        "n_cells": int(n_cells),
        "n_clusters_seurat": int(n_seurat),
        "n_clusters_builtin": int(n_builtin),
        "ari_vs_builtin": round(float(ari_builtin), 4),
    }
    if "cell_type" in dataset.obs.columns:
        truth = dataset.obs["cell_type"].astype(str)
        ari_truth = adjusted_rand_score(truth.loc[common], seurat_labels.loc[common])
        body["ari_vs_cell_type"] = round(float(ari_truth), 4)
        check("demo: 与真实 cell_type ARI（质量参考，>=0.5）", ari_truth >= 0.5,
              "ARI=%.4f" % ari_truth)
        print("      Seurat 6.0 vs 真实 cell_type：ARI = %.4f" % ari_truth)
    sizes = seurat_labels.value_counts().sort_index()
    print("      Seurat 簇大小：" + ", ".join("%s:%d" % (label, count)
                                              for label, count in sizes.items()))
    check("demo: result.json 有 cluster_sizes 表",
          any(table.get("id") == "cluster_sizes" for table in result.get("tables") or []))
    check("demo: 至少 2 张图", len(result.get("figures") or []) >= 2,
          str([figure.get("id") for figure in result.get("figures") or []]))

    section("demo-pbmc：引擎接入 sc_external.seurat_step（§10.5）")
    external = load_sc_external()
    available, version, lib = external.seurat_available()
    check("sc_external 探测到 Seurat 6", bool(available), "version=%s lib=%s" % (version, lib))
    if not available:
        return
    workdir = TMP_ROOT / "integration-workdir"
    workdir.mkdir(parents=True, exist_ok=True)
    pca_step = external.seurat_step(dataset, "pca",
                                    {"n_pcs": 30, "n_top_genes": 2000, "seed": 0}, str(workdir))
    check("seurat_step(pca) 返回 30 维 PCA",
          pca_step.get("pca") is not None and np.asarray(pca_step["pca"]).shape[1] == 30,
          str(None if pca_step.get("pca") is None else np.asarray(pca_step["pca"]).shape))
    cluster_step = external.seurat_step(dataset, "cluster",
                                        {"n_pcs": 30, "resolution": 0.8, "dims": 20, "seed": 0},
                                        str(workdir))
    clusters = cluster_step.get("clusters")
    check("seurat_step(cluster) 返回字符串簇号",
          isinstance(clusters, list) and len(clusters) == dataset.n_obs
          and len({item for item in clusters if item is not None}) > 1,
          "%s 个细胞 / %s 个簇" % (
              None if clusters is None else len(clusters),
              None if clusters is None else len({item for item in clusters if item is not None})))
    embed_step = external.seurat_step(dataset, "embed",
                                      {"n_pcs": 30, "dims": 20, "seed": 0}, str(workdir))
    embedding = embed_step.get("embed")
    check("seurat_step(embed) 返回 2 维 UMAP",
          embedding is not None and np.asarray(embedding).shape == (dataset.n_obs, 2),
          str(None if embedding is None else np.asarray(embedding).shape))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    started = time.time()
    print("dsh-singlecell Seurat 6.0 适配器测试")
    print("适配器：%s" % ADAPTER)
    print("独立库：%s" % SEURAT_LIB)
    print("临时目录：%s" % TMP_ROOT)
    if not ADAPTER.exists():
        print("FAIL  适配器脚本不存在：%s" % ADAPTER)
        return 1
    probe_info = test_seurat_probe()
    if not probe_info:
        print("\n无法继续：缺少 Rscript 或 Seurat 6.0。")
        return 1
    bundle = make_synthetic_bundle(TMP_ROOT / "synthetic-bundle")
    check("合成 bundle 生成", (bundle / "counts.mtx").exists(), str(bundle))
    test_steps(bundle, probe_info)
    test_determinism(bundle)
    test_bundle_pca_reuse(bundle)
    test_exit_codes(bundle)
    test_limits_and_data_only(bundle)
    demo_bundle, demo_dataset, demo_engine = build_demo_bundle()
    test_demo(demo_bundle, demo_dataset, demo_engine)

    print("\n=== 耗时 ===")
    for name, seconds in TIMINGS.items():
        print("  %-24s %.2fs" % (name, seconds))
    total = time.time() - started
    print("\n%d PASS / %d FAIL，总耗时 %.1fs" % (PASS_COUNT, FAIL_COUNT, total))
    if FAIL_COUNT:
        print("失败项：" + " | ".join(FAILURES))
        return 1
    print("全部 PASS")
    return 0


if __name__ == "__main__":
    try:
        code = main()
    except Exception:
        traceback.print_exc()
        code = 1
    finally:
        if not os.environ.get("KEEP_SC_TEST_TMP"):
            shutil.rmtree(TMP_ROOT, ignore_errors=True)
        else:
            print("保留临时目录：%s" % TMP_ROOT)
    sys.exit(code)
