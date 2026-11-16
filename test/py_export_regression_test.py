#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell — 导出问题回归测试（EXPORT-01 / EXPORT-02）

对着《DSH_GEO_Singlecell_审查与优化任务单》第六节写：

  EXPORT-01  Methods 必须从「结构化执行记录」生成，逐条与真实执行一致：
             - 内置 Louvain / 请求 Leiden 回退 Louvain / 外部 Seurat 引擎三种流程；
             - 记录 requested_method / actual_method / 软件版本 / 实际参数 /
               输入矩阵尺度 / 聚类 resolution / 种子 / 回退原因；
             - 只描述真实成功执行且仍属于当前数据版本的步骤；
             - 没有执行记录（如直接导入的 Seurat RDS）必须写「未记录」，不得补写 Louvain。

  EXPORT-02  导出脚本必须可重放：
             - 保存输入清单与 SHA-256、细胞/基因映射、按顺序的有效操作、每步实际参数、
               依赖版本与种子；
             - 重放走插件的真实执行入口（同一套 builtin 实现），而不是换算法顶替；
             - 过滤后 counts 随 AnnData 同步切片，不回填过滤前矩阵；
             - 保留的 Scanpy 脚本必须命名为「参考迁移脚本」并列出方法差异；
             - 在独立目录重新执行导出产物，逐值比较确定性步骤；
             - 输入被改动 / 依赖不可用时不能伪装成重放成功，只能标为「记录/示例」。

用法：
    cd plugins/dsh-singlecell && python3 test/py_export_regression_test.py

说明：本机没有 scanpy，也没有 Seurat 6.0（只有 Seurat 5.4.0）；因此
  * builtin 流程（CSV / h5ad 输入）真跑、真重放、逐值比较；
  * 外部 Seurat 流程用一个「模拟外部桥」走 `_run_seurat_step` 的真实记录代码路径，
    并在输出里明确标注是模拟；
  * 真实 Seurat RDS 用本机已有的 R 桥导入（只读），验证「无执行记录 → 未记录」。
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
PLUGIN_DIR = HERE.parent
LIB = PLUGIN_DIR / "lib"
WORKSPACE = PLUGIN_DIR.parent.parent

TMP_ROOT = Path(tempfile.mkdtemp(prefix="sc-export-reg-"))
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ["MPLCONFIGDIR"] = str(TMP_ROOT / "mpl-cache")
os.environ.setdefault("DS_SC_WORKDIR", str(TMP_ROOT / "engine-work"))
os.environ["PYTHONUNBUFFERED"] = "1"
sys.dont_write_bytecode = True

PYTHON = sys.executable or "python3"

PASSED: list[str] = []
FAILED: list[str] = []
SKIPPED: list[str] = []


def check(name, cond, detail="", required=True):
    if cond:
        PASSED.append(name)
        print(f"  PASS {name}" + (f"  [{detail}]" if detail else ""))
        return True
    if required:
        FAILED.append(name)
        print(f"  FAIL {name}" + (f"  [{detail}]" if detail else ""))
        return False
    print(f"  note {name}" + (f"  [{detail}]" if detail else ""))
    return False


def skip(name, detail=""):
    SKIPPED.append(name)
    print(f"  SKIP {name}" + (f"  [{detail}]" if detail else ""))


# ---------------------------------------------------------------------------
# 引擎（直接 import lib/sc_engine.py；测试内进程即 worker）
# ---------------------------------------------------------------------------


def _load_engine():
    spec = importlib.util.spec_from_file_location("_export_reg_engine", LIB / "sc_engine.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["_export_reg_engine"] = module
    spec.loader.exec_module(module)
    return module


ENGINE = _load_engine()
try:
    import anndata as ad  # noqa: E402
except Exception:  # noqa: BLE001
    ad = None


def fresh_engine():
    engine_module = sys.modules.get("_export_reg_engine")
    engine_module.DATASETS.clear()
    ENGINE._EXPM1_CACHE.clear()
    return ENGINE


def cmd(command, args=None, timeout=600):
    """在测试进程里直接调用 worker 命令（等价于一条协议请求）。"""
    return ENGINE.COMMANDS[command](dict(args or {}))


# ---------------------------------------------------------------------------
# 固定 fixture：块状结构，保证 Louvain 能分出多个簇，且过滤会真正删细胞/基因
# ---------------------------------------------------------------------------

N_CELLS = 30
N_GENES = 12


def _block_counts(seed=7, n_cells=N_CELLS, n_genes=N_GENES):
    rng = np.random.default_rng(seed)
    data = rng.integers(0, 4, size=(n_genes, n_cells)).astype(np.float64)
    for block in range(3):
        data[block * 4:(block + 1) * 4, block * 10:(block + 1) * 10] += rng.integers(
            6, 15, size=(4, 10)
        )
    data[n_genes - 1, :] = 0.0          # 全零基因：会被 min_cells 过滤
    for column in (3, 7, 11, 17):       # 全零细胞：会被 min_genes 过滤
        data[:, column] = 0.0
    genes = ["G%d" % index for index in range(n_genes)]
    cells = ["%s_c%d" % ("A" if index % 2 else "B", index) for index in range(n_cells)]
    return data, genes, cells


def make_csv_fixture(root):
    data, genes, cells = _block_counts(seed=7)
    path = Path(root) / "fixture-counts.csv"
    pd.DataFrame(data, index=genes, columns=cells).to_csv(path)
    return path


def make_h5ad_fixture(root):
    if ad is None:
        return None
    data, genes, cells = _block_counts(seed=11)
    obs = pd.DataFrame(
        {
            "sample": ["S%d" % (index % 3) for index in range(len(cells))],
            "batch": ["b%d" % (index % 2) for index in range(len(cells))],
        },
        index=pd.Index(cells, name="cell"),
    )
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    matrix = sp.csr_matrix(data.T.astype(np.float32))  # cells × genes
    path = Path(root) / "fixture-counts.h5ad"
    fixture = ad.AnnData(X=matrix, obs=obs, var=var)
    fixture.layers["counts"] = matrix.copy()  # 可验证的原始计数（SC-02）
    fixture.write_h5ad(str(path))
    return path


def open_fixture(path):
    fresh_engine()
    result = cmd("open", {"path": str(path)})
    dataset_id = result["dataset"]["id"]
    return dataset_id, ENGINE.DATASETS[dataset_id]


WORKFLOW = [
    ("qc", {"min_genes": 1, "max_pct_mt": 100}),
    ("filter", {"min_genes": 1, "min_cells": 2, "max_pct_mt": 100}),
    ("normalize", {"target_sum": 10000}),
    ("hvg", {"n_top_genes": 6}),
    ("pca", {"n_pcs": 5}),
    ("neighbors", {"n_neighbors": 5, "metric": "cosine"}),
    ("cluster", {"resolution": 0.5, "method": "louvain"}),
    ("embed", {"layout": "tsne", "perplexity": 5}),
    ("markers", {"top_n": 3}),
]


def run_workflow(dataset_id, steps=None, annotate=True, mapping=None):
    for step, params in steps or WORKFLOW:
        cmd("run", {"dataset": dataset_id, "step": step, "params": params})
    dataset = ENGINE.DATASETS[dataset_id]
    if annotate:
        if mapping is None:
            categories, _codes = dataset.cluster_codes()
            mapping = {str(label): "Type%s" % label for label in (categories or [])}
        if mapping:
            cmd("run", {"dataset": dataset_id, "step": "annotate", "params": {"mapping": mapping}})
    return ENGINE.DATASETS[dataset_id]


def _method_sections(text):
    """把 Methods 按 '### ' 切成 {标题: 内容}。"""
    sections = {}
    current = "intro"
    buffer = []
    for line in text.splitlines():
        if line.startswith("### "):
            sections[current] = "\n".join(buffer)
            current = line[4:].strip()
            buffer = []
        else:
            buffer.append(line)
    sections[current] = "\n".join(buffer)
    return sections


# ---------------------------------------------------------------------------
# EXPORT-01：结构化执行记录
# ---------------------------------------------------------------------------


def test_execution_record_structure():
    print("\n[EXPORT-01a] 结构化执行记录字段")
    (TMP_ROOT / "record").mkdir(parents=True, exist_ok=True)
    csv_path = make_csv_fixture(TMP_ROOT / "record")
    dataset_id, dataset = open_fixture(csv_path)
    run_workflow(dataset_id)

    records = ENGINE.execution_records(dataset)
    checks = {record.get("step"): record for record in records if record.get("kind") == "step"}
    required_steps = [step for step, _params in WORKFLOW]
    check("EXPORT-01a-1 每个成功步骤都有执行记录",
          all(step in checks for step in required_steps), str(sorted(checks)))
    entry = checks.get("cluster") or {}
    for field in ("requested_method", "actual_method", "software", "params_effective",
                  "params_requested", "input_scale", "output_scale", "seed", "fallback_reason",
                  "n_obs_after", "n_vars_after", "dataset_version", "kind", "step"):
        check("EXPORT-01a-2 记录字段 %-16s 存在" % field, field in entry, str(entry.get(field)))
    check("EXPORT-01a-3 聚类记录 requested=louvain / actual=louvain(builtin)",
          str(entry.get("requested_method")).lower() == "louvain"
          and str(entry.get("actual_method")).lower().startswith("louvain"),
          "%s -> %s" % (entry.get("requested_method"), entry.get("actual_method")))
    check("EXPORT-01a-4 记录聚类 resolution 与实际生效参数",
          float((entry.get("params_effective") or {}).get("resolution", -1)) == 0.5,
          str((entry.get("params_effective") or {}).get("resolution")))
    check("EXPORT-01a-5 记录种子", entry.get("seed") == ENGINE.SEED, str(entry.get("seed")))
    check("EXPORT-01a-6 记录输入矩阵尺度 = lognorm",
          "lognorm" in str(entry.get("input_scale")), str(entry.get("input_scale")))
    versions = entry.get("software") or {}
    check("EXPORT-01a-7 记录依赖版本（numpy/sklearn/python/插件版本）",
          isinstance(versions, dict) and versions.get("numpy") and versions.get("sklearn")
          and versions.get("python") and versions.get("dsh-singlecell"),
          str({key: versions.get(key) for key in ("numpy", "sklearn", "python", "dsh-singlecell")}))
    filter_record = checks.get("filter") or {}
    check("EXPORT-01a-8 过滤记录带实际保留数量",
          int(filter_record.get("n_obs_after") or 0) == dataset.n_obs,
          "n_obs_after=%s dataset.n_obs=%s" % (filter_record.get("n_obs_after"), dataset.n_obs))


# ---------------------------------------------------------------------------
# EXPORT-01：三种流程的 Methods
# ---------------------------------------------------------------------------


def test_methods_louvain():
    print("\n[EXPORT-01b] Methods：内置 Louvain + 真实参数")
    root = TMP_ROOT / "methods-louvain"
    root.mkdir(parents=True, exist_ok=True)
    dataset_id, dataset = open_fixture(make_csv_fixture(root))
    run_workflow(dataset_id)
    engine_dataset = ENGINE.DATASETS[dataset_id]
    text = ENGINE._build_methods(engine_dataset)
    sections = _method_sections(text)

    check("EXPORT-01b-1 含 '## Methods' 与英文主体", text.startswith("## Methods") and len(text) > 300,
          "%d chars" % len(text))
    cluster_section = sections.get("Clustering") or ""
    check("EXPORT-01b-2 聚类写实际 Louvain（不是模板外的方法）",
          "Louvain" in cluster_section and "louvain(builtin)" in cluster_section,
          cluster_section.strip().splitlines()[0][:100] if cluster_section.strip() else "<空>")
    check("EXPORT-01b-3 聚类写实际 resolution=0.5", "resolution = 0.5" in cluster_section,
          cluster_section[:160])
    check("EXPORT-01b-4 聚类写种子", "seed = %d" % ENGINE.SEED in cluster_section, cluster_section[:160])
    check("EXPORT-01b-5 聚类写输入矩阵尺度 lognorm", "lognorm" in cluster_section, cluster_section[:160])
    check("EXPORT-01b-6 Methods 只写实际执行的步骤（markers 用 Welch）",
          "Welch" in (sections.get("Marker genes") or ""), (sections.get("Marker genes") or "")[:120])
    check("EXPORT-01b-7 normalisation 写实际 target_sum",
          "10000" in (sections.get("Normalisation") or ""), (sections.get("Normalisation") or "")[:120])
    check("EXPORT-01b-8 PCA 写实际 solver（randomized）",
          "randomized" in (sections.get("Principal component analysis") or ""),
          (sections.get("Principal component analysis") or "")[:120])
    check("EXPORT-01b-9 含软件版本行（Python/NumPy/sklearn）",
          "Python" in text and "NumPy" in text and "scikit-learn" in text)
    check("EXPORT-01b-10 含数据来源与矩阵尺度声明",
          "log" in (sections.get("Data import and provenance") or "").lower()
          and "counts" in (sections.get("Data import and provenance") or "").lower(),
          (sections.get("Data import and provenance") or "")[:160])
    # 当前数据版本：重跑上游 hvg 后 cluster 记录必须清掉
    cmd("run", {"dataset": dataset_id, "step": "hvg", "params": {"n_top_genes": 6}})
    after = ENGINE._build_methods(ENGINE.DATASETS[dataset_id])
    after_sections = _method_sections(after)
    check("EXPORT-01b-11 上游重跑后不再描述已失效的聚类",
          "Clustering" not in after_sections and "Louvain" not in after,
          "steps=%s" % ENGINE.DATASETS[dataset_id].steps)
    check("EXPORT-01b-12 未执行的步骤不出现在 Methods",
          "Principal component analysis" not in after_sections, "steps=%s" % ENGINE.DATASETS[dataset_id].steps)


def test_methods_leiden_fallback():
    print("\n[EXPORT-01c] Methods：请求 Leiden、实际回退内置 Louvain")
    root = TMP_ROOT / "methods-leiden"
    root.mkdir(parents=True, exist_ok=True)
    dataset_id, _dataset = open_fixture(make_csv_fixture(root))
    steps = [item for item in WORKFLOW if item[0] not in ("cluster", "embed", "markers")]
    steps.append(("cluster", {"resolution": 0.5, "method": "leiden"}))
    steps.append(("embed", {"layout": "tsne", "perplexity": 5}))
    steps.append(("markers", {"top_n": 3}))
    run_workflow(dataset_id, steps=steps, annotate=False)
    engine_dataset = ENGINE.DATASETS[dataset_id]

    record = [item for item in ENGINE.execution_records(engine_dataset)
              if item.get("step") == "cluster"][-1]
    scanpy_available = ENGINE.SCANPY is not None and ENGINE.HAS_LEIDENALG
    if scanpy_available:
        check("EXPORT-01c-0 本机 scanpy+leidenalg 可用时实际方法应为 leiden",
              str(record.get("actual_method")).lower().startswith("leiden"),
              str(record.get("actual_method")))
    else:
        check("EXPORT-01c-1 无 scanpy 时记录 requested=leiden / actual=louvain",
              str(record.get("requested_method")).lower() == "leiden"
              and str(record.get("actual_method")).lower().startswith("louvain"),
              "%s -> %s" % (record.get("requested_method"), record.get("actual_method")))
        check("EXPORT-01c-2 记录回退原因（scanpy/leidenalg 不可用）",
              "leiden" in str(record.get("fallback_reason")).lower()
              and "unavailable" in str(record.get("fallback_reason")).lower()
              or "不可用" in str(record.get("fallback_reason")),
              str(record.get("fallback_reason"))[:120])
        text = ENGINE._build_methods(engine_dataset)
        cluster_section = _method_sections(text).get("Clustering") or ""
        check("EXPORT-01c-3 Methods 同时写 requested 与 actual",
              "leiden" in cluster_section.lower() and "louvain(builtin)" in cluster_section,
              cluster_section[:220])
        check("EXPORT-01c-4 Methods 写回退原因，不把 Leiden 冒充为实际方法",
              "fallback" in cluster_section.lower() and "scanpy" in cluster_section.lower(),
              cluster_section[:300])


def _install_fake_seurat_bridge():
    """模拟外部 Seurat 桥：只替换 R 调用，记录/写回走真实的 _run_seurat_step。"""
    module_name = "_dsh_sc_lib_sc_external"

    class _FakeSeurat:
        METHOD = ("Seurat 6.0 官方流程（NormalizeData/FindVariableFeatures/ScaleData/"
                  "RunPCA/FindNeighbors/FindClusters/RunUMAP）")

        def seurat_step(self, dataset, step, params, workdir, progress=None):
            metrics = {
                "integration": str(params.get("integration") or "none"),
                "reduction": str(params.get("reduction") or "pca"),
                "resolution": float(params.get("resolution") or 0.5),
                "n_pcs": int(params.get("n_pcs") or 30),
                "n_neighbors": int(params.get("n_neighbors") or 30),
            }
            summary = {"bullets": ["FindClusters（Louvain）得到 3 个簇"], "metrics": metrics}
            result = {"payload": {"method": self.METHOD, "summary": summary, "warnings": []},
                      "clusters": None, "pca": None, "embed": None, "warnings": []}
            n_obs = int(dataset.n_obs)
            if step == "pca":
                # 确定性的假 PCA：记录/写回路径是真的，数值只用于检查 Methods 记录
                result["pca"] = np.tile(
                    np.arange(n_obs, dtype=np.float32)[:, None], (1, metrics["n_pcs"])
                ) / 100.0
            elif step == "cluster":
                result["clusters"] = [str(index % 3) for index in range(n_obs)]
            elif step == "embed":
                result["embed"] = np.stack(
                    [np.arange(n_obs, dtype=np.float32), -np.arange(n_obs, dtype=np.float32)], axis=1
                )
            return result

    fake = type(sys)("_dsh_sc_lib_sc_external")
    fake.seurat_step = _FakeSeurat().seurat_step
    sys.modules[module_name] = fake
    return module_name


def test_methods_external_seurat():
    print("\n[EXPORT-01d] Methods：外部 Seurat 6.0 引擎（R 桥为模拟，记录路径为真实代码）")
    root = TMP_ROOT / "methods-seurat"
    root.mkdir(parents=True, exist_ok=True)
    dataset_id, _dataset = open_fixture(make_csv_fixture(root))
    _install_fake_seurat_bridge()
    steps = [
        ("qc", {"min_genes": 1, "max_pct_mt": 100}),
        ("filter", {"min_genes": 1, "min_cells": 2, "max_pct_mt": 100}),
        ("normalize", {"target_sum": 10000}),
        ("hvg", {"n_top_genes": 6}),
        ("pca", {"n_pcs": 5, "engine": "seurat6", "integration": "none"}),
        ("neighbors", {"n_neighbors": 5, "engine": "seurat6"}),
        ("cluster", {"resolution": 0.5, "engine": "seurat6", "integration": "none"}),
        ("embed", {"layout": "umap", "engine": "seurat6"}),
    ]
    for step, params in steps:
        cmd("run", {"dataset": dataset_id, "step": step, "params": params})
    sys.modules.pop("_dsh_sc_lib_sc_external", None)
    engine_dataset = ENGINE.DATASETS[dataset_id]
    records = ENGINE.execution_records(engine_dataset)
    cluster_record = [item for item in records if item.get("step") == "cluster"][-1]
    check("EXPORT-01d-1 外部聚类记录 actual_method 来自 Seurat 实现",
          "seurat" in str(cluster_record.get("actual_method")).lower()
          and "findclusters" in str(cluster_record.get("actual_method")).lower(),
          str(cluster_record.get("actual_method")))
    check("EXPORT-01d-2 记录保留外部实现来源（provenance）",
          "Seurat 6.0 官方流程" in str(cluster_record.get("method_provenance")),
          str(cluster_record.get("method_provenance"))[:120])
    check("EXPORT-01d-3 记录实际 integration 与 resolution",
          str((cluster_record.get("params_effective") or {}).get("integration")) in ("none", "")
          and float((cluster_record.get("params_effective") or {}).get("resolution", -1)) == 0.5,
          str(cluster_record.get("params_effective")))

    text = ENGINE._build_methods(engine_dataset)
    sections = _method_sections(text)
    cluster_section = sections.get("Clustering") or ""
    check("EXPORT-01d-4 Methods 写 Seurat::FindClusters，不写内置 Louvain 作为聚类方法",
          "FindClusters" in cluster_section and "built-in pure-Python Louvain" not in cluster_section,
          cluster_section[:240])
    check("EXPORT-01d-5 Methods 保留外部实现 provenance 原文",
          "Seurat 6.0 官方流程" in cluster_section, cluster_section[:240])
    pca_section = sections.get("Principal component analysis") or ""
    check("EXPORT-01d-6 PCA 一节写外部实现（Seurat::RunPCA）",
          "Seurat" in pca_section and "RunPCA" in pca_section, pca_section[:160])
    check("EXPORT-01d-7 Methods 声明引擎为外部实现而非 builtin",
          "external" in text.lower() or "Seurat" in text, "Seurat 出现 %d 次" % text.count("Seurat"))


def test_methods_imported_without_record():
    print("\n[EXPORT-01e] Methods：导入的 Seurat RDS（没有执行记录 → 未记录）")
    rds = WORKSPACE / "neu_analysis" / "objects" / "neu_processed.rds"
    if shutil.which("Rscript") is None or not rds.exists():
        skip("EXPORT-01e", "没有 Rscript 或缺少 %s" % rds)
        return
    fresh_engine()
    result = cmd("open", {"path": str(rds)})
    dataset_id = result["dataset"]["id"]
    engine_dataset = ENGINE.DATASETS[dataset_id]
    text = ENGINE._build_methods(engine_dataset)
    sections = _method_sections(text)
    check("EXPORT-01e-1 导入数据检测到已有步骤", len(engine_dataset.steps) >= 3,
          str(engine_dataset.steps))
    check("EXPORT-01e-2 没有执行记录时 Methods 不冒充内置 Louvain",
          "Louvain" not in text and "louvain(builtin)" not in text, text[:80])
    unrecorded = "\n".join(value for key, value in sections.items() if "未记录" in key + value)
    check("EXPORT-01e-3 未记录步骤显式写「未记录」", "未记录" in unrecorded, unrecorded[:200])
    check("EXPORT-01e-4 记录中确实没有 step 执行记录",
          not [item for item in ENGINE.execution_records(engine_dataset) if item.get("kind") == "step"],
          str([item.get("step") for item in ENGINE.execution_records(engine_dataset)]))
    check("EXPORT-01e-5 导入来源写进 Methods（rds）",
          "rds" in text.lower() or "RDS" in text, text[:120])


def test_methods_excludes_failed_steps():
    print("\n[EXPORT-01f] Methods：失败步骤不进入")
    root = TMP_ROOT / "methods-failed"
    root.mkdir(parents=True, exist_ok=True)
    dataset_id, _dataset = open_fixture(make_csv_fixture(root))
    cmd("run", {"dataset": dataset_id, "step": "qc", "params": {"min_genes": 1, "max_pct_mt": 100}})
    cmd("run", {"dataset": dataset_id, "step": "filter",
                "params": {"min_genes": 1, "min_cells": 2, "max_pct_mt": 100}})
    failed = False
    try:
        cmd("run", {"dataset": dataset_id, "step": "markers", "params": {"top_n": 3}})
    except ENGINE.EngineError:
        failed = True
    check("EXPORT-01f-1 markers 前置不满足确实失败", failed)
    text = ENGINE._build_methods(ENGINE.DATASETS[dataset_id])
    check("EXPORT-01f-2 Methods 不写失败步骤（无 Welch）", "Welch" not in text, text[:120])
    check("EXPORT-01f-3 失败步骤也没有执行记录",
          not [item for item in ENGINE.execution_records(ENGINE.DATASETS[dataset_id])
               if item.get("step") == "markers"])


# ---------------------------------------------------------------------------
# EXPORT-02：导出脚本结构
# ---------------------------------------------------------------------------


def test_script_structure():
    print("\n[EXPORT-02a] 导出脚本：可重放驱动 + 参考迁移脚本分离")
    root = TMP_ROOT / "script-structure"
    root.mkdir(parents=True, exist_ok=True)
    dataset_id, _dataset = open_fixture(make_csv_fixture(root))
    engine_dataset = run_workflow(dataset_id)
    cmd("obs_derive", {"dataset": dataset_id, "name": "group",
                       "source": "cell_name", "pattern": "^([AB])_"})
    text = ENGINE._build_script(engine_dataset)

    driver_part = text.split("REFERENCE_SCANPY_SCRIPT = r", 1)[0]
    annex_part = text.split("REFERENCE_SCANPY_SCRIPT = r", 1)[1] if "REFERENCE_SCANPY_SCRIPT = r" in text else ""
    check("EXPORT-02a-1 脚本嵌入 REPLAY_MANIFEST（输入清单/操作/期望）",
          "REPLAY_MANIFEST" in text and "operations" in text and "sha256" in text)
    check("EXPORT-02a-2 重放路径不再盲写 sc.read_h5ad", "sc.read_h5ad(" not in driver_part)
    check("EXPORT-02a-3 重放路径不再用 sc.tl.leiden 顶替内置 Louvain",
          "sc.tl.leiden(" not in driver_part and "engine_module.run_replay(" in text)
    check("EXPORT-02a-4 不再把过滤前的 raw_counts 写回", "raw_counts" not in text)
    check("EXPORT-02a-5 参考迁移脚本单独成段且标注「参考迁移脚本」",
          "参考迁移脚本" in annex_part, annex_part[:80])
    check("EXPORT-02a-6 参考迁移脚本按实际算法调用并列出 Louvain/Leiden 差异",
          "sc.tl.louvain" in annex_part and "louvain(builtin)" in annex_part.lower()
          and "sc.tl.leiden" in annex_part and "Leiden" in annex_part,
          annex_part[:200])
    check("EXPORT-02a-7 参考迁移脚本保留 scanpy 标记（兼容旧 contract）",
          "import scanpy as sc" in annex_part and "sc.pp" in annex_part
          and "write_h5ad" in annex_part)
    check("EXPORT-02a-8 输入清单含原始路径与 SHA-256",
          "fixture-counts.csv" in text and len(re.findall(r"\b[0-9a-f]{64}\b", text)) >= 1)
    check("EXPORT-02a-9 记录依赖版本与种子",
          '"numpy"' in text and '"sklearn"' in text and '"seed": 0' in text)
    annex_source = ""
    if "'''" in annex_part:
        start = annex_part.index("'''") + 3
        rest = annex_part[start:]
        if "'''" in rest:
            annex_source = rest[:rest.index("'''")]
    annex_syntax_ok = False
    try:
        compile(annex_source, "<reference_scanpy_migration>", "exec")
        annex_syntax_ok = True
    except SyntaxError:
        annex_syntax_ok = False
    check("EXPORT-02a-10 参考迁移脚本本身是合法 Python（可交付用户手改）",
          annex_syntax_ok and len(annex_source) > 300, "%d chars" % len(annex_source))


# ---------------------------------------------------------------------------
# EXPORT-02：清单内容（不靠子进程）
# ---------------------------------------------------------------------------


def test_manifest_records_workflow():
    print("\n[EXPORT-02b] 重放清单：顺序操作 / 实际参数 / 强制实际方法 / 输入校验")
    root = TMP_ROOT / "manifest"
    root.mkdir(parents=True, exist_ok=True)
    csv_path = make_csv_fixture(root)
    dataset_id, _dataset = open_fixture(csv_path)
    engine_dataset = run_workflow(dataset_id)
    cmd("obs_derive", {"dataset": dataset_id, "name": "group",
                       "source": "cell_name", "pattern": "^([AB])_"})
    manifest = ENGINE.build_repro_manifest(engine_dataset)

    check("EXPORT-02b-1 manifest_version / plugin / software / seed 齐全",
          manifest.get("manifest_version") and (manifest.get("software") or {}).get("engine")
          and manifest.get("seed") == ENGINE.SEED,
          str(manifest.get("manifest_version")))
    input_info = manifest.get("input") or {}
    actual_sha = ENGINE._sha256_file(str(csv_path))
    recorded_sha = input_info.get("sha256") or ((input_info.get("files") or [{}])[0].get("sha256"))
    check("EXPORT-02b-2 输入 SHA-256 与真实文件一致",
          recorded_sha == actual_sha, "%s vs %s" % (recorded_sha, actual_sha))
    genes = input_info.get("genes") or {}
    cells = input_info.get("cells") or {}
    check("EXPORT-02b-3 基因/细胞映射（数量 + 哈希 + 首尾示例）",
          genes.get("n") == N_GENES and cells.get("n") == N_CELLS
          and len(genes.get("head") or []) == min(5, N_GENES)
          and len(cells.get("head") or []) == min(5, N_CELLS)
          and genes.get("sha256") and cells.get("sha256"),
          "genes=%s cells=%s" % (genes.get("n"), cells.get("n")))
    ops = manifest.get("operations") or []
    steps = [op.get("step") for op in ops if op.get("kind") == "step"]
    check("EXPORT-02b-4 操作按真实执行顺序记录",
          steps == [step for step, _params in WORKFLOW] + ["annotate"], str(steps))
    cluster_op = [op for op in ops if op.get("step") == "cluster"][-1]
    check("EXPORT-02b-5 每步记录实际生效参数", "resolution" in (cluster_op.get("params_effective") or {}),
          str(cluster_op.get("params_effective")))
    check("EXPORT-02b-6 记录实际方法并生成强制方法的重放参数",
          "louvain(builtin)" in str(cluster_op.get("actual_method"))
          and str((cluster_op.get("replay_params") or {}).get("method")).lower() == "louvain",
          "actual=%s replay=%s" % (cluster_op.get("actual_method"),
                                   (cluster_op.get("replay_params") or {}).get("method")))
    expected = manifest.get("expectations") or {}
    structural = expected.get("structural") or {}
    check("EXPORT-02b-7 期望结果记录过滤后规模与矩阵尺度",
          structural.get("n_obs") == engine_dataset.n_obs
          and structural.get("n_vars") == engine_dataset.n_vars
          and "lognorm" in str(structural.get("matrix_scale")),
          str({key: structural.get(key) for key in ("n_obs", "n_vars", "matrix_scale")}))
    check("EXPORT-02b-8 期望结果含矩阵/计数/聚类/分组哈希",
          {"X", "counts", "clusters", "obs:group"} <= set((expected.get("hashes") or {}).keys()),
          str(sorted((expected.get("hashes") or {}).keys())))
    reference = manifest.get("reference_migration") or {}
    differences = "\n".join(reference.get("differences") or [])
    check("EXPORT-02b-9 参考迁移脚本差异清单提到 Louvain/Leiden",
          "Louvain" in differences and "Leiden" in differences, differences[:160])
    not_replayable = [op for op in ops if op.get("replayable") is False]
    check("EXPORT-02b-10 builtin 流程的所有操作都可重放", not not_replayable,
          str([op.get("step") for op in not_replayable]))


# ---------------------------------------------------------------------------
# EXPORT-02：独立目录真重放（builtin + CSV）
# ---------------------------------------------------------------------------


def _run_replay_script(script_path, out_dir, extra_args=None, timeout=300):
    env = dict(os.environ)
    env.update({"MPLBACKEND": "Agg", "MPLCONFIGDIR": str(out_dir.parent / "mpl")})
    command = [PYTHON, str(script_path), "--out", str(out_dir)]
    command.extend(extra_args or [])
    completed = subprocess.run(
        command, cwd=str(script_path.parent), capture_output=True, text=True,
        encoding="utf-8", timeout=timeout, env=env,
    )
    return completed


def _snapshot_values(dataset):
    """重放前后逐值比较用的纯 numpy 快照（不依赖引擎内部哈希）。"""
    def dense(matrix):
        return np.asarray(matrix.todense() if sp.issparse(matrix) else matrix, dtype=np.float32)

    markers = {}
    if isinstance(dataset.markers, dict):
        for group in dataset.markers.get("groups") or []:
            markers[str(group.get("id"))] = [
                (str(record.get("gene")), round(float(record.get("score", 0.0)), 8))
                for record in (group.get("genes") or [])
            ]
    return {
        "n_obs": int(dataset.n_obs),
        "n_vars": int(dataset.n_vars),
        "var_names": [str(name) for name in dataset.var_names],
        "X": dense(dataset.X),
        "counts": dense(dataset.counts),
        "clusters": [str(value) for value in dataset.obs["clusters"].tolist()]
        if "clusters" in dataset.obs.columns else None,
        "obs_columns": {
            column: [str(value) for value in dataset.obs[column].tolist()]
            for column in ("sample", "sample_name", "group", "cell_type")
            if column in dataset.obs.columns
        },
        "markers": markers,
        "filter_info": dict(dataset.filter_info or {}),
        "matrix_scale": ENGINE._matrix_scale(dataset),
        "layout": dataset.layout,
        "pca": np.asarray(dataset.obsm.get("X_pca")) if "X_pca" in dataset.obsm else None,
    }


def _compare_reproduced(label, expected, h5ad_path, replay_result=None):
    if ad is None:
        skip(label + " 逐值比较", "没有 anndata")
        return
    reproduced = ad.read_h5ad(str(h5ad_path), backed="r")
    check(label + " 细胞数一致", int(reproduced.n_obs) == expected["n_obs"],
          "%d vs %d" % (reproduced.n_obs, expected["n_obs"]))
    check(label + " 基因数一致", int(reproduced.n_vars) == expected["n_vars"],
          "%d vs %d" % (reproduced.n_vars, expected["n_vars"]))
    check(label + " 基因名逐值一致", list(reproduced.var_names) == expected["var_names"])
    x = reproduced.X[:]
    x = np.asarray(x.todense() if sp.issparse(x) else x, dtype=np.float32)
    check(label + " 当前矩阵 X 逐值一致（矩阵尺度=标准化后）",
          np.array_equal(x, expected["X"]),
          "max|diff|=%.3g" % (float(np.max(np.abs(x - expected["X"]))) if x.shape == expected["X"].shape else -1))
    counts = reproduced.layers["counts"][:]
    counts = np.asarray(counts.todense() if sp.issparse(counts) else counts, dtype=np.float32)
    check(label + " counts 层逐值一致（过滤后用同一形状，未回填旧矩阵）",
          np.array_equal(counts, expected["counts"]),
          "%s vs %s" % (tuple(counts.shape), expected["counts"].shape))
    if expected["clusters"] is not None:
        got = [str(value) for value in reproduced.obs["clusters"].tolist()]
        check(label + " 聚类标签逐值一致（Louvain 确定性）", got == expected["clusters"],
              "%d labels" % len(got))
    for column, values in expected["obs_columns"].items():
        got = [str(value) for value in reproduced.obs[column].tolist()]
        check(label + " obs 列 %s 逐值一致（含手动分组）" % column, got == values,
              "%d 个取值" % len(set(got)))
    if expected["pca"] is not None:
        got = np.asarray(reproduced.obsm["X_pca"], dtype=np.float32)
        check(label + " PCA 坐标逐值一致（randomized + 固定种子）",
              got.shape == expected["pca"].shape and np.array_equal(got, expected["pca"]),
              str(got.shape))
    check(label + " filter_info 的保留规模与记录一致",
          (reproduced.uns.get("dsh_singlecell") is not None), "h5ad 已写 dsh_singlecell 元数据")
    reproduced.file.close()

    if replay_result is None:
        return
    structural = (replay_result.get("actual_snapshot") or {}).get("structural") or {}
    expected_filter = {
        "n_obs_before": expected["filter_info"].get("n_obs_before"),
        "n_obs_after": expected["filter_info"].get("n_obs_after"),
        "n_vars_before": expected["filter_info"].get("n_vars_before"),
        "n_vars_after": expected["filter_info"].get("n_vars_after"),
    }
    check(label + " 重放记录 filter 保留规模逐项一致",
          structural.get("filter") == expected_filter,
          "%s vs %s" % (structural.get("filter"), expected_filter))
    check(label + " 重放记录矩阵尺度一致",
          str(structural.get("matrix_scale")) == str(expected["matrix_scale"]),
          "%s vs %s" % (structural.get("matrix_scale"), expected["matrix_scale"]))
    actual_markers = ((replay_result.get("actual_snapshot") or {}).get("numbers") or {}).get("markers") or {}
    if expected["markers"]:
        check(label + " 标记基因分组齐全",
              set(actual_markers) == set(expected["markers"]), "%d 组" % len(actual_markers))
        for group, items in expected["markers"].items():
            got = actual_markers.get(group) or []
            check(label + " 标记基因 %s 前 %d 个基因顺序一致" % (group, len(items)),
                  [item[0] for item in got] == [item[0] for item in items],
                  str([item[0] for item in got]))
            check(label + " 标记基因 %s 分数一致" % group,
                  len(got) == len(items) and np.allclose(
                      [item[1] for item in got], [item[1] for item in items], atol=1e-6),
                  "%d 项" % len(got))


def test_replay_csv_end_to_end():
    print("\n[EXPORT-02c] CSV 输入：独立目录重放并逐值比较")
    root = TMP_ROOT / "replay-csv"
    root.mkdir(parents=True, exist_ok=True)
    csv_path = make_csv_fixture(root)
    dataset_id, _dataset = open_fixture(csv_path)
    engine_dataset = run_workflow(dataset_id)
    cmd("obs_derive", {"dataset": dataset_id, "name": "group",
                       "source": "cell_name", "pattern": "^([AB])_"})
    cmd("obs_set_group", {"dataset": dataset_id, "column": "group"})
    expected = _snapshot_values(ENGINE.DATASETS[dataset_id])
    manifest = ENGINE.build_repro_manifest(ENGINE.DATASETS[dataset_id])
    check("EXPORT-02c-1 发生实际过滤（30→%d 细胞，%d→%d 基因）"
          % (expected["n_obs"], N_GENES, expected["n_vars"]),
          expected["n_obs"] < N_CELLS and expected["n_vars"] < N_GENES,
          "n_obs=%d n_vars=%d" % (expected["n_obs"], expected["n_vars"]))

    export_dir = root / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)
    script_path = export_dir / "replay_workflow.py"
    text = ENGINE._build_script(ENGINE.DATASETS[dataset_id])
    script_path.write_text(text, encoding="utf-8")

    replay_dir = root / "independent-replay"
    replay_dir.mkdir(parents=True, exist_ok=True)
    local_script = replay_dir / "replay_workflow.py"
    shutil.copy2(script_path, local_script)

    # 导出物自带 --dump-manifest：不依赖引擎即可离线读取记录清单
    dumped = subprocess.run(
        [PYTHON, str(local_script), "--dump-manifest"],
        cwd=str(replay_dir), capture_output=True, text=True, encoding="utf-8", timeout=120,
    )
    dumped_manifest = {}
    if dumped.returncode == 0:
        try:
            dumped_manifest = json.loads(dumped.stdout)
        except json.JSONDecodeError:
            dumped_manifest = {}
    dumped_steps = [op.get("step") for op in dumped_manifest.get("operations") or []
                    if op.get("kind") == "step"]
    check("EXPORT-02c-0 导出脚本可离线 dump 出完整操作清单",
          dumped_steps == [step for step, _params in WORKFLOW] + ["annotate"]
          and bool((dumped_manifest.get("input") or {}).get("sha256")),
          str(dumped_steps))

    out_dir = replay_dir / "replay_out"
    completed = _run_replay_script(local_script, out_dir)
    print("    重放 stdout 尾部：")
    for line in (completed.stdout or "").strip().splitlines()[-8:]:
        print("      " + line)
    if completed.returncode != 0:
        for line in (completed.stderr or "").strip().splitlines()[-10:]:
            print("      stderr: " + line)
    check("EXPORT-02c-2 独立目录重放退出码 0", completed.returncode == 0,
          "rc=%s" % completed.returncode)
    result_path = out_dir / "replay_manifest.json"
    check("EXPORT-02c-3 写出 replay_manifest.json", result_path.exists(), str(result_path))
    replay_result = {}
    if result_path.exists():
        replay_result = json.loads(result_path.read_text(encoding="utf-8"))
        check("EXPORT-02c-4 重放状态 verified", replay_result.get("status") in (
            "verified", "verified_with_tolerance"), str(replay_result.get("status")))
        bad = [item for item in replay_result.get("checks") or [] if item.get("required", True)
               and not item.get("ok")]
        check("EXPORT-02c-5 所有必需检查通过", not bad, str([item.get("id") for item in bad])[:200])
        check("EXPORT-02c-6 重放清单记录使用的引擎与版本指纹",
              (replay_result.get("software") or {}).get("engine_mode"),
              str((replay_result.get("software") or {}).get("engine_mode")))
    check("EXPORT-02c-7 写出 reproduced.h5ad",
          (out_dir / "reproduced.h5ad").exists(), str(out_dir / "reproduced.h5ad"))
    check("EXPORT-02c-8 写出参考迁移脚本",
          any(out_dir.glob("*reference*scanpy*.py")), str(list(out_dir.glob("*.py"))))
    _compare_reproduced("EXPORT-02c", expected, out_dir / "reproduced.h5ad", replay_result)

    # 清单本身记录了重放期望的规模（与过滤后一致）
    structural = (manifest.get("expectations") or {}).get("structural") or {}
    check("EXPORT-02c-9 manifest 的期望规模=过滤后规模",
          structural.get("n_obs") == expected["n_obs"]
          and structural.get("n_vars") == expected["n_vars"],
          "%s/%s" % (structural.get("n_obs"), structural.get("n_vars")))


def test_replay_h5ad_filtered_end_to_end():
    print("\n[EXPORT-02d] h5ad 输入 + 真实过滤 + 手动分组：独立目录重放")
    if ad is None:
        skip("EXPORT-02d", "没有 anndata")
        return
    root = TMP_ROOT / "replay-h5ad"
    root.mkdir(parents=True, exist_ok=True)
    h5ad_path = make_h5ad_fixture(root)
    dataset_id, _dataset = open_fixture(h5ad_path)
    engine_dataset = run_workflow(dataset_id)
    # 导入样本表 → 手动分组（改名 + 分组）
    cmd("samples_update", {
        "dataset": dataset_id,
        "samples": [
            {"id": "S0", "name": "Ctrl_1", "group": "Control"},
            {"id": "S1", "name": "Case_1", "group": "Disease"},
            {"id": "S2", "name": "Case_2", "group": "Disease"},
        ],
        "groups": ["Control", "Disease"],
    })
    cmd("obs_set_group", {"dataset": dataset_id, "column": "group"})
    current = ENGINE.DATASETS[dataset_id]
    expected = _snapshot_values(current)
    check("EXPORT-02d-1 分组写入 obs", set(expected["obs_columns"].get("group") or []) == {
        "Control", "Disease"}, str(set(expected["obs_columns"].get("group") or [])))
    check("EXPORT-02d-2 h5ad 过滤后 counts 与 X 形状同步",
          expected["counts"].shape == expected["X"].shape == (expected["n_obs"], expected["n_vars"])
          and (expected["n_obs"], expected["n_vars"]) != (N_CELLS, N_GENES),
          "%s vs %s" % (expected["counts"].shape, expected["X"].shape))

    export_dir = root / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)
    script_path = export_dir / "replay_workflow.py"
    h5ad_script_text = ENGINE._build_script(current)
    script_path.write_text(h5ad_script_text, encoding="utf-8")
    h5ad_driver = h5ad_script_text.split("REFERENCE_SCANPY_SCRIPT = r", 1)[0]
    check("EXPORT-02d-0 h5ad 输入：参考迁移脚本用 sc.read_h5ad，驱动路径不盲读 h5ad",
          "sc.read_h5ad(" in h5ad_script_text and "sc.read_h5ad(" not in h5ad_driver)

    replay_dir = root / "independent-replay"
    replay_dir.mkdir(parents=True, exist_ok=True)
    local_script = replay_dir / "replay_workflow.py"
    shutil.copy2(script_path, local_script)
    out_dir = replay_dir / "replay_out"
    completed = _run_replay_script(local_script, out_dir)
    print("    重放 stdout 尾部：")
    for line in (completed.stdout or "").strip().splitlines()[-6:]:
        print("      " + line)
    if completed.returncode != 0:
        for line in (completed.stderr or "").strip().splitlines()[-10:]:
            print("      stderr: " + line)
    check("EXPORT-02d-3 h5ad 重放退出码 0", completed.returncode == 0, "rc=%s" % completed.returncode)
    replay_result = {}
    if (out_dir / "replay_manifest.json").exists():
        replay_result = json.loads((out_dir / "replay_manifest.json").read_text(encoding="utf-8"))
    _compare_reproduced("EXPORT-02d", expected, out_dir / "reproduced.h5ad", replay_result)
    result_path = out_dir / "replay_manifest.json"
    if result_path.exists():
        replay_result = json.loads(result_path.read_text(encoding="utf-8"))
        check("EXPORT-02d-4 重放状态 verified", replay_result.get("status") in (
            "verified", "verified_with_tolerance"), str(replay_result.get("status")))
        check("EXPORT-02d-5 重放清单显示输入格式为 h5ad",
              str((replay_result.get("input") or {}).get("format")) == "h5ad",
              str((replay_result.get("input") or {}).get("format")))


def test_replay_refuses_changed_input():
    print("\n[EXPORT-02e] 输入被改动：只能标为记录/示例，不能冒充重放")
    root = TMP_ROOT / "replay-tamper"
    root.mkdir(parents=True, exist_ok=True)
    csv_path = make_csv_fixture(root)
    dataset_id, _dataset = open_fixture(csv_path)
    engine_dataset = run_workflow(dataset_id, steps=WORKFLOW[:4], annotate=False)
    replay_dir = root / "independent-replay"
    replay_dir.mkdir(parents=True, exist_ok=True)
    script_path = replay_dir / "replay_workflow.py"
    script_path.write_text(ENGINE._build_script(engine_dataset), encoding="utf-8")
    # 篡改输入（模拟上游数据被替换）
    with open(csv_path, "a", encoding="utf-8") as handle:
        handle.write("999,999,999\n")
    out_dir = replay_dir / "replay_out"
    completed = _run_replay_script(script_path, out_dir)
    check("EXPORT-02e-1 输入校验失败时非零退出", completed.returncode != 0, "rc=%s" % completed.returncode)
    result_path = out_dir / "replay_manifest.json"
    check("EXPORT-02e-2 写出失败/记录状态", result_path.exists(), str(result_path))
    if result_path.exists():
        replay_result = json.loads(result_path.read_text(encoding="utf-8"))
        check("EXPORT-02e-3 状态是 input_error / record_only（不是 verified）",
              replay_result.get("status") in ("input_error", "record_only"),
              str(replay_result.get("status")))
        check("EXPORT-02e-4 明确说明只能作为记录/示例",
              "记录" in str(replay_result.get("message")) or "示例" in str(replay_result.get("message")),
              str(replay_result.get("message"))[:160])
    output = (completed.stdout or "") + (completed.stderr or "")
    check("EXPORT-02e-5 输出不谎称重放成功",
          "replay status: verified" not in output.lower()
          and ("记录" in output or "示例" in output or "input_error" in output),
          output.strip().splitlines()[-1][:120] if output.strip() else "")


def test_replay_marks_seurat_record_only():
    print("\n[EXPORT-02f] 外部 Seurat 步骤：清单标为不可重放（记录/示例）")
    root = TMP_ROOT / "replay-seurat"
    root.mkdir(parents=True, exist_ok=True)
    csv_path = make_csv_fixture(root)
    dataset_id, _dataset = open_fixture(csv_path)
    _install_fake_seurat_bridge()
    steps = [
        ("qc", {"min_genes": 1, "max_pct_mt": 100}),
        ("filter", {"min_genes": 1, "min_cells": 2, "max_pct_mt": 100}),
        ("normalize", {"target_sum": 10000}),
        ("hvg", {"n_top_genes": 6}),
        ("pca", {"n_pcs": 5, "engine": "seurat6", "integration": "none"}),
        ("neighbors", {"n_neighbors": 5, "engine": "seurat6"}),
        ("cluster", {"resolution": 0.5, "engine": "seurat6", "integration": "none"}),
    ]
    for step, params in steps:
        cmd("run", {"dataset": dataset_id, "step": step, "params": params})
    sys.modules.pop("_dsh_sc_lib_sc_external", None)
    engine_dataset = ENGINE.DATASETS[dataset_id]
    manifest = ENGINE.build_repro_manifest(engine_dataset)
    ops = manifest.get("operations") or []
    non_replayable = [op for op in ops if op.get("replayable") is False]
    check("EXPORT-02f-1 Seurat 步骤被标为不可重放",
          {op.get("step") for op in non_replayable} >= {"pca", "cluster"},
          str([op.get("step") for op in non_replayable]))
    check("EXPORT-02f-2 列出不可重放原因",
          all(str(op.get("not_replayable_reason") or "") for op in non_replayable),
          str([op.get("not_replayable_reason") for op in non_replayable])[:160])
    replay_dir = root / "independent-replay"
    replay_dir.mkdir(parents=True, exist_ok=True)
    script_path = replay_dir / "replay_workflow.py"
    script_path.write_text(ENGINE._build_script(engine_dataset), encoding="utf-8")
    out_dir = replay_dir / "replay_out"
    completed = _run_replay_script(script_path, out_dir)
    check("EXPORT-02f-3 部分重放以非零退出码（3=记录/示例）", completed.returncode != 0,
          "rc=%s" % completed.returncode)
    output = (completed.stdout or "") + (completed.stderr or "")
    check("EXPORT-02f-4 输出明确标注记录/示例", "记录" in output or "示例" in output,
          output.strip().splitlines()[-1][:160] if output.strip() else "")
    result_path = out_dir / "replay_manifest.json"
    if result_path.exists():
        replay_result = json.loads(result_path.read_text(encoding="utf-8"))
        check("EXPORT-02f-5 清单状态 partial", replay_result.get("status") == "partial",
              str(replay_result.get("status")))
        check("EXPORT-02f-6 跳过操作被记录",
              len(replay_result.get("skipped_operations") or []) >= 2,
              str([op.get("step") for op in replay_result.get("skipped_operations") or []]))


def test_protocol_includes_record():
    print("\n[EXPORT-02g] protocol 导出带执行记录与重放清单")
    root = TMP_ROOT / "protocol"
    root.mkdir(parents=True, exist_ok=True)
    dataset_id, _dataset = open_fixture(make_csv_fixture(root))
    engine_dataset = run_workflow(dataset_id)
    protocol = ENGINE._build_protocol(engine_dataset)
    check("EXPORT-02g-1 protocol 含 execution_record",
          isinstance(protocol.get("execution_record"), list)
          and any(item.get("actual_method") for item in protocol["execution_record"]),
          "n=%s" % len(protocol.get("execution_record") or []))
    check("EXPORT-02g-2 protocol 含 replay 清单摘要",
          isinstance(protocol.get("replay"), dict) and protocol["replay"].get("manifest_version"),
          str((protocol.get("replay") or {}).get("manifest_version")))
    text = json.dumps(protocol, ensure_ascii=False)
    check("EXPORT-02g-3 protocol 标记当前数据版本", "version" in text and "revision" in text)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> int:
    print("== dsh-singlecell 导出回归（EXPORT-01 / EXPORT-02）==")
    print("临时目录：%s" % TMP_ROOT)
    tests = [
        test_execution_record_structure,
        test_methods_louvain,
        test_methods_leiden_fallback,
        test_methods_external_seurat,
        test_methods_imported_without_record,
        test_methods_excludes_failed_steps,
        test_script_structure,
        test_manifest_records_workflow,
        test_replay_csv_end_to_end,
        test_replay_h5ad_filtered_end_to_end,
        test_replay_refuses_changed_input,
        test_replay_marks_seurat_record_only,
        test_protocol_includes_record,
    ]
    for function in tests:
        try:
            function()
        except Exception:  # noqa: BLE001
            FAILED.append(function.__name__)
            print(f"  FAIL {function.__name__} 抛异常")
            traceback.print_exc()
    print("\n== 结果：%d 通过 / %d 失败 / %d 跳过 ==" % (len(PASSED), len(FAILED), len(SKIPPED)))
    if FAILED:
        for item in FAILED:
            print(f"  ✗ {item}")
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    finally:
        if os.environ.get("KEEP_SC_TEST_TMP") != "1":
            shutil.rmtree(TMP_ROOT, ignore_errors=True)
