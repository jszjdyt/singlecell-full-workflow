#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell 官方 Python 工具适配器测试（CONTRACT §10）。

覆盖：pySCENIC / CellPhoneDB / scVelo / Palantir / CellOracle / RegVelo。

用法：
  python3 test/py_external_py_test.py                     # 全量：生成/复用 demo bundle，逐个工具跑
  python3 test/py_external_py_test.py --tools pyscenic,celloracle
  python3 test/py_external_py_test.py --force-bundle      # 强制重建 demo bundle
  python3 test/py_external_py_test.py --bundle /path/dir  # 复用 Host（sc_external）导出的 bundle
  python3 test/py_external_py_test.py --scvelo-synthetic  # 额外跑 scVelo 正路径冒烟（合成剪接层）

判定：
  PASS = 适配器退出码 0 且 result.json/图/obs_extra 通过校验；
  PASS(prereq) = 预期缺前置的工具（如 scVelo 无 spliced/unspliced）正确给出退出码 3 + 中文原因；
  SKIP = 依赖未安装（附原因）；
  FAIL = 其它任何不符合契约的情况。
退出码：无 FAIL → 0；有 FAIL → 1。

说明：本测试会以 /opt/miniconda3/bin/python3（conda base，插件 externalPython）
运行适配器；若当前解释器不是它，会自动 execv 切换（demo bundle 生成需要 scanpy）。
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import time
import traceback
from pathlib import Path

CONDA_PY = "/opt/miniconda3/bin/python3"
PLUGIN_ROOT = Path(__file__).resolve().parents[1]
EXT_DIR = PLUGIN_ROOT / "lib" / "external"
TEST_DIR = PLUGIN_ROOT / "test"
CACHE = Path(os.path.expanduser("~/.dsh/sc-plugin-pydata/pyext-test"))
DEFAULT_BUNDLE = CACHE / "bundle-demo"
LAYERED_BUNDLE = CACHE / "bundle-demo-layered"
OUT_ROOT = CACHE / "outputs"
SCENIC_DB_CACHE = Path(os.path.expanduser("~/.dsh/sc-plugin-pydata/cistarget"))
CELLORACLE_BASE_GRN = Path(os.path.expanduser(
    "~/celloracle_data/promoter_base_GRN/"
    "hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"
))
CPDB_DB = Path(os.path.expanduser("~/.dsh/sc-plugin-pydata/cellphonedb/cellphonedb.zip"))

# ---------------------------------------------------------------------------
# 工具清单（与 lib/sc_external.py TOOLS 对齐）
# ---------------------------------------------------------------------------

TOOLS = [
    {
        "id": "pyscenic",
        "title": "pySCENIC 官方转录调控",
        "script": "pyscenic.py",
        "modules": ["pyscenic", "ctxcore", "arboreto"],
        "params": {"n_workers": 4, "n_top_genes": 1200},
        "timeout": 3600,
        "min_regulons": 1,
    },
    {
        "id": "cellphonedb",
        "title": "CellPhoneDB 官方通讯",
        "script": "cellphonedb.py",
        "modules": ["cellphonedb"],
        "params": {"n_perms": 100, "n_jobs": 4},
        "timeout": 1800,
    },
    {
        "id": "scvelo",
        "title": "scVelo 官方 RNA velocity",
        "script": "scvelo.py",
        "modules": ["scvelo"],
        "params": {"mode": "stochastic"},
        "timeout": 900,
        "expect_prereq": True,  # demo bundle（§10.2）没有 spliced/unspliced 层
    },
    {
        "id": "palantir",
        "title": "Palantir 官方轨迹",
        "script": "palantir.py",
        "modules": ["palantir"],
        "params": {"n_components": 10, "n_jobs": 4, "num_waypoints": 500},
        "timeout": 1800,
    },
    {
        "id": "celloracle",
        "title": "CellOracle 官方虚拟敲除",
        "script": "celloracle.py",
        "modules": ["celloracle"],
        "params": {"target": "IRF8", "mode": "knockout", "alpha": 10, "n_jobs": 4},
        "timeout": 2400,
    },
    {
        "id": "regvelo",
        "title": "RegVelo 官方调控速度",
        "script": "regvelo.py",
        "modules": ["regvelo"],
        "params": {},
        "timeout": 300,
        "expect_unavailable": True,  # torch<2.6 无 macOS arm64 + py3.13 wheel，装不上
    },
]


# ---------------------------------------------------------------------------
# 基础工具
# ---------------------------------------------------------------------------


def color(text, code):
    if not sys.stdout.isatty():
        return text
    return "\033[%sm%s\033[0m" % (code, text)


def info(message):
    print(color("[test]", "36") + " " + str(message), flush=True)


def ok(message):
    print(color("[PASS]", "32") + " " + str(message), flush=True)


def skip(message):
    print(color("[SKIP]", "33") + " " + str(message), flush=True)


def fail(message):
    print(color("[FAIL]", "31") + " " + str(message), flush=True)


def maybe_reexec_conda():
    """demo bundle 生成需要 scanpy/leidenalg，切到插件的 externalPython（conda base）。"""
    if os.environ.get("DS_SC_PYEXT_REEXEC") == "1":
        return
    if not Path(CONDA_PY).exists():
        return
    if os.path.realpath(sys.executable) == os.path.realpath(CONDA_PY):
        return
    info("当前解释器 %s 不是 externalPython，切换到 %s 执行测试" % (sys.executable, CONDA_PY))
    env = dict(os.environ)
    env["DS_SC_PYEXT_REEXEC"] = "1"
    os.execve(CONDA_PY, [CONDA_PY, str(Path(__file__).resolve())] + sys.argv[1:], env)


def run_command(command, timeout=None, env=None):
    started = time.time()
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, env=env or os.environ.copy()
        )
        return {
            "returncode": completed.returncode,
            "stdout": completed.stdout or "",
            "stderr": completed.stderr or "",
            "elapsed": time.time() - started,
            "timeout": False,
        }
    except subprocess.TimeoutExpired:
        return {"returncode": -9, "stdout": "", "stderr": "timeout", "elapsed": time.time() - started, "timeout": True}
    except Exception as error:
        return {"returncode": -1, "stdout": "", "stderr": str(error), "elapsed": time.time() - started, "timeout": False}


def probe_modules(modules):
    """用 externalPython 一次性探测若干模块，返回 {module: version|None}。"""
    code = (
        "import importlib, json\n"
        "out = {}\n"
        "for name in %r:\n"
        "    try:\n"
        "        m = importlib.import_module(name)\n"
        "        try:\n"
        "            from importlib.metadata import version as _v\n"
        "            out[name] = str(_v(name))\n"
        "        except Exception:\n"
        "            out[name] = str(getattr(m, '__version__', 'OK'))\n"
        "    except Exception as error:\n"
        "        out[name] = None\n"
        "print(json.dumps(out))\n" % (list(modules),)
    )
    result = run_command([CONDA_PY, "-c", code], timeout=300)
    payload = {}
    for line in (result["stdout"] or "").splitlines():
        line = line.strip()
        if line.startswith("{"):
            try:
                payload = json.loads(line)
            except Exception:
                pass
    if not payload:
        return {module: None for module in modules}
    return payload


# ---------------------------------------------------------------------------
# demo bundle 生成（§10.2）
# ---------------------------------------------------------------------------


def _write_obs_csv(path, obs, barcodes):
    def clean(value):
        text = str(value)
        return text.replace(",", " ").replace("\n", " ").replace("\r", " ")

    columns = ["cell"] + list(obs.columns)
    lines = ["\t".join(columns)] if False else [",".join(columns)]
    for barcode, (_, row) in zip(barcodes, obs.iterrows()):
        lines.append(",".join([clean(barcode)] + [clean(row[column]) for column in obs.columns]))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_bundle(dest, force=False, seed=0):
    """用 test/make_sample.py 的同一份 demo 数据生成 §10.2 bundle。"""
    dest = Path(dest)
    meta_path = dest / "meta.json"
    if meta_path.exists() and not force:
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            info("复用已有 demo bundle：%s（%d × %d，X_kind=%s）"
                 % (dest, meta.get("n_obs", -1), meta.get("n_vars", -1), meta.get("X_kind")))
            return meta
        except Exception:
            pass
    info("生成 demo bundle（test/make_sample.py → §10.2，含 lognorm/clusters/embed/pca）…")
    if str(TEST_DIR) not in sys.path:
        sys.path.insert(0, str(TEST_DIR))
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp
    from scipy.io import mmwrite

    import make_sample

    data = make_sample.generate_demo(n_cells=2400, n_genes=1200, seed=seed)
    counts = data["counts"]
    obs = data["obs"].copy()
    genes = [str(item) for item in data["var_names"]]
    keep = ~np.asarray(data["is_low_quality"], dtype=bool)  # 对应 qc/filter 后 2208 细胞
    counts = counts[keep]
    obs = obs.iloc[np.where(keep)[0]].copy()
    barcodes = [str(item) for item in obs.index]

    # normalize：log1p(CP10K)（与内置引擎 log 标准化一致）
    row_sum = np.asarray(counts.sum(axis=1)).ravel()
    row_sum[row_sum <= 0] = 1.0
    data_mtx = counts.multiply((1e4 / row_sum)[:, None]).tocsr()
    data_mtx.data = np.log1p(data_mtx.data)
    data_mtx = data_mtx.astype("float32")

    # 聚类/嵌入：忠实复刻工作台流程 qc→filter→normalize→hvg→pca→neighbors→cluster→embed
    adata = ad.AnnData(X=data_mtx.copy(), obs=obs.copy(),
                       var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    sc.pp.highly_variable_genes(adata, n_top_genes=min(500, adata.n_vars - 1))
    sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1, adata.n_obs - 1),
              use_highly_variable=True, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca", random_state=seed)
    clause = None
    try:
        sc.tl.leiden(adata, resolution=0.5, random_state=seed, key_added="clusters",
                     flavor="igraph", n_iterations=2, directed=False)
        clause = "leiden"
    except Exception as error:
        info("leiden 聚类失败（%s），回退用 cell_type 当 clusters" % error)
        adata.obs["clusters"] = adata.obs["cell_type"].astype(str)
        clause = "cell_type"
    sc.tl.umap(adata, random_state=seed)
    embed = np.asarray(adata.obsm["X_umap"], dtype="float64")
    pca = np.asarray(adata.obsm["X_pca"], dtype="float64")
    obs_out = adata.obs[["cell_type", "batch", "clusters"]].copy()
    obs_out["clusters"] = obs_out["clusters"].astype(str)

    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "features.tsv").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (dest / "barcodes.tsv").write_text("\n".join(barcodes) + "\n", encoding="utf-8")
    mmwrite(str(dest / "counts.mtx"), sp.coo_matrix(counts.T))
    mmwrite(str(dest / "data.mtx"), sp.coo_matrix(data_mtx.T))
    _write_obs_csv(dest / "obs.csv", obs_out, barcodes)
    with open(dest / "embed.csv", "w", encoding="utf-8") as handle:
        handle.write("cell,dim1,dim2\n")
        for index, barcode in enumerate(barcodes):
            handle.write("%s,%.6f,%.6f\n" % (barcode, embed[index, 0], embed[index, 1]))
    with open(dest / "pca.csv", "w", encoding="utf-8") as handle:
        handle.write("cell," + ",".join("PC%d" % (i + 1) for i in range(pca.shape[1])) + "\n")
        for index, barcode in enumerate(barcodes):
            handle.write("%s,%s\n" % (barcode, ",".join("%.6f" % value for value in pca[index])))
    meta = {
        "dataset": "demo-pyext",
        "name": "demo-pbmc",
        "n_obs": len(barcodes),
        "n_vars": len(genes),
        "steps": ["qc", "filter", "normalize", "hvg", "pca", "neighbors", "cluster", "embed"],
        "X_kind": "lognorm",
        "engine": "builtin",
        "clusters_method": clause,
        "n_clusters": int(obs_out["clusters"].nunique()),
    }
    (dest / "meta.json").write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    info("bundle 生成完成：%s（%d 细胞 × %d 基因，%d 簇）"
         % (dest, meta["n_obs"], meta["n_vars"], meta["n_clusters"]))
    return meta


def build_layered_bundle(source, dest, seed=0, force=False):
    """在 demo bundle 基础上合成 spliced/unspliced 层（仅用于 scVelo 正路径冒烟）。"""
    dest = Path(dest)
    if (dest / "meta.json").exists() and not force:
        info("复用已有合成层 bundle：%s" % dest)
        return
    info("生成合成剪接层 bundle（spliced=counts、unspliced=Poisson(0.3×counts)，仅冒烟）…")
    import numpy as np
    import scipy.sparse as sp
    from scipy.io import mmread, mmwrite

    source = Path(source)
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(source, dest)
    counts = mmread(str(source / "counts.mtx")).tocsr()
    rng = np.random.default_rng(seed)
    density = counts.copy()
    density.data = np.clip(density.data, 0, None)
    unspliced = counts.copy()
    unspliced.data = rng.poisson(np.clip(unspliced.data, 0, 50) * 0.3).astype("float32")
    unspliced.eliminate_zeros()
    mmwrite(str(dest / "spliced.mtx"), counts.astype("float32"))
    mmwrite(str(dest / "unspliced.mtx"), unspliced.tocsr().astype("float32"))
    meta_path = dest / "meta.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["layers"] = ["spliced", "unspliced"]
    meta["layers_note"] = "合成层：仅用于 scVelo 代码路径冒烟，不代表真实 RNA velocity 生物学结论"
    meta_path.write_text(json.dumps(meta, ensure_ascii=False), encoding="utf-8")
    info("合成层 bundle 完成：%s" % dest)


# ---------------------------------------------------------------------------
# 结果校验（§10.3）
# ---------------------------------------------------------------------------


def _png_info(path):
    """返回 (width, height, dpi, 左上角是否白底)。"""
    try:
        from PIL import Image

        with Image.open(path) as image:
            width, height = image.size
            dpi = image.info.get("dpi")
            if isinstance(dpi, tuple):
                dpi = float(dpi[0])
            elif dpi is not None:
                dpi = float(dpi)
            rgb = image.convert("RGB")
            corner = [rgb.getpixel((x, y)) for x in (0, 1, width - 1) for y in (0, 1, height - 1)]
            white = all(min(pixel) >= 240 for pixel in corner)
            return width, height, dpi, white
    except Exception:
        with open(path, "rb") as handle:
            head = handle.read(24)
        if len(head) >= 24 and head[:8] == b"\x89PNG\r\n\x1a\n":
            width, height = struct.unpack(">II", head[16:24])
            return width, height, None, None
        return 0, 0, None, None


def validate_result(out_dir, bundle_dir, tool, elapsed):
    """校验 result.json / 图 / obs_extra；返回 (ok, 摘要dict, 错误列表)。"""
    errors = []
    out_dir = Path(out_dir)
    result_path = out_dir / "result.json"
    if not result_path.exists():
        return False, {}, ["没有生成 result.json"]
    try:
        result = json.loads(result_path.read_text(encoding="utf-8"))
    except Exception as error:
        return False, {}, ["result.json 不是合法 JSON：%s" % error]
    for key in ("title", "method", "summary", "tables", "figures", "warnings", "notes"):
        if key not in result:
            errors.append("result.json 缺少键 %s" % key)
    summary = result.get("summary") or {}
    if not isinstance(summary, dict):
        errors.append("summary 不是 dict")
        summary = {}
    if not isinstance(summary.get("bullets"), list) or not summary.get("bullets"):
        errors.append("summary.bullets 为空")
    if not isinstance(summary.get("metrics"), dict) or not summary.get("metrics"):
        errors.append("summary.metrics 为空")
    figures = result.get("figures") or []
    if not figures:
        errors.append("figures 为空（至少应有 1 张图）")
    figure_notes = []
    for figure in figures:
        path = Path(str(figure.get("file") or ""))
        if not path.is_absolute():
            path = (out_dir / path).resolve()
        if not path.exists():
            errors.append("图不存在：%s" % path)
            continue
        width, height, dpi, white = _png_info(str(path))
        if width <= 0 or height <= 0:
            errors.append("图不是合法 PNG：%s" % path.name)
            continue
        if dpi is not None and dpi < 139:
            errors.append("图 dpi=%s < 140：%s" % (dpi, path.name))
        if white is False:
            errors.append("图不是白底（角落像素非白）：%s" % path.name)
        if figure.get("width") and int(figure["width"]) != int(width):
            errors.append("图 %s 的 width 元数据与实际不一致" % path.name)
        figure_notes.append("%s %dx%d@%s" % (path.name, width, height, dpi))
    tables = result.get("tables") or []
    for table in tables:
        if not isinstance(table, dict):
            errors.append("tables 里有非 dict 条目")
            continue
        for key in ("id", "title", "columns", "rows"):
            if key not in table:
                errors.append("表格 %s 缺少键 %s" % (table.get("id"), key))
        if not isinstance(table.get("rows"), list):
            errors.append("表格 %s 的 rows 不是 list" % table.get("id"))
        for column in table.get("columns") or []:
            if not isinstance(column, dict) or "key" not in column:
                errors.append("表格 %s 的 column 定义不合法" % table.get("id"))
    # obs_extra 校验
    obs_file = result.get("obs_file")
    obs_notes = []
    if obs_file:
        path = Path(str(obs_file))
        if not path.is_absolute():
            path = (out_dir / path).resolve()
        if not path.exists():
            errors.append("obs_file 不存在：%s" % path)
        else:
            try:
                import pandas as pd

                frame = pd.read_csv(path, index_col=0)
                barcodes = (Path(bundle_dir) / "barcodes.tsv").read_text(encoding="utf-8").splitlines()
                barcodes = [item for item in barcodes if item]
                unknown = [item for item in frame.index.astype(str) if item not in set(barcodes)]
                if unknown:
                    errors.append("obs_extra.csv 有 %d 个细胞名不在 barcode 中（示例：%s）"
                                  % (len(unknown), unknown[:2]))
                if frame.shape[1] == 0:
                    errors.append("obs_extra.csv 没有任何新增列")
                numeric_ok = 0
                for column in frame.columns:
                    values = pd.to_numeric(frame[column], errors="coerce")
                    if values.notna().sum() >= max(1, int(0.5 * len(values))):
                        numeric_ok += 1
                if numeric_ok == 0:
                    errors.append("obs_extra.csv 没有数值列")
                obs_notes.append("%d 行 × %d 列" % frame.shape)
            except Exception as error:
                errors.append("obs_extra.csv 校验失败：%s" % error)
    else:
        # CellPhoneDB 是簇级分析，不产出 obs_extra
        if tool["id"] not in ("cellphonedb", "regvelo"):
            errors.append("obs_file 缺失（该工具应输出 obs_extra.csv）")
    digest = {
        "tool": tool["id"],
        "seconds": round(elapsed, 2),
        "n_tables": len(tables),
        "table_rows": [len(table.get("rows") or []) for table in tables],
        "n_figures": len(figures),
        "figures": figure_notes,
        "obs": obs_notes,
        "metrics": summary.get("metrics", {}),
        "warnings": len(result.get("warnings") or []),
        "notes": len(result.get("notes") or []),
    }
    return (not errors), digest, errors


# ---------------------------------------------------------------------------
# 单个工具执行
# ---------------------------------------------------------------------------


def adapter_env():
    env = os.environ.copy()
    env.setdefault("DS_SC_SCENIC_DB", str(SCENIC_DB_CACHE.parent / "cistarget"))
    if SCENIC_DB_CACHE.exists():
        env["DS_SC_SCENIC_ANNOTATIONS_CACHE"] = str(SCENIC_DB_CACHE)
    if CPDB_DB.exists():
        env["DS_SC_CELLPHONEDB_DB"] = str(CPDB_DB)
    if CELLORACLE_BASE_GRN.exists():
        env["DS_SC_CELLORACLE_BASE_GRN"] = str(CELLORACLE_BASE_GRN)
    env.setdefault("OPENBLAS_NUM_THREADS", "2")
    env.setdefault("OMP_NUM_THREADS", "2")
    env.setdefault("MKL_NUM_THREADS", "2")
    return env


def run_tool(tool, bundle_dir, out_dir, versions, timeout_scale=1.0):
    script = EXT_DIR / tool["script"]
    if not script.exists():
        fail("%s: 适配器脚本不存在 %s" % (tool["id"], script))
        return {"status": "FAIL", "reason": "适配器脚本不存在"}
    missing = [module for module in tool["modules"] if not versions.get(module)]
    if missing:
        reason = "依赖未安装：%s" % ", ".join(missing)
        if tool.get("expect_unavailable"):
            # 即使依赖装不上，也要验证适配器"给出中文退出码 3"的契约
            out_dir.mkdir(parents=True, exist_ok=True)
            command = [
                CONDA_PY, str(script), "--bundle", str(bundle_dir),
                "--out", str(out_dir), "--params", json.dumps(tool["params"], ensure_ascii=False),
            ]
            result = run_command(command, timeout=300, env=adapter_env())
            message = (result["stderr"] or "").strip()
            if result["returncode"] == 3 and ("[prerequisite]" in message or "无法安装" in message):
                ok("%s: 依赖未安装（%s）；适配器正确给出中文退出码 3：%s"
                   % (tool["id"], ", ".join(missing), message.splitlines()[-1][:180]))
                return {"status": "PASS(prereq)", "reason": message.splitlines()[-1][:300],
                        "seconds": round(result["elapsed"], 2)}
            fail("%s: 依赖未安装（%s），但适配器没有给出退出码 3 + 中文原因（退出码 %s）"
                 % (tool["id"], ", ".join(missing), result["returncode"]))
            return {"status": "FAIL", "reason": "适配器未按契约报缺前置", "seconds": round(result["elapsed"], 2)}
        skip("%s: %s" % (tool["id"], reason))
        return {"status": "SKIP", "reason": reason}

    out_dir.mkdir(parents=True, exist_ok=True)
    command = [
        CONDA_PY, str(script),
        "--bundle", str(bundle_dir),
        "--out", str(out_dir),
        "--params", json.dumps(tool["params"], ensure_ascii=False),
    ]
    timeout = int(float(tool["timeout"]) * float(timeout_scale))
    info("运行 %s：%s" % (tool["id"], " ".join(command[-2:])))
    result = run_command(command, timeout=timeout, env=adapter_env())
    stdout, stderr, returncode = result["stdout"], result["stderr"], result["returncode"]
    log_path = out_dir / "adapter_stdout_stderr.log"
    log_path.write_text(
        "===== returncode: %s elapsed: %.1fs =====\n----- stdout -----\n%s\n----- stderr -----\n%s\n"
        % (returncode, result["elapsed"], stdout, stderr),
        encoding="utf-8",
    )
    if result["timeout"]:
        fail("%s: 超时（>%ds）" % (tool["id"], timeout))
        return {"status": "FAIL", "reason": "timeout", "seconds": round(result["elapsed"], 2)}
    if returncode == 3:
        reason = (stderr or "").strip().splitlines()[-1] if (stderr or "").strip() else ""
        if tool.get("expect_prereq") or tool.get("expect_unavailable"):
            ok("%s: 正确给出缺前置（退出码 3）：%s" % (tool["id"], reason[:180]))
            return {"status": "PASS(prereq)", "reason": reason, "seconds": round(result["elapsed"], 2)}
        fail("%s: 意外退出码 3（缺前置）：%s" % (tool["id"], reason[:300]))
        return {"status": "FAIL", "reason": reason, "seconds": round(result["elapsed"], 2)}
    if returncode != 0:
        tail = (stderr or "").strip()[-800:]
        fail("%s: 运行失败（退出码 %s）：%s" % (tool["id"], returncode, tail[-300:]))
        return {"status": "FAIL", "reason": tail, "seconds": round(result["elapsed"], 2)}

    passed, digest, errors = validate_result(out_dir, bundle_dir, tool, result["elapsed"])
    if not passed:
        for error in errors:
            print("        - " + str(error))
        fail("%s: result.json/图/obs_extra 校验未通过（%d 项）" % (tool["id"], len(errors)))
        return {"status": "FAIL", "reason": "; ".join(errors), "seconds": round(result["elapsed"], 2)}
    ok("%s: PASS（%.1fs；%d 表 / %d 图 / %s）"
       % (tool["id"], result["elapsed"], digest["n_tables"], digest["n_figures"],
          "、".join(digest["obs"]) if digest["obs"] else "无 obs_extra"))
    return {"status": "PASS", "digest": digest, "seconds": round(result["elapsed"], 2)}


def run_scvelo_synthetic(bundle_source, timeout_scale=1.0):
    """scVelo 正路径冒烟：合成 spliced/unspliced 层后真跑 velocity。"""
    tool = dict(next(item for item in TOOLS if item["id"] == "scvelo"))
    tool["params"] = dict(tool["params"])
    build_layered_bundle(bundle_source, LAYERED_BUNDLE)
    out_dir = OUT_ROOT / "scvelo_synthetic"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    result = run_tool(tool, LAYERED_BUNDLE, out_dir, {"scvelo": "installed"}, timeout_scale=timeout_scale)
    if result["status"] == "PASS":
        result["digest"]["note"] = "合成剪接层：仅验证代码路径，不代表真实 RNA velocity 生物学结论"
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description="dsh-singlecell 官方 Python 工具适配器测试（CONTRACT §10）")
    parser.add_argument("--tools", default="", help="逗号分隔的工具 id；默认全部")
    parser.add_argument("--bundle", default="", help="复用外部 bundle 目录（默认自动生成 demo bundle）")
    parser.add_argument("--force-bundle", action="store_true", help="强制重建 demo bundle")
    parser.add_argument("--scvelo-synthetic", action="store_true", help="额外跑 scVelo 正路径冒烟（合成层）")
    parser.add_argument("--timeout-scale", type=float, default=1.0, help="超时倍率（调试用）")
    parser.add_argument("--only-bundle", action="store_true", help="只生成/检查 bundle，不跑工具")
    args = parser.parse_args(argv)

    maybe_reexec_conda()
    info("externalPython = %s" % CONDA_PY)
    info("插件根目录 = %s" % PLUGIN_ROOT)
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    # 1) demo bundle
    external_bundle = Path(os.path.expanduser(args.bundle)) if args.bundle else None
    if external_bundle:
        if not (external_bundle / "barcodes.tsv").exists():
            fail("--bundle 指定的目录不是合法 bundle：%s" % external_bundle)
            return 1
        bundle_dir = external_bundle
        bundle_meta = json.loads((external_bundle / "meta.json").read_text(encoding="utf-8")) if (
            external_bundle / "meta.json").exists() else {}
        info("使用外部 bundle：%s" % bundle_dir)
    else:
        bundle_dir = DEFAULT_BUNDLE
        bundle_meta = build_bundle(bundle_dir, force=args.force_bundle)
    info("bundle：%d 细胞 × %d 基因，X_kind=%s，簇数=%s"
         % (bundle_meta.get("n_obs", -1), bundle_meta.get("n_vars", -1),
            bundle_meta.get("X_kind", "?"), bundle_meta.get("n_clusters", "?")))
    if args.only_bundle:
        ok("bundle 生成/复用完成")
        return 0

    # 2) 依赖探测
    all_modules = sorted({module for tool in TOOLS for module in tool["modules"]})
    versions = probe_modules(all_modules)
    info("依赖探测：")
    for module in all_modules:
        version = versions.get(module)
        if version:
            print("        %-14s %s" % (module, version))
        else:
            print("        %-14s %s" % (module, color("未安装", "33")))
    for label, path in [("SCENIC motif 库", SCENIC_DB_CACHE), ("CellPhoneDB 数据库", CPDB_DB),
                        ("CellOracle base GRN", CELLORACLE_BASE_GRN)]:
        print("        %-14s %s" % (label, path if path.exists() else color("未找到：%s" % path, "33")))

    # 3) 逐个跑
    selected = [item.strip() for item in args.tools.split(",") if item.strip()]
    tools = [tool for tool in TOOLS if not selected or tool["id"] in selected]
    unknown = set(selected) - {tool["id"] for tool in TOOLS}
    if unknown:
        fail("未知工具 id：%s" % ", ".join(sorted(unknown)))
        return 1
    results = {}
    started = time.time()
    for tool in tools:
        out_dir = OUT_ROOT / tool["id"]
        if out_dir.exists():
            shutil.rmtree(out_dir)
        try:
            results[tool["id"]] = run_tool(tool, bundle_dir, out_dir, versions, timeout_scale=args.timeout_scale)
        except Exception as error:
            traceback.print_exc()
            fail("%s: 测试框架异常：%s" % (tool["id"], error))
            results[tool["id"]] = {"status": "FAIL", "reason": str(error)}
    if args.scvelo_synthetic and (not selected or "scvelo" in selected):
        try:
            results["scvelo_synthetic"] = run_scvelo_synthetic(bundle_dir, timeout_scale=args.timeout_scale)
        except Exception as error:
            traceback.print_exc()
            fail("scvelo_synthetic: 测试框架异常：%s" % error)
            results["scvelo_synthetic"] = {"status": "FAIL", "reason": str(error)}

    # 4) 汇总
    total = time.time() - started
    passed = [key for key, value in results.items() if value["status"].startswith("PASS")]
    skipped = [key for key, value in results.items() if value["status"] == "SKIP"]
    failed = [key for key, value in results.items() if value["status"] == "FAIL"]
    print("")
    print("=" * 72)
    print("测试汇总（总计 %.1fs）" % total)
    for key in results:
        value = results[key]
        status = value["status"]
        seconds = value.get("seconds")
        line = "  %-18s %-13s %s" % (key, status, ("%.1fs" % seconds) if seconds else "")
        if value.get("reason") and status in ("FAIL", "SKIP"):
            line += "  " + str(value["reason"])[:160]
        if status == "FAIL":
            print(color(line, "31"))
        elif status == "SKIP":
            print(color(line, "33"))
        else:
            print(color(line, "32"))
    print("  PASS: %d    SKIP: %d    FAIL: %d" % (len(passed), len(skipped), len(failed)))
    print("=" * 72)
    if failed:
        fail("存在 FAIL：%s" % ", ".join(failed))
        return 1
    ok("全部通过（无 FAIL；SKIP 均为本机安装限制，已在报告说明）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
