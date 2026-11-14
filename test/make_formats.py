#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench — CONTRACT §8.5 导入格式 fixture 生成器（sc-io 所有）

在临时目录生成全格式导入层测试所需的全部 fixture，并在每个 fixture 目录里写一份
truth.json（n_obs / n_vars / 前几个基因 / 首细胞某基因值 / obs 期望列）。

生成清单（§8.2 / §8.5）：
  h5ad、h5ad(layers+counts)、h5ad.gz、10x_h5 v3(feature_type 过滤)、10x_h5 v2(genes/gene_names)、
  10x_mtx 三件套目录(v3 features.tsv)、10x_mtx v2(genes.tsv)、单个 matrix.mtx.gz(有 sidecar)、
  裸 matrix.mtx(无 sidecar)、loom、dense csv、tsv.gz、分号 txt、无表头 tsv、方向探针 csv、
  GEO series_matrix.txt.gz、zip(10x 三件套)、tar.gz(csv)、symlink tar.gz、
  多候选 zip、路径穿越 zip、空 zip、坏文件（坏 h5ad/h5/rds、空 csv、不存在路径）、
  Seurat rds、SCE rds、非 Seurat matrix rds、单对象 RData、多对象 RData、混合对象 RData。

用法：
    /Library/Frameworks/Python.framework/Versions/3.14/bin/python3 test/make_formats.py
    ... --out /tmp/my-formats
    ... --serve            # 生成后起本地 http.server 供手测 URL 导入

作为库使用（test/py_import_test.py 就是这么用的）：
    import make_formats
    manifest = make_formats.generate(root)     # 返回 manifest dict
    with make_formats.serve(root) as base_url: # 本地 http.server
        ...
"""

from __future__ import annotations

import argparse
import contextlib
import functools
import gzip
import http.server
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# 常量：所有 fixture 共用同一套细胞/基因/计数，truth 可交叉校验
# ---------------------------------------------------------------------------

GENES = ["GENEA", "GENEB", "GENEC", "GENED", "GENEE"]
CELLS = [
    "AAACCTGAGAAACCTA-1",
    "AAACCTGAGAAACCTC-1",
    "AAACCTGAGAAACCTG-1",
    "AAACCTGAGAAACCTT-1",
    "AAACCTGAGAAACCTA-2",
    "AAACCTGAGAAACCTC-2",
]

GENES_B = ["GENEX", "GENEY", "GENEZ"]
CELLS_B = ["TTTGGGAAACCCAAA-1", "TTTGGGAAACCCCCC-1", "TTTGGGAAACCCGGG-1", "TTTGGGAAACCCTTT-1"]

RSCRIPT_CANDIDATES = [
    "/usr/local/bin/Rscript",
    "/Library/Frameworks/R.framework/Resources/bin/Rscript",
]


def counts_cells_genes(cells=None, genes=None) -> np.ndarray:
    """值刻意做成 10*(i+1)+(j+1)：任何转置/错位都会被 truth 校验抓到。"""
    cells = CELLS if cells is None else cells
    genes = GENES if genes is None else genes
    return np.array(
        [[float(10 * (i + 1) + (j + 1)) for j in range(len(genes))] for i in range(len(cells))],
        dtype=np.float32,
    )


COUNTS_CG = counts_cells_genes()             # cells × genes（6×5）
COUNTS_GC = COUNTS_CG.T.copy()               # genes × cells（5×6）
COUNTS_B_CG = counts_cells_genes(CELLS_B, GENES_B)  # 4×3


def find_rscript():
    for candidate in RSCRIPT_CANDIDATES:
        if Path(candidate).exists():
            return candidate
    found = shutil.which("Rscript")
    return found


# ---------------------------------------------------------------------------
# truth.json 工具
# ---------------------------------------------------------------------------


def _truth_payload(
    fid: str,
    fmt: str,
    *,
    n_obs: int,
    n_vars: int,
    genes,
    cell_name=None,
    gene_name=None,
    value=None,
    counts_value=None,
    obs_columns=None,
    obs_columns_any=None,
    notes="",
):
    payload = {
        "id": fid,
        "format": fmt,
        "n_obs": int(n_obs),
        "n_vars": int(n_vars),
        "genes": [str(g) for g in list(genes)[:3]],
        "cell_name": cell_name,
        "gene_name": gene_name,
        "value": None if value is None else float(value),
        "counts_value": None if counts_value is None else float(counts_value),
        "obs_columns": list(obs_columns or []),
        "obs_columns_any": list(obs_columns_any or []),
        "notes": notes,
    }
    return payload


def _write_truth(fixture_dir: Path, payload: dict) -> str:
    fixture_dir.mkdir(parents=True, exist_ok=True)
    path = fixture_dir / "truth.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return str(path)


def _mk(root: Path, fid: str) -> Path:
    path = root / fid
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# 各格式 writer
# ---------------------------------------------------------------------------


def write_h5ad(path: Path, *, layers: bool = False):
    """AnnData：layers=True 时 X 已标准化、layers['counts'] 存原始计数。"""
    import anndata as ad

    obs = pd.DataFrame(
        {"batch": ["ctrl", "ctrl", "ctrl", "trt", "trt", "trt"],
         "score": [1.5, 2.5, 3.5, 4.5, 5.5, 6.5]},
        index=pd.Index(CELLS, name="cell"),
    )
    var = pd.DataFrame({"gene_id": ["ENSG%05d" % (i + 1) for i in range(len(GENES))]},
                       index=pd.Index(GENES, name="gene"))
    if layers:
        x = sp.csr_matrix(COUNTS_CG * 0.5)
        adata = ad.AnnData(X=x, obs=obs, var=var)
        adata.layers["counts"] = sp.csr_matrix(COUNTS_CG)
    else:
        adata = ad.AnnData(X=sp.csr_matrix(COUNTS_CG), obs=obs, var=var)
    adata.write_h5ad(str(path))


def write_10x_h5_v3(path: Path):
    """CellRanger v3：/matrix/{data,indices,indptr,shape,barcodes,features/*}，shape=genes×cells。

    features 里放 5 个 Gene Expression + 2 个 Antibody Capture，用于验证 feature_type 过滤。
    """
    import h5py

    genes = GENES + ["AB_CD3", "AB_CD19"]
    values = np.zeros((len(genes), len(CELLS)), dtype=np.float32)
    values[: len(GENES), :] = COUNTS_GC
    values[len(GENES):, :] = np.arange(len(CELLS), dtype=np.float32)[None, :] + 900.0
    spm = sp.csc_matrix(values)  # CellRanger 约定：genes×cells 的 CSC

    def sarr(items):
        return np.array([str(x).encode("utf-8") for x in items], dtype=f"S{max(len(str(x).encode('utf-8')) for x in items)}")

    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("barcodes", data=sarr(CELLS))
        matrix.create_dataset("data", data=spm.data.astype(np.int32))
        matrix.create_dataset("indices", data=spm.indices.astype(np.int32))
        matrix.create_dataset("indptr", data=spm.indptr.astype(np.int32))
        matrix.create_dataset("shape", data=np.array(spm.shape, dtype=np.int64))
        features = matrix.create_group("features")
        features.create_dataset("id", data=sarr(["ENSG%05d" % (i + 1) for i in range(len(genes))]))
        features.create_dataset("name", data=sarr(genes))
        features.create_dataset("feature_type",
                                data=sarr(["Gene Expression"] * len(GENES) + ["Antibody Capture"] * 2))
        features.create_dataset("genome", data=sarr(["GRCh38"] * len(genes)))


def write_10x_h5_v2(path: Path):
    """CellRanger v2：/matrix/{data,indices,indptr,shape,barcodes,genes/{gene_ids,gene_names}}。"""
    import h5py

    spm = sp.csc_matrix(COUNTS_GC.astype(np.float32))  # v2 同样是 CSC

    def sarr(items):
        return np.array([str(x).encode("utf-8") for x in items], dtype=f"S{max(len(str(x).encode('utf-8')) for x in items)}")

    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("barcodes", data=sarr(CELLS))
        matrix.create_dataset("data", data=spm.data.astype(np.int32))
        matrix.create_dataset("indices", data=spm.indices.astype(np.int32))
        matrix.create_dataset("indptr", data=spm.indptr.astype(np.int32))
        matrix.create_dataset("shape", data=np.array(spm.shape, dtype=np.int64))
        genes = matrix.create_group("genes")
        genes.create_dataset("gene_ids", data=sarr(["ENSG%05d" % (i + 1) for i in range(len(GENES))]))
        genes.create_dataset("gene_names", data=sarr(GENES))


def write_mtx_dir(path: Path, *, v2: bool = False, gz: bool = True):
    """10x mtx 三件套。v2 用 genes.tsv（2 列），v3 用 features.tsv（3 列，第 2 列基因名）。"""
    path.mkdir(parents=True, exist_ok=True)
    spm = sp.csr_matrix(COUNTS_GC.astype(np.float64))

    raw_mtx = path / "matrix.mtx"
    scipy.io.mmwrite(str(raw_mtx), spm)
    if gz:
        with open(raw_mtx, "rb") as src, gzip.open(path / "matrix.mtx.gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
        raw_mtx.unlink()

    barcodes = "\n".join(CELLS) + "\n"
    if gz:
        with gzip.open(path / "barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
            handle.write(barcodes)
    else:
        (path / "barcodes.tsv").write_text(barcodes, encoding="utf-8")

    if v2:
        lines = "\n".join("ENSG%05d\t%s" % (i + 1, gene) for i, gene in enumerate(GENES)) + "\n"
        name = "genes.tsv"
    else:
        lines = "\n".join("ENSG%05d\t%s\tGene Expression" % (i + 1, gene) for i, gene in enumerate(GENES)) + "\n"
        name = "features.tsv"
    target = path / (name + ".gz" if gz else name)
    if gz:
        with gzip.open(target, "wt", encoding="utf-8") as handle:
            handle.write(lines)
    else:
        target.write_text(lines, encoding="utf-8")


def write_loom(path: Path):
    """loom：/matrix 是 genes×cells；/row_attrs/Gene；/col_attrs/* 全部转 obs。"""
    import h5py

    def sarr(items):
        return np.array([str(x).encode("utf-8") for x in items], dtype=f"S{max(len(str(x).encode('utf-8')) for x in items)}")

    with h5py.File(path, "w") as handle:
        handle.create_dataset("matrix", data=COUNTS_GC.astype(np.float32))
        row_attrs = handle.create_group("row_attrs")
        row_attrs.create_dataset("Gene", data=sarr(GENES))
        col_attrs = handle.create_group("col_attrs")
        col_attrs.create_dataset("CellID", data=sarr(CELLS))
        col_attrs.create_dataset("batch", data=sarr(["ctrl", "ctrl", "ctrl", "trt", "trt", "trt"]))
        col_attrs.create_dataset("n_counts", data=np.arange(len(CELLS), dtype=np.int32) + 100)
        attrs = handle.create_group("attrs")
        attrs.create_dataset("title", data=sarr(["test loom"]))


def write_dense_csv(path: Path, frame: pd.DataFrame, *, sep: str = ",", header: bool = True, compression=None):
    frame.to_csv(path, sep=sep, header=header, compression=compression, encoding="utf-8")


def write_geo_series(path: Path):
    """GEO series matrix（gsm 为列、基因为行）+ !Sample_* 元数据。"""
    samples = ["GSM5000001", "GSM5000002", "GSM5000003", "GSM5000004"]
    titles = ["sample1", "sample2", "sample3", "sample4"]
    chars = ["cell type: T cell", "cell type: B cell", "cell type: T cell", "cell type: monocyte"]
    sample_counts = np.array(
        [[float(10 * (s + 1) + (g + 1)) for g in range(len(GENES))] for s in range(len(samples))],
        dtype=np.float32,
    )  # samples × genes

    def quoted(values):
        return "\t".join('"%s"' % v for v in values)

    lines = [
        '!Series_title\t"GSE99999 test series matrix"',
        '!Series_geo_accession\t"GSE99999"',
        '!Series_status\t"Public on Nov 06, 2026"',
        '!Series_submission_date\t"Nov 01, 2026"',
        '!Series_organism_ch1\t"Homo sapiens"',
        '!Series_summary\t"synthetic fixture for dsh-singlecell-workbench import tests"',
        "!Sample_title\t" + quoted(titles),
        "!Sample_geo_accession\t" + quoted(samples),
        "!Sample_organism_ch1\t" + quoted(["Homo sapiens"] * len(samples)),
        "!Sample_characteristics_ch1\t" + quoted(chars),
        "!Sample_library_strategy\t" + quoted(["RNA-Seq"] * len(samples)),
        "!series_matrix_table_begin",
        '"ID_REF"\t' + quoted(samples),
    ]
    for gene_index, gene in enumerate(GENES):
        lines.append('"%s"\t' % gene + "\t".join("%.1f" % sample_counts[s, gene_index] for s in range(len(samples))))
    lines.append("!series_matrix_table_end")
    text = "\n".join(lines) + "\n"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _zip_write(zf: zipfile.ZipFile, name: str, data):
    if not isinstance(data, (bytes, bytearray)):
        data = str(data).encode("utf-8")
    zf.writestr(name, data)


def _mtx_text(spm: sp.spmatrix) -> str:
    buffer = io.BytesIO()
    scipy.io.mmwrite(buffer, spm)
    return buffer.getvalue().decode("utf-8")


def _csv_text(frame: pd.DataFrame) -> str:
    return frame.to_csv(sep=",", header=True, encoding="utf-8")


def write_zip_tenx(path: Path):
    """zip 里是 10x 三件套目录（应被识别成 10x_mtx）。"""
    spm = sp.csr_matrix(COUNTS_GC.astype(np.float64))
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        _zip_write(zf, "bundle/matrix.mtx", _mtx_text(spm))
        _zip_write(zf, "bundle/barcodes.tsv", "\n".join(CELLS) + "\n")
        _zip_write(zf, "bundle/features.tsv",
                   "\n".join("ENSG%05d\t%s\tGene Expression" % (i + 1, g) for i, g in enumerate(GENES)) + "\n")
        _zip_write(zf, "bundle/README.txt", "10x bundle fixture\n")


def write_targz_csv(path: Path):
    frame = pd.DataFrame(COUNTS_CG, index=CELLS, columns=GENES)
    payload = _csv_text(frame).encode("utf-8")
    with tarfile.open(path, "w:gz") as tf:
        info = tarfile.TarInfo("export/expr.csv")
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))
        readme = b"tar.gz fixture\n"
        info2 = tarfile.TarInfo("export/README.md")
        info2.size = len(readme)
        tf.addfile(info2, io.BytesIO(readme))


def write_targz_symlink(path: Path):
    """tar.gz 里有正常 csv + 一个指向外部的符号链接（应跳过链接、仍能读到 csv）。"""
    frame = pd.DataFrame(COUNTS_CG, index=CELLS, columns=GENES)
    payload = _csv_text(frame).encode("utf-8")
    with tarfile.open(path, "w:gz") as tf:
        info = tarfile.TarInfo("data/expr.csv")
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))
        link = tarfile.TarInfo("data/evil_link.csv")
        link.type = tarfile.SYMTYPE
        link.linkname = "../outside-secret.csv"
        tf.addfile(link)


def write_zip_multi(path: Path):
    """多个候选文件 → ambiguous + candidates（member 相对路径）。"""
    frame_a = pd.DataFrame(COUNTS_CG, index=CELLS, columns=GENES)
    frame_b = pd.DataFrame(COUNTS_B_CG, index=CELLS_B, columns=GENES_B)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        _zip_write(zf, "a/expr_a.csv", _csv_text(frame_a))
        _zip_write(zf, "b/expr_b.csv", _csv_text(frame_b))
        _zip_write(zf, "readme.md", "two candidate matrices\n")


def write_zip_evil(path: Path):
    """路径穿越：../ 与绝对路径成员；解压必须拒绝，不能写到目标目录外。"""
    frame = pd.DataFrame(COUNTS_CG, index=CELLS, columns=GENES)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        _zip_write(zf, "good/expr.csv", _csv_text(frame))
        _zip_write(zf, "../sc_io_escape_marker.txt", "pwned by ../ member")
        _zip_write(zf, "/tmp/sc_io_escape_abs_marker.txt", "pwned by absolute member")


def write_zip_empty(path: Path):
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        _zip_write(zf, "notes/readme.txt", "no sc matrix inside\n")


def write_r_objects(root: Path, rscript: str) -> dict:
    """用 Rscript 生成 Seurat rds / SCE rds / matrix rds / 各种 RData。"""
    outdirs = {
        "seurat": _mk(root, "r_seurat"),
        "sce": _mk(root, "r_sce"),
        "matrix": _mk(root, "r_matrix"),
        "rdata_multi": _mk(root, "r_rdata_multi"),
        "rdata_single": _mk(root, "r_rdata_single"),
        "rdata_mixed": _mk(root, "r_rdata_mixed"),
    }
    r_code = r"""
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SingleCellExperiment)
})
genes <- c(%(genes)s)
cells <- c(%(cells)s)
genes_b <- c(%(genes_b)s)
cells_b <- c(%(cells_b)s)
counts_cg <- outer(seq_along(cells), seq_along(genes), function(i, j) 10L * i + j)
dimnames(counts_cg) <- list(cells, genes)
counts_gc <- t(counts_cg)
dimnames(counts_gc) <- list(genes, cells)
counts_b_gc <- outer(seq_along(genes_b), seq_along(cells_b), function(g, i) 10L * i + g)
dimnames(counts_b_gc) <- list(genes_b, cells_b)

seu <- CreateSeuratObject(counts = Matrix::Matrix(counts_gc, sparse = TRUE),
                          project = "GSE_TEST", min.cells = 0, min.features = 0)
seu$batch <- c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt")
saveRDS(seu, file.path("%(seurat_dir)s", "seurat.rds"))

sce <- SingleCellExperiment(
  assays = list(counts = Matrix::Matrix(counts_gc, sparse = TRUE)),
  colData = S4Vectors::DataFrame(batch = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"),
                                 score = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5))
)
saveRDS(sce, file.path("%(sce_dir)s", "sce.rds"))

saveRDS(counts_gc, file.path("%(matrix_dir)s", "matrix_genes_cells.rds"))

seu_b <- CreateSeuratObject(counts = Matrix::Matrix(counts_b_gc, sparse = TRUE),
                            project = "GSE_TEST_B", min.cells = 0, min.features = 0)
save(seu, seu_b, file = file.path("%(multi_dir)s", "multi_objects.RData"))
save(seu, file = file.path("%(single_dir)s", "single_seurat.RData"))
df <- as.data.frame(t(counts_b_gc))   # 4 cells x 3 genes，与 Seurat 不同维度
save(seu, df, file = file.path("%(mixed_dir)s", "mixed.RData"))
cat("R fixtures OK\n")
""" % {
        "genes": ", ".join('"%s"' % g for g in GENES),
        "cells": ", ".join('"%s"' % c for c in CELLS),
        "genes_b": ", ".join('"%s"' % g for g in GENES_B),
        "cells_b": ", ".join('"%s"' % c for c in CELLS_B),
        "seurat_dir": str(outdirs["seurat"]),
        "sce_dir": str(outdirs["sce"]),
        "matrix_dir": str(outdirs["matrix"]),
        "multi_dir": str(outdirs["rdata_multi"]),
        "single_dir": str(outdirs["rdata_single"]),
        "mixed_dir": str(outdirs["rdata_mixed"]),
    }
    script_path = root / "_make_r_objects.R"
    script_path.write_text(r_code, encoding="utf-8")
    result = subprocess.run(
        [rscript, "--vanilla", str(script_path)],
        capture_output=True, text=True, timeout=600,
        encoding="utf-8", errors="replace",
    )
    if result.returncode != 0:
        raise RuntimeError("R 生成 fixture 失败（%s）：%s\n%s" % (rscript, result.stdout, result.stderr))
    return {
        "seurat_dir": outdirs["seurat"],
        "sce_dir": outdirs["sce"],
        "matrix_dir": outdirs["matrix"],
        "rdata_multi_dir": outdirs["rdata_multi"],
        "rdata_single_dir": outdirs["rdata_single"],
        "rdata_mixed_dir": outdirs["rdata_mixed"],
    }


# ---------------------------------------------------------------------------
# 主生成逻辑
# ---------------------------------------------------------------------------


def generate(root, *, with_r=True, verbose=True) -> dict:
    """生成全部 fixture，返回 manifest dict。root 会被创建（不删除已有内容）。"""
    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)

    def say(message):
        if verbose:
            print("[make_formats] " + message)

    formats = []
    ambiguous = []
    errors = []
    url_tests = []
    transpose_probes = []

    def add_format(fid, fmt, path, truth, *, options=None, inspect_format=None, notes=""):
        fixture_dir = root / fid
        truth_path = _write_truth(fixture_dir, truth)
        entry = {
            "id": fid,
            "format": fmt,
            "path": str(path),
            "options": dict(options or {}),
            "truth": truth,
            "truth_path": truth_path,
            "inspect_format": inspect_format or fmt,
            "notes": notes,
        }
        formats.append(entry)
        return entry

    # ---- 1. h5ad（基础 / layers + counts / .h5ad.gz） ----
    d = _mk(root, "h5ad_simple")
    write_h5ad(d / "simple.h5ad")
    truth = _truth_payload("h5ad_simple", "h5ad", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           obs_columns=["batch", "score"])
    add_format("h5ad_simple", "h5ad", d / "simple.h5ad", truth, notes="AnnData 标准格式")

    d = _mk(root, "h5ad_layers")
    write_h5ad(d / "layers.h5ad", layers=True)
    truth = _truth_payload("h5ad_layers", "h5ad", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=5.5,
                           counts_value=11.0, obs_columns=["batch", "score"],
                           notes="X=标准化值，layers['counts'] 为原始计数")
    add_format("h5ad_layers", "h5ad", d / "layers.h5ad", truth, notes="layers['counts'] 分支")

    d = _mk(root, "h5ad_gz")
    write_h5ad(d / "simple.h5ad")
    with open(d / "simple.h5ad", "rb") as src, gzip.open(d / "simple.h5ad.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    (d / "simple.h5ad").unlink()
    truth = _truth_payload("h5ad_gz", "h5ad", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           obs_columns=["batch", "score"], notes="h5ad gzip 单文件")
    add_format("h5ad_gz", "h5ad", d / "simple.h5ad.gz", truth, notes="单文件 .h5ad.gz")

    # ---- 2. 10x_h5（v3 feature_type 过滤 / v2 genes） ----
    d = _mk(root, "tenx_h5_v3")
    write_10x_h5_v3(d / "filtered_feature_bc_matrix.h5")
    truth = _truth_payload("tenx_h5_v3", "10x_h5", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           notes="CellRanger v3，features 里 5 GE + 2 Antibody Capture，只保留 GE")
    add_format("tenx_h5_v3", "10x_h5", d / "filtered_feature_bc_matrix.h5", truth,
               notes="feature_type 过滤")

    d = _mk(root, "tenx_h5_v2")
    write_10x_h5_v2(d / "v2.h5")
    truth = _truth_payload("tenx_h5_v2", "10x_h5", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           notes="CellRanger v2，genes/gene_names")
    add_format("tenx_h5_v2", "10x_h5", d / "v2.h5", truth, notes="v2 gene_names 分支")

    # ---- 3. 10x_mtx（三件套目录 v3 / v2 / 单 .mtx.gz / 裸 .mtx） ----
    d = _mk(root, "tenx_mtx_dir")
    write_mtx_dir(d / "filtered_feature_bc_matrix", gz=True)
    truth = _truth_payload("tenx_mtx_dir", "10x_mtx", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("tenx_mtx_dir", "10x_mtx", d / "filtered_feature_bc_matrix", truth,
               notes="三件套目录 + v3 features.tsv（第 2 列基因名）")

    d = _mk(root, "tenx_mtx_v2_dir")
    write_mtx_dir(d / "genes_tsv_dir", v2=True, gz=True)
    truth = _truth_payload("tenx_mtx_v2_dir", "10x_mtx", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("tenx_mtx_v2_dir", "10x_mtx", d / "genes_tsv_dir", truth, notes="v2 genes.tsv")

    d = _mk(root, "tenx_mtx_single")
    write_mtx_dir(d, gz=True)
    truth = _truth_payload("tenx_mtx_single", "10x_mtx", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("tenx_mtx_single", "10x_mtx", d / "matrix.mtx.gz", truth,
               notes="只给 matrix.mtx.gz，同目录 sidecar 可用")

    d = _mk(root, "tenx_mtx_bare")
    scipy.io.mmwrite(str(d / "matrix.mtx"), sp.csr_matrix(COUNTS_GC.astype(np.float64)))
    truth = _truth_payload("tenx_mtx_bare", "10x_mtx", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           notes="无 sidecar，细胞/基因名自动生成")
    truth["names_generated"] = True
    add_format("tenx_mtx_bare", "10x_mtx", d / "matrix.mtx", truth, notes="无 sidecar 自动生成名")

    # ---- 4. loom ----
    d = _mk(root, "loom_simple")
    write_loom(d / "simple.loom")
    truth = _truth_payload("loom_simple", "loom", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           obs_columns=["CellID", "batch"])
    add_format("loom_simple", "loom", d / "simple.loom", truth, notes="genes×cells 转置 + col_attrs")

    # ---- 5. 文本矩阵（csv / tsv.gz / 分号 / 无表头 / 方向探针） ----
    d = _mk(root, "text_csv")
    write_dense_csv(d / "counts.csv", pd.DataFrame(COUNTS_CG, index=CELLS, columns=GENES))
    truth = _truth_payload("text_csv", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("text_csv", "matrix_text", d / "counts.csv", truth, notes="cells×genes，行名 barcode")

    d = _mk(root, "text_tsv_gz")
    write_dense_csv(d / "counts.tsv.gz", pd.DataFrame(COUNTS_GC, index=GENES, columns=CELLS),
                    sep="\t", compression="gzip")
    truth = _truth_payload("text_tsv_gz", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("text_tsv_gz", "matrix_text", d / "counts.tsv.gz", truth,
               notes="genes×cells，应自动转置")

    d = _mk(root, "text_semicolon")
    frame = pd.DataFrame(COUNTS_CG, index=CELLS, columns=GENES)
    (d / "counts_semicolon.txt").write_text(frame.to_csv(sep=";", header=True), encoding="utf-8")
    truth = _truth_payload("text_semicolon", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("text_semicolon", "matrix_text", d / "counts_semicolon.txt", truth,
               notes="分号分隔嗅探")

    d = _mk(root, "text_noheader")
    lines = ["\t".join([CELLS[i]] + ["%d" % int(COUNTS_CG[i, j]) for j in range(len(GENES))])
             for i in range(len(CELLS))]
    (d / "counts_noheader.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    truth = _truth_payload("text_noheader", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           notes="无表头，首列是行名；列名由实现自动生成")
    truth["names_generated"] = True
    add_format("text_noheader", "matrix_text", d / "counts_noheader.tsv", truth,
               notes="无表头 + 首列行名（列名自动生成）")

    d = _mk(root, "text_transpose_probe")
    probe_cg = np.array([[101.0 + 10 * i + j for j in range(3)] for i in range(5)], dtype=np.float32)
    probe = pd.DataFrame(probe_cg, index=["row_%d" % (i + 1) for i in range(5)],
                         columns=["col_%d" % (j + 1) for j in range(3)])
    write_dense_csv(d / "probe.csv", probe)
    truth = _truth_payload("text_transpose_probe", "matrix_text", n_obs=5, n_vars=3,
                           genes=["col_1", "col_2", "col_3"], notes="transpose 强制参数探针（5×3）")
    probe_entry = {
        "id": "text_transpose_probe",
        "format": "matrix_text",
        "path": str(d / "probe.csv"),
        "options": {},
        "truth": truth,
        "truth_path": _write_truth(d, truth),
        "shape_rows": 5,
        "shape_cols": 3,
    }
    transpose_probes.append(probe_entry)

    # ---- 6. GEO series matrix ----
    d = _mk(root, "geo_series")
    write_geo_series(d / "GSE99999_series_matrix.txt.gz")
    truth = _truth_payload("geo_series", "geo_series", n_obs=4, n_vars=5, genes=GENES,
                           cell_name="sample1", gene_name=GENES[0], value=11.0,
                           obs_columns_any=["geo_accession", "title", "characteristics"],
                           notes="4 sample × 5 gene；!Sample_* 元数据进 obs")
    add_format("geo_series", "geo_series", d / "GSE99999_series_matrix.txt.gz", truth,
               notes="GEO series matrix")

    # ---- 7. 压缩包 ----
    d = _mk(root, "archive_tenx_zip")
    write_zip_tenx(d / "bundle.zip")
    truth = _truth_payload("archive_tenx_zip", "10x_mtx", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("archive_tenx_zip", "10x_mtx", d / "bundle.zip", truth,
               inspect_format="archive", notes="zip 内 10x 三件套目录")

    d = _mk(root, "archive_csv_targz")
    write_targz_csv(d / "export.tar.gz")
    truth = _truth_payload("archive_csv_targz", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    add_format("archive_csv_targz", "matrix_text", d / "export.tar.gz", truth,
               inspect_format="archive", notes="tar.gz 内单个 csv")

    d = _mk(root, "archive_symlink_targz")
    write_targz_symlink(d / "with_symlink.tar.gz")
    truth = _truth_payload("archive_symlink_targz", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                           cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                           notes="含符号链接，应跳过链接只用正常 csv")
    add_format("archive_symlink_targz", "matrix_text", d / "with_symlink.tar.gz", truth,
               inspect_format="archive", notes="符号链接应被跳过")

    d = _mk(root, "archive_multi_zip")
    write_zip_multi(d / "multi.zip")
    multi_truth_a = _truth_payload("archive_multi_zip::expr_a", "matrix_text", n_obs=6, n_vars=5,
                                   genes=GENES, cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
    multi_truth_b = _truth_payload("archive_multi_zip::expr_b", "matrix_text", n_obs=4, n_vars=3,
                                   genes=GENES_B, cell_name=CELLS_B[0], gene_name=GENES_B[0], value=11.0)
    multi_dir = d
    _write_truth(multi_dir, multi_truth_a)
    ambiguous.append({
        "id": "archive_multi_zip",
        "path": str(d / "multi.zip"),
        "candidates_field": "member",
        "inspect_format": "archive",
        "retry": [
            {"options": {"member": "a/expr_a.csv"}, "truth": multi_truth_a},
            {"options": {"member": "b/expr_b.csv"}, "truth": multi_truth_b},
        ],
        "notes": "两个候选 csv → ambiguous，用 member 指定相对路径",
    })

    # ---- 8. R 对象 ----
    rscript = find_rscript() if with_r else None
    if rscript:
        say("用 %s 生成 R 对象 fixture…" % rscript)
        rdirs = write_r_objects(root, rscript)
        truth = _truth_payload("r_seurat", "rds", n_obs=6, n_vars=5, genes=GENES,
                               cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                               obs_columns=["orig.ident", "batch"], notes="Seurat v5 counts 层")
        add_format("r_seurat", "rds", rdirs["seurat_dir"] / "seurat.rds", truth, notes="Seurat rds")

        truth = _truth_payload("r_sce", "rds", n_obs=6, n_vars=5, genes=GENES,
                               cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                               obs_columns=["batch", "score"], notes="SingleCellExperiment counts + colData")
        add_format("r_sce", "rds", rdirs["sce_dir"] / "sce.rds", truth, notes="SCE rds")

        truth = _truth_payload("r_matrix", "rds", n_obs=6, n_vars=5, genes=GENES,
                               cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                               notes="非 Seurat：genes×cells 纯 matrix")
        add_format("r_matrix", "rds", rdirs["matrix_dir"] / "matrix_genes_cells.rds", truth,
                   notes="非 Seurat matrix rds")

        truth = _truth_payload("r_rdata_single", "rdata", n_obs=6, n_vars=5, genes=GENES,
                               cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                               obs_columns=["orig.ident", "batch"], notes="单对象 RData")
        add_format("r_rdata_single", "rdata", rdirs["rdata_single_dir"] / "single_seurat.RData", truth,
                   notes="单对象 RData")

        truth = _truth_payload("r_rdata_mixed", "rdata", n_obs=6, n_vars=5, genes=GENES,
                               cell_name=CELLS[0], gene_name=GENES[0], value=11.0,
                               obs_columns=["orig.ident"],
                               notes="Seurat + data.frame 混合，应优先 Seurat（不被 4×3 df 干扰）")
        add_format("r_rdata_mixed", "rdata", rdirs["rdata_mixed_dir"] / "mixed.RData", truth,
                   notes="混合对象，应优先 Seurat")

        multi_truth_a = _truth_payload("r_rdata_multi::seu", "rdata", n_obs=6, n_vars=5,
                                       genes=GENES, cell_name=CELLS[0], gene_name=GENES[0], value=11.0)
        multi_truth_b = _truth_payload("r_rdata_multi::seu_b", "rdata", n_obs=4, n_vars=3,
                                       genes=GENES_B, cell_name=CELLS_B[0], gene_name=GENES_B[0], value=11.0)
        _write_truth(rdirs["rdata_multi_dir"], multi_truth_a)
        ambiguous.append({
            "id": "r_rdata_multi",
            "path": str(rdirs["rdata_multi_dir"] / "multi_objects.RData"),
            "candidates_field": "object",
            "inspect_format": "rdata",
            "retry": [
                {"options": {"object": "seu_b"}, "truth": multi_truth_b},
            ],
            "notes": "两个 Seurat 对象 → ambiguous，用 object 指定",
        })
    else:
        say("警告：没找到 Rscript，跳过 R 对象 fixture（rds/SCE/RData）")

    # ---- 9. 坏文件 / 安全测试 ----
    d = _mk(root, "err_broken_h5ad")
    (d / "broken.h5ad").write_bytes(b"this is definitely not an HDF5 file\n" * 20)
    errors.append({"id": "err_broken_h5ad", "path": str(d / "broken.h5ad"),
                   "expect_codes": ["bad_request", "unsupported"],
                   "notes": "伪装成 .h5ad 的垃圾字节"})

    d = _mk(root, "err_broken_h5")
    (d / "broken.h5").write_bytes(b"\x89HDF\r\n\x1a\n" + b"\x00" * 64)
    errors.append({"id": "err_broken_h5", "path": str(d / "broken.h5"),
                   "expect_codes": ["bad_request", "unsupported"],
                   "notes": "HDF5 magic 但内容损坏"})

    d = _mk(root, "err_broken_rds")
    (d / "broken.rds").write_bytes(b"not an RDS serialization at all\n" * 4)
    errors.append({"id": "err_broken_rds", "path": str(d / "broken.rds"),
                   "expect_codes": ["bad_request", "unsupported"],
                   "notes": "垃圾 .rds，R 桥应给中文错误"})

    d = _mk(root, "err_empty_csv")
    (d / "empty.csv").write_bytes(b"")
    errors.append({"id": "err_empty_csv", "path": str(d / "empty.csv"),
                   "expect_codes": ["bad_request", "unsupported"],
                   "notes": "空 csv"})

    d = _mk(root, "err_evil_zip")
    write_zip_evil(d / "evil.zip")
    errors.append({"id": "err_evil_zip", "path": str(d / "evil.zip"),
                   "expect_codes": ["bad_request", "unsupported"],
                   "forbidden_markers": ["sc_io_escape_marker.txt", "sc_io_escape_abs_marker.txt"],
                   "notes": "路径穿越 zip（../ 与绝对路径成员）必须被拒绝"})

    d = _mk(root, "err_empty_zip")
    write_zip_empty(d / "empty.zip")
    errors.append({"id": "err_empty_zip", "path": str(d / "empty.zip"),
                   "expect_codes": ["bad_request", "unsupported"],
                   "notes": "压缩包里没有可识别矩阵"})

    errors.append({
        "id": "err_missing_path",
        "path": str(root / "does_not_exist" / "nope.csv"),
        "expect_codes": ["not_found", "bad_request"],
        "notes": "路径不存在",
    })

    # ---- 10. URL 测试（由 py_import_test 起本地 http.server 提供） ----
    url_tests.append({
        "id": "url_csv",
        "rel": "text_csv/counts.csv",
        "format": "matrix_text",
        "truth": _truth_payload("url_csv", "matrix_text", n_obs=6, n_vars=5, genes=GENES,
                                cell_name=CELLS[0], gene_name=GENES[0], value=11.0),
        "notes": "本地 http.server 上的 CSV 下载",
    })
    url_tests.append({
        "id": "url_zip",
        "rel": "archive_tenx_zip/bundle.zip",
        "format": "10x_mtx",
        "truth": _truth_payload("url_zip", "10x_mtx", n_obs=6, n_vars=5, genes=GENES,
                                cell_name=CELLS[0], gene_name=GENES[0], value=11.0),
        "notes": "本地 http.server 上的 zip 下载后按 10x_mtx 读",
    })
    url_tests.append({
        "id": "url_geo_gz",
        "rel": "geo_series/GSE99999_series_matrix.txt.gz",
        "format": "geo_series",
        "truth": _truth_payload("url_geo_gz", "geo_series", n_obs=4, n_vars=5, genes=GENES,
                                cell_name="sample1", gene_name=GENES[0], value=11.0),
        "notes": "本地 http.server 上的 .gz 下载",
    })

    # ---- 11. 真实 Seurat rds（只读验收；不在 fixture 内生成） ----
    rds_env = os.environ.get("DSH_SC_TEST_RDS_DIR", "").strip()
    real_rds_dir = Path(rds_env) if rds_env else (Path.home() / ".dsh" / "sc-test-rds")
    real_rds = []
    if real_rds_dir.is_dir():
        real_rds = sorted(str(p) for p in real_rds_dir.glob("*.rds"))
    manifest = {
        "root": str(root),
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "formats": formats,
        "ambiguous": ambiguous,
        "errors": errors,
        "url_tests": url_tests,
        "transpose_probes": transpose_probes,
        "real_rds": {"dir": str(real_rds_dir), "files": real_rds},
        "notes": "CONTRACT §8.5 fixture；truth.json 在每个 fixture 目录内",
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    say("fixture root = %s" % root)
    say("formats=%d ambiguous=%d errors=%d url=%d transpose=%d real_rds=%d" % (
        len(formats), len(ambiguous), len(errors), len(url_tests), len(transpose_probes), len(real_rds)))
    return manifest


# ---------------------------------------------------------------------------
# 本地 http.server 工具（py_import_test 也用它模拟 URL 下载，不联网）
# ---------------------------------------------------------------------------


class _QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, fmt, *args):  # noqa: A003 - 静默测试服务
        pass


@contextlib.contextmanager
def serve(root, port=0):
    """在 127.0.0.1 起一个只读静态服务，yield base_url；退出时关闭。"""
    root = Path(root).expanduser().resolve()
    handler = functools.partial(_QuietHandler, directory=str(root))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", int(port)), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        yield "http://127.0.0.1:%d" % httpd.server_address[1]
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join(timeout=5)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="生成 dsh-singlecell-workbench 全格式导入 fixture（CONTRACT §8.5）")
    parser.add_argument("--out", default=None, help="fixture 根目录（默认临时目录）")
    parser.add_argument("--no-r", action="store_true", help="跳过需要 Rscript 的 rds/RData fixture")
    parser.add_argument("--serve", action="store_true", help="生成后起本地 http.server 并打印 URL（Ctrl-C 退出）")
    args = parser.parse_args(argv)

    root = Path(args.out).expanduser().resolve() if args.out else Path(tempfile.mkdtemp(prefix="dsh-sc-formats-"))
    manifest = generate(root, with_r=not args.no_r)
    print("\nmanifest: %s" % (root / "manifest.json"))
    for entry in manifest["formats"]:
        print("  [format] %-24s %-12s n_obs=%-3d n_vars=%-3d %s" % (
            entry["id"], entry["format"], entry["truth"]["n_obs"], entry["truth"]["n_vars"],
            entry["path"]))
    for entry in manifest["ambiguous"]:
        print("  [ambiguous] %-24s candidates_field=%s" % (entry["id"], entry["candidates_field"]))
    for entry in manifest["errors"]:
        print("  [error] %-24s %s" % (entry["id"], entry["path"]))
    if args.serve:
        with serve(root) as base_url:
            print("\n本地 http.server: %s" % base_url)
            print("例如: %s/text_csv/counts.csv" % base_url)
            print("Ctrl-C 退出…")
            try:
                while True:
                    time.sleep(3600)
            except KeyboardInterrupt:
                pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
