#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell — 全格式导入层（CONTRACT §8）

把 GEO 上实际能下到的单细胞数据格式统一读成
  X: scipy.sparse.csr_matrix (cells × genes, float32)
  obs: pandas.DataFrame（行 = 细胞）
  var_names: list[str]（列 = 基因）

支持的格式（`SUPPORTED_FORMATS`）：
  h5ad / 10x_h5 / 10x_mtx / loom / matrix_text / geo_series / rds / rdata / archive / gz
  + URL 下载（http/https/ftp）与 GEO GSE 号（先列 supplement 文件再选）

只依赖标准库 + numpy/pandas/scipy/h5py/anndata（anndata 可选）+ 本机 Rscript（R 对象用）。

对外 API（键名固定，Host/引擎按它接线）：
  inspect_source(source, *, workdir, options=None) -> dict
  load_source(source, *, workdir, options=None, progress=None) -> dict
"""

from __future__ import annotations

import bz2
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

try:  # anndata 可选（h5ad 需要）
    import anndata as ad
except Exception:  # pragma: no cover
    ad = None

try:  # h5py 基本必备（10x h5 / loom）；缺失时只影响这两种格式
    import h5py
except Exception:  # pragma: no cover
    h5py = None


# ---------------------------------------------------------------------------
# 错误与常量
# ---------------------------------------------------------------------------


class EngineError(Exception):
    """可预期错误；code 与 CONTRACT §2.1/§8 一致。"""

    CODES = {"bad_request", "not_found", "prerequisite", "unsupported", "internal", "ambiguous"}

    def __init__(self, code, message, candidates=None):
        self.code = code if code in self.CODES else "internal"
        self.message = str(message)
        self.candidates = candidates or None
        super().__init__(self.message)


SUPPORTED_FORMATS = [
    "h5ad", "10x_h5", "10x_mtx", "loom", "matrix_text",
    "geo_series", "rds", "rdata", "archive", "auto",
]

TEXT_SUFFIXES = (".csv", ".tsv", ".txt", ".tab", ".matrix", ".counts")
DATA_SUFFIXES = (
    ".h5ad", ".h5", ".loom", ".mtx", ".csv", ".tsv", ".txt", ".tab",
    ".rds", ".rdata", ".rda", ".zip", ".tar", ".gz", ".bz2",
)
MAX_DOWNLOAD_BYTES = 8 * 1024 * 1024 * 1024
MAX_UNPACK_BYTES = 20 * 1024 * 1024 * 1024
R_TIMEOUT = 900
_HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"


def _err(code, message, candidates=None):
    return EngineError(code, message, candidates)


def _rscript():
    for candidate in (
        "/usr/local/bin/Rscript",
        "/Library/Frameworks/R.framework/Resources/bin/Rscript",
    ):
        if Path(candidate).exists():
            return candidate
    found = shutil.which("Rscript")
    if found:
        return found
    if os.name == "nt":  # Windows：常见安装位置
        import glob as _glob
        for pattern in (
            r"C:\\Program Files\\R\\R-*\\bin\\Rscript.exe",
            r"C:\\Program Files (x86)\\R\\R-*\\bin\\Rscript.exe",
            r"C:\\R\\R-*\\bin\\Rscript.exe",
        ):
            matches = sorted(_glob.glob(pattern), reverse=True)
            if matches:
                return matches[0]
    return None


def _as_csr(matrix):
    if matrix is None:
        return None
    if sp.issparse(matrix):
        result = sp.csr_matrix(matrix, dtype=np.float32)
    else:
        if isinstance(matrix, pd.DataFrame):
            matrix = matrix.to_numpy(dtype=np.float32, na_value=0.0)
        result = sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    if result.nnz:
        data = result.data
        bad = ~np.isfinite(data)
        if bad.any():
            result.data[bad] = 0.0
    return result


def _decode(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        return [_decode(item) for item in value]
    if isinstance(value, (list, tuple)):
        return [_decode(item) for item in value]
    return value


def _string_list(values, length, prefix="gene"):
    items = []
    if values is not None:
        for item in _decode(values):
            text = str(item).strip()
            items.append(text if text and text.lower() != "nan" else "")
    if len(items) != length:
        items = [item or "" for item in items[:length]]
        while len(items) < length:
            items.append("")
    for index, item in enumerate(items):
        if not item:
            items[index] = "%s_%d" % (prefix, index + 1)
    return items


def _payload(fmt, source, matrix, obs=None, var_names=None, counts=None, var=None,
             warnings=None, notes=None, meta=None, name=None):
    X = _as_csr(matrix)
    if X is None:
        raise _err("bad_request", "没有解析出表达矩阵（%s）" % source)
    n_obs, n_vars = X.shape
    if n_obs == 0 or n_vars == 0:
        raise _err("bad_request", "矩阵维度为 0（%d 细胞 × %d 基因）：%s" % (n_obs, n_vars, source))
    names = _string_list(var_names, n_vars, prefix="gene")
    if isinstance(obs, pd.DataFrame):
        obs_frame = obs.copy()
        if obs_frame.shape[0] != n_obs:
            if obs_frame.shape[0] > n_obs:
                obs_frame = obs_frame.iloc[:n_obs].copy()
            else:
                obs_frame = pd.DataFrame(index=pd.Index(["cell_%d" % (i + 1) for i in range(n_obs)], name="cell"))
    else:
        obs_frame = pd.DataFrame(index=pd.Index(["cell_%d" % (i + 1) for i in range(n_obs)], name="cell"))
    if obs_frame.index.has_duplicates:
        obs_frame.index = pd.Index(
            ["%s_%d" % (str(item), dup) if dup else str(item) for dup, item in
             zip(obs_frame.index.duplicated(keep=False).astype(int), obs_frame.index)],
            name=obs_frame.index.name or "cell",
        )
    counts_csr = _as_csr(counts) if counts is not None else None
    if counts_csr is not None and counts_csr.shape != X.shape:
        counts_csr = None
    warnings = list(warnings or [])
    notes = list(notes or [])
    meta = dict(meta or {})
    meta.setdefault("species", None)
    meta.setdefault("n_obs", int(n_obs))
    meta.setdefault("n_vars", int(n_vars))
    meta.setdefault("has_counts", counts_csr is not None)
    default_name = Path(str(source).split("?")[0].rstrip("/")).name or "dataset"
    return {
        "format": fmt,
        "source": str(source),
        "name": name or meta.pop("name", None) or default_name,
        "X": X,
        "counts": counts_csr,
        "obs": obs_frame,
        "var_names": names,
        "var": var if isinstance(var, pd.DataFrame) and len(var) == n_vars else None,
        "warnings": warnings,
        "notes": notes,
        "meta": meta,
    }


def _workdirs(workdir):
    root = Path(workdir or tempfile.gettempdir())
    dirs = {
        "root": root,
        "downloads": root / "downloads",
        "unpack": root / "unpack",
        "tmp": root / "tmp",
    }
    for path in dirs.values():
        try:
            path.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass
    return dirs


# ---------------------------------------------------------------------------
# URL / 压缩包
# ---------------------------------------------------------------------------


def _looks_like_url(source):
    return bool(re.match(r"^(https?|ftp)://", str(source).strip(), re.I))


def _looks_like_gse(source):
    return bool(re.fullmatch(r"GSE\d+", str(source).strip(), re.I))


def _gse_suppl_dir(gse):
    gse = gse.upper()
    digits = gse[3:]
    prefix = "GSE" + digits[:-3] + "nnn" if len(digits) > 3 else "GSEnnn"
    return "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/" % (prefix, gse)


def _http_get(url, timeout=60):
    request = urllib.request.Request(url, headers={"User-Agent": "dsh-singlecell/0.3"})
    return urllib.request.urlopen(request, timeout=timeout)


def _gse_candidates(gse):
    url = _gse_suppl_dir(gse)
    try:
        with _http_get(url, timeout=45) as response:
            html = response.read().decode("utf-8", errors="replace")
    except Exception as error:
        raise _err("bad_request", "无法读取 GEO 补充文件列表（%s）：%s" % (url, error))
    names = []
    for match in re.finditer(r'href="([^"]+)"', html):
        name = urllib.parse.unquote(match.group(1))
        if name.startswith(("..", "/", "http")) or name.endswith("/"):
            continue
        names.append(name)
    if not names:
        raise _err("bad_request", "GEO 目录里没有找到补充文件：%s" % url)
    candidates = [{"member": name, "name": name, "url": url + urllib.parse.quote(name)} for name in names]
    return {"format": "geo_suppl", "source": gse, "candidates": candidates,
            "notes": ["GEO 补充文件目录：%s" % url, "用 member 参数指定要下载的文件"]}


def _filename_from_response(url, response):
    disposition = response.headers.get("Content-Disposition") or ""
    match = re.search(r"filename\*\s*=\s*[^']*''([^;]+)", disposition, re.I)
    if not match:
        match = re.search(r'filename\s*=\s*"?([^";]+)"?', disposition, re.I)
    if match:
        name = urllib.parse.unquote(match.group(1).strip().strip('"'))
        if name:
            return Path(name).name
    parsed = urllib.parse.urlparse(url)
    name = Path(urllib.parse.unquote(parsed.path)).name
    return name or "download"


def _download(url, workdir, options=None, progress=None):
    options = options or {}
    dirs = _workdirs(workdir)
    max_bytes = int(options.get("max_bytes") or MAX_DOWNLOAD_BYTES)
    digest = hashlib.sha1(str(url).encode("utf-8")).hexdigest()[:10]
    if progress:
        progress("正在下载：%s" % url, 0)
    try:
        with _http_get(url, timeout=120) as response:
            filename = _filename_from_response(url, response)
            target = dirs["downloads"] / ("%s-%s" % (digest, filename))
            total = 0
            with open(target, "wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > max_bytes:
                        handle.close()
                        try:
                            target.unlink()
                        except Exception:
                            pass
                        raise _err("bad_request", "下载超过上限（%d GB）：%s" % (max_bytes // (1024 ** 3), url))
                    handle.write(chunk)
                    if progress and total % (20 * 1024 * 1024) < 1024 * 1024:
                        progress("已下载 %.0f MB" % (total / 1024 / 1024), min(95, int(total / max(1, max_bytes) * 100)))
    except EngineError:
        raise
    except Exception as error:
        raise _err("bad_request", "下载失败（%s）：%s" % (url, error))
    if progress:
        progress("下载完成：%s" % target.name, 100)
    return target


def _safe_member_path(root, member):
    root = Path(root).resolve()
    target = (root / member).resolve()
    if root != target and root not in target.parents:
        raise _err("bad_request", "压缩包内路径非法（疑似路径穿越）：%s" % member)
    return target


def _extract_archive(path, workdir, progress=None):
    path = Path(path)
    dirs = _workdirs(workdir)
    digest = hashlib.sha1(str(path).encode("utf-8") + str(path.stat().st_mtime_ns).encode()).hexdigest()[:10]
    target = dirs["unpack"] / digest
    if target.exists():
        return target
    target.mkdir(parents=True, exist_ok=True)
    total = 0
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                name = info.filename
                if name.endswith("/"):
                    continue
                if name.startswith(("/", "..")) or ".." in Path(name).parts:
                    continue
                if (info.external_attr >> 16) & 0o170000 == 0o120000:
                    continue
                total += int(info.file_size or 0)
                if total > MAX_UNPACK_BYTES:
                    raise _err("bad_request", "压缩包解压后超过上限（20 GB）")
                destination = _safe_member_path(target, name)
                destination.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(info) as source, open(destination, "wb") as out:
                    shutil.copyfileobj(source, out, 1024 * 1024)
                if progress:
                    progress("解压：%s" % name)
    else:
        try:
            archive = tarfile.open(path)
        except Exception as error:
            raise _err("bad_request", "无法识别的压缩包：%s（%s）" % (path.name, error))
        with archive:
            for member in archive.getmembers():
                if not member.isfile():
                    continue
                name = member.name.lstrip("./")
                if name.startswith(("/", "..")) or ".." in Path(name).parts:
                    continue
                total += int(member.size or 0)
                if total > MAX_UNPACK_BYTES:
                    raise _err("bad_request", "压缩包解压后超过上限（20 GB）")
                destination = _safe_member_path(target, name)
                destination.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                if source is None:
                    continue
                with source, open(destination, "wb") as out:
                    shutil.copyfileobj(source, out, 1024 * 1024)
                if progress:
                    progress("解压：%s" % name)
    return target


def _find_10x_dir(root):
    """在目录树里找 10x 三件套所在目录（matrix.mtx(.gz) + barcodes/features）。"""
    for current, _dirs, files in os.walk(root):
        names = {name.lower() for name in files}
        if any("matrix.mtx" in name for name in names):
            return Path(current)
    return None


def _archive_candidates(root):
    root = Path(root)
    found = []
    for current, _dirs, files in os.walk(root):
        for name in files:
            suffix = Path(name).suffix.lower()
            lowered = name.lower()
            if lowered.endswith((".tar.gz", ".tgz", ".tar.bz2")):
                suffix = ".tar.gz"
            if suffix not in (".h5ad", ".h5", ".loom", ".mtx", ".csv", ".tsv", ".txt", ".rds", ".rdata", ".rda", ".zip", ".tar", ".tar.gz", ".gz"):
                continue
            if lowered in ("barcodes.tsv", "features.tsv", "genes.tsv") or lowered.startswith(("barcodes.tsv", "features.tsv", "genes.tsv")):
                continue
            if lowered.startswith("matrix.mtx"):
                found.append({"member": str(Path(current).joinpath(name).relative_to(root)),
                              "name": name, "size": Path(current, name).stat().st_size,
                              "format": "10x_mtx"})
                continue
            found.append({"member": str(Path(current).joinpath(name).relative_to(root)),
                          "name": name, "size": Path(current, name).stat().st_size,
                          "format": None})
    for entry in found:
        if not entry.get("format"):
            entry["format"] = _format_by_suffix(entry["name"])
    # 有 10x mtx 时优先它，并把同目录视为一个候选
    mtx = [entry for entry in found if entry.get("format") == "10x_mtx"]
    if mtx:
        first = mtx[0]
        directory = str(Path(first["member"]).parent)
        return [{"member": directory or ".", "name": Path(directory or ".").name or first["name"],
                 "format": "10x_mtx", "size": sum(entry["size"] for entry in mtx)}]
    found.sort(key=lambda entry: (0 if entry.get("format") in ("h5ad", "h5", "loom") else 1, -entry.get("size", 0)))
    for entry in found:
        entry.pop("_", None)
    return found


def _format_by_suffix(name):
    lowered = str(name).lower()
    if lowered.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar", ".zip")):
        return "archive"
    if lowered.endswith((".h5ad", ".h5ad.gz")):
        return "h5ad"
    if lowered.endswith((".loom", ".loom.gz")):
        return "loom"
    if lowered.endswith((".mtx", ".mtx.gz")):
        return "10x_mtx"
    if lowered.endswith((".rds", ".rds.gz", ".rds.bz2")):
        return "rds"
    if lowered.endswith((".rdata", ".rda", ".rdata.gz", ".rda.gz")):
        return "rdata"
    if lowered.endswith(".h5") or lowered.endswith(".hdf5"):
        return "10x_h5"
    if any(lowered.endswith(suffix) for suffix in TEXT_SUFFIXES) or any(
        (suffix + ".gz" if suffix else suffix) in lowered for suffix in TEXT_SUFFIXES
    ):
        return "matrix_text"
    return None


# ---------------------------------------------------------------------------
# HDF5 系列（h5ad / 10x h5 / loom）
# ---------------------------------------------------------------------------


def _sniff_h5(path):
    if h5py is None:
        return None
    try:
        with h5py.File(path, "r") as handle:
            keys = set(handle.keys())
            if "matrix" in keys:
                node = handle["matrix"]
                if isinstance(node, h5py.Group) and "data" in node:
                    return "10x_h5"
                if isinstance(node, h5py.Dataset):
                    return "loom"
            if "X" in keys and "obs" in keys:
                return "h5ad"
            if "adata" in keys and isinstance(handle["adata"], h5py.Group):
                # 普通 .h5 里嵌了一个 h5ad 结构（GEO 上偶见）；当 h5ad 处理
                group = handle["adata"]
                if "X" in group.keys() and "obs" in group.keys():
                    return "h5ad"
            if "row_attrs" in keys and "matrix" in keys:
                return "loom"
    except Exception:
        return None
    return None


def _h5_shape(path):
    """inspect 用：尽量便宜地拿到维度（不加载矩阵）。"""
    if h5py is None:
        return None
    try:
        with h5py.File(path, "r") as handle:
            keys = set(handle.keys())
            if "matrix" in keys:
                node = handle["matrix"]
                if isinstance(node, h5py.Group) and "shape" in node:
                    shape = [int(value) for value in node["shape"][:]]
                    n_obs, n_vars = shape[1], shape[0]
                    # 10x v3：有 feature_type 时只算 Gene Expression
                    if "features" in node and isinstance(node["features"], h5py.Group):
                        features = node["features"]
                        if "feature_type" in features:
                            types = [str(item) for item in _decode(features["feature_type"][:])]
                            keep = [item for item in types if item.strip().lower() in ("gene expression", "")]
                            if keep and len(keep) < len(types):
                                n_vars = len(keep)
                    return {"n_obs": int(n_obs), "n_vars": int(n_vars)}
                if isinstance(node, h5py.Dataset):
                    # loom：/matrix 是 genes × cells
                    rows, cols = node.shape
                    return {"n_obs": int(cols), "n_vars": int(rows)}
            if "adata" in keys and isinstance(handle["adata"], h5py.Group):
                group = handle["adata"]
                if "X" in group.keys():
                    node = group["X"]
                    if isinstance(node, h5py.Dataset):
                        rows, cols = node.shape
                        return {"n_obs": int(rows), "n_vars": int(cols)}
                    if isinstance(node, h5py.Group):
                        shape = node.attrs.get("shape")
                        if shape is not None and len(shape) == 2:
                            return {"n_obs": int(shape[0]), "n_vars": int(shape[1])}
            if "X" in keys:
                node = handle["X"]
                if isinstance(node, h5py.Dataset):
                    rows, cols = node.shape
                    return {"n_obs": int(rows), "n_vars": int(cols)}
                if isinstance(node, h5py.Group):
                    shape = node.attrs.get("shape")
                    if shape is not None and len(shape) == 2:
                        return {"n_obs": int(shape[0]), "n_vars": int(shape[1])}
                    if "obs" in keys and isinstance(handle["obs"], h5py.Group) and "_index" in handle["obs"]:
                        n_obs = int(handle["obs"]["_index"].shape[0])
                        n_vars = None
                        if "var" in keys and isinstance(handle["var"], h5py.Group) and "_index" in handle["var"]:
                            n_vars = int(handle["var"]["_index"].shape[0])
                        return {"n_obs": n_obs, "n_vars": n_vars}
    except Exception:
        return None
    return None


def _load_10x_h5(path, options=None):
    if h5py is None:
        raise _err("unsupported", "缺少 h5py，无法读取 10x .h5")
    warnings = []
    with h5py.File(path, "r") as handle:
        if "matrix" not in handle or not isinstance(handle["matrix"], h5py.Group):
            raise _err("bad_request", "不是 10x CellRanger 的 .h5 文件（缺少 matrix 组）：%s" % Path(path).name)
        group = handle["matrix"]
        shape = tuple(int(value) for value in group["shape"][:])
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        # 10x 存的是 genes × cells 的 CSC
        matrix = sp.csc_matrix((data, indices, indptr), shape=shape)
        X = matrix.T.tocsr()
        barcodes = _string_list(
            group["barcodes"][:] if "barcodes" in group else None,
            X.shape[0], prefix="cell",
        )
        names = None
        if "features" in group:
            features = group["features"]
            if isinstance(features, h5py.Group):
                if "name" in features:
                    names = _decode(features["name"][:])
                elif "id" in features:
                    names = _decode(features["id"][:])
                if "feature_type" in features:
                    types = [str(item) for item in _decode(features["feature_type"][:])]
                    keep = [index for index, item in enumerate(types) if item.strip().lower() in ("gene expression", "")]
                    if keep and len(keep) < len(types):
                        X = X[:, keep]
                        if names is not None:
                            names = [names[index] for index in keep]
                        warnings.append("该 10x 文件含多种 feature 类型，已只保留 Gene Expression（%d/%d）" % (len(keep), len(types)))
            else:
                names = _decode(features[:])
        elif "genes" in group:
            genes = group["genes"]
            if isinstance(genes, h5py.Group):
                names = _decode(genes["gene_names"][:] if "gene_names" in genes else genes["id"][:])
            else:
                names = _decode(genes[:])
        elif "gene_names" in group:
            names = _decode(group["gene_names"][:])
    var_names = _string_list(names, X.shape[1], prefix="gene")
    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell"))
    notes = ["10x CellRanger HDF5：矩阵方向 genes×cells 已转置为 cells×genes"]
    meta = {"species": None, "assays": ["RNA"]}
    return _payload("10x_h5", path, X, obs=obs, var_names=var_names, counts=X,
                    warnings=warnings, notes=notes, meta=meta)


def _load_loom(path, options=None):
    if h5py is None:
        raise _err("unsupported", "缺少 h5py，无法读取 .loom")
    with h5py.File(path, "r") as handle:
        if "matrix" not in handle:
            raise _err("bad_request", "不是合法的 loom 文件（缺少 /matrix）：%s" % Path(path).name)
        matrix = handle["matrix"]
        dense = np.asarray(matrix[:])
        if dense.ndim != 2:
            raise _err("bad_request", "loom 的 /matrix 不是二维：%s" % Path(path).name)
        # loom 约定 matrix 是 genes × cells
        X = sp.csr_matrix(dense.T)
        genes = None
        if "row_attrs" in handle:
            attrs = handle["row_attrs"]
            for key in ("Gene", "gene", "genes", "var_names", "Accession"):
                if key in attrs:
                    genes = _decode(attrs[key][:])
                    break
            if "Chromosome" in attrs:
                pass
        if genes is None:
            genes = ["gene_%d" % (index + 1) for index in range(X.shape[1])]
        obs_data = {}
        if "col_attrs" in handle:
            for key, node in handle["col_attrs"].items():
                if not isinstance(node, h5py.Dataset):
                    continue
                values = _decode(node[:])
                if len(values) == X.shape[0]:
                    obs_data[str(key)] = values
        obs = pd.DataFrame(obs_data, index=pd.Index([str(item) for item in obs_data.get("CellID", ["cell_%d" % (i + 1) for i in range(X.shape[0])])], name="cell")) if obs_data else pd.DataFrame(index=pd.Index(["cell_%d" % (i + 1) for i in range(X.shape[0])], name="cell"))
        if "CellID" in obs.columns:
            obs.index = pd.Index([str(item) for item in obs["CellID"].tolist()], name="cell")
    var_names = _string_list(genes, X.shape[1], prefix="gene")
    notes = ["loom：矩阵方向 genes×cells 已转置为 cells×genes；/col_attrs 已转成 obs 列"]
    return _payload("loom", path, X, obs=obs, var_names=var_names, counts=X,
                    notes=notes, warnings=[], meta={"assays": ["RNA"]})


def _read_h5_adata_group(path):
    """普通 .h5 里内嵌的 `adata` 组：复制成临时 .h5ad 再让 anndata 读。"""
    if h5py is None or ad is None:
        return None
    try:
        with h5py.File(path, "r") as source:
            if "adata" not in source or not isinstance(source["adata"], h5py.Group):
                return None
            group = source["adata"]
            if "X" not in group.keys() or "obs" not in group.keys():
                return None
            tmp_dir = Path(tempfile.mkdtemp(prefix="sc-io-adata-"))
            temp_path = tmp_dir / "adata.h5ad"
            with h5py.File(temp_path, "w") as target:
                for key in group.keys():
                    source.copy(group[key], target, key)
                for key, value in group.attrs.items():
                    target.attrs[key] = value
        return ad.read_h5ad(str(temp_path))
    except Exception:
        return None


def _load_h5ad(path, options=None):
    if ad is None:
        raise _err("unsupported", "缺少 anndata，无法读取 .h5ad（可在 Host 的 Python 环境安装 anndata）")
    try:
        data = ad.read_h5ad(str(path))
    except Exception as error:
        data = _read_h5_adata_group(path)
        if data is None:
            raise _err("bad_request", "无法解析 .h5ad（%s）：%s" % (Path(path).name, error))
    if data.X is None:
        raise _err("bad_request", "h5ad 里没有 X 表达矩阵：%s" % Path(path).name)
    warnings = []
    counts = None
    layers = getattr(data, "layers", None) or {}
    if "counts" in layers:
        counts = layers["counts"]
    elif getattr(data, "raw", None) is not None and data.raw.X is not None:
        counts = data.raw.X
        warnings.append("未找到 layers['counts']，已用 adata.raw.X 作为原始计数")
    matrix = data.X
    try:
        dtype = matrix.dtype
        integer_like = np.issubdtype(dtype, np.integer)
    except Exception:
        integer_like = False
    if not integer_like and sp.issparse(matrix) and matrix.nnz:
        sample = matrix.data[: min(1000, matrix.nnz)]
        if not np.all(np.abs(sample - np.round(sample)) < 1e-6) and float(np.max(sample)) < 30:
            warnings.append("X 看起来已 log 标准化；内置流程的 normalize 会基于 counts 重算")
    obs = data.obs.copy() if data.obs is not None else None
    var = data.var.copy() if data.var is not None else None
    meta = {"assays": ["RNA"]}
    payload = _payload("h5ad", path, matrix, obs=obs, var_names=list(map(str, data.var_names)),
                       counts=counts, var=var, warnings=warnings, meta=meta)
    obsm = {}
    try:
        if "X_pca" in data.obsm:
            obsm["X_pca"] = np.asarray(data.obsm["X_pca"], dtype=np.float32)
        for key in ("X_umap", "X_tsne", "X_UMAP", "X_tSNE"):
            if key in data.obsm and np.asarray(data.obsm[key]).ndim == 2:
                obsm["X_emb"] = np.asarray(data.obsm[key], dtype=np.float32)[:, :2]
                break
    except Exception:
        obsm = {}
    var_features = []
    if "highly_variable" in getattr(data, "var", pd.DataFrame()).columns:
        try:
            mask = np.asarray(data.var["highly_variable"].to_numpy(), dtype=bool)
            var_features = [str(gene) for gene, keep in zip(data.var_names, mask) if keep]
        except Exception:
            var_features = []
    import_state = {
        "class": "AnnData",
        "reductions": [str(key) for key in list(getattr(data, "obsm", {}).keys())],
        "rank_genes": "rank_genes_groups" in getattr(data, "uns", {}),
        "variable_features": var_features,
        "embed_kind": "umap" if any("umap" in str(key).lower() for key in getattr(data, "obsm", {}).keys()) else "tsne",
    }
    payload["obsm"] = obsm
    payload["var_features"] = var_features
    payload["import_state"] = import_state
    payload["embed_kind"] = import_state["embed_kind"]
    return payload


# ---------------------------------------------------------------------------
# 10x mtx / 单 mtx
# ---------------------------------------------------------------------------


def _read_mtx(path):
    path = Path(path)
    try:
        if path.name.lower().endswith(".gz"):
            with gzip.open(path, "rb") as handle:
                return mmread(handle)
        return mmread(str(path))
    except Exception as error:
        raise _err("bad_request", "无法解析 MatrixMarket 文件（%s）：%s" % (path.name, error))


def _sidecar(directory, prefix, patterns):
    directory = Path(directory)
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def _read_lines_file(path):
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return [line.rstrip("\n").rstrip("\r") for line in handle]
    return [line.rstrip("\n").rstrip("\r") for line in path.read_text(encoding="utf-8", errors="replace").splitlines()]


def _load_10x_mtx(path, options=None):
    options = options or {}
    warnings = []
    path = Path(path)
    prefix = options.get("trio_prefix")
    barcodes_file = None
    features_file = None
    if path.is_dir() and prefix:
        # GEO 三联体：同一目录里按 GSM 前缀配对 matrix/barcodes/features
        directory = path
        mtx = None
        for item in sorted(directory.iterdir()):
            if not item.is_file():
                continue
            name = item.name.lower()
            if str(prefix).lower() not in name:
                continue
            if "matrix.mtx" in name:
                mtx = item
            elif "barcodes" in name:
                barcodes_file = item
            elif "features" in name or "genes" in name:
                features_file = item
        if mtx is None or barcodes_file is None or features_file is None:
            raise _err("not_found", "目录里没有找到前缀为 %s 的 10X 三联体：%s" % (prefix, directory))
    elif path.is_dir():
        directory = path
        mtx = _sidecar(directory, "matrix", ("matrix.mtx.gz", "matrix.mtx", "matrix.mtx.bz2"))
        if mtx is None:
            raise _err("not_found", "目录里没有 matrix.mtx(.gz)：%s" % directory)
    else:
        directory = path.parent
        mtx = path
    matrix = _read_mtx(mtx)
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    if barcodes_file is None:
        barcodes_file = _sidecar(directory, "barcodes", ("barcodes.tsv.gz", "barcodes.tsv"))
    if features_file is None:
        features_file = _sidecar(directory, "features", ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
    barcodes = _read_lines_file(barcodes_file) if barcodes_file else []
    features = _read_lines_file(features_file) if features_file else []
    feature_names = []
    keep_rows = None
    feature_types = []
    if features:
        for line in features:
            parts = line.split("\t")
            feature_names.append(parts[1] if len(parts) > 1 and parts[1].strip() else parts[0].strip())
            if len(parts) > 2:
                feature_types.append(parts[2].strip().lower())
    rows, cols = matrix.shape
    if barcodes and len(barcodes) == rows:
        X = sp.csr_matrix(matrix)
        cells = barcodes
    elif barcodes and len(barcodes) == cols:
        X = sp.csr_matrix(matrix.T)
        cells = barcodes
    elif feature_names and len(feature_names) == rows:
        X = sp.csr_matrix(matrix.T)
        cells = barcodes or ["cell_%d" % (index + 1) for index in range(cols)]
    else:
        X = sp.csr_matrix(matrix.T)
        cells = barcodes or ["cell_%d" % (index + 1) for index in range(X.shape[0])]
    if feature_types and len(feature_types) == len(feature_names):
        keep = [index for index, item in enumerate(feature_types) if item in ("gene expression", "")]
        if keep and len(keep) < len(feature_types):
            X = X[:, keep]
            feature_names = [feature_names[index] for index in keep]
            warnings.append("已只保留 Gene Expression 行（%d/%d）" % (len(keep), len(feature_types)))
    if not barcodes:
        warnings.append("没找到 barcodes.tsv(.gz)，细胞名自动生成")
    if not feature_names:
        warnings.append("没找到 features.tsv/genes.tsv，基因名自动生成")
    var_names = _string_list(feature_names, X.shape[1], prefix="gene")
    obs = pd.DataFrame(index=pd.Index(_string_list(cells, X.shape[0], prefix="cell"), name="cell"))
    notes = ["10x Matrix Market 三件套：matrix.mtx + barcodes/features"]
    return _payload("10x_mtx", path, X, obs=obs, var_names=var_names, counts=X,
                    warnings=warnings, notes=notes, meta={"assays": ["RNA"]})


# ---------------------------------------------------------------------------
# 文本矩阵 / GEO series matrix
# ---------------------------------------------------------------------------


def _open_text(path):
    path = Path(path)
    lowered = path.name.lower()
    if lowered.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if lowered.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _read_head(path, size=4 * 1024 * 1024):
    with _open_text(path) as handle:
        return handle.read(size)


def _looks_like_geo(text):
    return "!series_matrix_table_begin" in text


def _split_delimiter(sample_lines):
    candidates = ["\t", ",", ";", " "]
    best = "\t"
    best_score = -1
    for delimiter in candidates:
        counts = [len(line.split(delimiter)) for line in sample_lines[:8] if line.strip()]
        if not counts:
            continue
        score = int(np.median(counts))
        if score > best_score:
            best_score = score
            best = delimiter
    return best


def _is_numeric_text(text):
    try:
        float(str(text).strip().replace(",", ""))
        return True
    except Exception:
        return False


def _barcode_score(labels):
    if not labels:
        return 0.0
    hits = 0
    for label in labels[:200]:
        text = str(label)
        if re.search(r"[ACGT]{10,}", text) or re.search(r"-\d+$", text):
            hits += 1
    return hits / min(len(labels), 200)


def _gene_score(labels):
    if not labels:
        return 0.0
    hits = 0
    checked = 0
    for label in labels[:200]:
        text = str(label).strip()
        if not text:
            continue
        checked += 1
        # 基因符号通常以大写字母开头（MT-ND1 / AL627309.1 / CD3D）；细胞名多为小写/数字混合
        if re.fullmatch(r"[A-Z][A-Z0-9_.\-]*", text) and len(text) <= 30:
            hits += 1
    return hits / checked if checked else 0.0


def _parse_geo_series(text, source):
    table_lines = []
    sample_lines = []
    in_table = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if stripped.startswith("!series_matrix_table_end"):
            break
        if in_table:
            table_lines.append(line)
        elif stripped.startswith("!Sample_"):
            sample_lines.append(line)
    if not table_lines:
        raise _err("bad_request", "GEO series matrix 里没有找到表格段：%s" % source)
    rows = list(csv.reader(io.StringIO("\n".join(table_lines)), delimiter="\t", quotechar='"'))
    rows = [row for row in rows if row and any(cell.strip() for cell in row)]
    if len(rows) < 2:
        raise _err("bad_request", "GEO series matrix 表格为空：%s" % source)
    header = [cell.strip() for cell in rows[0]]
    genes = []
    values = []
    for row in rows[1:]:
        if not row:
            continue
        genes.append(row[0].strip())
        values.append([float(cell) if cell.strip() not in ("", "NA", "null", "NaN") else 0.0
                       for cell in row[1:]] * 1)
    matrix = np.asarray(values, dtype=np.float64)
    cells = header[1:1 + matrix.shape[1]]
    obs_data = {}
    title = []
    counter = 0
    for line in sample_lines:
        parts = next(csv.reader([line], delimiter="\t", quotechar='"'))
        if len(parts) < 2:
            continue
        key = parts[0].strip()
        values_row = [item.strip() for item in parts[1:]]
        if key == "!Sample_title":
            title = values_row
            continue
        if key == "!Sample_characteristics_ch1":
            counter += 1
            label = "characteristics_ch1"
            if values_row and ":" in values_row[0]:
                label = values_row[0].split(":", 1)[0].strip() or label
            obs_data[label if label not in obs_data else "%s_%d" % (label, counter)] = values_row
        elif key.startswith("!Sample_"):
            label = key[len("!Sample_"):]
            if label in ("geo_accession", "status", "submission_date", "last_update_date", "type",
                         "channel_count", "contact_name", "contact_email", "contact_institute",
                         "contact_address", "contact_city", "contact_state", "contact_zip/postal_code",
                         "contact_country", "supplementary_file", "data_row_count", "platform_id",
                         "platform_taxid", "sample_taxid", "sample_attributes", "relation",
                         "series_id", "title", "source_name_ch1", "organism_ch1", "molecule_ch1",
                         "extract_protocol_ch1", "library_strategy", "library_source", "library_selection",
                         "instrument_model", "description", "description.1", "data_processing"):
                if label in obs_data:
                    continue
                obs_data[label] = values_row
    obs_index = pd.Index(_string_list(title if title else cells, len(cells), prefix="sample"), name="sample")
    obs = pd.DataFrame(index=obs_index)
    for key, values_row in obs_data.items():
        if len(values_row) == len(cells):
            obs[key] = values_row
    warnings = ["这是 GEO series matrix（样本级汇总矩阵），不是单细胞计数矩阵；请确认表中的列确实是样本"]
    notes = ["从 !series_matrix_table_begin 段解析；!Sample_* 行已转成 obs 列"]
    X = sp.csr_matrix(matrix.T)
    return _payload("geo_series", source, X, obs=obs, var_names=_string_list(genes, X.shape[1], prefix="gene"),
                    warnings=warnings, notes=notes, meta={"assays": ["microarray/RNA"]})


def _parse_dense_text(text, path, options):
    """文本矩阵解析；任何底层异常都包装成可读的 EngineError。"""
    try:
        return _parse_dense_text_inner(text, path, options)
    except EngineError:
        raise
    except Exception as error:  # noqa: BLE001
        raise _err("bad_request", "无法解析文本矩阵（%s）：%s" % (Path(path).name, error))


def _parse_dense_text_inner(text, path, options):
    options = options or {}
    lines = [line for line in text.splitlines()[:200] if line.strip()]
    if not lines:
        raise _err("bad_request", "文本文件是空的：%s" % Path(path).name)
    delimiter = options.get("delimiter") or _split_delimiter(lines)
    header_cells = next(csv.reader([lines[0]], delimiter=delimiter, quotechar='"'))
    # 首行「除首列外全是数值」⇒ 这是数据行，无表头（常见：barcode + 计数）
    data_like_first_row = bool(header_cells[1:]) and all(
        _is_numeric_text(cell) for cell in header_cells[1:] if str(cell).strip()
    )
    header = (not data_like_first_row) and not all(
        _is_numeric_text(cell) for cell in header_cells if str(cell).strip()
    )
    first_col_numeric = all(_is_numeric_text(row.split(delimiter)[0]) for row in lines[1:21] if row.strip())
    read_kwargs = dict(
        sep=(r"\s+" if delimiter == " " else delimiter),
        header=0 if header else None,
        index_col=0 if not first_col_numeric else None,
        quoting=csv.QUOTE_MINIMAL,
    )
    if delimiter == " ":
        read_kwargs["engine"] = "python"
    else:
        read_kwargs["engine"] = "c"
        read_kwargs["low_memory"] = False
    frame = pd.read_csv(path, **read_kwargs)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    # 丢掉全空的列/行（例如文本末尾多余的分隔符）
    if frame.shape[1] > 1:
        frame = frame.loc[:, ~frame.isna().all(axis=0)]
    if frame.shape[0] > 1:
        frame = frame.loc[~frame.isna().all(axis=1), :]
    if frame.isna().all(axis=None):
        raise _err("bad_request", "文本矩阵里没有可解析的数值：%s" % Path(path).name)
    frame = frame.fillna(0.0)
    row_labels = [str(item) for item in frame.index]
    col_labels = [str(item) for item in frame.columns]
    transpose = options.get("transpose")
    notes = []
    if transpose is None:
        row_score = _barcode_score(row_labels)
        col_score = _barcode_score(col_labels)
        row_gene = _gene_score(row_labels)
        col_gene = _gene_score(col_labels)
        if row_score >= 0.5 and col_score < 0.5:
            transpose = False
        elif col_score >= 0.5 and row_score < 0.5:
            transpose = True
        elif row_gene >= 0.8 and col_gene < 0.8:
            transpose = True
        elif col_gene >= 0.8 and row_gene < 0.8:
            transpose = False
        else:
            transpose = bool(frame.shape[0] < frame.shape[1])
            notes.append("矩阵方向无法确定，按形状启发式处理（%d 行 × %d 列 → %s）；如不对请用 transpose 参数"
                         % (frame.shape[0], frame.shape[1], "转置" if transpose else "保持"))
    else:
        transpose = bool(transpose)
        notes.append("按用户指定的 transpose=%s 处理" % transpose)
    if transpose:
        frame = frame.T
        notes.append("已转置为 cells × genes")
    # 无表头/无行名时用 gene_N / cell_N 命名（避免列名变成 "0","1"... 这种位置数字）
    if not header:
        frame.columns = ["gene_%d" % (index + 1) for index in range(frame.shape[1])]
    if read_kwargs.get("index_col") is None:
        frame.index = ["cell_%d" % (index + 1) for index in range(frame.shape[0])]
    values = frame.to_numpy(dtype=np.float64)
    integer_like = np.all(np.abs(values - np.round(values)) < 1e-6) if values.size else False
    X = sp.csr_matrix(values.astype(np.float32))
    cells = _string_list([str(item) for item in frame.index], X.shape[0], prefix="cell")
    genes = _string_list([str(item) for item in frame.columns], X.shape[1], prefix="gene")
    obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
    warnings = []
    if not header:
        warnings.append("文本没有表头，列名自动生成")
    if X.shape[0] * X.shape[1] > 50_000_000:
        warnings.append("矩阵较大（%d × %d），已整体读入内存；如太慢可先用 scanpy 转 h5ad" % X.shape)
    return _payload("matrix_text", path, X, obs=obs, var_names=genes,
                    counts=X if integer_like else None,
                    warnings=warnings, notes=notes,
                    meta={"assays": ["RNA"], "has_counts": bool(integer_like)})


def _load_text(path, options=None):
    options = options or {}
    try:
        head = _read_head(path)
    except Exception as error:
        raise _err("bad_request", "无法读取文本文件（%s）：%s" % (Path(path).name, error))
    if _looks_like_geo(head):
        full = _read_head(path, size=10 ** 9) if Path(path).stat().st_size < 1024 ** 3 else None
        if full is None:
            raise _err("bad_request", "GEO series matrix 太大（>1 GB），请先裁剪：%s" % Path(path).name)
        return _parse_geo_series(full, path)
    return _parse_dense_text(head, path, options)


# ---------------------------------------------------------------------------
# R 对象（rds / RData）
# ---------------------------------------------------------------------------


def _load_r(path, workdir, options=None):
    options = options or {}
    rscript = _rscript()
    if not rscript:
        raise _err("unsupported", "没有找到 Rscript，无法读取 .rds/.RData。请安装 R（含 Seurat/SeuratObject/SingleCellExperiment）或先用 R 另存为 h5ad/10x mtx")
    bridge = Path(__file__).resolve().parent / "r_bridge.R"
    if not bridge.exists():
        raise _err("unsupported", "缺少 R 桥脚本：%s" % bridge)
    dirs = _workdirs(workdir)
    target = Path(tempfile.mkdtemp(prefix="r-", dir=str(dirs["tmp"])))
    command = [rscript, "--vanilla", str(bridge), "--in", str(Path(path).resolve()), "--out", str(target)]
    if options.get("object"):
        command += ["--object", str(options["object"])]
    if options.get("inspect"):
        command += ["--inspect"]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=R_TIMEOUT)
    except subprocess.TimeoutExpired:
        raise _err("bad_request", "R 读取超时（>%d 秒）：%s" % (R_TIMEOUT, Path(path).name))
    except Exception as error:
        raise _err("bad_request", "无法调用 Rscript：%s" % error)
    meta_path = target / "meta.json"
    meta = {}
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            meta = {}
    if completed.returncode != 0 or meta.get("status") == "error":
        detail = (meta.get("error") or completed.stderr or completed.stdout or "").strip()
        if meta_path.exists():
            try:
                detail = detail or meta_path.read_text(encoding="utf-8", errors="replace")[:800]
            except Exception:
                pass
        raise _err("bad_request", "R 读取失败（%s）：%s" % (Path(path).name, detail[-1200:] or "未知错误"))
    if meta.get("status") == "ambiguous":
        candidates = []
        for entry in meta.get("candidates") or []:
            if isinstance(entry, dict):
                candidates.append({
                    "object": entry.get("object") or entry.get("name"),
                    "class": entry.get("class"),
                    "label": entry.get("label") or entry.get("object"),
                    "n_obs": entry.get("n_obs") or (entry.get("dim") or [None, None])[0],
                    "n_vars": entry.get("n_vars") or (entry.get("dim") or [None, None])[1],
                })
        raise _err("ambiguous",
                   "R 文件里有多个可用的单细胞对象，请用 object 参数指定：%s" % "、".join(str(item.get("object") or item.get("name")) for item in candidates),
                   candidates=candidates)
    if options.get("inspect") or meta.get("inspect"):
        return {"format": "rdata" if Path(path).suffix.lower() in (".rdata", ".rda", ".rdata.gz", ".rda.gz") else "rds",
                "source": str(path), "class": meta.get("class"), "n_obs": meta.get("n_obs"),
                "n_vars": meta.get("n_vars"), "candidates": None,
                "notes": _decode(meta.get("notes") or []), "warnings": _decode(meta.get("warnings") or [])}
    matrix_path = target / "matrix.mtx"
    if not matrix_path.exists():
        raise _err("internal", "R 桥没有写出 matrix.mtx（%s）" % Path(path).name)
    matrix = _read_mtx(matrix_path)
    if sp.issparse(matrix):
        X = sp.csr_matrix(matrix)
    else:
        X = sp.csr_matrix(np.asarray(matrix))
    features_path = target / "features.tsv"
    barcodes_path = target / "barcodes.tsv"
    var_names = _read_lines_file(features_path) if features_path.exists() else None
    cells = _read_lines_file(barcodes_path) if barcodes_path.exists() else None
    obs = None
    obs_path = target / "obs.csv"
    if obs_path.exists():
        try:
            obs = pd.read_csv(obs_path, index_col=0)
            obs.index = pd.Index([str(item) for item in obs.index], name="cell")
        except Exception:
            obs = None
    if cells:
        obs = obs if obs is not None else pd.DataFrame(index=pd.Index(cells, name="cell"))
        if obs.shape[0] == len(cells):
            obs.index = pd.Index(cells, name="cell")
    warnings = [str(item) for item in _decode(meta.get("warnings") or [])]
    notes = [str(item) for item in _decode(meta.get("notes") or [])]
    cls = meta.get("class") or "R object"
    notes.append("R 对象：%s；assay=%s；矩阵方向 = %s" % (cls, meta.get("assay") or "-", meta.get("orientation") or "cells_genes"))
    has_counts = bool(meta.get("has_counts"))
    X_shaped = X
    if (meta.get("orientation") or "cells_genes") == "genes_cells":
        X_shaped = X.T.tocsr()
        notes.append("已按 meta.orientation=genes_cells 转置为 cells × genes")
    fmt = "rdata" if str(path).lower().endswith((".rdata", ".rda", ".rdata.gz", ".rda.gz")) else "rds"
    counts = X_shaped if (has_counts and np.all(np.abs(X_shaped.data - np.round(X_shaped.data)) < 1e-6)) else None
    payload = _payload(fmt, path, X_shaped, obs=obs, var_names=var_names, counts=counts,
                       warnings=warnings, notes=notes,
                       meta={"assays": [str(item) for item in _decode(meta.get("assays") or ["RNA"])],
                             "cls": cls, "has_counts": counts is not None})
    # 已有分析进度：pca/embed/variable features/import_state
    def _read_embedding(file_name):
        target_file = target / file_name
        if not target_file.exists():
            return None
        try:
            frame = pd.read_csv(target_file)
            expected = len(cells) if cells else (len(obs) if isinstance(obs, pd.DataFrame) else 0)
            if expected and frame.shape[0] != expected:
                return None
            return frame.iloc[:, 1:].to_numpy(dtype=np.float32)
        except Exception:
            return None

    obsm = {}
    pca_matrix = _read_embedding("pca.csv")
    if pca_matrix is not None:
        obsm["X_pca"] = pca_matrix
    embed_matrix = _read_embedding("embed.csv")
    if embed_matrix is not None and embed_matrix.shape[1] >= 2:
        obsm["X_emb"] = embed_matrix[:, :2]
    vf_path = target / "variable_features.txt"
    var_features = _read_lines_file(vf_path) if vf_path.exists() else []
    payload["obsm"] = obsm
    payload["var_features"] = [str(item).strip() for item in var_features if str(item).strip()]
    payload["import_state"] = meta.get("import_state") or {}
    payload["embed_kind"] = (payload["import_state"] or {}).get("embed_kind") or "umap"
    return payload


# ---------------------------------------------------------------------------
# 压缩包 / URL
# ---------------------------------------------------------------------------


def _load_archive(path, workdir, options=None, progress=None):
    options = options or {}
    root = _extract_archive(path, workdir, progress=progress)
    member = options.get("member")
    notes = ["已从压缩包 %s 解压" % Path(path).name]
    if member:
        target = _safe_member_path(root, member)
        if not target.exists():
            raise _err("bad_request", "压缩包里没有该成员：%s" % member)
        if target.is_dir():
            found = _find_10x_dir(target)
            if found is not None:
                target = found
            else:
                candidates = _archive_candidates(target)
                if not candidates:
                    raise _err("bad_request", "目录里没有可识别的数据文件：%s" % member)
                target = target / candidates[0]["member"]
        result = _load_path(target, workdir, options)
        result["notes"] = notes + list(result.get("notes") or [])
        return result
    found = _find_10x_dir(root)
    if found is not None:
        result = _load_10x_mtx(found)
        result["notes"] = notes + list(result.get("notes") or [])
        return result
    candidates = _archive_candidates(root)
    if not candidates:
        raise _err("bad_request", "压缩包里没有可识别的数据文件：%s" % Path(path).name)
    if len(candidates) > 1:
        raise _err("ambiguous",
                   "压缩包里有多个候选数据文件，请用 member 参数指定：%s" % "、".join(entry["member"] for entry in candidates[:12]),
                   candidates=[{key: entry.get(key) for key in ("member", "name", "size", "format")} for entry in candidates])
    result = _load_path(root / candidates[0]["member"], workdir, options)
    result["notes"] = notes + list(result.get("notes") or [])
    return result


# ---------------------------------------------------------------------------
# 格式识别与加载（统一入口）
# ---------------------------------------------------------------------------


def _resolve_source(source, workdir, options=None, progress=None):
    options = options or {}
    source = str(source).strip()
    if _looks_like_gse(source):
        listing = _gse_candidates(source)
        member = options.get("member")
        if not member:
            raise _err("ambiguous", listing["notes"][1], candidates=listing["candidates"])
        match = next((item for item in listing["candidates"] if item["member"] == member), None)
        if match is None:
            raise _err("bad_request", "GEO 目录里没有该文件：%s" % member, candidates=listing["candidates"])
        downloaded = _download(match["url"], workdir, options, progress)
        return downloaded, ["GEO %s 的补充文件：%s" % (source.upper(), member)]
    if _looks_like_url(source):
        downloaded = _download(source, workdir, options, progress)
        return downloaded, ["已从 URL 下载：%s" % source]
    path = Path(os.path.expanduser(source))
    if not path.exists():
        raise _err("not_found", "文件或目录不存在：%s" % source)
    return path, []


def _looks_archive(path):
    name = Path(path).name.lower()
    if name.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar")):
        return True
    if name.endswith(".zip"):
        return True
    try:
        return zipfile.is_zipfile(path)
    except Exception:
        return False


def _decompress_single(path, workdir):
    path = Path(path)
    name = path.name
    lowered = name.lower()
    for suffix in (".gz", ".bz2"):
        if lowered.endswith(suffix):
            inner_name = name[: -len(suffix)]
            inner_name = re.sub(r"\.(txt|tsv|csv)$", r"\1", inner_name)
            dirs = _workdirs(workdir)
            digest = hashlib.sha1(str(path).encode() + str(path.stat().st_mtime_ns).encode()).hexdigest()[:10]
            target = dirs["tmp"] / ("%s-%s" % (digest, inner_name))
            if not target.exists():
                opener = gzip.open if suffix == ".gz" else bz2.open
                with opener(path, "rb") as source, open(target, "wb") as out:
                    shutil.copyfileobj(source, out, 1024 * 1024)
            return target
    return path


def _detect_format(path, options=None):
    options = options or {}
    forced = options.get("format")
    if forced and str(forced).lower() not in ("auto", "自动"):
        return str(forced).lower()
    path = Path(path)
    if path.is_dir():
        if _find_10x_dir(path) is not None:
            return "10x_mtx"
        data_files = [item for item in path.iterdir() if item.is_file() and not item.name.startswith(".")]
        if len(data_files) == 1:
            return _detect_format(data_files[0], options)
        for priority in (".h5ad", ".h5", ".loom", ".rds", ".rdata", ".rda", ".csv", ".tsv", ".txt"):
            matches = [item for item in data_files if item.name.lower().endswith(priority)]
            if matches:
                return _detect_format(matches[0], options)
        raise _err("bad_request", "目录里没有可识别的数据文件（也没有 matrix.mtx 三件套）：%s" % path)
    lowered = path.name.lower()
    if _looks_archive(path):
        return "archive"
    if lowered.endswith((".h5ad", ".h5ad.gz")):
        return "h5ad"
    if lowered.endswith((".rds", ".rds.gz", ".rds.bz2")):
        return "rds"
    if lowered.endswith((".rdata", ".rda", ".rdata.gz", ".rda.gz")):
        return "rdata"
    if lowered.endswith(".mtx") or lowered.endswith(".mtx.gz"):
        return "10x_mtx"
    try:
        with open(path, "rb") as handle:
            head = handle.read(8)
    except Exception as error:
        raise _err("bad_request", "无法读取文件（%s）：%s" % (path.name, error))
    if head[:4] == b"\x1f\x8b" or head[:3] == b"BZh":
        inner = _decompress_single(path, Path(path).parent)
        return _detect_format(inner, options)
    if head == _HDF5_MAGIC:
        sniffed = _sniff_h5(path)
        if sniffed:
            return sniffed
        raise _err("unsupported", "HDF5 文件类型无法识别（既不是 h5ad，也不是 10x .h5 或 loom）：%s" % path.name)
    if lowered.endswith((".csv", ".tsv", ".txt", ".tab", ".csv.gz", ".tsv.gz", ".txt.gz", ".tab.gz")):
        return "matrix_text"
    if lowered.endswith((".loom", ".loom.gz")):
        return "loom"
    if lowered.endswith((".h5", ".hdf5")):
        sniffed = _sniff_h5(path)
        if sniffed:
            return sniffed
        raise _err("unsupported", "HDF5 文件里没有找到单细胞矩阵（期望 h5ad/10x/loom）：%s" % path.name)
    # 靠内容嗅探兜底
    try:
        head_text = path.read_text(encoding="utf-8", errors="replace")[:200000]
    except Exception:
        head_text = ""
    if "!series_matrix_table_begin" in head_text:
        return "geo_series"
    if head_text and any(delimiter in head_text for delimiter in (",", "\t")):
        return "matrix_text"
    if head[:2] == b"X\n" or head[:4] == b"RDX2" or head[:4] == b"RDX3":
        return "rdata" if head[:3] == b"RDX" else "rds"
    raise _err("unsupported", "无法识别的文件格式：%s（请用 format 参数指定）" % path.name)


def _load_path(path, workdir, options=None, progress=None):
    options = dict(options or {})
    path = Path(path)
    if path.is_dir():
        found = _find_10x_dir(path)
        if found is not None:
            return _load_10x_mtx(found, options)
        try:
            entries = [item for item in path.iterdir() if item.is_file() and not item.name.startswith(".")]
        except Exception as error:
            raise _err("bad_request", "无法读取目录（%s）：%s" % (path, error))
        if not entries:
            raise _err("bad_request", "目录是空的：%s" % path)
        if len(entries) == 1:
            return _load_path(entries[0], workdir, options, progress)
        candidates = []
        for entry in entries:
            fmt = _format_by_suffix(entry.name)
            if fmt is None:
                continue
            candidates.append({"member": entry.name, "name": entry.name,
                               "size": entry.stat().st_size, "format": fmt})
        if len(candidates) == 1:
            return _load_path(path / candidates[0]["member"], workdir, options, progress)
        if not candidates:
            raise _err("bad_request", "目录里没有可识别的数据文件：%s" % path)
        raise _err("ambiguous", "目录里有多个数据文件，请用 member 参数指定（传相对路径）",
                   candidates=candidates)
    name = path.name.lower()
    if name.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar")) or name.endswith(".zip"):
        return _load_archive(path, workdir, options, progress)
    fmt = _detect_format(path, options)
    if fmt == "h5ad":
        candidate = _decompress_single(path, workdir) if name.endswith(".gz") else path
        return _load_h5ad(candidate, options)
    if fmt == "10x_h5":
        return _load_10x_h5(path, options)
    if fmt == "loom":
        candidate = _decompress_single(path, workdir) if name.endswith(".gz") else path
        return _load_loom(candidate, options)
    if fmt in ("10x_mtx", "mtx"):
        return _load_10x_mtx(path, options)
    if fmt == "matrix_text":
        if name.endswith((".gz", ".bz2")):
            return _load_text(path, options)
        return _load_text(path, options)
    if fmt == "geo_series":
        return _load_text(path, options)
    if fmt in ("rds", "rdata"):
        candidate = path
        if name.endswith((".gz", ".bz2")):
            candidate = _decompress_single(path, workdir)
        return _load_r(candidate, workdir, options)
    if fmt == "archive":
        return _load_archive(path, workdir, options, progress)
    # 未知格式：尝试按 HDF5/文本最后兜底
    sniffed = _sniff_h5(path)
    if sniffed:
        return _load_path(path, workdir, {"format": sniffed, **options}, progress)
    if "!series_matrix_table_begin" in _read_head(path, 200000):
        return _load_text(path, options)
    raise _err("unsupported", "不支持的格式：%s" % path.name)


# ---------------------------------------------------------------------------
# inspect
# ---------------------------------------------------------------------------


def _source_size(path):
    path = Path(path)
    try:
        if path.is_dir():
            return int(sum(item.stat().st_size for item in path.rglob("*") if item.is_file()))
        return int(path.stat().st_size)
    except Exception:
        return 0


def _dir_candidates(path):
    entries = [item for item in Path(path).iterdir() if item.is_file() and not item.name.startswith(".")]
    candidates = []
    for entry in entries:
        fmt = _format_by_suffix(entry.name)
        if fmt:
            candidates.append({"member": entry.name, "name": entry.name,
                               "size": int(entry.stat().st_size), "format": fmt})
    return candidates


def _inspect_10x_dir(directory):
    directory = Path(directory)
    mtx = _sidecar(directory, "matrix", ("matrix.mtx.gz", "matrix.mtx", "matrix.mtx.bz2"))
    n_obs = None
    n_vars = None
    if mtx is not None:
        try:
            if mtx.name.lower().endswith(".gz"):
                handle = gzip.open(mtx, "rt", encoding="utf-8", errors="replace")
            else:
                handle = open(mtx, "rt", encoding="utf-8", errors="replace")
            with handle:
                for line in handle:
                    if line.startswith("%"):
                        continue
                    parts = line.split()
                    if len(parts) == 3:
                        n_vars = int(parts[0])
                        n_obs = int(parts[1])
                    break
        except Exception:
            pass
    return {"format": "10x_mtx", "n_obs": n_obs, "n_vars": n_vars,
            "files": [item.name for item in sorted(directory.iterdir()) if item.is_file()][:12]}


def inspect_source(source, *, workdir, options=None):
    options = dict(options or {})
    resolved, notes = _resolve_source(source, workdir, options)
    path = Path(resolved)
    size = _source_size(path)
    if path.is_dir():
        found = _find_10x_dir(path)
        if found is not None:
            info = _inspect_10x_dir(found)
            info.update({"source": str(path), "notes": notes + ["目录里找到 10x 三件套"], "size": size})
            return info
        candidates = _dir_candidates(path)
        if len(candidates) == 1:
            inner = inspect_source(path / candidates[0]["member"], workdir=workdir, options=options)
            inner["source"] = str(path)
            inner["size"] = size
            return inner
        if candidates:
            return {"format": "directory", "source": str(path), "notes": notes,
                    "candidates": candidates, "size": size}
        raise _err("bad_request", "目录里没有可识别的数据文件：%s" % path)
    if _looks_archive(path):
        root = _extract_archive(path, workdir)
        found = _find_10x_dir(root)
        if found is not None:
            info = _inspect_10x_dir(found)
            info.update({"source": str(path), "size": size,
                         "notes": notes + ["压缩包里找到 10x 三件套：%s" % found.name]})
            return info
        candidates = _archive_candidates(root)
        if len(candidates) == 1:
            inner = inspect_source(root / candidates[0]["member"], workdir=workdir, options=options)
            inner["source"] = str(path)
            inner["size"] = size
            inner["notes"] = notes + ["压缩包内文件：%s" % candidates[0]["member"]] + list(inner.get("notes") or [])
            return inner
        return {"format": "archive", "source": str(path), "notes": notes,
                "candidates": candidates, "size": size}
    name = path.name.lower()
    if name.endswith((".rds", ".rdata", ".rda", ".rds.gz", ".rdata.gz", ".rda.gz")):
        candidate = _decompress_single(path, workdir) if name.endswith((".gz", ".bz2")) else path
        fmt = "rdata" if name.endswith((".rdata", ".rda", ".rdata.gz", ".rda.gz")) else "rds"
        try:
            info = _load_r(candidate, workdir, {**options, "inspect": True})
            info["source"] = str(path)
            info["size"] = size
            info["notes"] = notes + list(info.get("notes") or [])
            return info
        except EngineError as error:
            if error.code == "ambiguous":
                return {"format": fmt, "source": str(path), "size": size,
                        "notes": notes + ["R 文件里有多个可用的单细胞对象，请用 object 参数指定"],
                        "candidates": error.candidates or []}
            raise
    fmt = _detect_format(path, options)
    result = {"format": fmt, "source": str(path), "size": size, "notes": notes}
    if fmt in ("10x_h5", "h5ad", "loom"):
        sniff_target = path
        if name.endswith((".gz", ".bz2")) and fmt in ("h5ad", "loom"):
            try:
                sniff_target = _decompress_single(path, workdir)
            except Exception:
                sniff_target = path
        shape = _h5_shape(sniff_target) or {}
        for key, value in shape.items():
            if value is not None:
                result[key] = value
    elif fmt == "10x_mtx":
        info = _inspect_10x_dir(path if path.is_dir() else path.parent)
        for key, value in info.items():
            if key != "format" and key not in result:
                result[key] = value
    elif fmt in ("matrix_text", "geo_series"):
        head = _read_head(path, 2 * 1024 * 1024)
        if "!series_matrix_table_begin" in head:
            result["format"] = "geo_series"
        first = next((line for line in head.splitlines() if line.strip()), "")
        delimiter = _split_delimiter([first])
        result["notes"] = result["notes"] + ["分隔符：%s" % ("TAB" if delimiter == "\t" else delimiter)]
    return result


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------


def load_source(source, *, workdir, options=None, progress=None):
    options = dict(options or {})
    resolved, notes = _resolve_source(source, workdir, options, progress)
    if progress:
        progress("正在识别格式：%s" % Path(resolved).name, 20)
    result = _load_path(Path(resolved), workdir, options, progress)
    result["notes"] = list(notes) + list(result.get("notes") or [])
    if options.get("name"):
        result["name"] = str(options["name"])
    if progress:
        progress("解析完成：%d 细胞 × %d 基因" % (result["meta"]["n_obs"], result["meta"]["n_vars"]), 100)
    return result


if __name__ == "__main__":  # 手动调试用：python3 lib/sc_io.py <path|url>
    import argparse

    parser = argparse.ArgumentParser(description="dsh-singlecell import layer debug CLI")
    parser.add_argument("source")
    parser.add_argument("--workdir", default=tempfile.mkdtemp(prefix="sc-io-"))
    parser.add_argument("--format", default=None)
    parser.add_argument("--member", default=None)
    parser.add_argument("--object", default=None)
    parser.add_argument("--inspect", action="store_true")
    args = parser.parse_args()
    opts = {key: value for key, value in
            {"format": args.format, "member": args.member, "object": args.object}.items() if value}
    if args.inspect:
        print(json.dumps(inspect_source(args.source, workdir=args.workdir, options=opts),
                         ensure_ascii=False, indent=2))
    else:
        out = load_source(args.source, workdir=args.workdir, options=opts,
                          progress=lambda message, pct=None: print("[%s%%] %s" % (pct if pct is not None else "-", message)))
        print(json.dumps({key: (value if key not in ("X", "counts", "obs", "var") else "…")
                          for key, value in out.items()}, ensure_ascii=False, indent=2, default=str))
        print("shape:", out["X"].shape, "format:", out["format"], "name:", out["name"])
