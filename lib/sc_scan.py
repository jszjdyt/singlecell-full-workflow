#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""扫描 GEO 下载的数据根目录，自动识别每个样本（对齐「单细胞分析器」的输入逻辑）。

识别规则：
  1. 子目录级：10X 目录（barcodes+features+matrix）、单样本目录（h5/h5ad/rds/RData/矩阵）
  2. GEO 三联体：同目录下按 GSM 前缀配对的 matrix/barcodes/features（可多样本）
  3. 散放文件：.h5(10X) / .h5ad / .rds / .RData / .csv/.tsv/.txt 矩阵
输出统一样本表格：{id, name, raw_name, format, format_label, path, detail, trio_prefix}

作者参考逻辑来自用户自研「单细胞分析器」（backend.scan_data_root），保持行为一致：
  - 10X 注释边车文件（barcodes/features/genes）绝不当作表达矩阵；
  - GEO 三联体按 GSM 前缀配对成样本；
  - 同名样本自动加后缀去重；未识别到样本时给可读错误。
"""

from __future__ import annotations

import gzip
import os
import re
from pathlib import Path

MTX_RE = re.compile(r"(^|[_\-.])(matrix)\.mtx(\.gz)?$", re.I)
BC_RE = re.compile(r"barcodes(\.tsv)?(\.gz)?$", re.I)
FEAT_RE = re.compile(r"(features|genes)(\.tsv)?(\.gz)?$", re.I)
TRIO_PREFIX_RE = re.compile(r"^(.*?)[_\-.]?(matrix|barcodes|features|genes)(\.tsv|\.mtx)?(\.gz)?$", re.I)
MATRIX_EXTS = {".csv", ".tsv", ".txt"}
HDF5_SIG = b"\x89HDF\r\n\x1a\n"

FORMAT_LABELS = {
    "10x_dir": "10X 目录",
    "10x_h5": "10X h5",
    "h5ad": "h5ad (Python)",
    "geo_mtx": "GEO 三联体",
    "expr_matrix": "表达矩阵",
    "seurat_rds": "Seurat RDS",
    "rdata": "RData (R)",
    "loom": "Loom",
    "geo_series": "GEO series matrix",
}


class ScanError(Exception):
    def __init__(self, code, message):
        self.code = str(code)
        self.message = str(message)
        super().__init__(self.message)


# ---------------------------------------------------------------------------
# 基础判断
# ---------------------------------------------------------------------------


def _first_line(path):
    try:
        if str(path).lower().endswith((".gz", ".bz2")):
            opener = gzip.open if str(path).lower().endswith(".gz") else __import__("bz2").open
            with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.strip():
                        return line
        else:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.strip():
                        return line
    except Exception:
        return None
    return None


def _is_hdf5(path):
    try:
        with open(path, "rb") as handle:
            return handle.read(8) == HDF5_SIG
    except Exception:
        return False


def _detect_delimiter(line):
    if line and line.count("\t") > line.count(","):
        return "\t"
    return ","


def _looks_like_matrix(line):
    if not line:
        return False
    line = line.strip()
    if line.startswith("%%MatrixMarket"):
        return False
    sep = _detect_delimiter(line)
    fields = [field for field in line.split(sep) if field.strip()]
    return len(fields) >= 2


def _matrix_stem(path):
    name = path.name
    for suffix in (".csv.gz", ".tsv.gz", ".txt.gz"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _split_trio_prefix(name):
    match = TRIO_PREFIX_RE.match(name)
    if match:
        return match.group(1).rstrip("_-.")
    return name


def _collect_geo_trios(files):
    """按 GSM 前缀把 matrix/barcodes/features 配对成样本：[(prefix, [files])]。"""
    by_prefix = {}
    for item in files:
        path = Path(item)
        if not path.is_file():
            continue
        if MTX_RE.search(path.name) or BC_RE.search(path.name) or FEAT_RE.search(path.name):
            by_prefix.setdefault(_split_trio_prefix(path.name), []).append(path)
    bundles = []
    for prefix, group in by_prefix.items():
        kinds = set()
        for path in group:
            if MTX_RE.search(path.name):
                kinds.add("mtx")
            elif BC_RE.search(path.name):
                kinds.add("bc")
            elif FEAT_RE.search(path.name):
                kinds.add("feat")
        if {"mtx", "bc", "feat"} <= kinds:
            bundles.append((prefix, group))
    return bundles


def _classify_dir(directory):
    """单个样本目录：返回 (format, detail, data_path)。"""
    directory = Path(directory)
    if not directory.is_dir():
        return None, None, None
    try:
        files = sorted(item for item in directory.iterdir() if item.is_file() and not item.name.startswith("."))
    except Exception:
        return None, None, None

    trio_files = {"mtx": [], "bc": [], "feat": []}
    for path in files:
        if MTX_RE.search(path.name):
            trio_files["mtx"].append(path)
        elif BC_RE.search(path.name):
            trio_files["bc"].append(path)
        elif FEAT_RE.search(path.name):
            trio_files["feat"].append(path)
    if trio_files["mtx"] and trio_files["bc"] and trio_files["feat"]:
        prefixes = {_split_trio_prefix(path.name) for group in trio_files.values() for path in group}
        if len(prefixes) == 1:
            prefix = prefixes.pop()
            label = FORMAT_LABELS["10x_dir"] if not prefix else FORMAT_LABELS["geo_mtx"]
            detail = "%d 个文件（barcodes/features/matrix%s）" % (
                len(files), "" if not prefix else "，前缀 %s" % prefix)
            return ("10x_dir" if not prefix else "geo_mtx"), detail, directory

    for path in files:
        lowered = path.name.lower()
        if lowered.endswith(".h5") and _is_hdf5(path):
            return "10x_h5", path.name, path
        if lowered.endswith(".h5ad"):
            return "h5ad", path.name, path
        if lowered.endswith(".loom"):
            return "loom", path.name, path
    for path in files:
        lowered = path.name.lower()
        if lowered.endswith(".rds"):
            return "seurat_rds", path.name, path
        if lowered.endswith((".rdata", ".rda")):
            return "rdata", path.name, path
    for path in files:
        lowered = path.name.lower()
        if path.suffix.lower() in MATRIX_EXTS or lowered.endswith((".csv.gz", ".tsv.gz", ".txt.gz")):
            line = _first_line(path)
            if line and "!series_matrix_table_begin" in line:
                continue
            if _looks_like_matrix(line):
                sep = _detect_delimiter(line)
                return "expr_matrix", "%s（%s分隔）" % (path.name, "TAB" if sep == "\t" else "逗号"), path
    return None, None, None


# ---------------------------------------------------------------------------
# 扫描主逻辑
# ---------------------------------------------------------------------------


def scan_data_root(root_dir, recursive=False):
    root = Path(str(root_dir)).expanduser()
    if not root.is_dir():
        return {"samples": [], "errors": ["目录不存在：%s" % root], "root": str(root)}

    samples = []
    errors = []

    def push(name, fmt, data_path, detail, trio_prefix=None):
        samples.append({
            "id": "",
            "name": name,
            "raw_name": name,
            "format": fmt,
            "format_label": FORMAT_LABELS.get(fmt, fmt),
            "path": str(data_path),
            "detail": detail,
            "trio_prefix": trio_prefix,
        })

    def is_sidecar(name):
        return bool(BC_RE.search(name) or FEAT_RE.search(name) or MTX_RE.search(name))

    def push_loose_files(directory, trio_prefixes):
        try:
            entries = sorted(item for item in Path(directory).iterdir()
                             if item.is_file() and not item.name.startswith("."))
        except Exception:
            return
        for path in entries:
            lowered = path.name.lower()
            if _split_trio_prefix(path.name) in trio_prefixes:
                continue
            if lowered.endswith(".rds"):
                push(_matrix_stem(path), "seurat_rds", path, path.name)
            elif lowered.endswith(".h5ad"):
                push(_matrix_stem(path), "h5ad", path, path.name)
            elif lowered.endswith(".loom"):
                push(_matrix_stem(path), "loom", path, path.name)
            elif lowered.endswith((".rdata", ".rda")):
                push(_matrix_stem(path), "rdata", path, path.name)
            elif lowered.endswith(".h5") and _is_hdf5(path):
                push(_matrix_stem(path), "10x_h5", path, path.name)
            elif path.suffix.lower() in MATRIX_EXTS or lowered.endswith((".csv.gz", ".tsv.gz", ".txt.gz")):
                if is_sidecar(path.name):
                    continue
                line = _first_line(path)
                if line and "!series_matrix_table_begin" in line:
                    push(_matrix_stem(path), "geo_series", path, path.name)
                    continue
                if line and _looks_like_matrix(line):
                    sep = _detect_delimiter(line)
                    push(_matrix_stem(path), "expr_matrix", path,
                         "%s（%s分隔）" % (path.name, "TAB" if sep == "\t" else "逗号"))

    if recursive:
        subdirs = []
        for current, dirnames, _filenames in os.walk(root):
            dirnames[:] = sorted(item for item in dirnames if not item.startswith("."))
            current_path = Path(current)
            if current_path != root:
                subdirs.append(current_path)
    else:
        try:
            subdirs = sorted(item for item in root.iterdir() if item.is_dir() and not item.name.startswith("."))
        except Exception as error:
            return {"samples": [], "errors": ["无法读取目录：%s" % error], "root": str(root)}

    for subdir in subdirs:
        fmt, detail, data_path = _classify_dir(subdir)
        if fmt:
            name = subdir.name
            prefix = None
            if fmt == "geo_mtx":
                try:
                    bundles = _collect_geo_trios(subdir.iterdir())
                except Exception:
                    bundles = []
                if len(bundles) == 1:
                    name = bundles[0][0]
                    prefix = name
            push(name, fmt, data_path, detail, trio_prefix=prefix)
            continue
        try:
            bundles = sorted(_collect_geo_trios(subdir.iterdir()), key=lambda item: item[0])
        except Exception:
            bundles = []
        if bundles:
            for prefix, files in bundles:
                push(prefix, "geo_mtx", subdir, "%d 个文件（前缀 %s）" % (len(files), prefix), trio_prefix=prefix)
            push_loose_files(subdir, {prefix for prefix, _files in bundles})
            continue
        # 常见 GEO 结构：root/sampleA/filtered_feature_bc_matrix/…（非递归时也向下看一层）
        nested = 0
        try:
            children = sorted(item for item in subdir.iterdir() if item.is_dir() and not item.name.startswith("."))
        except Exception:
            children = []
        for child in children:
            child_fmt, child_detail, child_path = _classify_dir(child)
            if not child_fmt:
                continue
            nested += 1
            prefix = None
            if child_fmt == "geo_mtx":
                try:
                    child_bundles = _collect_geo_trios(child.iterdir())
                except Exception:
                    child_bundles = []
                if len(child_bundles) == 1:
                    prefix = child_bundles[0][0]
            push(subdir.name if nested == 1 else child.name, child_fmt, child_path,
                 child_detail or "", trio_prefix=prefix)
        if nested:
            push_loose_files(subdir, set())
            continue
        push_loose_files(subdir, set())

    try:
        root_files = [item for item in root.iterdir() if item.is_file() and not item.name.startswith(".")]
    except Exception:
        root_files = []
    bundles = sorted(_collect_geo_trios(root_files), key=lambda item: item[0])
    for prefix, files in bundles:
        push(prefix, "geo_mtx", root, "%d 个文件（前缀 %s）" % (len(files), prefix), trio_prefix=prefix)
    push_loose_files(root, {prefix for prefix, _files in bundles})

    # 去重 + 稳定排序 + 唯一默认名
    seen = set()
    unique = []
    for sample in sorted(samples, key=lambda item: (item["path"], item["name"], item["format"])):
        key = (sample["name"], sample["path"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(sample)
    used = set()
    for index, sample in enumerate(unique, 1):
        base = sample["name"] or ("sample%d" % index)
        name = base
        suffix = 2
        while name in used:
            name = "%s_%d" % (base, suffix)
            suffix += 1
        sample["id"] = "s%d" % index
        sample["name"] = name
        used.add(name)

    if not unique:
        errors.append(
            "未识别到任何单细胞数据样本。支持：10X 目录、h5/h5ad、"
            "GEO 三联体（matrix/barcodes/features）、csv/tsv 矩阵、rds、RData。"
        )
    return {"samples": unique, "errors": errors, "root": str(root)}


# ---------------------------------------------------------------------------
# GEO 查询（输入 GSE → 标题 + GSM 说明，辅助分组）
# ---------------------------------------------------------------------------


def _ssl_context():
    """Framework Python 常缺 CA 证书：优先 certifi，其次 macOS 系统证书，最后才降级。"""
    import ssl

    try:
        import certifi  # type: ignore
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        pass
    for candidate in ("/etc/ssl/cert.pem", "/etc/ssl/certs/ca-certificates.crt"):
        try:
            import os
            if os.path.exists(candidate):
                return ssl.create_default_context(cafile=candidate)
        except Exception:
            continue
    context = ssl.create_default_context()
    try:
        context.check_hostname = False
        context.verify_mode = ssl.CERT_NONE
    except Exception:
        pass
    return context


def geo_lookup(gse):
    import urllib.parse
    import urllib.request
    import xml.etree.ElementTree as ET

    gse = str(gse or "").strip().upper()
    if not re.fullmatch(r"GSE\d{3,}", gse):
        raise ScanError("bad_request", "GSE 编号格式不正确（如 GSE246662）")
    try:
        search_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="
            + urllib.parse.quote(gse + "[ACCN]")
        )
        with urllib.request.urlopen(search_url, timeout=20, context=_ssl_context()) as response:
            search_root = ET.fromstring(response.read())
        ids = [node.text for node in search_root.iter("Id")]
        if not ids:
            raise ScanError("not_found", "GEO 未找到 %s" % gse)
        summary_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=" + ids[0]
        )
        with urllib.request.urlopen(summary_url, timeout=20, context=_ssl_context()) as response:
            root = ET.fromstring(response.read())
        title = ""
        samples = []
        for doc in root.iter("DocSum"):
            for item in doc.findall("Item"):
                if item.get("Name") == "title":
                    title = (item.text or "").strip()
                if item.get("Name") == "Samples":
                    for sample in item.findall("Item"):
                        accession = sample_title = None
                        for sub in sample.findall("Item"):
                            if sub.get("Name") == "Accession":
                                accession = sub.text
                            if sub.get("Name") == "Title":
                                sample_title = sub.text
                        if accession:
                            samples.append({"gsm": accession.strip(), "title": (sample_title or "").strip()})
        return {"gse": gse, "title": title, "samples": samples}
    except ScanError:
        raise
    except Exception as error:  # noqa: BLE001
        raise ScanError("bad_gateway", "GEO 查询失败：%s" % error)
