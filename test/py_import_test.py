#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell — 全格式导入层测试（CONTRACT §8.5，sc-io 所有）

对 `lib/sc_io.py` 的 `inspect_source` / `load_source` 跑全套格式回归：
  1. §8.2 每个格式 load 后与 truth.json 校验（n_obs/n_vars/前几个基因/首细胞某基因值/obs 列）；
  2. load_source 返回键名、X 类型（csr float32, cells×genes）、obs 行数、meta 固定键；
  3. inspect_source 格式判断 + size/n_obs/n_vars 正确性（不要求加载完整矩阵）；
  4. ambiguous 分支：多对象 RData / 多候选 zip → candidates → 用 object/member 重试成功；
  5. 坏文件 / 不存在路径 / 路径穿越 zip / 空压缩包 → 可读 Error 且不越界写文件；
  6. 本地 http.server URL 下载（含 max_bytes 上限）、transpose 强制参数；
  7. 真实 Seurat rds（neu_analysis/objects/*.rds，只读）逐个 load 并打印 n_obs/n_vars/obs 列数。

运行（必须用绝对路径 python3）：
    /Library/Frameworks/Python.framework/Versions/3.14/bin/python3 test/py_import_test.py

全部通过退出码 0，任何失败退出码 1 并打印可复现细节。
环境变量 DS_SC_KEEP_FIXTURES=1 时保留 /tmp 临时目录便于排查。
"""

from __future__ import annotations

import importlib.util
import json
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

TEST_DIR = Path(__file__).resolve().parent
PLUGIN_DIR = TEST_DIR.parent
LIB_DIR = PLUGIN_DIR / "lib"
CONTRACT_CODES = {"bad_request", "not_found", "prerequisite", "unsupported", "internal", "ambiguous"}
REQUIRED_KEYS = {"format", "source", "name", "X", "counts", "obs", "var_names", "var", "warnings", "notes", "meta"}
META_KEYS = {"n_obs", "n_vars", "species", "assays", "has_counts"}
CONTRACT_FORMATS = {"h5ad", "10x_h5", "10x_mtx", "loom", "matrix_text", "geo_series",
                    "rds", "rdata", "archive", "auto"}


# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError("无法加载模块：%s" % path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _engine_error(exc):
    """把 sc_io.EngineError（或 Host 转换后的异常）归一成 (code, message, candidates)。"""
    code = getattr(exc, "code", None)
    if code not in CONTRACT_CODES:
        code = None
    message = getattr(exc, "message", None) or str(exc) or ""
    candidates = getattr(exc, "candidates", None)
    if candidates is None:
        details = getattr(exc, "details", None)
        if isinstance(details, dict):
            candidates = details.get("candidates")
    return code, str(message), candidates


class Reporter:
    def __init__(self):
        self.total = 0
        self.failed = []
        self.section = ""
        self.section_counts = {}
        self.started = time.time()

    def start(self, title):
        self.section = title
        self.section_counts[title] = [0, 0]
        print("\n" + "=" * 78)
        print("== %s" % title)
        print("=" * 78)

    def check(self, name, ok, detail=""):
        self.total += 1
        self.section_counts[self.section][0] += 1
        if ok:
            print("  PASS  %s%s" % (name, ("  [" + detail + "]") if detail else ""))
        else:
            self.section_counts[self.section][1] += 1
            self.failed.append((self.section, name, detail))
            print("  FAIL  %s" % name)
            if detail:
                for line in str(detail).splitlines():
                    print("        " + line)
        return ok

    def check_many(self, name, errors, detail=""):
        if errors:
            return self.check(name, False, ("; ".join(errors)) if not detail else detail + "; " + "; ".join(errors))
        return self.check(name, True, detail)

    def summary(self):
        elapsed = time.time() - self.started
        print("\n" + "=" * 78)
        for title, (count, failed) in self.section_counts.items():
            print("  %-42s %3d 项, 失败 %d" % (title, count, failed))
        print("-" * 78)
        print("  合计: %d 项检查，失败 %d 项，用时 %.1fs" % (self.total, len(self.failed), elapsed))
        if not self.failed:
            print("  结果: ALL PASS")
        else:
            print("  结果: FAIL（%d 项）" % len(self.failed))
            for section, name, detail in self.failed:
                first = str(detail).splitlines()[0] if detail else ""
                print("    - [%s] %s :: %s" % (section, name, first))
        print("=" * 78)
        return len(self.failed) == 0


# ---------------------------------------------------------------------------
# load 结果校验（对齐 CONTRACT §8.1）
# ---------------------------------------------------------------------------


def _locate(result, truth):
    obs = result.get("obs")
    var_names = list(result.get("var_names") or [])
    cell_idx = gene_idx = None
    if truth.get("cell_name") is not None and isinstance(obs, pd.DataFrame):
        try:
            cell_idx = list(obs.index).index(truth["cell_name"])
        except Exception:
            cell_idx = None
    if truth.get("gene_name") is not None:
        try:
            gene_idx = var_names.index(truth["gene_name"])
        except Exception:
            gene_idx = None
    if (cell_idx is None or gene_idx is None) and truth.get("names_generated"):
        cell_idx, gene_idx = 0, 0
    return cell_idx, gene_idx


def check_loaded(result, truth, *, need_format=True):
    """逐条校验 load_source 返回 dict；返回错误字符串列表（空 = 全对）。"""
    errors = []
    if not isinstance(result, dict):
        return ["load_source 返回值不是 dict：%r" % type(result)]
    missing = REQUIRED_KEYS - set(result)
    if missing:
        errors.append("缺少返回键 %s" % sorted(missing))
    n_obs, n_vars = int(truth["n_obs"]), int(truth["n_vars"])

    if need_format and truth.get("format") and result.get("format") != truth["format"]:
        errors.append("format=%r，应 %r" % (result.get("format"), truth["format"]))

    X = result.get("X")
    if not isinstance(X, sp.csr_matrix):
        errors.append("X 不是 scipy.sparse.csr_matrix（%s）" % type(X).__name__)
    else:
        if X.dtype != np.float32:
            errors.append("X.dtype=%s，应 float32" % X.dtype)
        if tuple(X.shape) != (n_obs, n_vars):
            errors.append("X.shape=%s，应 (%d, %d)" % (tuple(X.shape), n_obs, n_vars))
        if X.nnz and not np.all(np.isfinite(X.data)):
            errors.append("X 含 NaN/Inf")

    obs = result.get("obs")
    if not isinstance(obs, pd.DataFrame):
        errors.append("obs 不是 pandas.DataFrame（%s）" % type(obs).__name__)
    else:
        if len(obs) != n_obs:
            errors.append("obs 行数=%d，应 n_obs=%d" % (len(obs), n_obs))
        if obs.index.has_duplicates:
            errors.append("obs 索引有重复")
        if not all(isinstance(item, str) for item in obs.index[:200]):
            errors.append("obs 索引含非字符串")

    var_names = result.get("var_names")
    if not isinstance(var_names, list):
        errors.append("var_names 不是 list")
    else:
        if len(var_names) != n_vars:
            errors.append("var_names 长度=%d，应 n_vars=%d" % (len(var_names), n_vars))
        if any(not isinstance(item, str) for item in var_names[:200]):
            errors.append("var_names 含非字符串")
        if truth.get("genes") and not truth.get("names_generated"):
            prefix = list(var_names[: len(truth["genes"])])
            if prefix != list(truth["genes"]):
                errors.append("前几个基因名=%r，应 %r" % (prefix, truth["genes"]))

    if not isinstance(result.get("warnings"), list):
        errors.append("warnings 不是 list")
    if not isinstance(result.get("notes"), list):
        errors.append("notes 不是 list")

    meta = result.get("meta")
    if not isinstance(meta, dict):
        errors.append("meta 不是 dict")
    else:
        missing_meta = META_KEYS - set(meta)
        if missing_meta:
            errors.append("meta 缺少固定键 %s" % sorted(missing_meta))
        if meta.get("n_obs") != n_obs or meta.get("n_vars") != n_vars:
            errors.append("meta.n_obs/n_vars=%r/%r，应 %d/%d" % (
                meta.get("n_obs"), meta.get("n_vars"), n_obs, n_vars))
        species = meta.get("species")
        if species is not None and not isinstance(species, str):
            errors.append("meta.species 类型=%s，应 str|None" % type(species).__name__)
        if "assays" in meta and not isinstance(meta["assays"], (list, tuple)):
            errors.append("meta.assays 不是 list")
        if "has_counts" in meta and not isinstance(meta["has_counts"], bool):
            errors.append("meta.has_counts 不是 bool")

    cell_idx, gene_idx = _locate(result, truth)
    if isinstance(X, sp.csr_matrix) and tuple(X.shape) == (n_obs, n_vars):
        if truth.get("value") is not None:
            if cell_idx is None or gene_idx is None:
                if not truth.get("names_generated"):
                    errors.append("在结果里找不到 truth 细胞/基因名 %r / %r" % (
                        truth.get("cell_name"), truth.get("gene_name")))
            else:
                value = float(X[cell_idx, gene_idx])
                if not np.isclose(value, float(truth["value"]), atol=1e-3, rtol=1e-4):
                    errors.append("X[%r, %r]=%s，应 %s" % (
                        truth.get("cell_name"), truth.get("gene_name"), value, truth["value"]))
        if truth.get("counts_value") is not None:
            counts = result.get("counts")
            if counts is None:
                errors.append("counts 为 None，但 truth 期望存在独立 counts 层")
            elif not sp.issparse(counts):
                errors.append("counts 不是稀疏矩阵（%s）" % type(counts).__name__)
            elif counts.dtype != np.float32:
                errors.append("counts.dtype=%s，应 float32" % counts.dtype)
            elif tuple(counts.shape) != (n_obs, n_vars):
                errors.append("counts.shape=%s，应 (%d, %d)" % (tuple(counts.shape), n_obs, n_vars))
            elif cell_idx is not None and gene_idx is not None:
                counts_value = float(counts[cell_idx, gene_idx])
                if not np.isclose(counts_value, float(truth["counts_value"]), atol=1e-3, rtol=1e-4):
                    errors.append("counts[%r, %r]=%s，应 %s" % (
                        truth.get("cell_name"), truth.get("gene_name"),
                        counts_value, truth["counts_value"]))

    if isinstance(obs, pd.DataFrame):
        for column in truth.get("obs_columns") or []:
            if column not in obs.columns:
                errors.append("obs 缺列 %r（现有 %s）" % (column, list(obs.columns)[:10]))
        any_names = truth.get("obs_columns_any") or []
        if any_names:
            lowered = [str(column).lower() for column in obs.columns]
            if not any(any(token.lower() in column for column in lowered) for token in any_names):
                errors.append("obs 列里没有 %s 之一（现有 %s）" % (any_names, list(obs.columns)[:10]))

    if truth.get("names_generated"):
        warnings = result.get("warnings") or []
        if not any("生成" in str(item) for item in warnings):
            errors.append("名字是自动生成的，但 warnings 没有说明（warnings=%r）" % (warnings,))
    return errors


def check_inspect(result, entry):
    """inspect_source 结果校验（§8.1：format 必对；size:int；n_obs/n_vars 若给必须对）。"""
    errors = []
    truth = entry["truth"]
    if not isinstance(result, dict):
        return ["inspect_source 返回值不是 dict：%r" % type(result)]
    expected = {entry["format"]}
    if entry.get("inspect_format"):
        expected.add(entry["inspect_format"])
    if result.get("format") not in expected:
        errors.append("format=%r，应在 %r 中" % (result.get("format"), sorted(expected)))
    if not isinstance(result.get("notes"), list):
        errors.append("notes 不是 list")
    if "size" not in result:
        errors.append("缺少 size 键（§8.1）")
    elif not isinstance(result["size"], int) or isinstance(result["size"], bool) or result["size"] < 0:
        errors.append("size=%r 不是非负 int" % (result.get("size"),))
    for key in ("n_obs", "n_vars"):
        if result.get(key) is not None:
            if result[key] != truth[key]:
                errors.append("inspect.%s=%r，应 %r" % (key, result[key], truth[key]))
    candidates = result.get("candidates")
    if candidates is not None and not isinstance(candidates, list):
        errors.append("candidates 不是 list")
    return errors


# ---------------------------------------------------------------------------
# 各测试段
# ---------------------------------------------------------------------------


def test_contract_surface(io, manifest, work, report):
    report.start("0. 契约表面（模块常量 / 函数签名 / options）")
    try:
        supported = set(io.SUPPORTED_FORMATS)
        missing = CONTRACT_FORMATS - supported
        report.check("SUPPORTED_FORMATS 含 §8.1 全部格式", not missing,
                     "" if not missing else "缺 %s（实际 %s）" % (sorted(missing), sorted(supported)))
    except Exception as error:  # pragma: no cover
        report.check("SUPPORTED_FORMATS 存在", False, "%s: %s" % (type(error).__name__, error))

    report.check("inspect_source / load_source 可调用",
                 callable(getattr(io, "inspect_source", None)) and callable(getattr(io, "load_source", None)))

    by_id = {entry["id"]: entry for entry in manifest["formats"]}
    entry = by_id["text_csv"]

    # options.names：load_source 结果 name 应被覆盖
    result = io.load_source(entry["path"], workdir=str(work), options={"name": "my-dataset"})
    report.check("options['name'] 覆盖结果 name", result.get("name") == "my-dataset",
                 "name=%r" % (result.get("name"),))

    # options.format='auto' / 强制格式
    result = io.load_source(entry["path"], workdir=str(work), options={"format": "matrix_text"})
    report.check_many("强制 options['format']=matrix_text", check_loaded(result, entry["truth"]))
    result = io.load_source(entry["path"], workdir=str(work), options={"format": "auto"})
    report.check_many("options['format']=auto 自动识别", check_loaded(result, entry["truth"]))

    # options.delimiter
    semi = by_id["text_semicolon"]
    result = io.load_source(semi["path"], workdir=str(work), options={"delimiter": ";"})
    report.check_many("options['delimiter'] 强制分隔符", check_loaded(result, semi["truth"]))

    # progress 回调（1/2 参兼容）
    calls = []

    def progress(*args, **kwargs):
        calls.append(args)

    result = io.load_source(entry["path"], workdir=str(work), progress=progress)
    report.check("progress 回调可用", bool(calls) and not check_loaded(result, entry["truth"]),
                 "progress 调用 %d 次" % len(calls))

    # workdir 为 Path 对象
    result = io.load_source(entry["path"], workdir=work)
    report.check_many("workdir 支持 pathlib.Path", check_loaded(result, entry["truth"]))


def test_formats(io, manifest, work, report):
    report.start("1. §8.2 逐格式 load（truth.json 校验）")
    for entry in manifest["formats"]:
        fid = entry["id"]
        truth = entry["truth"]
        try:
            result = io.load_source(entry["path"], workdir=str(work), options=entry.get("options") or {})
        except BaseException as error:
            code, message, _ = _engine_error(error)
            if code is None:
                detail = "抛了非 EngineError：%s: %s\n%s" % (
                    type(error).__name__, error, traceback.format_exc(limit=3))
            else:
                detail = "EngineError(%s): %s" % (code, message)
            report.check("load %-22s %-11s" % (fid, truth["format"]), False, detail)
            continue
        errors = check_loaded(result, truth)
        report.check_many(
            "load %-22s %-11s" % (fid, truth["format"]),
            errors,
            "%d×%d，%s" % (truth["n_obs"], truth["n_vars"], entry["path"]),
        )


def test_inspect(io, manifest, work, report):
    report.start("2. inspect_source（格式判断 / size / n_obs·n_vars 一致性）")
    for entry in manifest["formats"]:
        try:
            result = io.inspect_source(entry["path"], workdir=str(work), options=entry.get("options") or {})
        except BaseException as error:
            code, message, _ = _engine_error(error)
            detail = ("EngineError(%s): %s" % (code, message)) if code else \
                "抛了非 EngineError：%s: %s" % (type(error).__name__, error)
            report.check("inspect %-22s" % entry["id"], False, detail)
            continue
        errors = check_inspect(result, entry)
        dims = ""
        if result.get("n_obs") is not None:
            dims = "%s×%s" % (result.get("n_obs"), result.get("n_vars"))
        report.check_many("inspect %-22s %-11s" % (entry["id"], entry["format"]), errors,
                          (dims + " " if dims else "") + "size=%s" % (result.get("size"),))

    # 不存在的路径：inspect 应给 not_found 的 EngineError，而不是崩
    missing = str(Path(manifest["root"]) / "does_not_exist" / "nope.csv")
    try:
        io.inspect_source(missing, workdir=str(work))
        report.check("inspect 不存在路径给 EngineError", False, "未抛错")
    except BaseException as error:
        code, message, _ = _engine_error(error)
        report.check("inspect 不存在路径给 EngineError", code == "not_found" and bool(message),
                     "code=%r msg=%s" % (code, message[:100]))


def test_ambiguous(io, manifest, work, report):
    report.start("3. ambiguous 分支（多对象 RData / 多候选 zip → candidates → 重试）")
    for entry in manifest["ambiguous"]:
        fid = entry["id"]
        field = entry["candidates_field"]
        errors = []
        candidates = None
        try:
            io.load_source(entry["path"], workdir=str(work))
            errors.append("未抛 ambiguous 异常（应提示选择 %s）" % field)
        except BaseException as error:
            code, message, candidates = _engine_error(error)
            if code != "ambiguous":
                errors.append("code=%r，应 ambiguous（msg=%s）" % (code, message[:120]))
            if not candidates:
                errors.append("异常没有 candidates")
            else:
                for candidate in candidates:
                    if not isinstance(candidate, dict):
                        errors.append("candidate 不是 dict：%r" % (candidate,))
                        continue
                    if not candidate.get(field):
                        errors.append("candidate 缺 %r 定位字段：%r" % (field, candidate))
                    if field == "object":
                        if not candidate.get("class"):
                            errors.append("candidate 缺 class：%r" % (candidate,))
                        dim_info = candidate.get("dim") or (
                            [candidate.get("n_obs"), candidate.get("n_vars")]
                            if candidate.get("n_obs") and candidate.get("n_vars") else None)
                        if not dim_info:
                            errors.append("candidate 缺 dim 信息：%r" % (candidate,))
                    else:
                        if not isinstance(candidate.get("member"), str):
                            errors.append("candidate.member 不是 str：%r" % (candidate,))
        report.check_many("load   %-18s 抛 ambiguous + candidates" % fid, errors,
                          "candidates=%d" % (len(candidates or [])))

        # inspect 也要能拿到候选（返回 dict.candidates 或 EngineError ambiguous）
        inspect_errors = []
        inspect_candidates = None
        try:
            info = io.inspect_source(entry["path"], workdir=str(work))
            inspect_candidates = info.get("candidates")
            if not isinstance(info, dict):
                inspect_errors.append("inspect 返回不是 dict")
            if not inspect_candidates:
                inspect_errors.append("inspect 未返回 candidates")
        except BaseException as error:
            code, message, inspect_candidates = _engine_error(error)
            if code != "ambiguous" or not inspect_candidates:
                inspect_errors.append("inspect 抛了 %r 且无 candidates（%s）" % (code, message[:100]))
        report.check_many("inspect %-17s 暴露 candidates" % fid, inspect_errors,
                          "candidates=%d" % (len(inspect_candidates or [])))

        for retry in entry["retry"]:
            options = dict(retry["options"])
            label = ",".join("%s=%s" % (key, value) for key, value in options.items())
            try:
                result = io.load_source(entry["path"], workdir=str(work), options=options)
            except BaseException as error:
                code, message, _ = _engine_error(error)
                report.check("retry %s %s" % (fid, label), False,
                             "重试仍失败：%s(%s): %s" % (type(error).__name__, code, message))
                continue
            report.check_many("retry %s %s" % (fid, label), check_loaded(result, retry["truth"]),
                              "%d×%d" % (retry["truth"]["n_obs"], retry["truth"]["n_vars"]))


def test_errors(io, manifest, work, report):
    report.start("4. 错误分支（坏文件 / 不存在路径 / 路径穿越 / 空压缩包）")
    tempdir = Path(tempfile.gettempdir())
    for entry in manifest["errors"]:
        fid = entry["id"]
        allowed = set(entry["expect_codes"])
        errors = []
        raised = False
        try:
            result = io.load_source(entry["path"], workdir=str(work), options=entry.get("options") or {})
        except BaseException as error:
            raised = True
            code, message, _ = _engine_error(error)
            if code is None:
                errors.append("抛了非 EngineError 原始异常：%s: %s" % (type(error).__name__, error))
            else:
                if code not in allowed:
                    errors.append("code=%r 不在 %r" % (code, allowed))
                if not message.strip():
                    errors.append("错误消息为空")
                if len(message.strip()) < 4:
                    errors.append("错误消息过短：%r" % message)
        else:
            # 契约允许「路径穿越成员被跳过、正常成员继续加载」：只要没有越界写文件就接受
            if fid == "err_evil_zip" and result.get("format"):
                errors.append("")  # 占位，安全校验见下
                errors.pop()
            else:
                errors.append("未抛错（返回 format=%r）" % (result or {}).get("format"))
        # 安全不变量：任何 forbidden marker 都不能被写到目标目录之外
        for marker in entry.get("forbidden_markers") or []:
            hits = list(work.rglob(marker))
            if hits:
                errors.append("路径穿越写出了文件：%s" % hits[:3])
            for base in (tempdir, work.parent, Path(manifest["root"])):
                candidate = base / marker
                if candidate.exists():
                    errors.append("路径穿越写出了文件：%s" % candidate)
        detail = "已抛 EngineError" if raised else "已安全降级并继续加载（恶意成员被跳过）"
        report.check_many("error  %-18s %s" % (fid, entry["expect_codes"]), errors, detail)

    # archive 的 member 参数本身也不能路径穿越
    evil = next((entry for entry in manifest["errors"] if entry["id"] == "err_evil_zip"), None)
    if evil is not None:
        try:
            io.load_source(evil["path"], workdir=str(work), options={"member": "../../etc/passwd"})
            report.check("archive member=../../ 路径穿越被拒绝", False, "未报错")
        except BaseException as error:
            code, message, _ = _engine_error(error)
            report.check("archive member=../../ 路径穿越被拒绝",
                         code == "bad_request" and bool(message),
                         "code=%r msg=%s" % (code, message[:120]))


def test_transpose(io, manifest, work, report):
    report.start("5. matrix_text transpose 强制参数")
    for probe in manifest["transpose_probes"]:
        try:
            result_true = io.load_source(probe["path"], workdir=str(work), options={"transpose": True})
            result_false = io.load_source(probe["path"], workdir=str(work), options={"transpose": False})
        except BaseException as error:
            report.check("transpose 探针 %s" % probe["id"], False,
                         "%s: %s" % (type(error).__name__, error))
            continue
        errors = []
        shape_true = tuple(result_true["X"].shape)
        shape_false = tuple(result_false["X"].shape)
        rows, cols = probe["shape_rows"], probe["shape_cols"]
        expected_shapes = {(rows, cols), (cols, rows)}
        if {shape_true, shape_false} != expected_shapes:
            errors.append("transpose=True/False 形状=%r/%r，应各为 %r 之一且不同" % (
                shape_true, shape_false, [list(item) for item in expected_shapes]))
        try:
            dense_true = result_true["X"].toarray()
            dense_false = result_false["X"].toarray()
            if dense_true.shape == dense_false.shape:
                errors.append("transpose=True/False 形状相同=%r，未生效" % (dense_true.shape,))
            elif not np.allclose(dense_true, dense_false.T):
                errors.append("transpose=True/False 不是互为转置")
        except Exception as error:
            errors.append("转置内容校验失败：%s" % error)
        report.check_many("transpose 探针 %s" % probe["id"], errors,
                          "True=%r False=%r" % (shape_true, shape_false))


def test_urls(io, manifest, work, report, make_formats):
    report.start("6. URL 下载（本地 http.server，不联网）")
    with make_formats.serve(manifest["root"]) as base_url:
        for entry in manifest["url_tests"]:
            url = base_url + "/" + entry["rel"]
            try:
                result = io.load_source(url, workdir=str(work))
            except BaseException as error:
                code, message, _ = _engine_error(error)
                report.check("URL load %-12s" % entry["id"], False,
                             "%s(%s): %s" % (type(error).__name__, code, message))
                continue
            errors = check_loaded(result, entry["truth"])
            source = str(result.get("source") or "")
            if not (Path(source).exists() if source else False):
                errors.append("source 不是已下载的本地文件：%r" % source)
            report.check_many("URL load %-12s %-11s" % (entry["id"], entry["format"]), errors,
                              url)

        url = base_url + "/" + manifest["url_tests"][0]["rel"]
        try:
            info = io.inspect_source(url, workdir=str(work))
            report.check_many("URL inspect %-10s" % manifest["url_tests"][0]["id"],
                              check_inspect(info, manifest["url_tests"][0]))
        except BaseException as error:
            code, message, _ = _engine_error(error)
            report.check("URL inspect", False, "%s(%s): %s" % (type(error).__name__, code, message))

        # max_bytes 上限：小上限应给可读 bad_request，而不是写出超大文件
        try:
            io.load_source(url, workdir=str(work), options={"max_bytes": 32})
            report.check("URL max_bytes 上限生效", False, "max_bytes=32 仍未报错")
        except BaseException as error:
            code, message, _ = _engine_error(error)
            report.check("URL max_bytes 上限生效", code in ("bad_request", "unsupported") and bool(message),
                         "code=%r msg=%s" % (code, message[:100]))


def test_real_rds(io, manifest, work, report):
    report.start("7. 真实 Seurat rds 只读验收（可选：DSH_SC_TEST_RDS_DIR/*.rds）")
    files = list(manifest.get("real_rds", {}).get("files") or [])
    if not files:
        # 开源仓库里没有用户私有 rds：设置 DSH_SC_TEST_RDS_DIR 指向含 *.rds 的目录才会执行本段
        report.check("真实 rds 存在（可选，跳过）", True,
                     "未设置 DSH_SC_TEST_RDS_DIR，跳过真实 rds 验收（不影响其它格式）")
        return
    report.check("真实 rds 数量", len(files) >= 5, "找到 %d 个：%s" % (len(files), [Path(p).name for p in files]))
    loaded_info = []
    for index, path in enumerate(files):
        try:
            stat_before = os.stat(path)
            result = io.load_source(path, workdir=str(work))
            stat_after = os.stat(path)
        except BaseException as error:
            code, message, _ = _engine_error(error)
            report.check("真实 rds %-22s" % Path(path).name, False,
                         "%s(%s): %s" % (type(error).__name__, code, message))
            continue
        errors = []
        truth = {
            "n_obs": int(result["X"].shape[0]) if isinstance(result.get("X"), sp.csr_matrix) else -1,
            "n_vars": int(result["X"].shape[1]) if isinstance(result.get("X"), sp.csr_matrix) else -1,
            "genes": [],  # 基因名不预设，但要求长度一致
        }
        x = result.get("X")
        if not isinstance(x, sp.csr_matrix) or x.dtype != np.float32:
            errors.append("X 不是 csr float32：%r" % (getattr(x, "dtype", None),))
        if result.get("format") != "rds":
            errors.append("format=%r，应 rds" % (result.get("format"),))
        if truth["n_obs"] <= 0 or truth["n_vars"] <= 0:
            errors.append("维度非法：%s×%s" % (truth["n_obs"], truth["n_vars"]))
        obs_frame = result.get("obs")
        obs_rows = len(obs_frame) if isinstance(obs_frame, pd.DataFrame) else 0
        if obs_rows != truth["n_obs"]:
            errors.append("obs 行数=%d，应 %d" % (obs_rows, truth["n_obs"]))
        if len(result.get("var_names") or []) != truth["n_vars"]:
            errors.append("var_names 长度=%d，应 %d" % (len(result.get("var_names") or []), truth["n_vars"]))
        if stat_before.st_size != stat_after.st_size or stat_before.st_mtime_ns != stat_after.st_mtime_ns:
            errors.append("真实 rds 被修改（size/mtime 变化）")
        n_obs, n_vars = truth["n_obs"], truth["n_vars"]
        n_obs_cols = int(result["obs"].shape[1]) if isinstance(result.get("obs"), pd.DataFrame) else 0
        loaded_info.append((Path(path).name, n_obs, n_vars, n_obs_cols, result.get("meta", {}).get("assays")))
        report.check_many("真实 rds %-22s" % Path(path).name, errors,
                          "%d 细胞 × %d 基因，obs 列 %d | class/assays=%r" % (
                              n_obs, n_vars, n_obs_cols, result.get("meta", {}).get("assays")))
        if index == 0:
            try:
                info = io.inspect_source(path, workdir=str(work))
                errors = []
                if info.get("format") != "rds":
                    errors.append("inspect.format=%r，应 rds" % (info.get("format"),))
                if info.get("n_obs") not in (None, n_obs) or info.get("n_vars") not in (None, n_vars):
                    errors.append("inspect 维度=%rx%s，应 %dx%d" % (
                        info.get("n_obs"), info.get("n_vars"), n_obs, n_vars))
                if not isinstance(info.get("size"), int):
                    errors.append("inspect.size 不是 int：%r" % (info.get("size"),))
                report.check_many("真实 rds inspect %-18s" % Path(path).name, errors)
            except BaseException as error:
                code, message, _ = _engine_error(error)
                report.check("真实 rds inspect", False, "%s(%s): %s" % (
                    type(error).__name__, code, message))
    print("\n  真实 rds 汇总（name → n_obs × n_vars，obs 列数，assays）：")
    for name, n_obs, n_vars, n_obs_cols, assays in loaded_info:
        print("    %-24s %5d 细胞 × %5d 基因，obs 列 %2d，assays=%r" % (name, n_obs, n_vars, n_obs_cols, assays))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    started = time.time()
    print("dsh-singlecell 全格式导入层测试（CONTRACT §8.5）")
    print("  python        : %s" % sys.version.split()[0])
    print("  sc_io         : %s" % (LIB_DIR / "sc_io.py"))
    if not (LIB_DIR / "sc_io.py").exists():
        print("FAIL: 缺少 lib/sc_io.py")
        return 1

    sc_io = _load_module("dsh_sc_io_under_test", LIB_DIR / "sc_io.py")
    make_formats = _load_module("dsh_make_formats_under_test", TEST_DIR / "make_formats.py")

    temp_root = Path(tempfile.mkdtemp(prefix="dsh-sc-import-test-"))
    fixture_root = temp_root / "fixtures"
    work = temp_root / "work"
    work.mkdir(parents=True, exist_ok=True)
    report = Reporter()
    keep = os.environ.get("DS_SC_KEEP_FIXTURES") == "1"
    ok = False
    try:
        print("  临时目录      : %s" % temp_root)
        print("\n生成 fixture（CONTRACT §8.5）…")
        manifest = make_formats.generate(fixture_root, verbose=False)
        print("  fixture root  : %s" % manifest["root"])
        print("  formats=%d ambiguous=%d errors=%d url=%d transpose=%d real_rds=%d" % (
            len(manifest["formats"]), len(manifest["ambiguous"]), len(manifest["errors"]),
            len(manifest["url_tests"]), len(manifest["transpose_probes"]),
            len(manifest.get("real_rds", {}).get("files") or [])))

        test_contract_surface(sc_io, manifest, work, report)
        test_formats(sc_io, manifest, work, report)
        test_inspect(sc_io, manifest, work, report)
        test_ambiguous(sc_io, manifest, work, report)
        test_errors(sc_io, manifest, work, report)
        test_transpose(sc_io, manifest, work, report)
        test_urls(sc_io, manifest, work, report, make_formats)
        test_real_rds(sc_io, manifest, work, report)
        ok = report.summary()
    except BaseException:
        print("\n测试框架异常：")
        traceback.print_exc()
        ok = False
    finally:
        print("\n总用时 %.1fs" % (time.time() - started))
        if keep:
            print("DS_SC_KEEP_FIXTURES=1，保留临时目录：%s" % temp_root)
        else:
            shutil.rmtree(temp_root, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
