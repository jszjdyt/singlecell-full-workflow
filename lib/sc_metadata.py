#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell — 样本 / 分组 / 元数据（CONTRACT §11）

用户诉求（GEO 工作流）：
  每个 GSM 一个 CSV → 多文件合并成一个数据集 → 每个文件=一个样本 →
  逐个改样本名、给样本指定大组（young/old 等）→ 大组写回 obs 供组间分析使用。

本模块只处理 obs / 样本映射，不碰矩阵计算；由 sc_engine 的命令层调用。
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import scipy.sparse as sp
except Exception:  # pragma: no cover
    sp = None


class MetaError(Exception):
    def __init__(self, code, message):
        self.code = str(code)
        self.message = str(message)
        super().__init__(self.message)


# ---------------------------------------------------------------------------
# 多样本合并导入
# ---------------------------------------------------------------------------


def _sample_name(path):
    name = Path(str(path)).name
    for suffix in (".gz", ".bz2", ".zip"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in (".csv", ".tsv", ".txt", ".mtx", ".h5ad", ".h5", ".loom", ".rds", ".rdata", ".rda"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name or "sample"


def _as_csr(matrix):
    if matrix is None:
        return None
    if sp.issparse(matrix):
        return sp.csr_matrix(matrix, dtype=np.float32)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float32))


def list_sources(path):
    """目录 → 目录下可当样本处理的数据文件（排序稳定）；文件 → [文件]。"""
    target = Path(str(path))
    if not target.exists():
        raise MetaError("not_found", "路径不存在：%s" % target)
    if target.is_file():
        return [str(target)]
    files = []
    for item in sorted(target.iterdir()):
        if not item.is_file() or item.name.startswith("."):
            continue
        lowered = item.name.lower()
        if lowered.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz",
                             ".h5ad", ".h5", ".loom", ".mtx", ".mtx.gz", ".rds", ".rdata", ".rda")):
            files.append(str(item))
    if not files:
        raise MetaError("bad_request", "目录里没有可识别的数据文件：%s" % target)
    return files


def merge_sources(loader, sources, workdir, options=None, progress=None):
    """把多个来源合并成一个数据集：基因取交集，样本列 = 文件名。"""
    options = dict(options or {})
    loaded_list = []
    sample_meta = []
    warnings = []
    for index, item in enumerate(sources):
        if isinstance(item, dict):
            source = str(item.get("source") or item.get("path") or "")
            item_options = dict(options)
            item_options.update(item.get("options") or {})
            display = str(item.get("name") or _sample_name(source))
            group = item.get("group")
        else:
            source = str(item)
            item_options = options
            display = _sample_name(source)
            group = None
        if progress:
            progress("读取样本 %d/%d：%s" % (index + 1, len(sources), Path(source).name))
        loaded = loader.load_source(source, workdir=workdir, options=item_options)
        loaded_list.append((source, loaded))
        sample_meta.append({"source": source, "display": display, "group": group})
    if not loaded_list:
        raise MetaError("bad_request", "没有可合并的数据文件")

    first = loaded_list[0][1]
    genes = [str(item) for item in first["var_names"]]
    for source, loaded in loaded_list[1:]:
        other = set(str(item) for item in loaded["var_names"])
        genes = [gene for gene in genes if gene in other]
    if not genes:
        raise MetaError("bad_request", "各文件之间没有共同的基因，无法合并")
    if len(genes) < len(first["var_names"]):
        warnings.append("合并后基因取交集：%d / %d" % (len(genes), len(first["var_names"])))

    blocks = []
    counts_blocks = []
    obs_frames = []
    all_have_counts = True
    cells = []
    for index, (source, loaded) in enumerate(loaded_list):
        index_map = {str(gene): position for position, gene in enumerate(loaded["var_names"])}
        columns = [index_map[gene] for gene in genes]
        matrix = _as_csr(loaded["X"])[:, columns]
        blocks.append(matrix)
        counts = loaded.get("counts")
        if counts is None:
            all_have_counts = False
        else:
            counts = _as_csr(counts)
            if counts.shape[1] != len(loaded["var_names"]):
                all_have_counts = False
            else:
                counts_blocks.append(counts[:, columns])
        sample = sample_meta[index]["display"]
        frame = loaded.get("obs")
        frame = frame.copy() if isinstance(frame, pd.DataFrame) else pd.DataFrame(index=pd.Index(["cell_%d" % (i + 1) for i in range(matrix.shape[0])]))
        frame.insert(0, "sample", sample)
        if frame.index.has_duplicates:
            frame.index = pd.Index(["%s_%d" % (str(item), dup) for dup, item in zip(range(len(frame.index)), frame.index)])
        obs_frames.append(frame)
        cells.extend([str(item) for item in frame.index])
    if len(set(cells)) != len(cells):
        # 跨样本重名：统一加样本前缀
        for (source, _loaded), frame in zip(loaded_list, obs_frames):
            sample = _sample_name(source)
            frame.index = pd.Index(["%s|%s" % (sample, str(item)) for item in frame.index])
        warnings.append("跨样本存在重名细胞，已用 <样本名>|<原名> 去重")
    X = sp.vstack(blocks, format="csr") if sp is not None else np.vstack(blocks)
    counts = sp.vstack(counts_blocks, format="csr") if (all_have_counts and counts_blocks and sp is not None) else None
    obs = pd.concat(obs_frames, axis=0)
    obs["sample"] = obs["sample"].astype("object")
    samples = []
    group_by_sample = {entry["display"]: entry.get("group") for entry in sample_meta}
    for sample, group in obs.groupby("sample", sort=False):
        samples.append({"id": str(sample), "name": str(sample),
                        "group": group_by_sample.get(str(sample)), "n_cells": int(group.shape[0])})
    if any(entry.get("group") for entry in sample_meta):
        group_lookup = {str(entry["display"]): entry.get("group") for entry in sample_meta}
        obs["group"] = [group_lookup.get(str(item)) for item in obs["sample"].tolist()]
        meta_group_column = "group"
    else:
        meta_group_column = None
    meta = {
        "samples": samples,
        "sample_column": "sample",
        "group_column": meta_group_column,
        "sources": [sample_meta_entry["source"] for sample_meta_entry in sample_meta],
    }
    return {
        "X": X,
        "counts": counts,
        "obs": obs,
        "var_names": genes,
        "warnings": warnings,
        "meta": meta,
        "n_obs": int(X.shape[0]),
        "n_vars": int(X.shape[1]),
    }


def _persist_column(dataset, name, values):
    """把列同时写到当前 obs 和 full_obs（qc/filter 步骤会从 full_obs 重建 obs，元数据不能丢）。"""
    series = pd.Series(list(values), index=dataset.obs.index)
    dataset.obs[name] = series
    full = getattr(dataset, "full_obs", None)
    if isinstance(full, pd.DataFrame):
        if name in full.columns:
            full.drop(columns=[name], inplace=True)
        full[name] = series.reindex(full.index)
    return series


# ---------------------------------------------------------------------------
# 样本表 / 分组
# ---------------------------------------------------------------------------


def _sample_table(dataset):
    obs = dataset.obs
    window = dataset.uns.get("samples") if isinstance(getattr(dataset, "uns", None), dict) else None
    if isinstance(window, list) and window:
        table = []
        for entry in window:
            sample_id = str(entry.get("id") or entry.get("name"))
            n_cells = 0
            if "sample" in obs.columns:
                n_cells = int((obs["sample"].astype(str) == sample_id).sum())
            table.append({
                "id": sample_id,
                "name": str(entry.get("name") or sample_id),
                "group": entry.get("group"),
                "n_cells": n_cells,
            })
        return table
    if "sample" in obs.columns:
        table = []
        for sample, group in obs.groupby("sample", sort=False):
            table.append({"id": str(sample), "name": str(sample), "group": None, "n_cells": int(group.shape[0])})
        return table
    return []


def _groups_of(dataset):
    column = dataset.uns.get("group_column") if isinstance(getattr(dataset, "uns", None), dict) else None
    if column and column in dataset.obs.columns:
        values = [str(item) for item in dataset.obs[column].dropna().unique().tolist()]
        return sorted(values)
    return []


def obs_list(dataset):
    columns = []
    for name in dataset.obs.columns:
        series = dataset.obs[name]
        kind = "numeric" if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series) else "categorical"
        entry = {"key": str(name), "kind": kind, "n_unique": int(series.nunique(dropna=True))}
        if kind == "categorical" and entry["n_unique"] <= 30:
            entry["values"] = [str(item) for item in series.dropna().unique().tolist()[:30]]
        columns.append(entry)
    return {
        "columns": columns,
        "samples": _sample_table(dataset),
        "groups": _groups_of(dataset),
        "group_column": dataset.uns.get("group_column") if isinstance(getattr(dataset, "uns", None), dict) else None,
        "sample_column": "sample" if "sample" in dataset.obs.columns else None,
    }


def samples_update(dataset, samples, groups=None, sample_column="sample", group_column="group"):
    if "sample" not in dataset.obs.columns and sample_column not in dataset.obs.columns:
        raise MetaError("prerequisite", "该数据集没有样本列：请用多文件/目录合并导入（每个文件=一个样本），或先 obs_derive 建列")
    table = {entry["id"]: dict(entry) for entry in _sample_table(dataset)}
    if not table:
        raise MetaError("prerequisite", "没有可编辑的样本（缺少 sample 列）")
    for entry in samples or []:
        sample_id = str(entry.get("id") or entry.get("name") or "")
        if sample_id not in table:
            # 允许用当前显示名匹配
            matched = [key for key, value in table.items() if value.get("name") == sample_id]
            if not matched:
                continue
            sample_id = matched[0]
        if entry.get("name"):
            table[sample_id]["name"] = str(entry["name"])
        if entry.get("group") is not None:
            table[sample_id]["group"] = str(entry["group"]) if str(entry["group"]).strip() else None
    sample_lookup = {key: {"name": value.get("name") or key, "group": value.get("group")} for key, value in table.items()}
    source = dataset.obs[sample_column].astype(str) if sample_column in dataset.obs.columns else dataset.obs["sample"].astype(str)
    _persist_column(dataset, "sample_name", [sample_lookup.get(item, {}).get("name") or item for item in source.tolist()])
    active_group = group_column if group_column else "group"
    values = [sample_lookup.get(item, {}).get("group") for item in source.tolist()]
    if any(value is not None for value in values):
        _persist_column(dataset, active_group, [value if value is not None else np.nan for value in values])
        dataset.uns["group_column"] = active_group
    dataset.uns["samples"] = [{"id": key, "name": value["name"], "group": value["group"]} for key, value in table.items()]
    if groups:
        dataset.uns["groups_declared"] = [str(item) for item in groups]
    return {
        "samples": _sample_table(dataset),
        "groups": _groups_of(dataset),
        "group_column": dataset.uns.get("group_column"),
        "updated_cells": int(dataset.obs.shape[0]),
    }


def samples_auto_group(dataset, pattern, sample_column="sample"):
    table = _sample_table(dataset)
    if not table:
        raise MetaError("prerequisite", "没有可编辑的样本（缺少 sample 列）：请用多文件/目录合并导入")
    try:
        regex = re.compile(str(pattern), re.IGNORECASE)
    except re.error as error:
        raise MetaError("bad_request", "分组规则不是合法正则：%s" % error)
    mapping = {}
    unmatched = []
    for entry in table:
        texts = [str(entry["id"]), str(entry["name"] or entry["id"])]
        match = None
        for text in texts:
            match = regex.search(text)
            if match:
                break
        if match:
            mapping[entry["id"]] = match.group(1) if match.groups() else match.group(0)
        else:
            unmatched.append(entry["id"])
    return {"mapping": mapping, "unmatched": unmatched, "pattern": str(pattern)}


def samples_import_map(dataset, path, sample_column=None, group_column=None):
    target = Path(str(path))
    if not target.exists():
        raise MetaError("not_found", "对照文件不存在：%s" % target)
    try:
        frame = pd.read_csv(target)
    except Exception as error:
        raise MetaError("bad_request", "无法读取对照 CSV：%s" % error)
    frame.columns = [str(item).strip() for item in frame.columns]
    sample_key = sample_column or ("sample" if "sample" in frame.columns else frame.columns[0])
    group_key = group_column or ("group" if "group" in frame.columns else (frame.columns[1] if len(frame.columns) > 1 else None))
    if group_key is None:
        raise MetaError("bad_request", "对照 CSV 至少需要两列：样本名、大组")
    mapping = {}
    for _, row in frame.iterrows():
        mapping[str(row[sample_key]).strip()] = str(row[group_key]).strip()
    table = _sample_table(dataset)
    if not table:
        raise MetaError("prerequisite", "没有可编辑的样本（缺少 sample 列）")
    updates = []
    matched = 0
    for entry in table:
        key = str(entry["id"])
        if key in mapping:
            updates.append({"id": key, "group": mapping[key]})
            matched += 1
    return {"updates": updates, "mapping": mapping, "matched": matched, "total": len(table), "sample_key": sample_key, "group_key": group_key}


def obs_set_group(dataset, column):
    column = str(column)
    if column not in dataset.obs.columns:
        raise MetaError("bad_request", "obs 里没有列 %s：%s" % (column, "、".join(map(str, dataset.obs.columns))))
    unique = int(dataset.obs[column].nunique(dropna=True))
    if unique < 2:
        raise MetaError("bad_request", "列 %s 只有一个取值，不能作为分组列" % column)
    dataset.uns["group_column"] = column
    return {"group_column": column, "groups": _groups_of(dataset), "n_groups": unique}


# ---------------------------------------------------------------------------
# 外部 metadata join / 细胞名派生
# ---------------------------------------------------------------------------


def _normalize_key(value):
    text = str(value).strip()
    text = re.sub(r"-\d+$", "", text)
    return text.lower()


def obs_join(dataset, path, key=None, cell_column=None):
    target = Path(str(path))
    if not target.exists():
        raise MetaError("not_found", "metadata 文件不存在：%s" % target)
    try:
        frame = pd.read_csv(target)
    except Exception as error:
        raise MetaError("bad_request", "无法读取 metadata CSV：%s" % error)
    if frame.empty:
        raise MetaError("bad_request", "metadata CSV 是空的")
    frame.columns = [str(item).strip() for item in frame.columns]
    cell_key = key or cell_column
    if cell_key:
        if cell_key not in frame.columns:
            raise MetaError("bad_request", "metadata 里没有列 %s（现有：%s）" % (cell_key, "、".join(frame.columns[:10])))
    else:
        best = None
        for column in frame.columns:
            values = frame[column].astype(str)
            hits = sum(1 for item in values.head(200) if _normalize_key(item) in {_normalize_key(c) for c in dataset.obs.index[:5000]})
            if best is None or hits > best[1]:
                best = (column, hits)
        cell_key = best[0] if best else frame.columns[0]
    lookup = {}
    for _, row in frame.iterrows():
        lookup[_normalize_key(row[cell_key])] = row
    matched = 0
    new_columns = [column for column in frame.columns if column != cell_key]
    data = {column: [np.nan] * dataset.n_obs for column in new_columns}
    for position, cell in enumerate(dataset.obs.index):
        row = lookup.get(_normalize_key(cell))
        if row is None:
            continue
        matched += 1
        for column in new_columns:
            data[column][position] = row[column]
    for column in new_columns:
        column_name = str(column)
        if column_name in dataset.obs.columns and column_name != cell_key:
            column_name = "%s_meta" % column_name
        values = data[str(column)]
        if pd.api.types.is_numeric_dtype(frame[column]):
            values = pd.to_numeric(pd.Series(values), errors="coerce").tolist()
        _persist_column(dataset, column_name, values)
    warnings = []
    if matched == 0:
        warnings.append("没有匹配到任何细胞：检查 metadata 的细胞名列与数据细胞名（barcode）是否一致")
    return {
        "joined_columns": [str(column) for column in new_columns],
        "key_column": str(cell_key),
        "matched": matched,
        "total": int(dataset.n_obs),
        "match_rate": round(matched / max(1, dataset.n_obs), 4),
        "warnings": warnings,
    }


def obs_derive(dataset, name, source="cell_name", pattern=None, group=1, value=None):
    column = str(name or "").strip()
    if not column:
        raise MetaError("bad_request", "需要提供新列名 name")
    if source == "cell_name":
        if not pattern:
            raise MetaError("bad_request", "source=cell_name 需要提供正则 pattern（第 1 个捕获组作为值）")
        try:
            regex = re.compile(str(pattern))
        except re.error as error:
            raise MetaError("bad_request", "正则不合法：%s" % error)
        values = []
        unmatched = 0
        for cell in dataset.obs.index:
            match = regex.search(str(cell))
            if not match:
                values.append(None)
                unmatched += 1
                continue
            try:
                values.append(match.group(int(group)))
            except Exception:
                values.append(match.group(0))
        _persist_column(dataset, column, values)
        if unmatched:
            return {"column": column, "unmatched": unmatched, "n_unique": int(dataset.obs[column].nunique(dropna=True)),
                    "warnings": ["有 %d 个细胞名未匹配正则（值为空）" % unmatched]}
        return {"column": column, "unmatched": 0, "n_unique": int(dataset.obs[column].nunique(dropna=True))}
    if source == "constant":
        text = str(value if value is not None else (pattern or ""))
        _persist_column(dataset, column, [text] * dataset.n_obs)
        return {"column": column, "n_unique": 1}
    raise MetaError("bad_request", "未知 source：%s（可用 cell_name / constant）" % source)
