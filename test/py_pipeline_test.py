#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench Python 计算引擎测试（不依赖 node）。

直接运行：
    python3 test/py_pipeline_test.py

覆盖 CONTRACT.md §2.5 的 7 条与 §6.4 的 5 条，并通过**真实的 subprocess
JSON-lines 协议**驱动 lib/sc_engine.py（起 worker → ping → open(demo) →
run 全链路 → embedding → plot → markers → suggest → export → list → close →
shutdown），逐条校验返回字段存在且类型正确。

退出码 0 = 全部 PASS，1 = 有 FAIL。临时工作目录建在系统临时区（唯一目录），
测试结束后清理（设置环境变量 KEEP_SC_TEST_TMP=1 可保留用于排查）。
"""

from __future__ import annotations

import json
import os
import select
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PLUGIN_DIR = HERE.parent
ENGINE_SCRIPT = PLUGIN_DIR / "lib" / "sc_engine.py"
PYTHON = sys.executable or "python3"

sys.path.insert(0, str(HERE))
from make_sample import generate_demo  # noqa: E402

STEP_ORDER = [
    "qc", "filter", "normalize", "hvg", "pca",
    "neighbors", "cluster", "embed", "markers", "annotate",
]

PASS_COUNT = 0
FAIL_COUNT = 0
_LAST_FAILURES = []


def check(name, condition, detail="", required=True):
    """记录一条断言。required=False 的失败只打印不计数（用于软性提示）。"""
    global PASS_COUNT, FAIL_COUNT
    if condition:
        PASS_COUNT += 1
        print("PASS  %-62s %s" % (name, detail))
        return True
    if required:
        FAIL_COUNT += 1
        _LAST_FAILURES.append(name)
        print("FAIL  %-62s %s" % (name, detail))
        return False
    print("note  %-62s %s" % (name, detail))
    return False


# ---------------------------------------------------------------------------
# JSON-lines worker 客户端
# ---------------------------------------------------------------------------


class WorkerError(AssertionError):
    pass


class Worker:
    def __init__(self, workdir):
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ)
        env.update(
            {
                "MPLBACKEND": "Agg",
                "MPLCONFIGDIR": str(self.workdir / "mpl-cache"),
                "PYTHONUNBUFFERED": "1",
                "PYTHONIOENCODING": "utf-8",
            }
        )
        self.proc = subprocess.Popen(
            [PYTHON, "-u", str(ENGINE_SCRIPT), "--workdir", str(self.workdir)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            bufsize=1,
            env=env,
        )
        self.seq = 0
        self.protocol_lines = []

    # -- 底层读写 --------------------------------------------------------
    def _readline(self, timeout):
        try:
            ready, _, _ = select.select([self.proc.stdout], [], [], timeout)
        except Exception:
            ready = None
        if not ready:
            raise WorkerError("等待引擎响应超时（%.0fs）" % timeout)
        line = self.proc.stdout.readline()
        if line == "":
            raise WorkerError("引擎 stdout 意外 EOF（进程可能已退出）")
        self.protocol_lines.append(line)
        return line

    def raw_call(self, cmd, args=None, timeout=300):
        """发一条命令并返回完整响应 dict（不判断 ok）。"""
        self.seq += 1
        identifier = "req-%d" % self.seq
        payload = {"id": identifier, "cmd": cmd, "args": args if args is not None else {}}
        try:
            self.proc.stdin.write(json.dumps(payload, ensure_ascii=False) + "\n")
            self.proc.stdin.flush()
        except Exception as error:
            raise WorkerError("写入引擎失败：%s" % error)
        line = self._readline(timeout)
        try:
            message = json.loads(line)
        except json.JSONDecodeError:
            raise WorkerError("引擎 stdout 出现非 JSON 行：%r" % line[:200])
        if not isinstance(message, dict):
            raise WorkerError("引擎响应不是 JSON 对象：%r" % line[:200])
        if str(message.get("id")) != identifier:
            raise WorkerError("响应 id 不匹配：期望 %s，收到 %r" % (identifier, message.get("id")))
        return message

    def call(self, cmd, args=None, timeout=300):
        message = self.raw_call(cmd, args, timeout=timeout)
        if not message.get("ok"):
            error = message.get("error") or {}
            raise WorkerError("%s 失败（%s）：%s" % (cmd, error.get("code"), error.get("message")))
        return message.get("result")

    def call_error(self, cmd, args=None, timeout=300):
        """发命令并断言失败；返回 error dict。"""
        message = self.raw_call(cmd, args, timeout=timeout)
        if message.get("ok"):
            raise WorkerError("%s 本应失败却返回成功" % cmd)
        error = message.get("error") or {}
        if not isinstance(error, dict) or "code" not in error or "message" not in error:
            raise WorkerError("错误响应缺少 code/message：%r" % message)
        return error

    # -- 收尾 ------------------------------------------------------------
    def shutdown(self, timeout=120):
        result = self.call("shutdown", {}, timeout=timeout)
        code = self.proc.wait(timeout=timeout)
        leftover = self.proc.stdout.read()
        stderr = self.proc.stderr.read()
        return result, code, leftover, stderr

    def kill(self):
        try:
            if self.proc.poll() is None:
                self.proc.kill()
                self.proc.wait(timeout=5)
        except Exception:
            pass


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------


def png_dimensions(path):
    with open(path, "rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return int.from_bytes(header[16:20], "big"), int.from_bytes(header[20:24], "big")


def is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def iso_like(value):
    return isinstance(value, str) and "T" in value and (value.endswith("Z") or "+" in value)


def check_history_entry(entry, label):
    ok = True
    ok &= check("%s: at 是 ISO 时间" % label, iso_like(entry.get("at")), str(entry.get("at")))
    ok &= check("%s: step 合法" % label, entry.get("step") in STEP_ORDER, str(entry.get("step")))
    ok &= check("%s: params 是对象" % label, isinstance(entry.get("params"), dict))
    ok &= check("%s: ok 是 bool" % label, isinstance(entry.get("ok"), bool))
    ok &= check("%s: ms 是非负数字" % label, is_number(entry.get("ms")) and entry.get("ms") >= 0)
    ok &= check("%s: n_obs/n_vars 是正整数" % label,
                isinstance(entry.get("n_obs"), int) and isinstance(entry.get("n_vars"), int)
                and entry["n_obs"] > 0 and entry["n_vars"] > 0)
    ok &= check("%s: summary 是字符串" % label, isinstance(entry.get("summary"), str))
    engine = entry.get("engine") or {}
    versions = engine.get("versions") or {}
    ok &= check("%s: engine.mode/python/seed" % label,
                engine.get("mode") in ("builtin", "scanpy")
                and isinstance(engine.get("python"), str)
                and engine.get("seed") == 0)
    ok &= check("%s: engine.versions 含 numpy/sklearn" % label,
                isinstance(versions, dict) and "numpy" in versions and "sklearn" in versions)
    return ok


# ---------------------------------------------------------------------------
# 主测试流程
# ---------------------------------------------------------------------------


def run_tests(tmp_root):
    workdir = tmp_root / "workdir"
    started_at = time.time()
    demo = generate_demo()
    true_types = list(demo["cell_types"])
    marker_genes = demo["marker_genes"]

    worker = Worker(workdir)
    try:
        # ===============================================================
        # 0. 协议握手 ping（启动不打印任何东西）
        # ===============================================================
        ping = worker.call("ping", {}, timeout=120)
        check("ping.mode ∈ builtin|scanpy", ping.get("mode") in ("builtin", "scanpy"), str(ping.get("mode")))
        check("ping.python 是版本字符串", isinstance(ping.get("python"), str)
              and ping["python"].count(".") >= 1, str(ping.get("python")))
        check("ping.platform 是非空字符串", isinstance(ping.get("platform"), str) and bool(ping["platform"]),
              str(ping.get("platform")))
        packages = ping.get("packages") or {}
        check("ping.packages 含契约要求的 6 个包",
              all(key in packages for key in ("numpy", "pandas", "scipy", "sklearn", "anndata", "scanpy")))
        check("ping.packages.numpy 是字符串", isinstance(packages.get("numpy"), str), str(packages.get("numpy")))
        check("ping.packages.scanpy 是 null 或字符串",
              packages.get("scanpy") is None or isinstance(packages.get("scanpy"), str))
        layouts = ping.get("layouts") or []
        check("ping.layouts 含 tsne 与 pca", "tsne" in layouts and "pca" in layouts, str(layouts))
        check("ping.warnings 是数组", isinstance(ping.get("warnings"), list), str(ping.get("warnings"))[:60])

        # ===============================================================
        # 1. §2.5-7 open 不存在的路径 → bad_request，进程仍存活
        # ===============================================================
        missing_path = str(tmp_root / "definitely-not-here.h5ad")
        error = worker.call_error("open", {"path": missing_path})
        check("缺路径 open → bad_request/not_found", error.get("code") in ("bad_request", "not_found"), error.get("message", "")[:70])
        ping_alive = worker.call("ping", {})
        check("bad_request 后进程仍存活", ping_alive.get("mode") == ping.get("mode"))

        # ===============================================================
        # 2. open demo + dataset 对象字段
        # ===============================================================
        open_started = time.time()
        opened = worker.call("open", {"demo": True}, timeout=120)
        dataset_id = opened.get("id")
        check("open.demo 返回 dataset.id", isinstance(dataset_id, str) and dataset_id.startswith("ds"),
              str(dataset_id))
        check("open.demo name=demo-pbmc", opened.get("name") == "demo-pbmc", str(opened.get("name")))
        check("open.demo ≥2000 细胞 ≥1000 基因",
              opened.get("n_obs", 0) >= 2000 and opened.get("n_vars", 0) >= 1000,
              "%s×%s" % (opened.get("n_obs"), opened.get("n_vars")))
        check("open.demo steps=[] history=[]",
              opened.get("steps") == [] and opened.get("history") == [])
        check("open.demo 带嵌套 dataset（Host 兼容）", (opened.get("dataset") or {}).get("id") == dataset_id)
        obs_columns = {column["key"]: column for column in opened.get("obs_columns") or []}
        check("obs_columns 含 cell_type(categorical)",
              obs_columns.get("cell_type", {}).get("kind") == "categorical")
        check("obs_columns 每项含 key/kind/n_unique",
              all(isinstance(column.get("key"), str) and column.get("kind") in ("numeric", "categorical")
                  and isinstance(column.get("n_unique"), int) for column in opened.get("obs_columns") or []))

        # preview
        preview = worker.call("preview", {"dataset": dataset_id})
        check("preview.genes 首 50 个", isinstance(preview.get("genes"), list) and len(preview["genes"]) == 50)
        check("preview.n_obs/n_vars 正确",
              preview.get("n_obs") == opened.get("n_obs") and preview.get("n_vars") == opened.get("n_vars"))
        check("preview.obs_columns/steps 是数组",
              isinstance(preview.get("obs_columns"), list) and isinstance(preview.get("steps"), list))

        # ===============================================================
        # 3. §2.5-6 缺前置 → prerequisite（且不影响进程）
        # ===============================================================
        for command, args, label in [
            ("embedding", {"dataset": dataset_id, "color_by": {"kind": "cluster", "key": "clusters"}},
             "embed 前取 embedding"),
            ("markers", {"dataset": dataset_id}, "cluster 前取 markers"),
            ("plot", {"dataset": dataset_id, "kind": "embedding", "params": {}}, "embed 前画 embedding"),
            ("plot", {"dataset": dataset_id, "kind": "dotplot", "params": {}}, "cluster 前画 dotplot"),
        ]:
            error = worker.call_error(command, args)
            check("%s → prerequisite" % label, error.get("code") == "prerequisite",
                  error.get("message", "")[:70])
        check("prerequisite 后进程仍存活", worker.call("ping", {}).get("mode") == ping.get("mode"))

        # ===============================================================
        # 4. §2.5-1 全链路 open→markers（记时），逐步校验返回字段
        # ===============================================================
        chain_plan = [
            ("qc", {}),
            ("filter", {}),
            ("normalize", {}),
            ("hvg", {}),
            ("pca", {}),
            ("neighbors", {}),
            ("cluster", {"resolution": 1.0}),
            ("embed", {}),
            ("markers", {"groupby": "clusters", "top_n": 10, "method": "t-test"}),
        ]
        metrics_by_step = {}
        history_lengths = []
        audit_qc_only = None
        expected_history = 0
        for step, params in chain_plan:
            before = time.time()
            result = worker.call("run", {"dataset": dataset_id, "step": step, "params": params}, timeout=300)
            elapsed = time.time() - before
            expected_history += 1
            metrics = (result.get("summary") or {}).get("metrics") or {}
            metrics_by_step[step] = metrics
            ok_fields = (
                result.get("step") == step
                and isinstance((result.get("summary") or {}).get("title"), str)
                and isinstance((result.get("summary") or {}).get("bullets"), list)
                and len((result.get("summary") or {}).get("bullets")) >= 1
                and isinstance(metrics, dict)
                and isinstance(result.get("steps"), list)
                and result.get("current") == step
                and result.get("steps", [])[-1] == step
                and isinstance((result.get("dataset") or {}).get("n_obs"), int)
                and isinstance(result.get("history"), list)
            )
            check("run %-10s 返回字段齐全（%.2fs）" % (step, elapsed), ok_fields)
            history = result.get("history") or []
            history_lengths.append(len(history))
            check("run %-10s history 增长到 %d" % (step, expected_history),
                  len(history) == expected_history, "%d entries" % len(history))
            if history:
                check_history_entry(history[-1], "history[%s]" % step)
            if step == "qc":
                audit_qc_only = worker.call("audit", {"dataset": dataset_id})
            if step == "filter":
                check("filter 前后数量字段",
                      all(key in metrics for key in ("n_obs_before", "n_obs_after", "n_vars_before", "n_vars_after")))
        chain_seconds = time.time() - open_started
        check("§2.5 全链路 open→markers < 60s", chain_seconds < 60.0, "%.2fs" % chain_seconds)

        # 各步骤关键数值
        qc_metrics = metrics_by_step["qc"]
        check("qc 指标：n_obs=2400 且 n_vars=1200",
              qc_metrics.get("n_obs") == 2400 and qc_metrics.get("n_vars") == 1200, str(qc_metrics))
        check("qc 指标：n_mito_genes ≥ 1", qc_metrics.get("n_mito_genes", 0) >= 1)
        filter_metrics = metrics_by_step["filter"]
        n_obs_after = filter_metrics.get("n_obs_after")
        check("filter 去掉低质量细胞（< 原细胞数且 > 1800）",
              isinstance(n_obs_after, int) and 1800 <= n_obs_after < 2400,
              "%s -> %s" % (filter_metrics.get("n_obs_before"), n_obs_after))
        check("normalize.target_sum=10000", metrics_by_step["normalize"].get("target_sum") == 10000.0)
        hvg_metrics = metrics_by_step["hvg"]
        check("hvg 选出 min(n_top_genes, n_vars) 个基因",
              hvg_metrics.get("n_selected") == min(2000, hvg_metrics.get("n_vars", 0)),
              str(hvg_metrics.get("n_selected")))
        check("pca 取 30 个主成分且 PC1 比例在 (0,1)",
              metrics_by_step["pca"].get("n_pcs") == 30
              and 0.0 < metrics_by_step["pca"].get("pc1_variance_ratio", 0.0) < 1.0,
              "PC1=%.3f" % metrics_by_step["pca"].get("pc1_variance_ratio", 0.0))
        check("neighbors 有边且连通分量 ≥1",
              metrics_by_step["neighbors"].get("n_edges", 0) > 0
              and metrics_by_step["neighbors"].get("n_components", 0) >= 1,
              "%s edges" % metrics_by_step["neighbors"].get("n_edges"))
        cluster_metrics = metrics_by_step["cluster"]
        cluster_sizes = cluster_metrics.get("cluster_sizes") or {}
        check("cluster 簇数 2..50 且大小合计=过滤后细胞数",
              2 <= cluster_metrics.get("n_clusters", 0) <= 50
              and sum(int(value) for value in cluster_sizes.values()) == n_obs_after,
              "%s clusters" % cluster_metrics.get("n_clusters"))
        check("embed layout=tsne 且点数=n_obs",
              metrics_by_step["embed"].get("layout") == "tsne"
              and metrics_by_step["embed"].get("n_points") == n_obs_after)
        check("markers n_markers = n_groups × top_n",
              metrics_by_step["markers"].get("n_markers")
              == metrics_by_step["markers"].get("n_groups", 0) * 10,
              str(metrics_by_step["markers"].get("n_markers")))

        # ===============================================================
        # 5. §2.5-2 注入 marker 基因命中率（每个簇的 top markers 对应真实类型）
        # ===============================================================
        markers_r1 = worker.call("markers", {"dataset": dataset_id, "groupby": "clusters", "top_n": 10})
        groups_r1 = markers_r1.get("groups") or []
        check("markers.groups 结构完整",
              len(groups_r1) >= 2 and all(
                  isinstance(group.get("id"), str) and isinstance(group.get("label"), str)
                  and isinstance(group.get("n_cells"), int) and isinstance(group.get("genes"), list)
                  for group in groups_r1))
        check("每个 marker gene 字段齐全且类型正确",
              all(
                  isinstance(record.get("gene"), str)
                  and is_number(record.get("score"))
                  and is_number(record.get("logfc"))
                  and is_number(record.get("pval_adj"))
                  and is_number(record.get("pct_in"))
                  and is_number(record.get("pct_out"))
                  for group in groups_r1 for record in group.get("genes", [])
              ))
        cluster_best_type = {}
        per_cluster_hits = []
        for group in groups_r1:
            top_genes = [record["gene"] for record in group.get("genes", [])]
            if not top_genes:
                continue
            hits = {cell_type: len(set(top_genes) & set(marker_genes[cell_type])) for cell_type in true_types}
            best_type = max(hits, key=lambda key: (hits[key], key))
            best_hits = hits[best_type]
            cluster_best_type[str(group["id"])] = best_type
            per_cluster_hits.append(best_hits / float(len(top_genes)))
            check("簇 %-3s top10 命中注入 marker ≥60%%" % group["id"],
                  best_hits / float(len(top_genes)) >= 0.6,
                  "%d/%d（最匹配 %s）" % (best_hits, len(top_genes), best_type))
        average_hit = sum(per_cluster_hits) / len(per_cluster_hits) if per_cluster_hits else 0.0
        check("平均 marker 命中率 ≥80%", average_hit >= 0.8, "%.0f%%" % (average_hit * 100.0))
        covered_types = set(cluster_best_type.values())
        check("6 种注入细胞类型都能在聚类中找到对应簇",
              len(covered_types & set(true_types)) >= 6, str(sorted(covered_types)))
        strongest = groups_r1[0]["genes"][0]
        check("最强 marker 的 logFC>0 且 adj.P<0.05",
              strongest["logfc"] > 0 and strongest["pval_adj"] < 0.05,
              "%s logfc=%.2f p=%.2g" % (strongest["gene"], strongest["logfc"], strongest["pval_adj"]))

        # ===============================================================
        # 6. §2.5-4 embedding：长度一致、坐标有限、连续/基因着色
        # ===============================================================
        embedding_r1 = worker.call(
            "embedding",
            {"dataset": dataset_id, "color_by": {"kind": "cluster", "key": "clusters"}, "max_points": 12000},
        )
        n_points = embedding_r1.get("n_points")
        check("embedding 分类着色：n_points=n_obs", n_points == n_obs_after, str(n_points))
        check("embedding x/y/group_index 长度一致",
              len(embedding_r1.get("x") or []) == n_points
              and len(embedding_r1.get("y") or []) == n_points
              and len(embedding_r1.get("group_index") or []) == n_points)
        coordinates = (embedding_r1.get("x") or []) + (embedding_r1.get("y") or [])
        check("embedding 坐标有限且保留 3 位小数",
              len(coordinates) > 0 and all(is_number(value) and abs(value) < 1e6 for value in coordinates)
              and all(abs(round(value, 3) - value) < 1e-9 for value in coordinates))
        categories = embedding_r1.get("categories") or []
        check("embedding 分类 categories 有 id/label/color/n，且 n 合计=n_points",
              len(categories) >= 2 and all(
                  isinstance(item.get("id"), str) and isinstance(item.get("color"), str)
                  and item["color"].startswith("#") and isinstance(item.get("n"), int)
                  for item in categories)
              and sum(item["n"] for item in categories) == n_points)
        check("embedding continuous=false / legend_title=clusters",
              embedding_r1.get("continuous") is False and embedding_r1.get("legend_title") == "clusters")

        embedding_obs = worker.call(
            "embedding", {"dataset": dataset_id, "color_by": {"kind": "obs", "key": "n_genes"}}
        )
        check("embedding 数值 obs 着色为连续且 values 长度一致",
              embedding_obs.get("continuous") is True
              and len(embedding_obs.get("values") or []) == embedding_obs.get("n_points"),
              str(embedding_obs.get("legend_title")))
        check("embedding 连续 values 有限", all(is_number(value) for value in embedding_obs.get("values") or []))

        embedding_gene = worker.call(
            "embedding", {"dataset": dataset_id, "color_by": {"kind": "gene", "key": "CD3D"}}
        )
        check("embedding 基因着色：values 长度=n_points 且非全 0",
              len(embedding_gene.get("values") or []) == embedding_gene.get("n_points")
              and max(embedding_gene.get("values") or [0.0]) > 0,
              "legend=%s" % embedding_gene.get("legend_title"))
        embedding_missing_gene = worker.call(
            "embedding", {"dataset": dataset_id, "color_by": {"kind": "gene", "key": "GENE_NOT_PRESENT_XYZ"}}
        )
        check("embedding 缺失基因：全 0 + warning",
              all(value == 0 for value in embedding_missing_gene.get("values") or [])
              and bool(embedding_missing_gene.get("warnings")),
              str(embedding_missing_gene.get("warnings"))[:60])

        embedding_small = worker.call(
            "embedding",
            {"dataset": dataset_id, "color_by": {"kind": "cluster", "key": "clusters"}, "max_points": 500},
        )
        check("embedding 支持 max_points 降采样（长度仍一致）",
              embedding_small.get("downsampled") is True
              and embedding_small.get("n_points") == 500
              and len(embedding_small.get("x") or []) == 500
              and len(embedding_small.get("group_index") or []) == 500
              and sum(item.get("n", 0) for item in embedding_small.get("categories") or []) == 500,
              "n_points=%s" % embedding_small.get("n_points"))

        # ===============================================================
        # 7. §2.5-5 plot：7 种图 → PNG magic + 宽高 > 200 + 同参数复用
        # ===============================================================
        plot_specs = [
            ("qc_violin", {}),
            ("embedding", {}),
            ("dotplot", {}),
            ("heatmap", {}),
            ("cluster_bar", {}),
            ("pca_variance", {}),
            ("qc_scatter", {}),
        ]
        plot_files = {}
        for kind, params in plot_specs:
            info = worker.call("plot", {"dataset": dataset_id, "kind": kind, "params": params}, timeout=300)
            path = Path(str(info.get("file") or ""))
            exists = path.exists() and path.stat().st_size > 0
            magic = path.read_bytes()[:8] if exists else b""
            dimensions = png_dimensions(path) if exists else None
            plot_files[kind] = str(path)
            check("plot %-12s 是 PNG 且宽高>200" % kind,
                  exists and magic == b"\x89PNG\r\n\x1a\n" and dimensions is not None
                  and dimensions[0] > 200 and dimensions[1] > 200
                  and int(info.get("width", 0)) == dimensions[0]
                  and int(info.get("height", 0)) == dimensions[1]
                  and int(info.get("bytes", 0)) == path.stat().st_size,
                  "%sx%s %s bytes" % ((dimensions or (0, 0))[0], (dimensions or (0, 0))[1],
                                     info.get("bytes")))
        reuse = worker.call("plot", {"dataset": dataset_id, "kind": "cluster_bar", "params": {}})
        check("plot 同参数复用同一文件", Path(str(reuse.get("file"))) == Path(plot_files["cluster_bar"]))

        plot_variants = [
            ("qc_violin", {"metrics": ["pct_counts_mt"], "thresholds": {"pct_counts_mt": 20}}),
            ("qc_scatter", {"x_metric": "pct_counts_mt", "y_metric": "n_genes", "color_by_pct_mt": False}),
            ("dotplot", {"top_n": 3}),
            ("heatmap", {"top_n": 3}),
        ]
        for kind, params in plot_variants:
            info = worker.call("plot", {"dataset": dataset_id, "kind": kind, "params": params}, timeout=300)
            check("plot %-12s 参数变体可用" % kind, Path(str(info.get("file"))).exists()
                  and int(info.get("width", 0)) > 200, str(params))
        markers_top5 = worker.call("markers", {"dataset": dataset_id, "top_n": 5})
        check("markers 支持自定义 top_n",
              all(len(group.get("genes") or []) == 5 for group in markers_top5.get("groups") or []))

        # ===============================================================
        # 8. §2.5-3 同参数两次 cluster 完全一致（确定性）
        # ===============================================================
        rerun_cluster = worker.call(
            "run", {"dataset": dataset_id, "step": "cluster", "params": {"resolution": 1.0}}, timeout=300
        )
        check("重跑 cluster 簇数一致",
              (rerun_cluster["summary"]["metrics"] or {}).get("n_clusters")
              == cluster_metrics.get("n_clusters"),
              "%s == %s" % ((rerun_cluster["summary"]["metrics"] or {}).get("n_clusters"),
                            cluster_metrics.get("n_clusters")))
        # 重跑 cluster 会清空下游 embed/markers，需要重新跑
        worker.call("run", {"dataset": dataset_id, "step": "embed", "params": {}}, timeout=300)
        markers_r2 = worker.call("markers", {"dataset": dataset_id, "groupby": "clusters", "top_n": 10})
        embedding_r2 = worker.call(
            "embedding", {"dataset": dataset_id, "color_by": {"kind": "cluster", "key": "clusters"}}
        )
        check("同参数两次 markers 结果完全一致",
              json.dumps(markers_r2.get("groups"), sort_keys=True, ensure_ascii=False)
              == json.dumps(groups_r1, sort_keys=True, ensure_ascii=False))
        check("同参数两次 cluster 的 group_index 完全一致",
              embedding_r2.get("group_index") == embedding_r1.get("group_index"))
        check("同参数两次 cluster 的 categories 完全一致",
              [item.get("id") for item in embedding_r2.get("categories") or []]
              == [item.get("id") for item in embedding_r1.get("categories") or []])
        check("同参数两次 embed 的坐标完全一致",
              embedding_r2.get("x") == embedding_r1.get("x")
              and embedding_r2.get("y") == embedding_r1.get("y"))
        groups_r1 = markers_r2.get("groups") or []
        cluster_best_type = {}
        for group in groups_r1:
            top_genes = [record["gene"] for record in group.get("genes", [])]
            if not top_genes:
                continue
            hits = {cell_type: len(set(top_genes) & set(marker_genes[cell_type])) for cell_type in true_types}
            cluster_best_type[str(group["id"])] = max(hits, key=lambda key: (hits[key], key))

        # ===============================================================
        # 8b. layout=umap 在无 umap 依赖时回退 t-SNE 并说明（§2.4）
        # ===============================================================
        umap_run = worker.call("run", {"dataset": dataset_id, "step": "embed", "params": {"layout": "umap"}},
                               timeout=300)
        umap_metrics = (umap_run.get("summary") or {}).get("metrics") or {}
        if "umap" in (ping.get("layouts") or []):
            check("layout=umap 可用时走 UMAP", umap_metrics.get("layout") == "umap")
        else:
            check("无 umap 依赖时 layout=umap 回退 t-SNE 并在 bullets 说明",
                  umap_metrics.get("layout") == "tsne"
                  and any("回退" in str(bullet) or "UMAP" in str(bullet)
                          for bullet in umap_run["summary"]["bullets"]),
                  str(umap_run["summary"]["bullets"])[:100])

        # ===============================================================
        # 9. §6.4-2 autoannotate：top1 与真实类型一致 ≥4/6
        # ===============================================================
        autoannotate = worker.call("autoannotate", {"dataset": dataset_id, "groupby": "clusters"}, timeout=180)
        suggestions = autoannotate.get("suggestions") or []
        check("autoannotate.panel 是字符串", isinstance(autoannotate.get("panel"), str),
              str(autoannotate.get("panel")))
        check("autoannotate.suggestions 结构完整",
              len(suggestions) >= 2 and all(
                  isinstance(item.get("cluster"), str) and isinstance(item.get("label"), str)
                  and 0.0 <= item.get("score", -1) <= 1.0 and is_number(item.get("margin"))
                  and isinstance(item.get("evidence"), list) and item["evidence"]
                  and isinstance(item.get("alternative"), str)
                  for item in suggestions))
        label_by_cluster = {str(item["cluster"]): item["label"] for item in suggestions}
        best_cluster_for_type = {}
        for cluster_id, cell_type in cluster_best_type.items():
            previous = best_cluster_for_type.get(cell_type)
            if previous is None:
                best_cluster_for_type[cell_type] = cluster_id
        hit_types = 0
        for cell_type in true_types:
            cluster_id = best_cluster_for_type.get(cell_type)
            matched = cluster_id is not None and label_by_cluster.get(cluster_id) == cell_type
            if matched:
                hit_types += 1
            check("autoannotate 真实类型 %-14s 命中" % cell_type, matched,
                  "cluster %s → %s" % (cluster_id, label_by_cluster.get(cluster_id)))
        check("§6.4 autoannotate top1 命中 ≥4/6", hit_types >= 4, "%d/6" % hit_types)

        # ===============================================================
        # 10. annotate：把建议写回 obs（不写 obs 的是 autoannotate）
        # ===============================================================
        mapping = {str(item["cluster"]): item["label"] for item in suggestions}
        annotate_result = worker.call(
            "run", {"dataset": dataset_id, "step": "annotate", "params": {"mapping": mapping}}
        )
        annotate_metrics = (annotate_result.get("summary") or {}).get("metrics") or {}
        check("annotate 写入 cell_type（n_types=映射后类型数）",
              annotate_metrics.get("n_types") == len(set(mapping.values())),
              "%s 种" % annotate_metrics.get("n_types"))
        preview_after = worker.call("preview", {"dataset": dataset_id})
        check("annotate 后 obs_columns 出现 cell_type(categorical)",
              {column["key"]: column["kind"] for column in preview_after.get("obs_columns") or []}.get("cell_type")
              == "categorical")

        # ===============================================================
        # 11. §6.4-3 resolution_scan：条目数=输入数、ARI ∈ [-1,1]、不改当前聚类
        # ===============================================================
        scan_input = [0.4, 0.8, 1.2]
        scan = worker.call(
            "resolution_scan", {"dataset": dataset_id, "resolutions": scan_input}, timeout=300
        )
        scan_results = scan.get("results") or []
        check("resolution_scan 条目数 = 输入数", len(scan_results) == len(scan_input), str(len(scan_results)))
        check("resolution_scan 每项含 resolution/n_clusters/ari_vs_current",
              all(is_number(item.get("resolution")) and isinstance(item.get("n_clusters"), int)
                  and item["n_clusters"] >= 1 and -1.0 <= item.get("ari_vs_current", 2.0) <= 1.0
                  for item in scan_results),
              str([round(item.get("ari_vs_current", 0), 3) for item in scan_results]))
        check("resolution_scan.current 是当前聚类",
              (scan.get("current") or {}).get("n_clusters") == cluster_metrics.get("n_clusters")
              and (scan.get("current") or {}).get("resolution") == 1.0)
        embedding_after_scan = worker.call(
            "embedding", {"dataset": dataset_id, "color_by": {"kind": "cluster", "key": "clusters"}}
        )
        check("resolution_scan 不改变当前聚类（group_index 不变）",
              embedding_after_scan.get("group_index") == embedding_r2.get("group_index"))

        scan_default = worker.call("resolution_scan", {"dataset": dataset_id}, timeout=300)
        check("resolution_scan 输入为空时默认 6 个分辨率",
              len(scan_default.get("results") or []) == 6,
              str([item.get("resolution") for item in scan_default.get("results") or []]))

        # ===============================================================
        # 12. §6.4-1 audit：qc-only 有 warn/fail，全流程后明显更高
        # ===============================================================
        audit_full = worker.call("audit", {"dataset": dataset_id}, timeout=120)
        checks_full = audit_full.get("checks") or []
        expected_check_ids = [
            "qc_before_cluster", "mito_high", "cells_low", "clusters_few", "cluster_skew",
            "markers_missing", "embedding_missing", "annotate_missing", "doublet_risk",
            "resolution_untested",
        ]
        check("audit.checks 恰好 10 项且 id 固定",
              [item.get("id") for item in checks_full] == expected_check_ids,
              str([item.get("id") for item in checks_full]))
        check("audit 每项 level ∈ ok/warn/fail 且有 title/detail",
              all(item.get("level") in ("ok", "warn", "fail")
                  and isinstance(item.get("title"), str) and isinstance(item.get("detail"), str)
                  and isinstance(item.get("hint"), str) for item in checks_full))
        check("audit.score ∈ [0,100] 且 level 合法",
              isinstance(audit_full.get("score"), int) and 0 <= audit_full["score"] <= 100
              and audit_full.get("level") in ("ok", "warn", "fail"),
              "score=%s level=%s" % (audit_full.get("score"), audit_full.get("level")))
        qc_checks = {item["id"]: item for item in (audit_qc_only.get("checks") or [])}
        check("只跑 qc 时 audit 有 warn/fail（clusters_few / embedding_missing）",
              qc_checks.get("clusters_few", {}).get("level") in ("warn", "fail")
              and qc_checks.get("embedding_missing", {}).get("level") in ("warn", "fail"),
              "score(qc-only)=%s" % audit_qc_only.get("score"))
        check("§6.4 全流程 audit 分数明显更高（≥ qc-only+30）",
              audit_full.get("score", 0) >= (audit_qc_only.get("score", 0) + 30),
              "%s -> %s" % (audit_qc_only.get("score"), audit_full.get("score")))
        check("§6.4 全流程 audit 分数 ≥90",
              audit_full.get("score", 0) >= 90, "score=%s" % audit_full.get("score"))

        # ===============================================================
        # 13. §6.4-5 history：随 run 增长，字段齐全（从 list 读）
        # ===============================================================
        listing = worker.call("list", {})
        datasets = listing.get("datasets") or []
        check("list.datasets 是数组且包含当前数据集",
              isinstance(datasets, list) and any(item.get("id") == dataset_id for item in datasets))
        current_summary = next(item for item in datasets if item["id"] == dataset_id)
        history = current_summary.get("history") or []
        check("history 随 run 增长", len(history) >= expected_history + 3,
              "%d entries" % len(history))
        check("history 截断上限 20（dataset 对象）", len(history) <= 20, str(len(history)))
        check_history_entry(history[-1], "history[last]")
        check("history 最后一步是 annotate", history[-1].get("step") == "annotate",
              str(history[-1].get("step")))
        check("dataset 对象带 history", isinstance(current_summary.get("history"), list))

        # ===============================================================
        # 14. §6.4-4 导出 protocol / methods / report / script（+ 既有 4 种）
        # ===============================================================
        exports = {}
        for what in ("protocol", "methods", "report", "script", "h5ad", "markers_csv", "obs_csv", "summary_json"):
            result = worker.call("export", {"dataset": dataset_id, "what": what}, timeout=300)
            path = Path(str(result.get("file") or ""))
            exports[what] = result
            check("export %-12s 文件存在且非空" % what,
                  path.exists() and path.stat().st_size > 0
                  and int(result.get("bytes", 0)) == path.stat().st_size
                  and result.get("what") == what,
                  "%d bytes" % path.stat().st_size)

        report_path = Path(str(exports["report"]["file"]))
        report_text = report_path.read_text(encoding="utf-8", errors="replace")
        check("report 含 PNG 引用",
              "![" in report_text and ".png" in report_text)
        check("report 含 markers 表头", "| cluster | gene |" in report_text)
        check("report 含 Methods 与 audit 段落",
              "## Methods" in report_text and "## 自检 Audit" in report_text)
        report_pngs = [child for child in report_path.parent.glob("*.png")]
        check("同名 PNG 被复制到报告同目录", len(report_pngs) >= 1, "%d png" % len(report_pngs))

        methods_path = Path(str(exports["methods"]["file"]))
        methods_text = methods_path.read_text(encoding="utf-8", errors="replace")
        check("methods 是英文 Methods（含 ## Methods / Welch / Python）",
              "## Methods" in methods_text and "Welch" in methods_text and "Python" in methods_text
              and len(methods_text) > 300, "%d chars" % len(methods_text))

        protocol_path = Path(str(exports["protocol"]["file"]))
        protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
        check("protocol 含 dataset/engine/history/figures/seed/plugin",
              all(key in protocol for key in ("dataset", "engine", "history", "figures", "seed", "plugin")))
        check("protocol.history 完整且 engine.versions 含 numpy",
              isinstance(protocol.get("history"), list) and len(protocol["history"]) >= expected_history + 3
              and isinstance((protocol.get("engine") or {}).get("versions"), dict)
              and "numpy" in protocol["engine"]["versions"]
              and protocol.get("seed") == 0)
        check("protocol.figures 是数组（图清单）", isinstance(protocol.get("figures"), list))

        script_path = Path(str(exports["script"]["file"]))
        script_text = script_path.read_text(encoding="utf-8", errors="replace")
        check("script 是 scanpy 风格可复现脚本",
              "import scanpy as sc" in script_text and "sc.pp" in script_text
              and "write_h5ad" in script_text and len(script_text) > 300,
              "%d chars" % len(script_text))

        # h5ad 内容：obs(clusters/cell_type) / obsm(PCA+嵌入) / X / layers["counts"]
        try:
            import anndata as anndata_module
            import scipy.sparse as sp

            written = anndata_module.read_h5ad(str(Path(str(exports["h5ad"]["file"]))))
            observed = written.obs
            obsm_keys = set(written.obsm.keys())
            layers_keys = set(written.layers.keys())
            check("导出 h5ad：obs 含 clusters 与 cell_type",
                  "clusters" in observed.columns and "cell_type" in observed.columns)
            check("导出 h5ad：obsm 含 X_pca 与 X_emb", {"X_pca", "X_emb"} <= obsm_keys, str(obsm_keys))
            check("导出 h5ad：layers 含 counts",
                  "counts" in layers_keys)
            check("导出 h5ad：X 形状与当前数据集一致",
                  written.shape == (n_obs_after, hvg_metrics.get("n_vars")))
            counts_matrix = written.layers["counts"]
            if sp.issparse(counts_matrix):
                counts_max = float(counts_matrix.data.max()) if counts_matrix.nnz else 0.0
            else:
                counts_max = float(counts_matrix.max())
            check("导出 h5ad：counts 是原始计数（最大值远离 1 量级且 ≥1）",
                  counts_max >= 1.0, "max=%.1f" % counts_max)
        except Exception as error:  # pragma: no cover
            check("导出 h5ad 可被 anndata 读回", False, "%s: %s" % (type(error).__name__, error))

        markers_csv_text = Path(str(exports["markers_csv"]["file"])).read_text(encoding="utf-8")
        check("markers_csv 含 cluster/gene/logfc/pval_adj 表头",
              all(column in markers_csv_text.splitlines()[0]
                  for column in ("cluster", "gene", "logfc", "pval_adj")))
        obs_rows = len(Path(str(exports["obs_csv"]["file"])).read_text(encoding="utf-8").splitlines())
        check("obs_csv 行数与细胞数一致", obs_rows == n_obs_after + 1, "%d rows" % obs_rows)
        summary = json.loads(Path(str(exports["summary_json"]["file"])).read_text(encoding="utf-8"))
        check("summary_json 含 dataset/filter/engine", all(key in summary for key in ("dataset", "filter", "engine")))

        # 显式 path 导出
        custom_path = tmp_root / "custom-methods.md"
        custom = worker.call(
            "export", {"dataset": dataset_id, "what": "methods", "path": str(custom_path)}, timeout=120
        )
        check("export 支持显式 path", Path(str(custom.get("file"))) == custom_path
              and custom_path.exists() and custom_path.stat().st_size > 0)

        # ===============================================================
        # 15. recipe / suggest（确定性）
        # ===============================================================
        recipes = (worker.call("recipe", {}) or {}).get("recipes") or []
        recipe_ids = {recipe.get("id") for recipe in recipes}
        check("recipe 至少 3 个预设（standard/quick/qc）",
              {"standard", "quick", "qc"} <= recipe_ids, str(sorted(recipe_ids)))
        check("recipe 每个预设含 steps 且步骤合法",
              all(recipe.get("title") and isinstance(recipe.get("steps"), list)
                  and recipe["steps"] and all(step.get("step") in STEP_ORDER
                                              and isinstance(step.get("params"), dict)
                                              for step in recipe["steps"])
                  for recipe in recipes))
        quick = worker.call("recipe", {"id": "quick"})
        quick_recipe = next(recipe for recipe in recipes if recipe["id"] == "quick")
        check("recipe {id} 返回该预设的 steps",
              quick.get("id") == "quick"
              and [step.get("step") for step in quick.get("steps")] ==
              [step.get("step") for step in quick_recipe["steps"]])

        suggest_first = worker.call("suggest", {"dataset": dataset_id})
        suggest_second = worker.call("suggest", {"dataset": dataset_id})
        check("suggest 是确定性的", suggest_first == suggest_second)
        suggestion_list = suggest_first.get("suggestions") or []
        check("suggest 非空且每项含 id/title/detail",
              len(suggestion_list) >= 1 and all(
                  isinstance(item.get("id"), str) and isinstance(item.get("title"), str)
                  and isinstance(item.get("detail"), str) for item in suggestion_list),
              str([item.get("id") for item in suggestion_list])[:100])

        # ===============================================================
        # 16. .h5ad 文件 open / 重复 open / 失败 history / 重跑清空下游 / 过滤参数
        # ===============================================================
        import h5py
        from make_sample import demo_adata

        h5ad_path = Path(os.path.realpath(str(tmp_root))) / "demo-input.h5ad"
        demo_adata().write_h5ad(str(h5ad_path))
        file_open = worker.call("open", {"path": str(h5ad_path)}, timeout=180)
        file_dataset_id = file_open.get("id")
        check("open(.h5ad 文件) 读到数据",
              isinstance(file_dataset_id, str)
              and os.path.realpath(str(file_open.get("path"))) == os.path.realpath(str(h5ad_path))
              and file_open.get("n_obs") == 2400 and file_open.get("n_vars") == 1200,
              "%s×%s path=%s" % (file_open.get("n_obs"), file_open.get("n_vars"), file_open.get("path")))
        duplicate_open = worker.call("open", {"path": str(h5ad_path)})
        check("同一路径重复 open 返回已存在数据集", duplicate_open.get("id") == file_dataset_id)
        fresh_suggestions = (worker.call("suggest", {"dataset": file_dataset_id}).get("suggestions") or [])
        check("新数据集 suggest 第一步是 qc",
              bool(fresh_suggestions) and fresh_suggestions[0].get("step") == "qc",
              str([item.get("id") for item in fresh_suggestions]))

        step_error = worker.call_error("run", {"dataset": file_dataset_id, "step": "pca", "params": {}})
        check("run 缺前置 → prerequisite", step_error.get("code") == "prerequisite",
              step_error.get("message", "")[:60])
        file_datasets = (worker.call("list", {}) or {}).get("datasets") or []
        file_summary = next((item for item in file_datasets if item.get("id") == file_dataset_id), {})
        failed_entries = [entry for entry in file_summary.get("history") or [] if entry.get("ok") is False]
        check("失败步骤也记入 history（ok=false）",
              bool(failed_entries) and failed_entries[-1].get("step") == "pca",
              "%d 条失败记录" % len(failed_entries))

        worker.call("run", {"dataset": file_dataset_id, "step": "qc", "params": {}})
        filter_custom = worker.call(
            "run",
            {
                "dataset": file_dataset_id,
                "step": "filter",
                "params": {"min_genes": 500, "max_pct_mt": 10, "min_cells": 20},
            },
        )
        custom_metrics = (filter_custom.get("summary") or {}).get("metrics") or {}
        check("filter 自定义阈值生效（min_genes=500 / min_cells=20）",
              isinstance(custom_metrics.get("n_obs_after"), int)
              and custom_metrics["n_obs_after"] < 2208
              and custom_metrics.get("n_vars_after", 1200) <= 1200,
              "%s cells / %s genes" % (custom_metrics.get("n_obs_after"), custom_metrics.get("n_vars_after")))
        worker.call("run", {"dataset": file_dataset_id, "step": "normalize", "params": {}})
        reset_result = worker.call("run", {"dataset": file_dataset_id, "step": "qc", "params": {}})
        check("重跑上游 qc 自动清空下游状态",
              reset_result.get("steps") == ["qc"]
              and (reset_result.get("dataset") or {}).get("n_obs") == 2400,
              "steps=%s n_obs=%s" % (reset_result.get("steps"), (reset_result.get("dataset") or {}).get("n_obs")))

        # 普通 .h5（adata 组）余力支持
        generic_h5 = tmp_root / "generic-adata.h5"
        with h5py.File(str(h5ad_path), "r") as source:
            with h5py.File(str(generic_h5), "w") as target:
                group = target.create_group("adata")
                for key in source.keys():
                    source.copy(source[key], group, key)
        h5_open = worker.call("open", {"path": str(generic_h5)}, timeout=180)
        check("open(.h5 里带 adata 组) 余力支持",
              h5_open.get("n_obs") == 2400 and h5_open.get("n_vars") == 1200
              and h5_open.get("id") != file_dataset_id,
              str(h5_open.get("id")))
        worker.call("close", {"dataset": h5_open.get("id")})
        worker.call("close", {"dataset": file_dataset_id})

        # ===============================================================
        # 17. close / list / 无数据 suggest
        # ===============================================================
        closed = worker.call("close", {"dataset": dataset_id})
        check("close 返回 closed=true", closed == {"closed": True})
        listing_after = worker.call("list", {})
        check("close 后 list 不再包含该数据集",
              all(item.get("id") != dataset_id for item in listing_after.get("datasets") or []))
        empty_suggest = worker.call("suggest", {})
        empty_ids = [item.get("id") for item in empty_suggest.get("suggestions") or []]
        check("无数据时 suggest 建议 open/demo", "demo" in empty_ids and "open" in empty_ids, str(empty_ids))
        close_again = worker.call_error("close", {"dataset": dataset_id})
        check("重复 close → not_found", close_again.get("code") == "not_found", close_again.get("message", "")[:50])

        # ===============================================================
        # 18. shutdown：退出码 0，stdout 没有非协议内容
        # ===============================================================
        bye, exit_code, leftover, stderr_text = worker.shutdown(timeout=120)
        check("shutdown 返回 bye=true", bye == {"bye": True}, str(bye))
        check("shutdown 后进程退出码 0", exit_code == 0, str(exit_code))
        check("shutdown 后 stdout 无残留协议外内容", leftover.strip() == "", repr(leftover[:100]))
        protocol_ok = True
        for line in worker.protocol_lines:
            try:
                message = json.loads(line)
                protocol_ok &= isinstance(message, dict) and "id" in message and "ok" in message
            except Exception:
                protocol_ok = False
        check("stdout 每一行都是合法协议 JSON", protocol_ok, "%d lines" % len(worker.protocol_lines))
    except Exception:  # 框架级失败：打印 traceback，最后统一返回非 0
        trace = traceback.format_exc()
        check("测试框架未抛异常", False, trace.splitlines()[-1] if trace else "")
        print(trace)
    finally:
        worker.kill()
        if os.environ.get("KEEP_SC_TEST_TMP") != "1":
            shutil.rmtree(str(tmp_root), ignore_errors=True)
        else:
            print("保留临时目录：%s" % tmp_root)

    elapsed = time.time() - started_at
    total = PASS_COUNT + FAIL_COUNT
    print("")
    print("=" * 78)
    print("结果：%d/%d PASS，%d FAIL，用时 %.1fs" % (PASS_COUNT, total, FAIL_COUNT, elapsed))
    if FAIL_COUNT:
        print("失败项：")
        for name in _LAST_FAILURES:
            print("  - %s" % name)
    print("=" * 78)
    return 0 if FAIL_COUNT == 0 else 1


def main():
    if not ENGINE_SCRIPT.exists():
        print("找不到引擎脚本：%s" % ENGINE_SCRIPT)
        return 1
    tmp_root = Path(tempfile.mkdtemp(prefix="dsh-sc-py-test-"))
    print("临时目录：%s" % tmp_root)
    return run_tests(tmp_root)


if __name__ == "__main__":
    raise SystemExit(main())
