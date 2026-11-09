# -*- coding: utf-8 -*-
"""扫描根目录 / 样本表合并 / 已有进度检测（sc_scan / sc_metadata / sc_state）单元测试。"""
import gzip
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))

import sc_scan  # noqa: E402
import sc_metadata  # noqa: E402
import sc_io  # noqa: E402
import sc_state  # noqa: E402


def build_fixtures(root: Path):
    rng = np.random.default_rng(7)
    genes = ["G%04d" % i for i in range(1, 41)]
    geo = root / "geo"
    geo.mkdir(parents=True)
    from scipy import sparse
    from scipy.io import mmwrite
    for gsm, n_cells in (("GSM001", 20), ("GSM002", 25)):
        cells = ["%s_c%d-1" % (gsm, i) for i in range(n_cells)]
        matrix = sparse.random(40, n_cells, density=0.3, random_state=1, format="csc",
                               data_rvs=lambda n: rng.integers(0, 5, n).astype(float))
        with gzip.open(geo / ("%s_matrix.mtx.gz" % gsm), "wb") as handle:
            mmwrite(handle, matrix.tocoo())
        (geo / ("%s_features.tsv.gz" % gsm)).write_bytes(gzip.compress(("\n".join(genes) + "\n").encode()))
        (geo / ("%s_barcodes.tsv.gz" % gsm)).write_bytes(gzip.compress(("\n".join(cells) + "\n").encode()))
    tenx = root / "sampleA" / "filtered_feature_bc_matrix"
    tenx.mkdir(parents=True)
    matrix = sparse.random(40, 15, density=0.3, random_state=2, format="csc")
    mmwrite(tenx / "matrix.mtx", matrix.tocoo())
    (tenx / "features.tsv").write_text("\n".join("%s\t%s\tGene Expression" % (g, g) for g in genes), encoding="utf-8")
    (tenx / "barcodes.tsv").write_text("\n".join("b%d" % i for i in range(15)), encoding="utf-8")
    return genes


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print("  [%s] %s%s" % (status, name, (" — " + detail) if detail else ""))
    return bool(ok)


def main():
    print("=" * 70)
    print("py_scan_test：根目录扫描 / 样本表合并 / 进度检测")
    passed = failed = 0
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        genes = build_fixtures(root)
        result = sc_scan.scan_data_root(root)
        samples = result["samples"]
        names = [item["name"] for item in samples]
        formats = [item["format"] for item in samples]
        for name, ok, detail in [
            ("扫描到 3 个样本（2 GEO 三联体 + 1 10X 目录）", len(samples) == 3, str(names)),
            ("GEO 三联体按 GSM 前缀配对", formats[:2] == ["geo_mtx", "geo_mtx"], str(formats)),
            ("嵌套 filtered_feature_bc_matrix 被识别", "sampleA" in names, str(names)),
            ("三联体返回 trio_prefix", samples[0].get("trio_prefix") == "GSM001", str(samples[0].get("trio_prefix"))),
            ("样本 id 唯一", len({item["id"] for item in samples}) == len(samples), ""),
        ]:
            passed += check(name, ok, detail); failed += (not ok)

        # 样本表合并：自定义名 + 大组 + trio_prefix
        geo = [item for item in samples if item["format"] == "geo_mtx"]
        sources = [{"source": item["path"], "name": "PT%d" % (index + 1),
                    "group": "PT" if index == 0 else "NT",
                    "options": {"trio_prefix": item["trio_prefix"]}} for index, item in enumerate(geo)]
        merged = sc_metadata.merge_sources(sc_io, sources, str(root / "work"))
        obs = merged["obs"]
        for name, ok, detail in [
            ("合并后细胞数 = 各样本之和", merged["X"].shape[0] == 45, str(merged["X"].shape)),
            ("sample 列使用自定义样本名", sorted(set(obs["sample"].tolist())) == ["PT1", "PT2"], str(set(obs["sample"].tolist()))),
            ("group 列来自样本表", sorted(set(obs["group"].dropna().tolist())) == ["NT", "PT"], str(set(obs["group"].dropna().tolist()))),
            ("group_column=group", merged["meta"].get("group_column") == "group", str(merged["meta"].get("group_column"))),
        ]:
            passed += check(name, ok, detail); failed += (not ok)

        # 已有进度检测：用轻量假 dataset
        class FakeDataset:
            def __init__(self, n_obs=10):
                import pandas as pd
                self.n_obs = n_obs
                self.obs = pd.DataFrame(index=["c%d" % i for i in range(n_obs)])
                self.var = pd.DataFrame(index=genes)
                self.var_names = np.asarray(genes)
                self.obsm = {}
                self.uns = {}
                self.layers = {}
                from scipy import sparse
                self.X = sparse.csr_matrix(np.abs(np.random.default_rng(0).normal(size=(n_obs, len(genes)))))

        dataset = FakeDataset()
        dataset.obs["seurat_clusters"] = [str(i % 3) for i in range(10)]
        dataset.obs["cell_type"] = ["T", "B"] * 5
        loaded = {"obsm": {"X_pca": np.zeros((10, 5), dtype=np.float32),
                           "X_emb": np.zeros((10, 2), dtype=np.float32)},
                  "import_state": {"reductions": ["pca", "umap"], "graphs": ["RNA_snn"], "variable_features": ["G0001"]},
                  "var_features": ["G0001"]}
        state = sc_state.apply(dataset, loaded)
        for name, ok, detail in [
            ("检测到聚类/嵌入/注释步骤", {"cluster", "embed", "annotate"} <= set(state["steps"]), str(state["steps"])),
            ("写入既有聚类列 clusters", "clusters" in dataset.obs.columns, ""),
            ("写入既有嵌入 X_emb", "X_emb" in dataset.obsm, ""),
            ("cell_type 映射且非默认", "cell_type" in dataset.obs.columns and dataset.cell_type_default is False, ""),
            ("mapped 报告 3 个簇", (state["mapped"] or {}).get("n_clusters") == 3, str((state["mapped"] or {}).get("n_clusters"))),
        ]:
            passed += check(name, ok, detail); failed += (not ok)

        # 原始数据不应冒充进度
        raw = FakeDataset()
        from scipy import sparse as _sp
        raw.X = _sp.csr_matrix(np.round(np.abs(np.random.default_rng(1).normal(size=(10, len(genes))))))
        state_raw = sc_state.detect(raw, {})
        ok_raw = check("原始 counts 不冒充已有进度", state_raw["steps"] == [], str(state_raw["steps"]))
        passed += ok_raw
        failed += (not ok_raw)

    print("-" * 70)
    print("结果：%d PASS / %d FAIL" % (passed, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
