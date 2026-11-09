#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell 合成 demo 数据生成器（既是脚本也是可 import 的模块）。

用途
----
1. `lib/sc_engine.py` 的 `open {demo:true}` 在内存里调用 `generate_demo()`，
   不落盘；demo 名 `demo-pbmc`。
2. 本文件命令行把同一份数据写成 `.h5ad`：
       python3 test/make_sample.py --out /tmp/demo.h5ad
3. `test/py_pipeline_test.py` 直接 import 本模块，拿到真实细胞类型
   （`obs["cell_type"]`）与每类注入的 marker 基因，用来校验聚类 / 标记基因 /
   自动注释的正确性。

设计
----
- 2400 细胞 × 1200 基因（可调），6 种真实 PBMC 细胞类型，每类 12–14 个真实
  marker 基因强表达；另有一个“低质量”亚群（线粒体占比高、检测基因少），
  保证 QC / 过滤步骤有意义。
- 所有 16 类内置参考 marker 面板（`BUILTIN_MARKER_PANEL`）的基因都进入表达
  矩阵（低质量亚群以外是低表达），所以 `autoannotate` 可以在 demo 上工作。
- 随机种子固定（默认 seed=0），同参数两次生成的数据**逐元素完全一致**。
- 生成 < 3 秒（默认参数在 Apple Silicon 上约 0.2–0.4 秒）。
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# 内置参考 marker 面板（autoannotate 用，≥14 类，每类 10–12 个真实基因）
# 键 = 细胞类型名；值 = 参考 marker 基因列表。
# ---------------------------------------------------------------------------
BUILTIN_MARKER_PANEL = {
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK", "TRAC", "TRBC1", "TRBC2", "LAT", "ZAP70", "LTB"],
    "CD4 T": ["CD4", "IL7R", "CCR7", "LEF1", "TCF7", "ICOS", "CD40LG", "TRAT1", "CD28", "MAL"],
    "CD8 T": ["CD8A", "CD8B", "GZMK", "CCL5", "GZMH", "CST7", "KLRD1", "EOMES", "RUNX3", "PRF1", "NKG7"],
    "NK cell": ["NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "KLRF1", "NCAM1", "FCGR3A", "TYROBP", "CST7", "CTSW", "GZMA"],
    "B cell": ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "CD37", "CD22", "TNFRSF13C", "CD74", "HLA-DRA"],
    "Plasma cell": ["MZB1", "JCHAIN", "SDC1", "XBP1", "IGHG1", "IGKC", "DERL3", "PRDM1", "TNFRSF17", "SSR4"],
    "Monocyte": ["LYZ", "CD14", "FCGR3A", "S100A8", "S100A9", "CST3", "FCN1", "CTSD", "HLA-DRA", "VCAN", "CD68", "CSF1R"],
    "Macrophage": ["CD68", "CD163", "MRC1", "MSR1", "C1QA", "C1QB", "APOE", "CTSB", "CTSL", "FABP5"],
    "Dendritic cell": ["FCER1A", "CST3", "CD1C", "CLEC10A", "HLA-DRA", "HLA-DPA1", "CD74", "ITGAX", "LAMP3", "CCR7", "IRF8", "CD1E"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9", "CXCR2", "MNDA", "FGR", "SLC25A37", "IFITM2", "LILRB2"],
    "Mast cell": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC", "CMA1", "GATA2", "SLC18A2", "HPGDS"],
    "Platelet": ["PPBP", "PF4", "ITGA2B", "GP9", "TUBB1", "SELP", "NRGN", "CLU", "GP1BA", "TREML1", "RGS18", "MYL9"],
    "Erythrocyte": ["HBB", "HBA1", "HBA2", "ALAS2", "AHSP", "GYPA", "SLC4A1", "CA1", "BPGM", "EPB42"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR", "EMCN", "CDH5", "ESAM", "RAMP2", "ENG", "PLVAP"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "FBLN1", "SPARC", "FN1", "MMP2"],
    "Epithelial": ["EPCAM", "KRT18", "KRT8", "KRT19", "KRT7", "CDH1", "CLDN4", "MUC1", "TACSTD2", "SFN"],
}

#: demo 的 6 种真实细胞类型（都来自内置面板）
DEMO_TYPES = ["T cell", "B cell", "NK cell", "Monocyte", "Dendritic cell", "Platelet"]

#: demo 每种类型额外注入的 2 个 marker（面板之外），凑到 12–14 个强 marker
_DEMO_EXTRA_MARKERS = {
    "T cell": ["CD6", "ITM2A"],
    "B cell": ["FCRLA", "VPREB3"],
    "NK cell": ["SH2D1B", "GZMM"],
    "Monocyte": ["CD300E", "LILRB1"],
    "Dendritic cell": ["FLT3", "BATF3"],
    "Platelet": ["ITGB3", "PROS1"],
}

#: 线粒体基因（MT- 前缀，QC 的 pct_counts_mt 依赖）
MITO_GENES = [
    "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8",
    "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4", "MT-CYB",
]

#: 背景基因里的真实常用基因（其余用 BG0001… 补齐）
#: 细胞周期背景基因（S 期 + G2M 期）：让 demo 数据也能演示 cell_cycle 分析
_CELL_CYCLE_BACKGROUND = [
    # S 期（DNA 复制）
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS", "RRM1", "RRM2",
    "GINS2", "CDC6", "CDT1", "POLD2", "POLE2", "UNG", "UHRF1", "CENPU", "CENPM", "FEN1",
    "MCM10", "NASP", "TMPO", "SLC25A5", "PRIM1", "POLA1", "POLA2", "DSCC1", "CHAF1B", "MRPL36",
    # G2M 期（有丝分裂）
    "CDK1", "CCNB1", "CCNB2", "CCNA2", "TOP2A", "AURKA", "AURKB", "BIRC5", "CDC20", "CDC25C",
    "PLK1", "KIF11", "KIF20B", "KIF23", "CENPF", "CENPE", "CENPA", "UBE2C", "TPX2", "TACC3",
    "ANLN", "ECT2", "PRC1", "RACGAP1", "MKI67", "CDCA8", "DLGAP5", "HJURP", "NUF2", "SPAG5",
]


_REAL_BACKGROUND = [
    "ACTB", "ACTG1", "GAPDH", "PGK1", "PPIA", "B2M", "TMSB10", "TMSB4X",
    "EEF1A1", "EEF1B2", "EEF2", "EIF4B", "HSP90AA1", "HSPA1A", "HSPA1B", "HSPB1",
    "FTH1", "FTL", "UBC", "UBA52", "UBB", "MALAT1", "XIST", "NEAT1",
    "COX4I1", "COX5A", "ATP5F1E", "NDUFA4", "SDHB", "CYCS", "CCT3", "CCT5",
    "CANX", "CALR", "PPIB", "SERPINA1", "SERPINB1", "GSTP1", "TXN", "TXNIP",
    "PRDX1", "PRDX2", "LDHA", "PKM", "ENO1", "TPI1", "ALDOA", "GPI",
    "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8", "RPL9", "RPL10",
    "RPL11", "RPL12", "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL19",
    "RPL21", "RPL22", "RPL23", "RPL24", "RPL26", "RPL27", "RPL28", "RPL29",
    "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL36", "RPL37", "RPL38",
    "RPL39", "RPL41", "RPS2", "RPS3", "RPS4X", "RPS5", "RPS6", "RPS7",
    "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15",
    "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24",
    "RPS25", "RPS26", "RPS27", "RPS28", "RPS29", "RPS30", "RPL7A", "RPL10A",
    "RPS3A", "RPS27A", "POLR2A", "SRP9", "SRP14", "YBX1", "YBX3", "HMGB1",
    "HMGB2", "H3-3A", "H3-3B", "H4C1", "H2AC1", "H1-4", "PTMA", "SET",
    "NPM1", "NCL", "HNRNPA1", "HNRNPA2B1", "SRSF1", "SRSF2", "SRSF3", "SRSF7",
    "DDX5", "DDX17", "EIF3A", "EIF3B", "EIF5B", "VIM", "CFL1", "GSN",
    "TAGLN2", "ANXA2", "ANXA1", "S100A11", "LGALS1", "LGALS3", "IFITM3", "IFITM1",
    "BST2", "CD52", "SPCS3", "SSR4P1", "ATP1A1", "ATP1B3", "GNAI2", "GNB2L1",
    "RACK1", "STAT1", "STAT3", "JUNB", "FOS", "FOSB", "JUN", "EGR1",
    "IER2", "DUSP1", "KLF2", "KLF6", "CEBPD", "NFKBIA", "ZFP36", "BTG1",
    "BTG2", "CXCR4", "CD69", "CD44", "CD53", "CD48", "CD46", "CD59",
    "CD47", "CD55", "CD63", "CD81", "CD82", "CD99", "CD151", "CD164",
]


def demo_marker_genes():
    """返回 demo 每个细胞类型强表达的 marker 基因（12–14 个/类）。"""
    out = {}
    for cell_type in DEMO_TYPES:
        genes = list(BUILTIN_MARKER_PANEL[cell_type])
        for gene in _DEMO_EXTRA_MARKERS.get(cell_type, []):
            if gene not in genes:
                genes.append(gene)
        out[cell_type] = genes
    return out


def _build_gene_names(n_genes, marker_map):
    """构造基因名列表：线粒体 + demo marker + 面板基因 + 真实背景 + 补齐。"""
    required = list(MITO_GENES)
    for cell_type in DEMO_TYPES:
        required.extend(marker_map[cell_type])
    for cell_type in sorted(BUILTIN_MARKER_PANEL):
        required.extend(BUILTIN_MARKER_PANEL[cell_type])
    required.extend(_CELL_CYCLE_BACKGROUND)
    required.extend(_REAL_BACKGROUND)

    names = []
    seen = set()
    for gene in required:
        if gene not in seen:
            seen.add(gene)
            names.append(gene)
    if n_genes < len(names):
        raise ValueError(
            f"n_genes={n_genes} 太小：需要至少 {len(names)} 个基因才能放下内置 marker 面板与 demo marker"
        )
    index = 1
    while len(names) < n_genes:
        gene = "BG%04d" % index
        index += 1
        if gene not in seen:
            seen.add(gene)
            names.append(gene)
    return names[:n_genes]


def generate_demo(n_cells=2400, n_genes=1200, seed=0, low_quality_frac=0.08):
    """生成 demo 数据（完全确定性）。

    返回 dict：
      counts        scipy.sparse.csr_matrix(float32)，原始计数 (n_cells × n_genes)
      var_names     list[str]，基因名（顺序与 counts 列一致）
      obs           pandas.DataFrame，index=细胞名，列 cell_type / batch
      marker_genes  dict[细胞类型 -> list[str]]，注入的强 marker
      cell_types    list[str]，6 种真实类型
      is_low_quality numpy bool 数组，哪些细胞是低质量亚群
      mito_genes    list[str]
    """
    rng = np.random.default_rng(int(seed))
    marker_map = demo_marker_genes()
    gene_names = _build_gene_names(int(n_genes), marker_map)
    n_total = len(gene_names)
    gene_index = {gene: i for i, gene in enumerate(gene_names)}

    n_lq = int(round(int(n_cells) * float(low_quality_frac)))
    n_lq = max(0, min(n_lq, int(n_cells) - len(DEMO_TYPES)))
    n_clean = int(n_cells) - n_lq

    # 每个细胞的真实类型（低质量亚群单独标出）
    clean_types = np.array([DEMO_TYPES[i % len(DEMO_TYPES)] for i in range(n_clean)], dtype=object)
    rng.shuffle(clean_types)
    cell_types = np.empty(int(n_cells), dtype=object)
    cell_types[:] = "Low quality"
    is_low_quality = np.zeros(int(n_cells), dtype=bool)
    if n_lq > 0:
        lq_positions = rng.choice(int(n_cells), size=n_lq, replace=False)
        is_low_quality[lq_positions] = True
    clean_mask = ~is_low_quality
    cell_types[clean_mask] = clean_types

    # 文库大小（细胞间差异）；低质量亚群整体偏低
    size_factor = np.clip(rng.lognormal(mean=0.0, sigma=0.22, size=int(n_cells)), 0.55, 1.9)
    size_factor[is_low_quality] *= 0.13

    # 背景表达权重（lognormal，长尾）
    base_weight = rng.lognormal(mean=0.0, sigma=0.75, size=n_total)
    base_weight = base_weight / base_weight.sum()
    lam = np.outer(size_factor * 950.0, base_weight)

    # 每种细胞类型给它自己的 marker 基因加强表达
    for cell_type in DEMO_TYPES:
        indices = np.array([gene_index[g] for g in marker_map[cell_type]], dtype=np.int64)
        rates = rng.uniform(3.5, 7.5, size=indices.size)
        mask = cell_types == cell_type
        if mask.any():
            lam[np.ix_(mask, indices)] += size_factor[mask, None] * rates[None, :]

    # 线粒体基因：正常细胞少量，低质量亚群大量（让 pct_counts_mt 过滤有意义）
    mito_indices = np.array([gene_index[g] for g in MITO_GENES], dtype=np.int64)
    normal_mask = ~is_low_quality
    if normal_mask.any():
        lam[np.ix_(normal_mask, mito_indices)] += size_factor[normal_mask, None] * 0.6
    if is_low_quality.any():
        mito_boost = rng.uniform(12.0, 18.0, size=(int(is_low_quality.sum()), len(MITO_GENES)))
        lam[np.ix_(is_low_quality, mito_indices)] += mito_boost

    counts = rng.poisson(lam).astype(np.int32)
    counts_csr = sp.csr_matrix(counts.astype(np.float32))
    counts_csr.sort_indices()

    batches = np.array(["batch1", "batch2"] * (int(n_cells) // 2 + 1), dtype=object)[: int(n_cells)]
    rng.shuffle(batches)
    cell_index = pd.Index(["cell%05d" % i for i in range(int(n_cells))], name="cell")
    obs = pd.DataFrame(
        {
            "cell_type": pd.Series(np.asarray(cell_types, dtype=object), index=cell_index, dtype=object),
            "batch": pd.Series(np.asarray(batches, dtype=object), index=cell_index, dtype=object),
        },
        index=cell_index,
    )
    return {
        "counts": counts_csr,
        "var_names": list(gene_names),
        "obs": obs,
        "marker_genes": {k: list(v) for k, v in marker_map.items()},
        "cell_types": list(DEMO_TYPES),
        "is_low_quality": is_low_quality,
        "mito_genes": list(MITO_GENES),
    }


def demo_adata(**kwargs):
    """生成 demo 的 AnnData（X=counts，layers['counts']=原始计数）。"""
    import anndata as ad

    data = generate_demo(**kwargs)
    obs = data["obs"].copy()
    var = pd.DataFrame(index=pd.Index(data["var_names"], name="gene"))
    adata = ad.AnnData(X=data["counts"], obs=obs, var=var)
    adata.layers["counts"] = data["counts"].copy()
    adata.uns["source"] = "dsh-singlecell demo generator (test/make_sample.py)"
    return adata


def _summarize(data, seconds):
    obs = data["obs"]
    type_counts = {str(k): int(v) for k, v in obs["cell_type"].value_counts().sort_index().items()}
    return {
        "n_cells": int(obs.shape[0]),
        "n_genes": int(len(data["var_names"])),
        "cell_types": type_counts,
        "n_low_quality": int(data["is_low_quality"].sum()),
        "n_marker_genes": int(sum(len(v) for v in data["marker_genes"].values())),
        "n_mito_genes": int(len(data["mito_genes"])),
        "seconds": round(float(seconds), 3),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="生成 dsh-singlecell 的合成 demo 单细胞数据")
    parser.add_argument("--out", default=None, help="输出 .h5ad 路径（省略则只在内存生成并打印摘要）")
    parser.add_argument("--cells", type=int, default=2400)
    parser.add_argument("--genes", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--low-quality-frac", type=float, default=0.08)
    parsed = parser.parse_args(argv)

    started = time.time()
    data = generate_demo(
        n_cells=parsed.cells, n_genes=parsed.genes, seed=parsed.seed,
        low_quality_frac=parsed.low_quality_frac,
    )
    info = _summarize(data, time.time() - started)
    if parsed.out:
        out_path = os.path.abspath(os.path.expanduser(parsed.out))
        parent = os.path.dirname(out_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        adata = demo_adata(
            n_cells=parsed.cells, n_genes=parsed.genes, seed=parsed.seed,
            low_quality_frac=parsed.low_quality_frac,
        )
        adata.write_h5ad(out_path)
        info["out"] = out_path
        info["bytes"] = int(os.path.getsize(out_path))
    print(json.dumps(info, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
