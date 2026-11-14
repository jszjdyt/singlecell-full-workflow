#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench — 内置生物学先验（sc-analysis-b 模块专用）

内容
----
1. ``LR_PAIRS``：人源配体–受体（ligand–receptor）对，约 170 对，覆盖趋化因子、
   细胞因子、Notch、Wnt、生长因子、ECM/黏附、免疫检查点、Semaphorin/Ephrin 等
   常见通讯家族。用于 ``lr_communication`` 的 CellChat/NATMI-lite 打分。
2. ``TF_LIST``：人源转录因子（TF）符号列表，约 400 个，覆盖核受体、homeobox、
   forkhead、ETS、bZIP、GATA/KLF、IRF/STAT、E2F、SOX 等主要家族。用于
   ``regulon_activity`` 的 SCENIC-lite regulon 构建与 ``virtual_screen`` 候选。
3. ``TF_FAMILY``：TF → 家族名的映射（仅对列表中的家族成员给出，供表格/解释使用）。

诚实说明
--------
- 这些先验是**人工整理的文献常用列表**，不是 CellChatDB / CellPhoneDB / AnimalTFDB /
  Lambert et al. 2018 等官方数据库的完整导出，也没有做物种直系同源校验。
- 配体–受体对只收录"配体基因 × 受体基因"的简单二元关系，**不建模复合物**
  （如 IL12 = IL12A+IL12B、TGFB 二聚体）与辅助受体（如 IL6 的 IL6ST 已单列），
  因此与官方数据库的条目口径不同。
- TF 列表以"序列特异性 DNA 结合转录因子"为主，包含少量经典共调控因子/核受体；
  不保证覆盖全部人源 TF（官方约 1600 个）。
- 本文件只提供先验，不实现任何算法；算法在 ``sc_analyses_dynamics.py``。
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# 1. 配体–受体对（ligand, receptor, family）
# ---------------------------------------------------------------------------

_LR_CHEMOKINE = [
    ("CCL2", "CCR2", "chemokine"),
    ("CCL3", "CCR1", "chemokine"),
    ("CCL3", "CCR5", "chemokine"),
    ("CCL4", "CCR1", "chemokine"),
    ("CCL4", "CCR5", "chemokine"),
    ("CCL5", "CCR1", "chemokine"),
    ("CCL5", "CCR5", "chemokine"),
    ("CCL7", "CCR2", "chemokine"),
    ("CCL8", "CCR2", "chemokine"),
    ("CCL8", "CCR5", "chemokine"),
    ("CCL11", "CCR3", "chemokine"),
    ("CCL13", "CCR2", "chemokine"),
    ("CCL13", "CCR3", "chemokine"),
    ("CCL14", "CCR1", "chemokine"),
    ("CCL17", "CCR4", "chemokine"),
    ("CCL19", "CCR7", "chemokine"),
    ("CCL20", "CCR6", "chemokine"),
    ("CCL21", "CCR7", "chemokine"),
    ("CCL22", "CCR4", "chemokine"),
    ("CCL25", "CCR9", "chemokine"),
    ("CCL28", "CCR10", "chemokine"),
    ("CXCL1", "CXCR2", "chemokine"),
    ("CXCL2", "CXCR2", "chemokine"),
    ("CXCL3", "CXCR2", "chemokine"),
    ("CXCL5", "CXCR2", "chemokine"),
    ("CXCL6", "CXCR2", "chemokine"),
    ("CXCL8", "CXCR1", "chemokine"),
    ("CXCL8", "CXCR2", "chemokine"),
    ("CXCL9", "CXCR3", "chemokine"),
    ("CXCL10", "CXCR3", "chemokine"),
    ("CXCL11", "CXCR3", "chemokine"),
    ("CXCL12", "CXCR4", "chemokine"),
    ("CXCL12", "ACKR3", "chemokine"),
    ("CXCL13", "CXCR5", "chemokine"),
    ("CXCL16", "CXCR6", "chemokine"),
    ("CX3CL1", "CX3CR1", "chemokine"),
    ("XCL1", "XCR1", "chemokine"),
    ("XCL2", "XCR1", "chemokine"),
]

_LR_CYTOKINE = [
    ("TNF", "TNFRSF1A", "cytokine"),
    ("TNF", "TNFRSF1B", "cytokine"),
    ("LTA", "TNFRSF1A", "cytokine"),
    ("LTA", "TNFRSF1B", "cytokine"),
    ("LTB", "LTBR", "cytokine"),
    ("IL1A", "IL1R1", "cytokine"),
    ("IL1B", "IL1R1", "cytokine"),
    ("IL1RN", "IL1R1", "cytokine"),
    ("IL2", "IL2RA", "cytokine"),
    ("IL2", "IL2RB", "cytokine"),
    ("IL4", "IL4R", "cytokine"),
    ("IL6", "IL6R", "cytokine"),
    ("IL6", "IL6ST", "cytokine"),
    ("IL7", "IL7R", "cytokine"),
    ("IL10", "IL10RA", "cytokine"),
    ("IL10", "IL10RB", "cytokine"),
    ("IL12A", "IL12RB1", "cytokine"),
    ("IL12B", "IL12RB1", "cytokine"),
    ("IL13", "IL13RA1", "cytokine"),
    ("IL15", "IL2RB", "cytokine"),
    ("IL16", "CD4", "cytokine"),
    ("IL17A", "IL17RA", "cytokine"),
    ("IL18", "IL18R1", "cytokine"),
    ("IL21", "IL21R", "cytokine"),
    ("IL23A", "IL23R", "cytokine"),
    ("IL27", "IL27RA", "cytokine"),
    ("IFNG", "IFNGR1", "cytokine"),
    ("IFNG", "IFNGR2", "cytokine"),
    ("IFNA1", "IFNAR1", "cytokine"),
    ("IFNA2", "IFNAR1", "cytokine"),
    ("IFNB1", "IFNAR1", "cytokine"),
    ("CSF1", "CSF1R", "cytokine"),
    ("CSF2", "CSF2RA", "cytokine"),
    ("CSF3", "CSF3R", "cytokine"),
    ("TSLP", "CRLF2", "cytokine"),
    ("OSM", "OSMR", "cytokine"),
    ("LIF", "LIFR", "cytokine"),
    ("KITLG", "KIT", "cytokine"),
    ("FLT3LG", "FLT3", "cytokine"),
    ("EPO", "EPOR", "cytokine"),
    ("THPO", "MPL", "cytokine"),
    ("MIF", "CD74", "cytokine"),
    ("MIF", "CXCR4", "cytokine"),
    ("TGFB1", "TGFBR1", "cytokine"),
    ("TGFB1", "TGFBR2", "cytokine"),
    ("TGFB2", "TGFBR1", "cytokine"),
    ("TGFB3", "TGFBR2", "cytokine"),
    ("INHBA", "ACVR2A", "cytokine"),
    ("NODAL", "ACVR2B", "cytokine"),
    ("BMP2", "BMPR1A", "cytokine"),
    ("BMP4", "BMPR1A", "cytokine"),
    ("BMP7", "BMPR1B", "cytokine"),
    ("GDF5", "BMPR1B", "cytokine"),
]

_LR_SIGNALING = [
    ("DLL1", "NOTCH1", "Notch"),
    ("DLL1", "NOTCH2", "Notch"),
    ("DLL4", "NOTCH1", "Notch"),
    ("JAG1", "NOTCH1", "Notch"),
    ("JAG1", "NOTCH2", "Notch"),
    ("JAG2", "NOTCH1", "Notch"),
    ("JAG2", "NOTCH2", "Notch"),
    ("WNT3A", "FZD1", "Wnt"),
    ("WNT5A", "FZD1", "Wnt"),
    ("WNT5A", "ROR2", "Wnt"),
    ("WNT7A", "FZD5", "Wnt"),
    ("WNT2B", "FZD4", "Wnt"),
    ("WNT11", "FZD4", "Wnt"),
    ("PDGFA", "PDGFRA", "growth factor"),
    ("PDGFB", "PDGFRB", "growth factor"),
    ("PDGFC", "PDGFRB", "growth factor"),
    ("PDGFD", "PDGFRB", "growth factor"),
    ("EGF", "EGFR", "growth factor"),
    ("TGFA", "EGFR", "growth factor"),
    ("HBEGF", "EGFR", "growth factor"),
    ("AREG", "EGFR", "growth factor"),
    ("EREG", "EGFR", "growth factor"),
    ("BTC", "EGFR", "growth factor"),
    ("NRG1", "ERBB3", "growth factor"),
    ("NRG1", "ERBB4", "growth factor"),
    ("NRG2", "ERBB4", "growth factor"),
    ("NRG3", "ERBB4", "growth factor"),
    ("VEGFA", "KDR", "growth factor"),
    ("VEGFA", "FLT1", "growth factor"),
    ("VEGFA", "NRP1", "growth factor"),
    ("VEGFB", "FLT1", "growth factor"),
    ("VEGFC", "FLT4", "growth factor"),
    ("PGF", "FLT1", "growth factor"),
    ("FGF1", "FGFR1", "growth factor"),
    ("FGF2", "FGFR1", "growth factor"),
    ("FGF7", "FGFR2", "growth factor"),
    ("FGF10", "FGFR2", "growth factor"),
    ("HGF", "MET", "growth factor"),
    ("MST1", "MST1R", "growth factor"),
    ("IGF1", "IGF1R", "growth factor"),
    ("IGF2", "IGF1R", "growth factor"),
    ("IGF2", "IGF2R", "growth factor"),
    ("INS", "INSR", "growth factor"),
    ("ANGPT1", "TEK", "growth factor"),
    ("ANGPT2", "TEK", "growth factor"),
    ("ANGPT4", "TEK", "growth factor"),
]

_LR_ECM = [
    ("SPP1", "CD44", "ECM"),
    ("SPP1", "ITGAV", "ECM"),
    ("SPP1", "ITGB1", "ECM"),
    ("FN1", "ITGAV", "ECM"),
    ("FN1", "ITGB1", "ECM"),
    ("FN1", "CD44", "ECM"),
    ("LAMA1", "ITGA6", "ECM"),
    ("LAMB1", "ITGB1", "ECM"),
    ("COL1A1", "ITGB1", "ECM"),
    ("COL1A1", "DDR1", "ECM"),
    ("COL4A1", "ITGB1", "ECM"),
    ("THBS1", "CD36", "ECM"),
    ("THBS1", "CD47", "ECM"),
    ("TNC", "ITGAV", "ECM"),
    ("ANXA1", "FPR1", "ECM"),
]

_LR_ADHESION = [
    ("ICAM1", "ITGAL", "adhesion"),
    ("ICAM1", "ITGB2", "adhesion"),
    ("VCAM1", "ITGA4", "adhesion"),
    ("MADCAM1", "ITGB7", "adhesion"),
    ("SELL", "SELPLG", "adhesion"),
    ("SELP", "SELPLG", "adhesion"),
    ("PECAM1", "PECAM1", "adhesion"),
    ("CDH1", "CDH1", "adhesion"),
    ("CDH2", "CDH2", "adhesion"),
    ("LGALS9", "HAVCR2", "adhesion"),
    ("LGALS9", "CD44", "adhesion"),
    ("LGALS9", "P4HB", "adhesion"),
]

_LR_IMMUNE = [
    ("CD274", "PDCD1", "immune checkpoint"),
    ("PDCD1LG2", "PDCD1", "immune checkpoint"),
    ("CTLA4", "CD80", "immune checkpoint"),
    ("CTLA4", "CD86", "immune checkpoint"),
    ("CD80", "CD28", "immune checkpoint"),
    ("CD86", "CD28", "immune checkpoint"),
    ("LAG3", "HLA-DRA", "immune checkpoint"),
    ("TIGIT", "PVR", "immune checkpoint"),
    ("TIGIT", "NECTIN2", "immune checkpoint"),
    ("CD226", "PVR", "immune checkpoint"),
    ("CD47", "SIRPA", "immune checkpoint"),
    ("CD40", "CD40LG", "immune checkpoint"),
    ("CD40LG", "CD40", "immune checkpoint"),
    ("FASLG", "FAS", "immune checkpoint"),
    ("TNFSF10", "TNFRSF10A", "immune checkpoint"),
    ("TNFSF10", "TNFRSF10B", "immune checkpoint"),
    ("TNFSF13B", "TNFRSF13B", "immune checkpoint"),
    ("TNFSF13B", "TNFRSF13C", "immune checkpoint"),
    ("TNFSF12", "TNFRSF12A", "immune checkpoint"),
    ("CD70", "CD27", "immune checkpoint"),
    ("CD27", "CD70", "immune checkpoint"),
    ("KLRK1", "MICA", "immune checkpoint"),
    ("KLRK1", "MICB", "immune checkpoint"),
    ("KLRK1", "ULBP1", "immune checkpoint"),
    ("KLRK1", "ULBP2", "immune checkpoint"),
    ("HLA-DRA", "CD4", "immune checkpoint"),
]

_LR_GUIDANCE = [
    ("SEMA4D", "PLXNB1", "guidance"),
    ("SEMA4D", "PLXNB2", "guidance"),
    ("SEMA4A", "PLXNB2", "guidance"),
    ("SEMA4C", "PLXNB2", "guidance"),
    ("SEMA3A", "NRP1", "guidance"),
    ("SEMA3A", "PLXNA1", "guidance"),
    ("SEMA7A", "PLXNC1", "guidance"),
    ("EFNA1", "EPHA2", "guidance"),
    ("EFNB1", "EPHB2", "guidance"),
    ("EFNB2", "EPHB4", "guidance"),
]

_LR_COAG = [
    ("F2", "F2R", "coagulation"),
    ("PLAU", "PLAUR", "coagulation"),
    ("C3", "C3AR1", "complement"),
    ("C5", "C5AR1", "complement"),
]

#: 全部配体–受体对：list[(ligand, receptor, family)]
LR_PAIRS = list(
    _LR_CHEMOKINE
    + _LR_CYTOKINE
    + _LR_SIGNALING
    + _LR_ECM
    + _LR_ADHESION
    + _LR_IMMUNE
    + _LR_GUIDANCE
    + _LR_COAG
)

#: 配体–受体家族（保持出现顺序去重）
LR_FAMILIES = []
for _pair in LR_PAIRS:
    if _pair[2] not in LR_FAMILIES:
        LR_FAMILIES.append(_pair[2])

# ---------------------------------------------------------------------------
# 2. 人源转录因子列表（按家族分组后拼接）
# ---------------------------------------------------------------------------

# 核受体（经典配体依赖转录因子，含孤儿受体）
_TF_NUCLEAR = [
    "AR", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "NR0B1", "NR0B2",
    "NR1D1", "NR1D2", "NR1H2", "NR1H3", "NR1H4", "NR1I2", "NR1I3",
    "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6",
    "NR3C1", "NR3C2", "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1",
    "PGR", "PPARA", "PPARD", "PPARG", "RARA", "RARB", "RARG",
    "RXRA", "RXRB", "RXRG", "THRA", "THRB", "VDR",
]

# Homeobox（含 HOX 簇、PAX、POU、NKX、LHX、SIX、DLX 等）
_TF_HOMEOBOX = [
    "BARX1", "BARX2", "CDX1", "CDX2", "CDX4", "CUX1", "CUX2", "DBX1", "DBX2",
    "DLX1", "DLX2", "DLX3", "DLX4", "DLX5", "DLX6", "EMX1", "EMX2", "EN1",
    "EN2", "EVX1", "EVX2", "GSX2", "HESX1", "HLX",
    "HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9",
    "HOXA10", "HOXA11", "HOXA13",
    "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8",
    "HOXB9", "HOXB13",
    "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXC10", "HOXC11",
    "HOXC12", "HOXC13",
    "HOXD1", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11",
    "HOXD12", "HOXD13",
    "ISL1", "ISL2", "LBX1", "LBX2", "LHX1", "LHX2", "LHX3", "LHX4", "LHX5",
    "LHX6", "LHX8", "LHX9", "LMX1A", "LMX1B", "MEIS1", "MEIS2", "MEIS3",
    "MEOX1", "MEOX2", "MNX1", "MSX1", "MSX2",
    "NKX2-1", "NKX2-2", "NKX2-3", "NKX2-5", "NKX2-8", "NKX3-1",
    "NKX6-1", "NKX6-2", "NKX6-3", "OTP",
    "PAX1", "PAX2", "PAX3", "PAX4", "PAX5", "PAX6", "PAX7", "PAX8", "PAX9",
    "PBX1", "PBX2", "PBX3", "PBX4", "PDX1", "PITX1", "PITX2", "PITX3",
    "POU1F1", "POU2F1", "POU2F2", "POU2F3", "POU3F1", "POU3F2", "POU3F3",
    "POU3F4", "POU4F1", "POU4F2", "POU4F3", "POU5F1", "POU6F1", "POU6F2",
    "PROP1", "PRRX1", "PRRX2", "PROX1", "PROX2", "RAX", "RAX2",
    "SIX1", "SIX2", "SIX3", "SIX4", "SIX5", "SIX6",
    "TLX1", "TLX2", "TLX3", "UNCX", "VAX1", "VAX2", "VENTX", "VSX1", "VSX2",
]

# Forkhead
_TF_FOX = [
    "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXB2", "FOXC1", "FOXC2", "FOXD1",
    "FOXD2", "FOXD3", "FOXD4", "FOXE1", "FOXE3", "FOXF1", "FOXF2", "FOXG1",
    "FOXH1", "FOXI1", "FOXI2", "FOXI3", "FOXJ1", "FOXJ2", "FOXJ3", "FOXK1",
    "FOXK2", "FOXL1", "FOXL2", "FOXM1", "FOXN1", "FOXN2", "FOXN3", "FOXN4",
    "FOXO1", "FOXO3", "FOXO4", "FOXO6", "FOXP1", "FOXP2", "FOXP3", "FOXP4",
    "FOXQ1", "FOXS1",
]

# HMG / SOX / TCF-LEF
_TF_HMG = [
    "HMGA1", "HMGA2", "HBP1", "LEF1", "TCF7", "TCF7L1", "TCF7L2",
    "SOX1", "SOX2", "SOX3", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8", "SOX9",
    "SOX10", "SOX11", "SOX12", "SOX13", "SOX14", "SOX15", "SOX17", "SOX18",
    "SOX21", "SOX30", "SRY", "TOX", "TOX2", "TOX3", "TOX4",
]

# ETS
_TF_ETS = [
    "ELF1", "ELF2", "ELF3", "ELF4", "ELF5", "ELK1", "ELK3", "ELK4", "ERF",
    "ERG", "ETS1", "ETS2", "ETV1", "ETV2", "ETV3", "ETV4", "ETV5", "ETV6",
    "ETV7", "FEV", "FLI1", "GABPA", "SPDEF", "SPIB", "SPI1",
]

# bZIP（AP-1 / CREB / CEBP / MAF / NF-E2 / NF-κB 等）
_TF_BZIP = [
    "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6", "ATF7",
    "BACH1", "BACH2", "BATF", "BATF2", "BATF3",
    "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG",
    "CREB1", "CREB3", "CREB3L1", "CREB3L2", "CREB3L3", "CREB3L4", "CREB5",
    "CREM", "DDIT3", "DBP", "FOS", "FOSB", "FOSL1", "FOSL2", "JDP2",
    "JUN", "JUNB", "JUND", "MAF", "MAFA", "MAFB", "MAFF", "MAFG", "MAFK",
    "NFE2", "NFE2L1", "NFE2L2", "NFE2L3", "NFIL3", "NFKB1", "NFKB2",
    "REL", "RELA", "RELB", "TEF", "XBP1",
]

# IRF / STAT
_TF_IRF_STAT = [
    "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
]

# GATA / KLF / SP / GFI / SREBF
_TF_GATA_KLF = [
    "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6", "GFI1", "GFI1B",
    "KLF1", "KLF2", "KLF3", "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9",
    "KLF10", "KLF11", "KLF12", "KLF13", "KLF14", "KLF15", "KLF16", "KLF17",
    "SP1", "SP2", "SP3", "SP4", "SP5", "SP6", "SP7", "SP8", "SP9",
    "SREBF1", "SREBF2", "ZFPM1", "ZFPM2",
]

# 造血 / 淋系 / 免疫
_TF_HEME = [
    "BCL11A", "BCL11B", "BCL3", "BCL6", "BHLHE40", "BHLHE41",
    "EBF1", "EBF3", "IKZF1", "IKZF2", "IKZF3", "IKZF4", "IKZF5",
    "LYL1", "MYB", "MYBL1", "MYBL2", "TAL1", "TAL2",
    "TCF3", "TCF4", "TCF12", "TCF21", "TCF15",
    "HES1", "HES2", "HES4", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL",
    "ID1", "ID2", "ID3", "ID4", "RUNX1", "RUNX2", "RUNX3",
]

# 细胞周期 / MYC / MITF / USF
_TF_CELLCYCLE = [
    "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8",
    "TFDP1", "TFDP2", "MYC", "MYCN", "MXD1", "MXD3", "MXD4", "MAX", "MNT",
    "MLX", "MLXIP", "MLXIPL", "TFE3", "TFEB", "TFEC", "MITF",
    "USF1", "USF2", "USF3", "TP53", "TP63", "TP73", "WT1",
    "EGR1", "EGR2", "EGR3", "EGR4",
]

# 其它经典转录因子（发育、应激、代谢、表观相关）
_TF_OTHER = [
    "AHR", "ARNT", "ARNTL", "ARNTL2", "ASCL1", "ASCL2", "ASCL3",
    "ATOH1", "ATOH7", "ATOH8", "BNC1", "BNC2", "CLOCK", "CTCF", "DEAF1",
    "DMRT1", "DMRT2", "DMRT3", "E4F1", "GMEB1", "GMEB2",
    "GRHL1", "GRHL2", "GRHL3", "HAND1", "HAND2", "HIC1", "HIC2",
    "HIF1A", "HIF3A", "HLF", "HNF1A", "HNF1B", "HNF4A", "HNF4G",
    "HSF1", "HSF2", "HSF4", "IRX1", "IRX2", "IRX3", "IRX4", "IRX5", "IRX6",
    "MAZ", "MECOM", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "NANOG",
    "NFAT5", "NFATC1", "NFATC2", "NFATC3", "NFATC4",
    "NFIA", "NFIB", "NFIC", "NFIX", "NPAS1", "NPAS2", "NPAS3", "NPAS4",
    "NRF1", "OSR1", "OSR2", "PATZ1", "PKNOX1", "PKNOX2",
    "PLAG1", "PLAGL1", "PLAGL2",
    "PRDM1", "PRDM2", "PRDM4", "PRDM5", "PRDM6", "PRDM8", "PRDM10",
    "PRDM14", "PRDM16", "PREB", "RBPJ", "REST",
    "RFX1", "RFX2", "RFX3", "RFX4", "RFX5", "RFX6", "RFX7", "RFX8",
    "SALL1", "SALL2", "SALL3", "SALL4", "SATB1", "SATB2",
    "SNAI1", "SNAI2", "SNAI3", "SRF", "TBR1",
    "TBX1", "TBX2", "TBX3", "TBX4", "TBX5", "TBX6", "TBX15", "TBX18",
    "TBX19", "TBX20", "TBX21", "TBX22",
    "TEAD1", "TEAD2", "TEAD3", "TEAD4",
    "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E", "TFAP4", "TFCP2",
    "TGIF1", "TGIF2", "THAP1", "THAP11", "TRPS1",
    "TSHZ1", "TSHZ2", "TSHZ3", "VEZF1", "YY1", "YY2",
    "ZBTB1", "ZBTB7A", "ZBTB7B", "ZBTB16", "ZBTB17", "ZBTB18", "ZBTB20",
    "ZBTB24", "ZBTB32", "ZBTB33", "ZBTB46", "ZFX", "ZFY",
    "ZIC1", "ZIC2", "ZIC3", "ZIC4", "ZIC5",
    "ZNF143", "ZNF148", "ZNF217", "ZNF263", "ZNF281", "ZNF331", "ZNF384",
    "ZNF410", "ZNF423", "ZNF462", "ZNF483", "ZNF521", "ZNF638", "ZNF687",
    "ZNF750", "ZNF76", "ZNF800",
]

#: 全部人源 TF 符号（去重、保持顺序）
TF_LIST = []
_seen_tf = set()
for _gene in (
    _TF_NUCLEAR
    + _TF_HOMEOBOX
    + _TF_FOX
    + _TF_HMG
    + _TF_ETS
    + _TF_BZIP
    + _TF_IRF_STAT
    + _TF_GATA_KLF
    + _TF_HEME
    + _TF_CELLCYCLE
    + _TF_OTHER
):
    if _gene not in _seen_tf:
        _seen_tf.add(_gene)
        TF_LIST.append(_gene)
del _seen_tf

#: TF 家族映射（只对已显式分家族的前 9 组给出；其它归入 "other"）
TF_FAMILY = {}
for _family, _genes in (
    ("nuclear receptor", _TF_NUCLEAR),
    ("homeobox", _TF_HOMEOBOX),
    ("forkhead", _TF_FOX),
    ("HMG/SOX/TCF", _TF_HMG),
    ("ETS", _TF_ETS),
    ("bZIP", _TF_BZIP),
    ("IRF/STAT", _TF_IRF_STAT),
    ("GATA/KLF/SP", _TF_GATA_KLF),
    ("hematopoietic", _TF_HEME),
    ("cell cycle/MYC", _TF_CELLCYCLE),
):
    for _gene in _genes:
        TF_FAMILY.setdefault(_gene, _family)


# ---------------------------------------------------------------------------
# 3. 便捷访问器（保持纯函数，方便测试）
# ---------------------------------------------------------------------------


def lr_pairs():
    """返回配体–受体对副本：list[(ligand, receptor, family)]。"""
    return [(str(ligand), str(receptor), str(family)) for ligand, receptor, family in LR_PAIRS]


def tf_list():
    """返回人源 TF 符号列表副本。"""
    return list(TF_LIST)


def lr_genes():
    """返回 LR 对涉及的全部基因（配体在前、受体在后，去重保序）。"""
    genes = []
    seen = set()
    for ligand, receptor, _family in LR_PAIRS:
        for gene in (ligand, receptor):
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return genes


def tf_family(gene):
    """返回 TF 家族名；未知基因返回 None。"""
    return TF_FAMILY.get(str(gene))


__all__ = [
    "LR_PAIRS",
    "LR_FAMILIES",
    "TF_LIST",
    "TF_FAMILY",
    "lr_pairs",
    "tf_list",
    "lr_genes",
    "tf_family",
]
