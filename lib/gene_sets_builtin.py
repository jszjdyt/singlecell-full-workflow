#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dsh-singlecell-workbench 内置基因集（sc-analysis-a 专用，无第三方依赖）。

内容：
  - ``CELL_CYCLE_S`` / ``CELL_CYCLE_G2M``：Tirosh/Regev 人细胞周期 S 与 G2M 基因集
    （Seurat ``cc.genes`` 常用版本；基因符号为人类 HGNC symbol）。
  - ``BUILTIN_GMT``：30 条常用通路/基因集（GO/KEGG/Hallmark 风格的精选小集），
    用于 ORA / GSEA / ssGSEA / AUCell / WGCNA 等分析的内置基因集来源。

这些都是**精选的小型参考集**（每条 10~60 个经典标志基因），不是 MSigDB 全量数据；
在报告/结果里必须如实标注"内置精选集"。用户可通过 ``gmt_file`` 参数传入自己的
GMT 文件（标准格式：``term<TAB>description<TAB>gene1<TAB>gene2...``）。
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Tirosh / Regev 细胞周期基因集（人类 symbol）
# ---------------------------------------------------------------------------

CELL_CYCLE_S = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
]

CELL_CYCLE_G2M = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


# ---------------------------------------------------------------------------
# 内置小 GMT：30 条通路（GO / KEGG / Hallmark 风格精选标志基因）
# order 即 UI 展示顺序（细胞周期在前，方便演示）
# ---------------------------------------------------------------------------

BUILTIN_GMT = {
    "Cell cycle (KEGG)": [
        "CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CCNA2", "CCNB1", "CCNB2",
        "CCND1", "CCNE1", "CCNE2", "CDC6", "CDC20", "CDC25A", "CDC25B",
        "CDC25C", "CDC45", "E2F1", "E2F2", "E2F3", "E2F4", "MCM2", "MCM3",
        "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "RB1", "RBL1", "RBL2",
        "TP53", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "SKP2", "TFDP1",
        "ORC1", "ORC2", "ORC6", "MAD2L1", "BUB1", "BUB3", "PLK1", "WEE1",
        "AURKA", "AURKB", "MKI67", "TOP2A",
    ],
    "Apoptosis": [
        "BAX", "BAK1", "BCL2", "BCL2L1", "BCL2L11", "BAD", "BID", "BIRC5",
        "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10", "APAF1",
        "CYCS", "DIABLO", "TP53", "TNF", "TNFSF10", "FAS", "FASLG",
        "TNFRSF10A", "TNFRSF10B", "XIAP", "AIFM1", "ENDOG", "PARP1",
        "PMAIP1", "BBC3", "MCL1", "DFFA", "DFFB", "CFLAR",
    ],
    "Interferon alpha response": [
        "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2",
        "IFIT3", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
        "OAS3", "OASL", "STAT1", "STAT2", "IRF7", "IRF9", "GBP1", "GBP2",
        "GBP4", "GBP5", "RSAD2", "XAF1", "BST2", "IFIH1", "DDX58", "HERC5",
        "SAMD9", "SAMD9L", "USP18",
    ],
    "Interferon gamma response": [
        "IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "STAT2", "IRF1",
        "IRF8", "CXCL9", "CXCL10", "CXCL11", "CCL5", "GBP1", "GBP2", "GBP3",
        "GBP4", "GBP5", "GBP6", "IDO1", "HLA-A", "HLA-B", "HLA-C", "HLA-E",
        "HLA-F", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "PSMB10", "SOCS1",
        "SOCS3", "CD74", "CIITA",
    ],
    "Hypoxia": [
        "HIF1A", "ARNT", "EPAS1", "VEGFA", "VEGFB", "LDHA", "LDHB", "PGK1",
        "PDK1", "PDK3", "BNIP3", "BNIP3L", "NDRG1", "DDIT4", "ANGPTL4",
        "SLC2A1", "SLC2A3", "HK2", "PFKFB3", "PFKFB4", "ADM", "MIF",
        "EGLN1", "EGLN2", "EGLN3", "CA9", "P4HA1", "P4HA2", "PGAM1",
        "ENO1", "TPI1", "ALDOA", "GAPDH", "VHL",
    ],
    "Epithelial-mesenchymal transition (EMT)": [
        "VIM", "CDH1", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "TWIST2",
        "ZEB1", "ZEB2", "TGFB1", "TGFBR1", "TGFBR2", "COL1A1", "COL1A2",
        "COL3A1", "COL5A1", "SPARC", "MMP2", "MMP9", "MMP14", "ACTA2",
        "ITGB1", "ITGA5", "CTNNB1", "LEF1", "JUN", "FOS", "SERPINE1",
        "THBS1", "TIMP1", "DDR2", "LUM", "FBN1",
    ],
    "Oxidative phosphorylation (OXPHOS)": [
        "NDUFA1", "NDUFA2", "NDUFA4", "NDUFA8", "NDUFA9", "NDUFB1",
        "NDUFB2", "NDUFB7", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFV1",
        "NDUFV2", "SDHA", "SDHB", "SDHC", "SDHD", "UQCRB", "UQCRC1",
        "UQCRC2", "UQCRFS1", "UQCRQ", "CYC1", "COX4I1", "COX5A", "COX5B",
        "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7B", "COX7C", "COX8A",
        "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5MC1", "ATP5MC2",
        "ATP5MC3", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5ME",
        "ATP5MF", "ATP5MG", "ATP5IF1",
    ],
    "Ribosome": [
        "RPS3", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10",
        "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS16", "RPS17",
        "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25",
        "RPS26", "RPS27", "RPS27A", "RPS28", "RPS29", "RPL3", "RPL4",
        "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10",
        "RPL10A", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15",
        "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23",
        "RPL23A", "RPL24", "RPL26", "RPL27", "RPL28", "RPL29", "RPL30",
        "RPL31", "RPL32", "RPL34", "RPL35", "RPL36", "RPL37", "RPL38",
        "RPL39", "RPL41", "RPLP0", "RPLP1", "RPLP2",
    ],
    "T cell receptor signaling": [
        "CD3D", "CD3E", "CD3G", "CD247", "LCK", "FYN", "ZAP70", "LAT",
        "LCP2", "ITK", "PLCG1", "PRKCQ", "PIK3CD", "PIK3R1", "GRB2", "SOS1",
        "HRAS", "NRAS", "RAF1", "MAP2K1", "MAPK1", "MAPK3", "JUN", "FOS",
        "NFATC1", "NFATC2", "NFKB1", "RELA", "CD28", "ICOS", "VAV1", "RAC1",
        "CDC42", "WAS", "PAK1", "PRKCD",
    ],
    "B cell receptor signaling": [
        "CD79A", "CD79B", "BLNK", "SYK", "BTK", "PLCG2", "PIK3CD", "VAV1",
        "GRB2", "SOS1", "HRAS", "MAP2K1", "MAPK1", "JUN", "NFATC1", "NFKB1",
        "RELA", "CD19", "MS4A1", "CD22", "PTPN6", "LYN", "HCK", "CARD11",
        "BCL10", "MALT1", "PRKCB", "CARD9",
    ],
    "NK cell cytotoxicity": [
        "NKG7", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1",
        "KLRD1", "KLRF1", "KLRK1", "NCR1", "NCR3", "FCGR3A", "NCAM1",
        "SH2D1A", "CD160", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR3DL1",
        "KLRB1", "IFNG", "TNF", "FASLG", "SPON2", "CCL3", "CCL4", "CCL5",
        "CST7", "CTSW", "PRDM1",
    ],
    "Antigen processing and presentation": [
        "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-DRA", "HLA-DRB1",
        "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "B2M", "TAP1",
        "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSMA2", "PSMA3",
        "PSMA4", "PSMA5", "PSMB1", "PSMB2", "CD74", "CIITA", "IFNG",
        "IFNGR1", "JAK1", "JAK2", "STAT1", "LGMN", "CTSB", "CTSD",
    ],
    "Inflammatory response": [
        "IL1A", "IL1B", "IL6", "IL18", "TNF", "CXCL1", "CXCL2", "CXCL3",
        "CXCL5", "CXCL6", "CXCL8", "CCL2", "CCL3", "CCL4", "CCL5", "CCL20",
        "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "IKBKG", "IL1R1",
        "IL6R", "IL6ST", "SOCS3", "PTGS2", "VCAM1", "ICAM1", "SELE",
        "TLR2", "TLR4", "MYD88", "IRAK1", "TRAF6", "NLRP3", "IL1RN",
    ],
    "Complement system": [
        "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C3AR1", "C4A",
        "C4B", "C5", "C5AR1", "C7", "C9", "CFB", "CFD", "CFH", "CFI",
        "CFHR1", "CR1", "CR2", "SERPING1", "CD55", "CD46", "CD59",
        "MASP1", "MASP2", "VTN", "CFP",
    ],
    "IL6-JAK-STAT3 signaling": [
        "IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "JAK3", "TYK2", "STAT1",
        "STAT3", "STAT5A", "STAT5B", "SOCS1", "SOCS3", "PTPN11", "GRB2",
        "SOS1", "HRAS", "MAP2K1", "MAPK1", "PIAS3", "IL10", "IL10RA",
        "IL10RB", "CEBPD", "MYC", "BCL2L1", "MCL1", "VEGFA", "OSMR",
        "IL11RA", "IL13RA1", "IL4R",
    ],
    "PI3K-AKT-mTOR signaling": [
        "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "AKT1", "AKT2",
        "AKT3", "MTOR", "RPTOR", "RICTOR", "TSC1", "TSC2", "RHEB",
        "RPS6KB1", "RPS6KB2", "EIF4EBP1", "PTEN", "PDPK1", "GSK3B", "FOXO1",
        "FOXO3", "RPS6", "RPS6KA1", "RAF1", "MAP2K1", "MAPK1", "IGF1R",
        "INSR", "IRS1", "IRS2", "PDK1", "PDPK1",
    ],
    "MYC targets": [
        "MYC", "MAX", "MAD1L1", "MXD1", "MNT", "TERT", "CCND1", "CDK4",
        "E2F1", "NPM1", "NCL", "RPL3", "RPL5", "RPS6", "ODC1", "LDHA",
        "NME1", "NME2", "PES1", "DKC1", "GNL3", "SLC7A5", "SRM", "TYMS",
        "PA2G4", "EIF4E", "TFDP1", "CDKN2B", "BAX", "TP53", "HSPD1",
        "PTMA", "NPM3", "RRP1",
    ],
    "p53 pathway": [
        "TP53", "MDM2", "MDM4", "CDKN1A", "GADD45A", "GADD45B", "GADD45G",
        "ATM", "ATR", "CHEK1", "CHEK2", "BAX", "BCL2", "BBC3", "PMAIP1",
        "FAS", "E2F1", "RB1", "CDKN2A", "SERPINE1", "PTGS2", "TNFRSF10B",
        "ZMAT3", "TIGAR", "CCNG1", "RPS27L", "EEF2K", "PLK3", "SFN",
        "TRIAP1", "DDB2", "XPC", "RRM2B", "TSC22D1",
    ],
    "DNA repair": [
        "RAD51", "RAD51B", "RAD51C", "RAD51D", "BRCA1", "BRCA2", "PALB2",
        "XRCC5", "XRCC6", "PRKDC", "LIG4", "XRCC1", "LIG1", "LIG3", "MSH2",
        "MSH3", "MSH6", "MLH1", "PMS2", "PCNA", "RFC1", "RFC2", "FEN1",
        "EXO1", "RPA1", "RPA2", "RPA3", "ATM", "ATR", "CHEK1", "CHEK2",
        "PARP1", "OGG1", "APEX1", "XPA", "XPC", "ERCC1", "ERCC2", "POLH",
        "POLI", "POLK", "REV3L",
    ],
    "Glycolysis": [
        "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP", "PFKFB1",
        "PFKFB2", "PFKFB3", "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH",
        "PGK1", "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "LDHA",
        "LDHB", "SLC2A1", "SLC2A3", "PDK1", "PDK2", "PDK3", "PDK4", "PGM1",
        "PGM2", "FBP1", "FBP2", "G6PD", "HIF1A",
    ],
    "Fatty acid oxidation": [
        "CPT1A", "CPT1B", "CPT1C", "CPT2", "SLC25A20", "ACADM", "ACADS",
        "ACADL", "ACADVL", "HADHA", "HADHB", "ECHS1", "ECI1", "ECI2",
        "ACAA1", "ACAA2", "ACAT1", "ACAT2", "HMGCS1", "HMGCS2", "HMGCR",
        "ACACA", "ACACB", "FASN", "SCD", "ESRRA", "PPARA", "PPARD",
        "PPARGC1A", "PPARGC1B", "UCP2", "UCP3", "CD36", "LPL", "FABP3",
        "FABP4", "FABP5", "SLC27A1",
    ],
    "Cholesterol homeostasis": [
        "HMGCR", "HMGCS1", "HMGCS2", "INSIG1", "INSIG2", "SREBF1", "SREBF2",
        "LDLR", "LRP1", "APOE", "APOC1", "APOC2", "APOC3", "APOA1", "APOA2",
        "ABCA1", "ABCG1", "NPC1", "NPC2", "SOAT1", "SCARB1", "PCSK9",
        "IDI1", "IDI2", "MVK", "MVD", "SQLE", "LSS", "DHCR7", "DHCR24",
        "CYP51A1", "MSMO1", "NSDHL", "EBP", "FDFT1", "LDLRAP1", "ANGPTL3",
    ],
    "Unfolded protein response": [
        "HSPA5", "HSPA8", "HSP90B1", "HSPA1A", "HSPA1B", "HSPH1", "ATF4",
        "ATF6", "ATF6B", "ERN1", "EIF2AK3", "XBP1", "DDIT3", "CEBPB",
        "CREB3L3", "HERPUD1", "SEL1L", "EDEM1", "EDEM2", "EDEM3", "DNAJB9",
        "DNAJC3", "SEC61A1", "SEC61B", "SEC62", "SEC63", "CALR", "CANX",
        "PDIA3", "PDIA4", "PDIA6", "DERL1", "DERL2", "SDF2L1", "SELENOS",
        "VAPB", "MBTPS1", "MBTPS2",
    ],
    "TGF-beta signaling": [
        "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMAD1",
        "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9",
        "ACVR1", "ACVR2A", "ACVR2B", "BMP2", "BMP4", "BMP6", "BMPR1A",
        "BMPR1B", "BMPR2", "ID1", "ID2", "ID3", "ID4", "JUN", "JUNB", "FOS",
        "SKIL", "SERPINE1", "CCN2", "COL1A1", "COL1A2", "SMURF1", "SMURF2",
        "INHBA", "INHBB", "LTBP1",
    ],
    "WNT signaling": [
        "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6",
        "WNT7A", "WNT7B", "WNT10A", "WNT10B", "WNT11", "FZD1", "FZD2",
        "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10",
        "LRP5", "LRP6", "DVL1", "DVL2", "DVL3", "APC", "AXIN1", "AXIN2",
        "CTNNB1", "TCF7", "TCF7L1", "TCF7L2", "LEF1", "GSK3B", "CSNK1A1",
        "BTRC", "SFRP1", "SFRP2", "DKK1", "DKK2", "DKK3", "MYC", "CCND1",
        "LGR5", "NKD1", "RNF43", "ZNRF3",
    ],
    "NOTCH signaling": [
        "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4",
        "JAG1", "JAG2", "HES1", "HES2", "HES4", "HES5", "HES6", "HES7",
        "HEY1", "HEY2", "HEYL", "LFNG", "RFNG", "MFNG", "PSEN1", "PSEN2",
        "NCSTN", "APH1A", "APH1B", "PSENEN", "MAML1", "MAML2", "MAML3",
        "RBPJ", "DTX1", "DTX2", "DTX3", "DTX4", "TLE1", "TLE2", "TLE3",
        "TLE4", "NUMB", "NUMBL",
    ],
    "MAPK signaling": [
        "HRAS", "KRAS", "NRAS", "RAF1", "BRAF", "MAP2K1", "MAP2K2",
        "MAP2K3", "MAP2K4", "MAP2K6", "MAPK1", "MAPK3", "MAPK4", "MAPK6",
        "MAPK7", "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13",
        "MAPK14", "MAPKAPK2", "MAPKAPK3", "ELK1", "JUN", "FOS", "MYC",
        "RPS6KA1", "RPS6KA2", "RPS6KA3", "MKNK1", "MKNK2", "DUSP1",
        "DUSP2", "DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
        "PTPN11", "SOS1", "SOS2", "GRB2", "SHC1", "SHC2",
    ],
    "Chemokine signaling": [
        "CCL2", "CCL3", "CCL3L1", "CCL4", "CCL4L2", "CCL5", "CCL7", "CCL8",
        "CCL13", "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22",
        "CCL24", "CCL26", "CCL28", "CXCL1", "CXCL2", "CXCL3", "CXCL5",
        "CXCL6", "CXCL8", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13",
        "CXCL16", "CXCL17", "CCR2", "CCR5", "CCR7", "CCR10", "CXCR3",
        "CXCR4", "CXCR5", "CXCR6", "ACKR1", "ACKR2", "ACKR3",
    ],
    "Platelet activation and coagulation": [
        "PF4", "PF4V1", "PPBP", "GNG11", "SPARC", "GPX1", "NRGN", "SELP",
        "PECAM1", "THBS1", "ITGA2B", "ITGB3", "GP1BA", "GP1BB", "GP6",
        "GP9", "VWF", "F2R", "F2", "FGA", "FGB", "FGG", "SERPINC1",
        "SERPIND1", "PROS1", "PROC", "ANXA5", "TMSB4X", "CLU", "CD9",
        "PEAR1", "MYL9", "TLN1", "VCL",
    ],
    "E2F targets": [
        "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8",
        "RB1", "RBL1", "RBL2", "TFDP1", "TFDP2", "MCM2", "MCM3", "MCM4",
        "MCM5", "MCM6", "MCM7", "CDC6", "CDT1", "PCNA", "TYMS", "DHFR",
        "CCNA2", "CCNE1", "CCNE2", "CDK2", "CDK4", "MKI67", "TK1", "SKP2",
        "GMNN", "POLE", "POLD1", "POLA1", "RRM1", "RRM2", "UNG", "FEN1",
        "CHAF1A", "CHAF1B", "NASP", "SLBP", "WDR76",
    ],
    "Heme and erythroid signature": [
        "HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBE1", "HBM",
        "HBQ1", "ALAS2", "AHSP", "ALAD", "ALAS1", "GATA1", "KLF1", "CA1",
        "CA2", "SLC4A1", "EPB42", "ANK1", "SPTA1", "SPTB", "EPOR", "TFRC",
        "GYPA", "GYPB", "GYPE", "RHAG", "RHD", "RHCE", "KLF3", "NFE2",
        "TAL1", "MAP2K3",
    ],
    "Plasma cell and immunoglobulin": [
        "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGLC7", "IGHG1",
        "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGHD",
        "IGHE", "JCHAIN", "MZB1", "SDC1", "XBP1", "PRDM1", "BHLHA15",
        "IRF4", "SLAMF7", "CD38", "MS4A1", "CD19", "CD27", "TNFRSF17",
        "MZB1", "SSR4", "SEC11C", "DERL3", "CRELD2", "FCRL5", "ITGB7",
    ],
    "T cell exhaustion": [
        "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "TOX2",
        "ENTPD1", "CD160", "CD244", "BTLA", "VSIR", "LAYN", "CXCL13",
        "TCF7", "SELL", "IL7R", "CCR7", "GZMK", "BATF", "IRF4", "PRDM1",
        "EOMES", "TBX21", "NR4A1", "NR4A2", "NR4A3", "CD69", "ITGAE",
        "CD39", "ENTPD1",
    ],
    "Myeloid activation": [
        "LYZ", "CD14", "FCGR3A", "FCN1", "S100A8", "S100A9", "S100A12",
        "VCAN", "CD68", "CD163", "MRC1", "MSR1", "MARCO", "TLR4", "TLR2",
        "MYD88", "NFKB1", "RELA", "IL1B", "IL6", "TNF", "CXCL8", "CCL2",
        "CCL3", "CCL4", "CCL5", "FCGR1A", "FCGR2A", "ITGAM", "ITGAX",
        "CSF1R", "CSF2RA", "CSF3R", "SPI1", "CEBPB", "CEBPA", "IRF8",
        "LILRB1", "LILRB2", "AIF1", "TYROBP", "CTSS", "CTSB", "CTSD",
    ],
    "Dendritic cell signature": [
        "CD1C", "CD1E", "CLEC9A", "CLEC10A", "CLEC4C", "FCER1A", "CD83",
        "CD86", "CD80", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
        "ITGAX", "ITGAM", "IRF4", "IRF8", "BATF3", "FLT3", "ZBTB46",
        "CCR7", "LAMP3", "IDO1", "CCL17", "CCL19", "CCL22", "XCR1",
        "CADM1", "LILRA4", "IL3RA", "GZMB", "TCF4", "SPIB",
    ],
    "G2M checkpoint": [
        "CDK1", "CCNB1", "CCNB2", "CCNA2", "CDC20", "CDC25A", "CDC25B",
        "CDC25C", "PLK1", "AURKA", "AURKB", "BUB1", "BUB1B", "BUB3",
        "MAD2L1", "MAD2L2", "NDC80", "NUF2", "TPX2", "TOP2A", "MKI67",
        "KIF11", "KIF20B", "KIF23", "KIF2C", "CENPE", "CENPF", "CENPA",
        "ANLN", "ECT2", "RACGAP1", "TACC3", "CKAP2", "CKAP2L", "GTSE1",
        "UBE2C", "BIRC5", "HMMR", "NEK2", "TTK", "PSRC1", "HJURP",
        "CDCA3", "CDCA8", "DLGAP5", "LBR", "TMPO", "SMC4", "RANGAP1",
        "NCAPD2", "ANP32E", "CBX5",
    ],
}

# 兼容别名：部分调用方使用小写命名
PATHWAY_GMT = BUILTIN_GMT
GENE_SETS = BUILTIN_GMT


def all_gene_sets():
    """返回内置基因集的浅拷贝；键为通路名，值为基因 symbol 列表。"""
    return {name: list(genes) for name, genes in BUILTIN_GMT.items()}


def cell_cycle_sets():
    """返回 (S 基因集, G2M 基因集)。"""
    return list(CELL_CYCLE_S), list(CELL_CYCLE_G2M)


__all__ = [
    "CELL_CYCLE_S",
    "CELL_CYCLE_G2M",
    "BUILTIN_GMT",
    "PATHWAY_GMT",
    "GENE_SETS",
    "all_gene_sets",
    "cell_cycle_sets",
]
