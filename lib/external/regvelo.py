#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""官方 RegVelo 适配器（CONTRACT §10.3）—— 可用才接，不可用给明确原因。

本机实测结论（2026-11-09，/opt/miniconda3 Python 3.13.11 / macOS arm64）：
  regvelo 0.4.2 硬性依赖 ``torch<2.6.0`` 与 ``scvi-tools<1.2.1``，而：
    * torch < 2.6 在 PyPI 上没有 macOS arm64 + CPython 3.13 的二进制 wheel
      （torch 2.6.0 才提供 cp313 macOS arm64 wheel）；
    * scvi-tools < 1.2.1 不支持 Python 3.13。
  因此 ``pip install regvelo`` 在本环境 **无法安装**（已用
  ``pip install --dry-run --only-binary=:all: regvelo`` 复现：
  "regvelo 0.4.2 depends on torch<2.6.0 ... no matching distribution"）。

所以本适配器：
  * 若 ``import regvelo`` 失败 → 退出码 3，stderr 给上面这段中文原因；
  * 若 RegVelo 已装但 bundle 缺 spliced/unspliced 层 → 退出码 3（RNA velocity 类工具需要层数据）；
  * 若 RegVelo 已装且层数据齐全 → 仍要求 params.prior（TF regulon 先验，通常是 pySCENIC
    / SCENIC+ 输出）；没有先验时不会假造，退出码 3 说明还缺什么。

官方实现：RegVelo (cistopic/regvelo) 用 scvi-tools 的变分推断联合建模 spliced/unspliced
与 TF regulon 先验，训练成本高（官方示例在 GPU 上跑数十分钟），不适合在 ~5000 细胞上限的
交互式工作台里默认开启。本机未能安装成功，故不提供未经验证的训练代码。
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_EXT_DIR = str(Path(__file__).resolve().parent)
sys.path.insert(0, _EXT_DIR)
import pyext_common as C  # noqa: E402

# 适配器文件名与第三方包同名（pyscenic.py / cellphonedb.py / palantir.py …），
# 运行脚本时 Python 会把本目录放进 sys.path[0]，若不移除，import palantir 会导入
# 适配器自身而不是官方包（自遮蔽）。pyext_common 已导入，随后移除本目录。
while _EXT_DIR in sys.path:
    sys.path.remove(_EXT_DIR)

INSTALL_REASON = (
    "RegVelo 无法安装：regvelo 0.4.2 依赖 torch<2.6.0 与 scvi-tools<1.2.1；"
    "torch<2.6 没有 macOS arm64 + Python 3.13 的 wheel（torch 2.6 才提供 cp313 macOS wheel），"
    "scvi-tools<1.2.1 也不支持 Python 3.13。"
    "如需 RegVelo，请另外创建 Python 3.11 的 conda 环境（conda create -n regvelo python=3.11; "
    "pip install regvelo），再把 externalPython 指向该环境。"
)


def _layer_path(bundle_dir, name):
    root = Path(bundle_dir)
    for candidate in (root / ("%s.mtx" % name), root / "layers" / ("%s.mtx" % name)):
        if candidate.is_file():
            return candidate
    return None


def run(bundle_dir, out, params):
    warnings = []
    notes = []
    try:
        import regvelo  # noqa: F401

        version = C.version_of(regvelo)
    except Exception as error:
        notes.append(INSTALL_REASON)
        raise C.PrerequisiteError(INSTALL_REASON + "（import 失败：%s）" % str(error)[:200])
    notes.append("RegVelo %s" % version)
    missing = []
    if _layer_path(bundle_dir, "spliced") is None:
        missing.append("spliced")
    if _layer_path(bundle_dir, "unspliced") is None:
        missing.append("unspliced")
    if missing:
        raise C.PrerequisiteError(
            "RegVelo 需要 spliced/unspliced 层（当前 bundle 缺少：%s，§10.2 只导出 counts/data）。"
            "请导入 velocyto / STARsolo / kb-python 的 spliced+unspliced 矩阵后再运行。" % "、".join(missing)
        )
    prior = params.get("prior") or params.get("regulon_prior") or params.get("tf_prior")
    if not prior:
        raise C.PrerequisiteError(
            "RegVelo 需要 TF regulon 先验（params.prior：pySCENIC/SCENIC+ 导出的 TF-target 表），"
            "当前未提供。RegVelo 官方训练成本高，交互式工作台建议先用 pySCENIC + scVelo。"
        )
    raise C.PrerequisiteError(
        "RegVelo 已安装但本适配器未启用官方训练流程（未在本机验证过完整训练，拒绝给未经验证的结果）。"
        "如需官方 RegVelo 结果，请在独立环境按官方 notebook 运行。"
    )


def main():
    return C.run_adapter(run)


if __name__ == "__main__":
    sys.exit(main())
