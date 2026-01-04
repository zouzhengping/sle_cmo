#!/usr/bin/env python3
"""
Stage 4C — normalization & HVGs for GSE137029 SLE cohort
输入:  data/qc/GSE137029_sle.qc_filtered.h5ad
输出:  data/qc/GSE137029_sle.hvg.h5ad   (包含 hvg 与 normalized counts)
"""

import scanpy as sc
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.qc_filtered.h5ad"
OUT_H5AD = ROOT / "data/qc" / "GSE137029_sle.hvg.h5ad"

print("[INFO] 读取过滤后对象:", IN_H5AD)
adata = sc.read_h5ad(IN_H5AD)

print("[INFO] 归一化 total-counts (1e4) ...")
sc.pp.normalize_total(adata, target_sum=1e4)

print("[INFO] log1p 转换 ...")
sc.pp.log1p(adata)

print("[INFO] 寻找 HVGs (highly variable genes) ...")
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3")

print("[INFO] HVG 数量:", adata.var["highly_variable"].sum())

print("[INFO] 保存 HVG 对象:", OUT_H5AD)
adata.write_h5ad(OUT_H5AD)

print("[DONE] Stage 4C 完成")
