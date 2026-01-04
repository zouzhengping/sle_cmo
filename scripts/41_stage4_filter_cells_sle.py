#!/usr/bin/env python3
"""
Stage 4B: Filter cells for GSE137029 (sle cohort)

输入:
  data/qc/GSE137029_sle.raw_qc.h5ad

输出:
  data/qc/GSE137029_sle.qc_filtered.h5ad
  data/qc/GSE137029_sle.qc_summary.tsv

过滤逻辑(可根据需要调整):
  - 保留 n_genes_by_counts 在 [200, 5000] 之间
  - 保留 pct_counts_mt < 15
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.raw_qc.h5ad"
OUT_H5AD = ROOT / "data/qc" / "GSE137029_sle.qc_filtered.h5ad"
OUT_QC   = ROOT / "data/qc" / "GSE137029_sle.qc_summary.tsv"

print("[INFO] ROOT:", ROOT)
print("[INFO] 读取原始 QC 对象:", IN_H5AD)

adata = sc.read_h5ad(IN_H5AD)
n_cells_before = adata.n_obs
n_genes = adata.n_vars
print(f"[INFO] 原始: n_cells={n_cells_before}, n_genes={n_genes}")

# 过滤条件
min_genes = 200
max_genes = 5000
max_pct_mt = 15.0

print("[INFO] 应用过滤阈值:")
print(f"       n_genes_by_counts >= {min_genes}")
print(f"       n_genes_by_counts <= {max_genes}")
print(f"       pct_counts_mt < {max_pct_mt}")

keep_mask = (
    (adata.obs["n_genes_by_counts"] >= min_genes) &
    (adata.obs["n_genes_by_counts"] <= max_genes) &
    (adata.obs["pct_counts_mt"] < max_pct_mt)
)

n_keep = int(keep_mask.sum())
print(f"[INFO] 通过过滤的细胞数: {n_keep} (占比 {n_keep / n_cells_before * 100:.2f}%)")

adata_filtered = adata[keep_mask].copy()

print("[INFO] 过滤后 AnnData 维度: n_cells={}, n_genes={}".format(
    adata_filtered.n_obs, adata_filtered.n_vars
))

# 保存过滤后的对象
adata_filtered.write_h5ad(OUT_H5AD)
print("[INFO] 已保存过滤后对象:", OUT_H5AD)

# 输出简单 QC 汇总
summary = {
    "n_cells_before": [n_cells_before],
    "n_cells_after": [adata_filtered.n_obs],
    "n_genes": [adata_filtered.n_vars],
    "min_genes_threshold": [min_genes],
    "max_genes_threshold": [max_genes],
    "max_pct_mt_threshold": [max_pct_mt],
}

df_sum = pd.DataFrame(summary)
df_sum.to_csv(OUT_QC, sep="\t", index=False)
print("[INFO] 已写入 QC 汇总:", OUT_QC)

print("[DONE] Stage 4B filtering 完成。")
