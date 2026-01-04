#!/usr/bin/env python3
"""
Stage 4D (part 1): T-cell scoring & subsetting for GSE137029 (sle cohort)

输入:
  data/qc/GSE137029_sle.hvg.h5ad

输出:
  data/qc/GSE137029_sle.tcells.raw.h5ad
  data/qc/GSE137029_sle.tcells_score_summary.tsv
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.hvg.h5ad"
OUT_T    = ROOT / "data/qc" / "GSE137029_sle.tcells.raw.h5ad"
OUT_SUM  = ROOT / "data/qc" / "GSE137029_sle.tcells_score_summary.tsv"

print("[INFO] 读取 HVG 对象:", IN_H5AD)
adata = sc.read_h5ad(IN_H5AD)
print("[INFO] n_cells={}, n_genes={}".format(adata.n_obs, adata.n_vars))

# 基础 T-cell marker
tcell_markers = [
    "CD3D", "CD3E", "CD3G",
    "TRAC", "TRBC1", "TRBC2",
    "CD247"
]

print("[INFO] T 细胞 marker 列表:")
print("       " + ", ".join(tcell_markers))

print("[INFO] 计算 tcell_score ...")
sc.tl.score_genes(adata, gene_list=tcell_markers, score_name="tcell_score")

desc = adata.obs["tcell_score"].describe()
print("[INFO] tcell_score 描述统计:")
print(desc.to_string())

# 简单阈值: 大于 0 视为 T cell（先宽一点）
threshold = 0.0
t_mask = adata.obs["tcell_score"] > threshold
n_tcells = int(t_mask.sum())

print(f"[INFO] 使用阈值 tcell_score > {threshold} 选择 T cells")
print(f"       T 细胞数: {n_tcells} / {adata.n_obs} "
      f"({n_tcells / adata.n_obs * 100:.2f}%)")

# 子集
adata_t = adata[t_mask].copy()
print("[INFO] T 细胞子集维度: n_cells={}, n_genes={}".format(
    adata_t.n_obs, adata_t.n_vars
))

# 保存子集
adata_t.write_h5ad(OUT_T)
print("[INFO] 已写入 T 细胞子集对象:", OUT_T)

# 保存打分汇总
df_sum = desc.to_frame(name="tcell_score")
df_sum.to_csv(OUT_SUM, sep="\t")
print("[INFO] 已写入 tcell_score 汇总:", OUT_SUM)

print("[DONE] Stage 4D part1 (T-cell scoring & subsetting) 完成。")
