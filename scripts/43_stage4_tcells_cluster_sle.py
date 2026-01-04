#!/usr/bin/env python3
"""
Stage 4D: Clustering & T-cell extraction for GSE137029 (sle cohort)

输入:
  data/qc/GSE137029_sle.hvg.h5ad

输出:
  data/qc/GSE137029_sle.clustered.h5ad      # 全部细胞, 带 PCA/Leiden 等信息
  data/qc/GSE137029_sle.tcells.h5ad         # 初步筛出的 T 细胞子集
"""

import scanpy as sc
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.hvg.h5ad"
OUT_FULL = ROOT / "data/qc" / "GSE137029_sle.clustered.h5ad"
OUT_T    = ROOT / "data/qc" / "GSE137029_sle.tcells.h5ad"

print("[INFO] 读取 HVG 对象:", IN_H5AD)
adata = sc.read_h5ad(IN_H5AD)

print("[INFO] 当前对象维度: n_cells={}, n_genes={}".format(adata.n_obs, adata.n_vars))

# ------------------------------
# 1. PCA + 邻接图 + Leiden 聚类
# ------------------------------
print("[INFO] 进行 PCA (基于 HVG)...")
# 如果前面用 highly_variable 标记了 HVGs, 这里直接用
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")

print("[INFO] 计算邻接图 (neighbors)...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

print("[INFO] 计算 UMAP (可选, 方便后续可视化)...")
sc.tl.umap(adata)

print("[INFO] Leiden 聚类 (resolution=0.5)...")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0_5")

print("[INFO] Leiden 聚类完成，示例前 10 个细胞聚类标签:")
print(adata.obs["leiden_0_5"].head(10).to_string())

# ------------------------------
# 2. T 细胞 marker 打分
# ------------------------------
# 基础 T 细胞 marker (CD3/TCR 等)
tcell_markers = [
    "CD3D", "CD3E", "CD3G",
    "TRAC", "TRBC1", "TRBC2",
    "CD247"
]

print("[INFO] 对 T 细胞 marker 进行打分, markers:")
print("       " + ", ".join(tcell_markers))

sc.tl.score_genes(adata, gene_list=tcell_markers, score_name="tcell_score")

print("[INFO] tcell_score 概览:")
print(adata.obs["tcell_score"].describe())

# 简单阈值: tcell_score > 0.1 视为 T 细胞 (可后续调整)
threshold = 0.1
t_mask = adata.obs["tcell_score"] > threshold
n_tcells = int(t_mask.sum())

print(f"[INFO] 以阈值 tcell_score > {threshold} 筛选 T cells:")
print(f"       T 细胞数: {n_tcells} / {adata.n_obs} "
      f"({n_tcells / adata.n_obs * 100:.2f}%)")

adata_t = adata[t_mask].copy()
print("[INFO] T 细胞子集维度: n_cells={}, n_genes={}".format(
    adata_t.n_obs, adata_t.n_vars
))

# ------------------------------
# 3. 保存结果
# ------------------------------
print("[INFO] 保存带聚类信息的完整对象:", OUT_FULL)
adata.write_h5ad(OUT_FULL)

print("[INFO] 保存 T 细胞子集对象:", OUT_T)
adata_t.write_h5ad(OUT_T)

print("[DONE] Stage 4D (Clustering + T-cell extraction) 完成。")
