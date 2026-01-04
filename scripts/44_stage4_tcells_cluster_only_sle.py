#!/usr/bin/env python3
"""
Stage 4D (part 2): Clustering within T-cell subset for GSE137029 (sle cohort)

输入:
  data/qc/GSE137029_sle.tcells.raw.h5ad

输出:
  data/qc/GSE137029_sle.tcells.clustered.h5ad

说明:
  - 不再依赖 leidenalg，而是使用 sklearn KMeans 进行简单聚类
  - 仍然保留 PCA / neighbors / UMAP，便于后续可视化与 TSCM 打分
"""

import scanpy as sc
from sklearn.cluster import KMeans
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.tcells.raw.h5ad"
OUT_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.h5ad"

print("[INFO] 读取 T 细胞子集对象:", IN_H5AD)
adata_t = sc.read_h5ad(IN_H5AD)
print("[INFO] T 细胞子集维度: n_cells={}, n_genes={}".format(
    adata_t.n_obs, adata_t.n_vars
))

# 若已有 PCA，就直接用；否则重新算
if "X_pca" in adata_t.obsm:
    print("[INFO] 检测到已有 PCA（X_pca），直接复用。")
else:
    use_hvg = "highly_variable" in adata_t.var.columns
    print(f"[INFO] 重新计算 PCA, use_highly_variable={use_hvg}")
    sc.pp.pca(
        adata_t,
        n_comps=50,
        use_highly_variable=use_hvg,
        svd_solver="arpack"
    )

print("[INFO] neighbors ...")
sc.pp.neighbors(adata_t, n_neighbors=15, n_pcs=50)

print("[INFO] UMAP ...（用于可视化）")
sc.tl.umap(adata_t)

# 使用 KMeans 聚类替代 Leiden
n_clusters = 20   # 先给一个中等数量，后续可以调整
print(f"[INFO] 使用 KMeans 聚类, n_clusters={n_clusters}")
X_pca = adata_t.obsm["X_pca"]

km = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
clusters = km.fit_predict(X_pca)

adata_t.obs["cluster_k20"] = clusters.astype(str)

print("[INFO] 聚类标签示例 (前 10 行):")
print(adata_t.obs["cluster_k20"].head(10).to_string())

print("[INFO] 保存聚类后的 T 细胞对象:", OUT_H5AD)
adata_t.write_h5ad(OUT_H5AD)

print("[DONE] Stage 4D part2 (T-cell clustering via KMeans) 完成。")
