#!/usr/bin/env python3
"""
Stage 4F: 在 SLE T 细胞中标记 TSCM 倾向细胞

输入:
  data/qc/GSE137029_sle.tcells.clustered.tscm.h5ad
    - 要求 obs 里已有:
        * cluster_k20   : KMeans 聚类编号 (字符串或整数)
        * 各种 score_*  : TSCM, naive, effector 等打分 (可选, 已有更好)

输出:
  data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad
    - 新增 obs 列:
        * tscm_label : "TSCM_high" 或 "T_other"
        * is_tscm_high : True/False
"""

from pathlib import Path
import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm.h5ad"
OUT_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm_labeled.h5ad"

# 手动指定 TSCM 倾向最高的一批 cluster
TSCM_CLUSTERS = ["13", "16", "17", "5"]  # 如需调整，改这里即可

print("[INFO] 读取带 TSCM 打分的 T 细胞对象:", IN_H5AD)
adata = sc.read_h5ad(IN_H5AD)
print(f"[INFO] n_cells={adata.n_obs}, n_genes={adata.n_vars}")

if "cluster_k20" not in adata.obs.columns:
    raise SystemExit("[ERROR] obs 中未发现 cluster_k20，请确认 Stage 4D(KMeans 聚类) 已完成。")

# 统一把 cluster_k20 转为字符串，方便比较
clusters_str = adata.obs["cluster_k20"].astype(str)
adata.obs["cluster_k20"] = clusters_str

is_tscm_high = clusters_str.isin(TSCM_CLUSTERS)
adata.obs["is_tscm_high"] = is_tscm_high
adata.obs["tscm_label"] = np.where(is_tscm_high, "TSCM_high", "T_other")

# 简单统计一下
n_tscm = int(is_tscm_high.sum())
n_other = int((~is_tscm_high).sum())
print("[INFO] TSCM_high 细胞数:", n_tscm)
print("[INFO] T_other    细胞数:", n_other)

print("[INFO] TSCM cluster 列表及细胞数:")
for cl in sorted(set(TSCM_CLUSTERS)):
    n_cl = int((clusters_str == cl).sum())
    print(f"  cluster_k20 = {cl}: {n_cl} cells")

# 保存新对象
adata.write_h5ad(OUT_H5AD)
print("[DONE] 已写出带 tscm_label 的对象:", OUT_H5AD)
