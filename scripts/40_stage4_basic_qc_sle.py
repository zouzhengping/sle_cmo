#!/usr/bin/env python3
"""
Stage 4: Basic QC for GSE137029 (sle cohort)

- 读取 10x MatrixMarket 格式的计数矩阵（mtx + barcodes + genes）
- 构建 AnnData 对象（cells x genes）
- 计算基础 QC 指标（每个细胞的 gene 数、UMI 数、线粒体比例等）
- 保存为 h5ad，供后续 T/TSCM 分析使用
"""

import scanpy as sc
import scipy.io
import scipy.sparse
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MTX = ROOT / "data/raw/GSE137029/sle/GSE137029_sle.mtx.gz"
BC  = ROOT / "data/raw/GSE137029/sle/GSE137029_sle.barcodes.txt.gz"
GN  = ROOT / "data/raw/GSE137029/sle/GSE137029_sle.genes.txt.gz"

OUT_DIR = ROOT / "data/qc"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_H5AD = OUT_DIR / "GSE137029_sle.raw_qc.h5ad"

print("[INFO] ROOT:", ROOT)
print("[INFO] MTX :", MTX)
print("[INFO] BC  :", BC)
print("[INFO] GN  :", GN)

if not MTX.exists() or not BC.exists() or not GN.exists():
    raise SystemExit("[ERROR] mtx / barcodes / genes 有缺失，请检查路径。")

print("[INFO] 读取稀疏矩阵 (MatrixMarket, genes x cells)...")
X = scipy.io.mmread(str(MTX))  # 形状: genes x cells

print("[INFO] 读取 barcodes（细胞）...")
barcodes = pd.read_csv(BC, header=None)[0].astype(str).tolist()

print("[INFO] 读取 genes（特征）...")
genes_df = pd.read_csv(GN, header=None, sep="\t")
# 兼容 1 列或 2 列的基因文件
if genes_df.shape[1] >= 2:
    gene_ids = genes_df.iloc[:, 0].astype(str).tolist()
    gene_names = genes_df.iloc[:, 1].astype(str).tolist()
else:
    gene_ids = genes_df.iloc[:, 0].astype(str).tolist()
    gene_names = gene_ids

print(f"[INFO] matrix shape (原始): {X.shape}")
print(f"[INFO] #genes in file : {len(gene_ids)}")
print(f"[INFO] #barcodes      : {len(barcodes)}")

# 判断矩阵格式：如果行数等于 barcodes 数，则为 cells x genes；如果行数等于 genes 数，则为 genes x cells
if X.shape[0] == len(barcodes) and X.shape[1] == len(gene_ids):
    print("[INFO] 检测到矩阵格式为 cells x genes，无需转置")
    X_final = X
elif X.shape[0] == len(gene_ids) and X.shape[1] == len(barcodes):
    print("[INFO] 检测到矩阵格式为 genes x cells，需要转置")
    X_final = X.T
else:
    print(f"[WARN] 矩阵维度 ({X.shape[0]}, {X.shape[1]}) 与预期不匹配")
    print(f"[WARN] 基因数: {len(gene_ids)}, 细胞数: {len(barcodes)}")
    # 尝试自动判断：如果行数更接近 barcodes 数，假设为 cells x genes
    if abs(X.shape[0] - len(barcodes)) < abs(X.shape[0] - len(gene_ids)):
        print("[INFO] 自动判断为 cells x genes 格式")
        X_final = X
    else:
        print("[INFO] 自动判断为 genes x cells 格式，进行转置")
        X_final = X.T

print("[INFO] 构建 AnnData（cells x genes）对象...")
# 确保矩阵是 CSR 格式（AnnData 保存 h5ad 时需要）
if scipy.sparse.issparse(X_final):
    if X_final.format != 'csr':
        print(f"[INFO] 转换稀疏矩阵格式: {X_final.format} -> csr")
        X_final = X_final.tocsr()
else:
    # 如果是密集矩阵，也转换为 CSR 以节省空间
    print("[INFO] 将密集矩阵转换为稀疏矩阵 (CSR)")
    X_final = scipy.sparse.csr_matrix(X_final)

adata = sc.AnnData(X=X_final)  # 使用正确格式的矩阵
adata.obs_names = barcodes
adata.var_names = gene_names
adata.var["gene_id"] = gene_ids

# 线粒体基因标记（MT- 开头）
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

print("[INFO] 计算基础 QC 指标（n_genes_by_counts, total_counts, pct_counts_mt）...")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)

print("[INFO] AnnData 维度: n_cells={}, n_genes={}".format(adata.n_obs, adata.n_vars))
print("[INFO] 示例 QC 指标（前 5 个细胞）:")
print(
    adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]]
    .head()
    .to_string()
)

print("[INFO] 保存 h5ad 到:", OUT_H5AD)
adata.write_h5ad(OUT_H5AD)

print("[DONE] Stage 4 basic QC 完成。")
