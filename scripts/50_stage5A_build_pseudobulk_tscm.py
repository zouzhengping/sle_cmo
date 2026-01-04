#!/usr/bin/env python3
"""
Stage 5A: build donor-level pseudo-bulk counts for candidate TSCM clusters

输入:
  data/qc/GSE137029_sle.tcells.clustered.tscm.h5ad
    - 要求 adata.obs 至少包含列:
        * cluster_k20   : KMeans 聚类编号 (字符串, 如 "17")
        * donor_id      : 供者 ID (如 "D001", "D002"...)
        * condition     : 组别 (如 "SLE", "Healthy")

输出:
  data/pseudobulk/tscm_candidate/cluster{CLUSTER}_counts.tsv
      - 行: 基因 (adata.var_names)
      - 列: 每个 donor 一列 (e.g. SLE_D001, HC_D005 ...)
  data/pseudobulk/tscm_candidate/metadata.tsv
      - 每一行对应一个 pseudo-bulk 样本 (donor+cluster)
      - 包含: cluster_k20, donor_id, condition, n_cells, counts_file
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm.h5ad"
OUT_DIR = ROOT / "data/pseudobulk" / "tscm_candidate"


def main():
    parser = argparse.ArgumentParser(
        description="Stage 5A: build donor-level pseudo-bulk for candidate TSCM clusters"
    )
    parser.add_argument(
        "--clusters",
        nargs="+",
        default=["17", "13", "5"],
        help="候选 TSCM cluster_k20 列表 (默认: 17 13 5)",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        default=30,
        help="每个 donor 在某个 cluster 中至少需要的细胞数 (默认: 30)",
    )
    args = parser.parse_args()

    clusters = [str(c) for c in args.clusters]

    print("[INFO] 读取 TSCM 打分后的 T 细胞对象:", IN_H5AD)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[INFO] n_cells={adata.n_obs}, n_genes={adata.n_vars}")

    # 基本检查
    for col in ["cluster_k20", "donor_id", "condition"]:
        if col not in adata.obs.columns:
            raise SystemExit(
                f"[ERROR] obs 中缺少必要列: {col}，请先在上游步骤中补充 donor/condition 信息。"
            )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    gene_names = adata.var_names.to_list()

    meta_records = []

    for cl in clusters:
        print(f"[INFO] 处理 cluster_k20 = {cl} ...")
        mask_cl = adata.obs["cluster_k20"] == cl
        n_cells_cl = int(mask_cl.sum())
        print(f"[INFO]   该 cluster 细胞数: {n_cells_cl}")

        if n_cells_cl == 0:
            print(f"[WARN]   cluster {cl} 没有任何细胞，跳过。")
            continue

        ad_cl = adata[mask_cl].copy()

        # 只保留 SLE / Healthy 两种
        valid_cond = ad_cl.obs["condition"].isin(["SLE", "Healthy"])
        ad_cl = ad_cl[valid_cond].copy()
        if ad_cl.n_obs == 0:
            print(f"[WARN]   cluster {cl} 中没有 SLE/Healthy 标注的细胞，跳过。")
            continue

        # 使用 raw counts (如果存在)，否则用当前 X
        if ad_cl.raw is not None:
            print("[INFO]   使用 raw.X 作为计数矩阵进行 pseudo-bulk。")
            X_all = ad_cl.raw.X
            var_names = ad_cl.raw.var_names
        else:
            print("[WARN]   未检测到 .raw，使用当前 X 作为 counts (可能是归一化后的)。")
            X_all = ad_cl.X
            var_names = ad_cl.var_names

        # 确保基因顺序一致
        if not np.array_equal(var_names, ad_cl.var_names):
            # 一般不会走到这里，简单对齐一下
            print("[WARN]   raw.var_names 与 var.names 不一致，尝试对齐。")
        # 统一转为列表
        genes = var_names.to_list()

        # 按 donor_id + condition 聚合
        group = ad_cl.obs.groupby(["condition", "donor_id"])

        # 为当前 cluster 构建一个 DataFrame: genes x pseudo-bulk-samples
        cols = []
        cols_arrays = []

        for (cond, donor), idx in group.indices.items():
            idx = np.array(idx)
            n_cells_donor = idx.size
            if n_cells_donor < args.min_cells:
                print(
                    f"[INFO]   跳过 {cl} - {cond} - {donor}: "
                    f"细胞数 {n_cells_donor} < min_cells={args.min_cells}"
                )
                continue

            label = f"{cond}_{donor}"
            print(
                f"[INFO]   构建 pseudo-bulk: cluster={cl}, label={label}, n_cells={n_cells_donor}"
            )

            # 取该 donor 在该 cluster 中的细胞
            if sparse.issparse(X_all):
                sub = X_all[idx, :]
                summed = np.asarray(sub.sum(axis=0)).ravel()
            else:
                sub = X_all[idx, :]
                summed = sub.sum(axis=0)

            cols.append(label)
            cols_arrays.append(summed.astype(float))

            # 记录 metadata
            meta_records.append(
                {
                    "cluster_k20": cl,
                    "condition": cond,
                    "donor_id": donor,
                    "label": label,
                    "n_cells": int(n_cells_donor),
                    "counts_file": f"cluster{cl}_counts.tsv",
                }
            )

        if not cols:
            print(f"[WARN]   cluster {cl} 没有任何 donor 满足 min_cells 要求，跳过输出。")
            continue

        # 拼成 genes x samples 的表
        mat = np.vstack(cols_arrays).T  # shape: n_genes x n_samples
        df = pd.DataFrame(mat, index=genes, columns=cols)

        out_counts = OUT_DIR / f"cluster{cl}_counts.tsv"
        df.to_csv(out_counts, sep="\t")
        print("[INFO]   写出 counts 矩阵:", out_counts)

    # 写 metadata
    if meta_records:
        df_meta = pd.DataFrame(meta_records)
        df_meta.to_csv(OUT_DIR / "metadata.tsv", sep="\t", index=False)
        print("[INFO] 写出 pseudo-bulk metadata:", OUT_DIR / "metadata.tsv")
    else:
        print("[WARN] 没有任何 pseudo-bulk 样本被生成，请检查 donor_id/condition/min_cells 设置。")

    print("[DONE] Stage 5A pseudo-bulk 构建完成。")


if __name__ == "__main__":
    main()
