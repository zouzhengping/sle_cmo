#!/usr/bin/env python3
"""
Stage 4G: 在 SLE T 细胞中进行差异分析 (TSCM_high vs T_other)

输入:
  data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad
    - 需要已经经过:
        * 归一化 + log1p
        * HVG / PCA / neighbors (不一定必须)
    - obs 至少包含:
        * tscm_label: "TSCM_high" or "T_other"

输出:
  results/tables/tscm_vs_other_de.full.tsv      # 所有基因的 DE 结果
  results/tables/tscm_vs_other_de.top50.tsv     # logFC 排前 50 个上调基因
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm_labeled.h5ad"
OUT_DIR = ROOT / "results" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("[INFO] 读取带 tscm_label 的 T 细胞对象:", IN_H5AD)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[INFO] n_cells={adata.n_obs}, n_genes={adata.n_vars}")

    if "tscm_label" not in adata.obs.columns:
        raise SystemExit("[ERROR] obs 中未发现 tscm_label，请先运行 Stage 4F (46 脚本)。")

    # 确认两组大小
    grp_counts = adata.obs["tscm_label"].value_counts()
    print("[INFO] tscm_label 细胞数:")
    print(grp_counts)

    # 每组最多保留 N 个细胞，减少计算压力
    MAX_CELLS_PER_GROUP = 50000  # 可按需要调整
    use_idx = []

    for label in ["TSCM_high", "T_other"]:
        mask = adata.obs["tscm_label"] == label
        idx = np.where(mask)[0]
        n = idx.size
        if n == 0:
            print(f"[WARN] 组 {label} 没有细胞，将导致 DE 失败。")
            continue
        if n > MAX_CELLS_PER_GROUP:
            print(f"[INFO] 组 {label} 细胞数 {n} > {MAX_CELLS_PER_GROUP}，随机抽样。")
            sel = np.random.choice(idx, size=MAX_CELLS_PER_GROUP, replace=False)
            use_idx.append(sel)
        else:
            use_idx.append(idx)

    if len(use_idx) == 0:
        raise SystemExit("[ERROR] 两组都没有可用细胞。")

    use_idx = np.concatenate(use_idx)
    adata_sub = adata[use_idx].copy()
    print(f"[INFO] 抽样后对象维度: n_cells={adata_sub.n_obs}, n_genes={adata_sub.n_vars}")

    # ========= 新增：按表达细胞比例过滤基因，去掉极低表达噪音 =========
    print("[INFO] 进行基因过滤：保留在任一组中表达细胞比例 >= 5% 的基因")
    X = adata_sub.X
    if sp.issparse(X):
        X_bin = X > 0
    else:
        X_bin = X > 0

    labels = adata_sub.obs["tscm_label"].values
    mask_tscm = labels == "TSCM_high"
    mask_other = labels == "T_other"

    expr_frac_tscm = np.asarray(X_bin[mask_tscm].mean(axis=0)).ravel()
    expr_frac_other = np.asarray(X_bin[mask_other].mean(axis=0)).ravel()

    MIN_FRAC = 0.05  # 至少 5% 细胞表达，可根据需要调到 0.1
    keep = (expr_frac_tscm >= MIN_FRAC) | (expr_frac_other >= MIN_FRAC)

    n_before = adata_sub.n_vars
    adata_sub = adata_sub[:, keep].copy()
    n_after = adata_sub.n_vars
    print(f"[INFO] 基因过滤: {n_before} -> {n_after} (保留表达比例 >= {MIN_FRAC} 的基因)")
    # ========= 基因过滤结束 =========

    # 使用当前 X (已 log1p) 做差异分析
    print("[INFO] 开始差异基因分析: TSCM_high vs T_other (wilcoxon)...")
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="tscm_label",
        groups=["TSCM_high"],
        reference="T_other",
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        n_genes=adata_sub.n_vars,
        use_raw=False,
    )

    # 提取结果
    res = adata_sub.uns["rank_genes_groups"]
    scores = pd.DataFrame({
        "gene":         res["names"]["TSCM_high"],
        "logfoldchange": res["logfoldchanges"]["TSCM_high"],
        "pvals":        res["pvals"]["TSCM_high"],
        "pvals_adj":    res["pvals_adj"]["TSCM_high"],
        "scores":       res["scores"]["TSCM_high"],
    })

    # 去掉缺失值/无穷大
    scores = scores.replace([np.inf, -np.inf], np.nan).dropna()

    # 保存全表
    out_full = OUT_DIR / "tscm_vs_other_de.full.tsv"
    scores.to_csv(out_full, sep="\t", index=False)
    print("[INFO] 已写出完整 DE 结果:", out_full)

    # 选取上调基因（logFC > 0），按 logFC 排前 50 个
    up = scores[scores["logfoldchange"] > 0].copy()
    up = up.sort_values("logfoldchange", ascending=False)
    top50 = up.head(50)

    out_top = OUT_DIR / "tscm_vs_other_de.top50.tsv"
    top50.to_csv(out_top, sep="\t", index=False)
    print("[INFO] 已写出 top50 上调 TSCM 基因:", out_top)

    print("[INFO] top 10 上调基因预览:")
    print(top50.head(10).to_string(index=False))

    print("[DONE] Stage 4G (TSCM vs other DE) 完成。")


if __name__ == "__main__":
    main()
