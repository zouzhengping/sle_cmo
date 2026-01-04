#!/usr/bin/env python3
"""
Stage 4E: TSCM / naive / memory / effector scoring in T-cell subset (GSE137029 SLE cohort)

输入:
  data/qc/GSE137029_sle.tcells.clustered.h5ad

输出:
  data/qc/GSE137029_sle.tcells.cluster_k20_scores.tsv   # 每个 cluster 的平均打分 + 细胞数
  data/qc/GSE137029_sle.tcells.clustered.tscm.h5ad      # 带所有 score 的 T 细胞对象
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.h5ad"
OUT_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm.h5ad"
OUT_TSV  = ROOT / "data/qc" / "GSE137029_sle.tcells.cluster_k20_scores.tsv"

print("[INFO] 读取聚类后的 T 细胞对象:", IN_H5AD)
adata = sc.read_h5ad(IN_H5AD)
print("[INFO] n_cells={}, n_genes={}".format(adata.n_obs, adata.n_vars))

if "cluster_k20" not in adata.obs.columns:
    raise SystemExit("[ERROR] 未找到 cluster_k20，请确认 44_stage4_tcells_cluster_only_sle.py 已成功运行。")

# --------------------------
# 1. 定义各类 marker 列表
# --------------------------

# 干性记忆性 T 细胞 (TSCM) 倾向 marker
tscm_markers = [
    "CCR7",   # 归巢受体
    "SELL",   # CD62L
    "IL7R",   # CD127
    "LEF1",
    "TCF7",
    "BCL2",
    "CXCR3",
    "CD27",
    "CD28",
    "FAS",    # CD95
]

# naive T 细胞 marker（CD95-/效应低）
naive_markers = [
    "CCR7",
    "SELL",
    "IL7R",
    "LEF1",
    "TCF7",
]

# 中央记忆 / 效应记忆 (粗略)
memory_markers = [
    "CCR7",
    "IL7R",
    "SELL",
    "CD27",
    "CD28",
]

# 效应 / 细胞毒性 marker
effector_markers = [
    "GZMB",
    "PRF1",
    "GNLY",
    "NKG7",
    "IFNG",
    "KLRD1",
    "KLRB1",
]

# 疲劳 / 耗竭相关 marker
exhaustion_markers = [
    "PDCD1",   # PD-1
    "LAG3",
    "CTLA4",
    "HAVCR2",  # TIM-3
    "TIGIT",
    "TOX",
]

# CD4 / CD8 倾向
cd4_markers = ["CD4", "CCR7", "IL7R"]
cd8_markers = ["CD8A", "CD8B", "GZMB", "PRF1"]

marker_panels = {
    "score_tscm": tscm_markers,
    "score_naive": naive_markers,
    "score_memory": memory_markers,
    "score_effector": effector_markers,
    "score_exhaustion": exhaustion_markers,
    "score_cd4like": cd4_markers,
    "score_cd8like": cd8_markers,
}

# --------------------------
# 2. 对每一组 marker 做打分
# --------------------------
var_names_upper = adata.var_names.str.upper()
adata.var["SYMBOL_UPPER"] = var_names_upper

def map_markers(panel):
    """把 marker 列表映射到实际存在于数据中的基因名（忽略大小写）"""
    found = []
    missing = []
    for g in panel:
        mask = var_names_upper == g.upper()
        if mask.any():
            # 取第一个匹配到的实际名
            found.append(adata.var_names[mask.argmax()])
        else:
            missing.append(g)
    return found, missing

all_missing = {}

for score_name, genes in marker_panels.items():
    found, missing = map_markers(genes)
    all_missing[score_name] = missing
    print(f"[INFO] {score_name}: {len(found)} markers used, {len(missing)} missing")
    if missing:
        print(f"       缺失 marker: {', '.join(missing)}")
    if len(found) == 0:
        print(f"[WARN] {score_name} 没有可用 marker，跳过该打分。")
        continue
    sc.tl.score_genes(adata, gene_list=found, score_name=score_name)

# --------------------------
# 3. 汇总到 cluster_k20 层面
# --------------------------
print("[INFO] 按 cluster_k20 汇总平均得分和细胞数...")
score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
group = adata.obs.groupby("cluster_k20")

summary = group[score_cols].mean()
summary["n_cells"] = group.size()

# 方便排序：增加一个 “tscm-ness” 简单指标 = score_tscm - score_effector
if "score_tscm" in summary.columns and "score_effector" in summary.columns:
    summary["tscm_minus_effector"] = summary["score_tscm"] - summary["score_effector"]

summary = summary.sort_values(
    by=[c for c in ["tscm_minus_effector", "score_tscm"] if c in summary.columns],
    ascending=False
)

print("[INFO] 聚类评分表（前几行）:")
print(summary.head().to_string())

summary.to_csv(OUT_TSV, sep="\t")
print("[INFO] 已写入 cluster 评分表:", OUT_TSV)

# --------------------------
# 4. 写回带 score 的 AnnData
# --------------------------
print("[INFO] 保存带所有 score 的 T 细胞对象:", OUT_H5AD)
adata.write_h5ad(OUT_H5AD)

print("[DONE] Stage 4E (TSCM scoring) 完成。")
