#!/usr/bin/env python3
"""
创建瘦身版TSCM signature (core gene set)
去掉核糖体/housekeeping基因，保留：
1. T细胞分化/归巢相关基因（CCR7, LEF1, TCF7, IL7R, CD27, SELL等）
2. 少量代谢/翻译相关的代表性基因
"""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# 读取完整DE结果
de_file = ROOT / "results/tables/tscm_vs_other_de.full.tsv"
df_de = pd.read_csv(de_file, sep="\t")

# 筛选显著上调基因
df_sig = df_de[(df_de["logfoldchange"] > 0.5) & (df_de["pvals_adj"] < 0.05)].copy()
df_sig = df_sig.sort_values("logfoldchange", ascending=False)

print(f"总显著上调基因数: {len(df_sig)}")

# 定义要排除的基因模式（核糖体/housekeeping）
exclude_patterns = [
    r'^RPS\d+',  # 核糖体小亚基蛋白
    r'^RPL\d+',  # 核糖体大亚基蛋白
    r'^RP[0-9]',  # 其他核糖体蛋白
    r'^EIF\d+',  # 翻译起始因子
    r'^EEF\d+',  # 翻译延伸因子
    r'^MRPS\d+',  # 线粒体核糖体小亚基
    r'^MRPL\d+',  # 线粒体核糖体大亚基
    r'^ACTB$',  # β-actin
    r'^GAPDH$',  # GAPDH
    r'^TUBB\d+',  # β-tubulin
    r'^UBA52$',  # 泛素
    r'^FAU$',  # 40S核糖体蛋白S30
    r'^PABPC\d+',  # Poly(A)结合蛋白
]

import re
exclude_mask = pd.Series([False] * len(df_sig), index=df_sig.index)
for pattern in exclude_patterns:
    exclude_mask |= df_sig["gene"].str.match(pattern, case=False, na=False)

df_filtered = df_sig[~exclude_mask].copy()
print(f"排除核糖体/housekeeping基因后: {len(df_filtered)}")

# 定义核心TSCM marker基因（必须保留）
core_markers = [
    "CCR7", "LEF1", "TCF7", "IL7R", "CD27", "SELL",  # T细胞分化/归巢
    "LRRN3", "TRABD2A", "SATB1", "FOXP1",  # T细胞转录因子
    "BCL2", "CD28", "CD45RA", "CD95",  # T细胞表面标志
]

# 确保核心marker都在
core_genes = []
for marker in core_markers:
    matches = df_filtered[df_filtered["gene"].str.upper() == marker.upper()]
    if len(matches) > 0:
        core_genes.append(matches.iloc[0])
    else:
        # 如果不在过滤后的列表中，从原始列表中找
        matches_orig = df_sig[df_sig["gene"].str.upper() == marker.upper()]
        if len(matches_orig) > 0:
            core_genes.append(matches_orig.iloc[0])
            print(f"  注意: {marker} 在原始列表中但被过滤，已添加回来")

# 选择top基因（按logFC排序）
n_top = 50 - len(core_genes)  # 目标50个基因
if n_top > 0:
    top_genes = df_filtered.head(n_top)
    core_gene_set = pd.concat([pd.DataFrame(core_genes), top_genes]).drop_duplicates(subset="gene")
else:
    core_gene_set = pd.DataFrame(core_genes)

# 添加少量代谢/翻译相关的代表性基因（选择logFC最高的几个）
metabolic_genes = df_sig[df_sig["gene"].str.contains(r'^(PDK|COX|UQCR|NDUF|ATP)', case=False, na=False)]
if len(metabolic_genes) > 0:
    # 选择top 3-5个代谢基因
    top_metabolic = metabolic_genes.head(5)
    core_gene_set = pd.concat([core_gene_set, top_metabolic]).drop_duplicates(subset="gene")

# 按logFC排序
core_gene_set = core_gene_set.sort_values("logfoldchange", ascending=False)

print(f"\n最终Core TSCM gene set: {len(core_gene_set)} 个基因")
print(f"\nTop 20 genes:")
print(core_gene_set[["gene", "logfoldchange", "pvals_adj"]].head(20).to_string(index=False))

# 保存core gene set
output_file = ROOT / "results/tables/tscm_core_signature.tsv"
core_gene_set.to_csv(output_file, sep="\t", index=False)
print(f"\n已保存: {output_file}")

# 保存基因列表（用于后续分析）
gene_list_file = ROOT / "results/tables/tscm_core_signature_genes.txt"
with open(gene_list_file, 'w') as f:
    for gene in core_gene_set["gene"]:
        f.write(f"{gene}\n")
print(f"已保存基因列表: {gene_list_file}")

# 统计信息
print(f"\n=== 统计信息 ===")
print(f"Core gene set大小: {len(core_gene_set)}")
print(f"Full signature大小: {len(df_sig)}")
print(f"排除的核糖体/housekeeping基因数: {exclude_mask.sum()}")

# 分类统计
t_cell_related = core_gene_set["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1', case=False, na=False).sum()
print(f"T细胞相关基因数: {t_cell_related}")
print(f"其他功能基因数: {len(core_gene_set) - t_cell_related}")



