#!/usr/bin/env python3
"""
创建core TSCM signature和full signature的对比分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path
from PIL import Image
import io

ROOT = Path(__file__).resolve().parents[1]

# 读取数据
core_file = ROOT / "results/tables/tscm_core_signature.tsv"
full_file = ROOT / "results/tables/tscm_vs_other_de.full.tsv"

df_core = pd.read_csv(core_file, sep="\t")
df_full = pd.read_csv(full_file, sep="\t")

# 筛选full signature中的显著基因
df_full_sig = df_full[(df_full["logfoldchange"] > 0.5) & (df_full["pvals_adj"] < 0.05)].copy()

print("=" * 60)
print("Core vs Full Signature 对比分析")
print("=" * 60)
print(f"\nCore signature: {len(df_core)} 个基因")
print(f"Full signature: {len(df_full_sig)} 个基因")
print(f"Core占比: {len(df_core)/len(df_full_sig)*100:.1f}%")

# 统计核糖体基因
ribosomal_patterns = [r'^RPS\d+', r'^RPL\d+', r'^RP[0-9]', r'^EIF\d+', r'^EEF\d+', r'^MRPS\d+', r'^MRPL\d+']
ribosomal_mask = pd.Series([False] * len(df_full_sig), index=df_full_sig.index)
for pattern in ribosomal_patterns:
    ribosomal_mask |= df_full_sig["gene"].str.match(pattern, case=False, na=False)

n_ribosomal = ribosomal_mask.sum()
print(f"\nFull signature中核糖体/housekeeping基因数: {n_ribosomal}")
print(f"核糖体基因占比: {n_ribosomal/len(df_full_sig)*100:.1f}%")

# 创建对比表格
comparison_data = {
    "Signature": ["Core TSCM", "Full TSCM"],
    "Gene Count": [len(df_core), len(df_full_sig)],
    "Ribosomal/Housekeeping": [0, n_ribosomal],
    "T-cell Related": [
        df_core["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1', case=False, na=False).sum(),
        df_full_sig["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1', case=False, na=False).sum()
    ],
    "Mean logFC": [df_core["logfoldchange"].mean(), df_full_sig["logfoldchange"].mean()],
    "Median logFC": [df_core["logfoldchange"].median(), df_full_sig["logfoldchange"].median()],
}

df_comparison = pd.DataFrame(comparison_data)
print("\n对比统计:")
print(df_comparison.to_string(index=False))

# 保存对比表格
out_file = ROOT / "results/tables/tscm_signature_comparison.tsv"
df_comparison.to_csv(out_file, sep="\t", index=False)
print(f"\n已保存对比表格: {out_file}")

# 创建可视化
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Core vs Full TSCM Signature Comparison', fontsize=14, fontweight='bold', x=0.02, ha='left')

# A: 基因数量对比
ax1 = axes[0, 0]
categories = ['Core', 'Full']
counts = [len(df_core), len(df_full_sig)]
colors = ['#2E86AB', '#A23B72']
bars = ax1.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
ax1.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
ax1.set_title('A. Gene Count Comparison', fontsize=12, fontweight='bold')
ax1.text(0.05, 0.95, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
for bar, count in zip(bars, counts):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
            f'{count}', ha='center', va='bottom', fontsize=11, fontweight='bold')
ax1.grid(axis='y', alpha=0.3, linestyle='--')

# B: logFC分布对比
ax2 = axes[0, 1]
ax2.hist(df_core["logfoldchange"], bins=30, alpha=0.6, color='#2E86AB', label='Core', edgecolor='black', linewidth=0.5)
ax2.hist(df_full_sig["logfoldchange"], bins=30, alpha=0.6, color='#A23B72', label='Full', edgecolor='black', linewidth=0.5)
ax2.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax2.set_title('B. LogFC Distribution', fontsize=12, fontweight='bold')
ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3, linestyle='--')

# C: Top基因对比（Venn图风格）
ax3 = axes[1, 0]
top_core = df_core.head(10)
top_full = df_full_sig.head(10)
common_genes = set(top_core["gene"]) & set(top_full["gene"])

# 创建堆叠条形图显示top 10基因
top_genes_combined = pd.concat([
    top_core[["gene", "logfoldchange"]].assign(Source='Core'),
    top_full[["gene", "logfoldchange"]].assign(Source='Full')
]).head(20)

# 简化：只显示前10个
top_10_core = top_core.head(10)
y_pos = np.arange(len(top_10_core))
ax3.barh(y_pos, top_10_core["logfoldchange"], color='#2E86AB', alpha=0.7, edgecolor='black', linewidth=0.5)
ax3.set_yticks(y_pos)
ax3.set_yticklabels(top_10_core["gene"], fontsize=9)
ax3.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
ax3.set_title('C. Top 10 Core Signature Genes', fontsize=12, fontweight='bold')
ax3.text(0.05, 0.95, 'C', transform=ax3.transAxes, fontsize=16, fontweight='bold', 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax3.grid(axis='x', alpha=0.3, linestyle='--')

# D: 功能分类对比
ax4 = axes[1, 1]
categories_func = ['T-cell\nRelated', 'Metabolic/\nTranslation', 'Other']
core_counts = [
    df_core["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1', case=False, na=False).sum(),
    df_core["gene"].str.contains(r'PDK|COX|UQCR|NDUF|ATP', case=False, na=False).sum(),
    len(df_core) - df_core["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1|PDK|COX|UQCR|NDUF|ATP', case=False, na=False).sum()
]
full_counts = [
    df_full_sig["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1', case=False, na=False).sum(),
    df_full_sig["gene"].str.contains(r'PDK|COX|UQCR|NDUF|ATP', case=False, na=False).sum(),
    len(df_full_sig) - df_full_sig["gene"].str.contains(r'CCR7|LEF1|TCF7|IL7R|CD27|SELL|LRRN3|TRABD2A|SATB1|FOXP1|PDK|COX|UQCR|NDUF|ATP', case=False, na=False).sum() - n_ribosomal
]

x = np.arange(len(categories_func))
width = 0.35
ax4.bar(x - width/2, core_counts, width, label='Core', color='#2E86AB', alpha=0.7, edgecolor='black', linewidth=0.5)
ax4.bar(x + width/2, full_counts, width, label='Full', color='#A23B72', alpha=0.7, edgecolor='black', linewidth=0.5)
ax4.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
ax4.set_title('D. Functional Category Comparison', fontsize=12, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels(categories_func, fontsize=10)
ax4.text(0.05, 0.95, 'D', transform=ax4.transAxes, fontsize=16, fontweight='bold', 
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax4.legend(fontsize=10)
ax4.grid(axis='y', alpha=0.3, linestyle='--')

plt.tight_layout(rect=[0, 0, 1, 0.98])

# 保存图片
out_fig = ROOT / "results/figures/plos_one_format_optimized/figure_supplement_signature_comparison.tiff"
buf = io.BytesIO()
plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
buf.seek(0)
img = Image.open(buf)
img = img.convert('RGB')
img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
buf.close()
print(f"\n已保存对比图: {out_fig}")

plt.close()

print("\n✅ 对比分析完成！")



