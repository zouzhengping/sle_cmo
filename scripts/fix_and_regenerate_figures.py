#!/usr/bin/env python3
"""
修复并重新生成所有图表，解决以下问题：
1. GSE88884图表排列过于紧密，字体重叠
2. GSE88884图表横坐标乱码（方框）- 字体问题
3. 火山图横坐标0点居中
4. 所有图表标题位置改为左边，添加Figure编号
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from scipy import stats
from scipy.stats import kruskal, mannwhitneyu, kendalltau
from PIL import Image
import io
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 设置字体 - 解决乱码问题
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10

# 输出目录
OUT_DIR = ROOT / "results/validation/sledai"
FIG_DIR = ROOT / "results/figures/plos_one_format_optimized"
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("修复并重新生成图表")
print("=" * 80)

# ============================================================================
# 1. 修复GSE88884 TSCM vs SLEDAI图表
# ============================================================================
print("\n[1] 修复GSE88884 TSCM vs SLEDAI图表...")

data_file = OUT_DIR / "GSE88884_tscm_vs_sledai.tsv"
if not data_file.exists():
    print(f"  ⚠ 未找到数据文件: {data_file}")
else:
    df = pd.read_csv(data_file, sep="\t", index_col=0)
    df = df[df['SLEDAI'].notna() & df['tscm_score'].notna()].copy()
    
    # 计算相关性
    pearson_r, pearson_p = stats.pearsonr(df['tscm_score'], df['SLEDAI'])
    spearman_r, spearman_p = stats.spearmanr(df['tscm_score'], df['SLEDAI'])
    
    # 创建图表 - 使用更大的间距
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Figure 4. Association between TSCM signature and SLEDAI disease activity', 
                 fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
    
    # 子图A: 散点图
    ax1 = axes[0]
    ax1.scatter(df['SLEDAI'], df['tscm_score'], alpha=0.6, s=30, color='steelblue', edgecolors='black', linewidths=0.5)
    
    # 添加回归线
    z = np.polyfit(df['SLEDAI'], df['tscm_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(df['SLEDAI'].min(), df['SLEDAI'].max(), 100)
    ax1.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2, label=f'Regression line')
    
    ax1.set_xlabel('SLEDAI Score', fontsize=12, fontweight='bold')
    ax1.set_ylabel('TSCM Signature Score', fontsize=12, fontweight='bold')
    ax1.text(0.05, 0.95, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax1.text(0.05, 0.85, f'Pearson r = {pearson_r:.3f}, p = {pearson_p:.4f}\nSpearman ρ = {spearman_r:.3f}, p = {spearman_p:.4f}', 
             transform=ax1.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.legend(fontsize=9)
    
    # 子图B: 箱线图
    ax2 = axes[1]
    
    # 定义分组
    df['SLEDAI_group'] = pd.cut(
        df['SLEDAI'],
        bins=[0, 4, 8, 12, 100],
        labels=['Low (0-4)', 'Moderate (4-8)', 'High (8-12)', 'Very high (>12)']
    )
    
    group_order = ['Low (0-4)', 'Moderate (4-8)', 'High (8-12)', 'Very high (>12)']
    data_for_plot = [df[df['SLEDAI_group'] == g]['tscm_score'].values 
                     for g in group_order if g in df['SLEDAI_group'].values]
    labels_for_plot = [g for g in group_order if g in df['SLEDAI_group'].values]
    
    bp = ax2.boxplot(data_for_plot, labels=labels_for_plot, patch_artist=True, 
                     widths=0.6, showmeans=True, meanline=True)
    colors = ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3']
    for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Kruskal-Wallis检验
    groups_data = [df[df['SLEDAI_group'] == g]['tscm_score'].values for g in labels_for_plot]
    h_stat, kw_p = kruskal(*groups_data)
    
    ax2.set_ylabel('TSCM Signature Score', fontsize=12, fontweight='bold')
    ax2.set_xlabel('SLEDAI Activity Group', fontsize=12, fontweight='bold')
    ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.text(0.05, 0.85, f'Kruskal-Wallis: H = {h_stat:.2f}, p = {kw_p:.4f}', 
             transform=ax2.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.tick_params(axis='x', rotation=0, labelsize=10)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为标题留出空间
    
    # 保存为TIFF（使用PIL进行压缩）
    out_fig = FIG_DIR / "GSE88884_tscm_vs_sledai.tiff"
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
    buf.seek(0)
    img = Image.open(buf)
    img = img.convert('RGB')
    img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
    buf.close()
    print(f"  ✓ 已保存: {out_fig}")
    plt.close()

# ============================================================================
# 2. 修复GSE88884分组详细图表
# ============================================================================
print("\n[2] 修复GSE88884分组详细图表...")

if not data_file.exists():
    print(f"  ⚠ 未找到数据文件: {data_file}")
else:
    df = pd.read_csv(data_file, sep="\t", index_col=0)
    df = df[df['SLEDAI'].notna() & df['tscm_score'].notna()].copy()
    
    # 定义分组
    df['group_standard'] = pd.cut(
        df['SLEDAI'],
        bins=[0, 4, 8, 12, 100],
        labels=['Low (0-4)', 'Moderate (4-8)', 'High (8-12)', 'Very high (>12)']
    )
    
    # 计算统计
    group_order = ['Low (0-4)', 'Moderate (4-8)', 'High (8-12)', 'Very high (>12)']
    group_medians_sledai = []
    group_means_tscm = []
    group_sems = []
    for group in group_order:
        if group in df['group_standard'].values:
            group_data = df[df['group_standard'] == group]
            group_medians_sledai.append(group_data['SLEDAI'].median())
            group_means_tscm.append(group_data['tscm_score'].mean())
            group_sems.append(group_data['tscm_score'].std() / np.sqrt(len(group_data)))
    
    # 趋势检验
    groups_data = [df[df['group_standard'] == g]['tscm_score'].values for g in group_order if g in df['group_standard'].values]
    all_scores = np.concatenate(groups_data)
    all_ranks = np.concatenate([[i] * len(g) for i, g in enumerate(groups_data)])
    tau, tau_p = kendalltau(all_ranks, all_scores)
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle('Figure 4C. TSCM signature score trend across SLEDAI activity groups', 
                 fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
    
    ax.errorbar(range(len(group_medians_sledai)), group_means_tscm, yerr=group_sems, 
                marker='o', markersize=10, linewidth=2.5, capsize=6, capthick=2.5,
                color='steelblue', markerfacecolor='white', markeredgecolor='steelblue', 
                markeredgewidth=2, elinewidth=2)
    
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, fontsize=11, fontweight='bold')
    ax.set_ylabel('TSCM Signature Score (Mean ± SEM)', fontsize=12, fontweight='bold')
    ax.set_xlabel('SLEDAI Activity Group', fontsize=12, fontweight='bold')
    ax.text(0.05, 0.95, 'C', transform=ax.transAxes, fontsize=16, fontweight='bold', 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(0.05, 0.85, f"Kendall's τ = {tau:.4f}, p = {tau_p:.4f}", 
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.grid(alpha=0.3, linestyle='--', axis='y')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存为TIFF（使用PIL进行压缩）
    out_fig = FIG_DIR / "GSE88884_tscm_by_sledai_groups_detailed.tiff"
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
    buf.seek(0)
    img = Image.open(buf)
    img = img.convert('RGB')
    img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
    buf.close()
    print(f"  ✓ 已保存: {out_fig}")
    plt.close()

# ============================================================================
# 3. 修复火山图 - 横坐标0点居中
# ============================================================================
print("\n[3] 修复火山图...")

de_file = ROOT / "results/tables/tscm_vs_other_de.full.tsv"
if not de_file.exists():
    print(f"  ⚠ 未找到DE结果文件: {de_file}")
else:
    df_de = pd.read_csv(de_file, sep="\t")
    
    # 筛选显著基因
    sig_up = df_de[(df_de["logfoldchange"] > 0.5) & (df_de["pvals_adj"] < 0.05)].copy()
    sig_down = df_de[(df_de["logfoldchange"] < -0.5) & (df_de["pvals_adj"] < 0.05)].copy()
    ns = df_de[(df_de["pvals_adj"] >= 0.05) | 
               ((df_de["pvals_adj"] < 0.05) & (df_de["logfoldchange"].abs() < 0.5))].copy()
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.suptitle('Figure 2A. Differential expression analysis: TSCM-high vs other T cells', 
                 fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
    
    # 绘制点
    ax.scatter(ns["logfoldchange"], -np.log10(ns["pvals_adj"] + 1e-300), 
               c="gray", alpha=0.5, s=15, label="Not significant", edgecolors='none')
    
    if len(sig_down) > 0:
        ax.scatter(sig_down["logfoldchange"], 
                  -np.log10(sig_down["pvals_adj"] + 1e-300),
                  c="blue", alpha=0.6, s=25, label="Downregulated", edgecolors='black', linewidths=0.5)
    
    if len(sig_up) > 0:
        ax.scatter(sig_up["logfoldchange"], 
                  -np.log10(sig_up["pvals_adj"] + 1e-300),
                  c="red", alpha=0.6, s=25, label="Upregulated", edgecolors='black', linewidths=0.5)
    
    # 标记关键TSCM marker
    tscm_markers = ["CCR7", "LEF1", "TCF7", "SELL", "IL7R", "CD27", "LRRN3"]
    for marker in tscm_markers:
        if marker in df_de["gene"].values:
            row = df_de[df_de["gene"] == marker].iloc[0]
            ax.scatter(row["logfoldchange"], 
                      -np.log10(row["pvals_adj"] + 1e-300),
                      c="orange", s=150, marker="*", edgecolors="black", linewidths=1.5, zorder=5)
            ax.annotate(marker, 
                       (row["logfoldchange"], -np.log10(row["pvals_adj"] + 1e-300)),
                       xytext=(8, 8), textcoords="offset points", fontsize=10, fontweight="bold",
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    
    # 计算x轴范围，使0点居中
    x_max = max(abs(df_de["logfoldchange"].min()), abs(df_de["logfoldchange"].max()))
    x_range = x_max * 1.1  # 增加10%的边距
    
    # 设置x轴范围，使0点居中
    ax.set_xlim(-x_range, x_range)
    
    # 添加阈值线
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1.5, alpha=0.7, label='p = 0.05')
    ax.axvline(0, color="black", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(-0.5, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    
    ax.set_xlabel("Log2 Fold Change", fontsize=12, fontweight="bold")
    ax.set_ylabel("-Log10 Adjusted P-value", fontsize=12, fontweight="bold")
    ax.text(0.05, 0.95, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold', 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存为TIFF（使用PIL进行压缩）
    out_fig = FIG_DIR / "tscm_volcano_plot.tiff"
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
    buf.seek(0)
    img = Image.open(buf)
    img = img.convert('RGB')
    img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
    buf.close()
    print(f"  ✓ 已保存: {out_fig}")
    plt.close()

print("\n" + "=" * 80)
print("✓ 所有图表修复完成！")
print("=" * 80)
print(f"\n输出目录: {FIG_DIR}")
print("\n生成的图表:")
for fig_file in sorted(FIG_DIR.glob("*.tiff")):
    size_mb = fig_file.stat().st_size / (1024 * 1024)
    print(f"  - {fig_file.name}: {size_mb:.2f} MB")

