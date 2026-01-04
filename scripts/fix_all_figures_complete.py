#!/usr/bin/env python3
"""
修复所有图片的标题、编号和位置
按照以下顺序：
Figure 1. Identification of TSCM cells in SLE T-cell populations
Figure 2. TSCM signature definition and functional enrichment
Figure 3. Distribution of TSCM signature scores in validation dataset
Figure 4. Association between TSCM signature and SLEDAI disease activity (合并A, B, C)
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
from scipy.stats import kruskal, kendalltau
import warnings
warnings.filterwarnings('ignore')
from PIL import Image
import io

try:
    import scanpy as sc
    sc.settings.verbosity = 3
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 设置字体 - 解决乱码问题
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10

# 输出目录
FIG_DIR = ROOT / "results/figures/plos_one_format_optimized"
OUT_DIR = ROOT / "results/validation/sledai"
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("修复所有图片的标题、编号和位置")
print("=" * 80)

# ============================================================================
# Figure 1: Identification of TSCM cells in SLE T-cell populations
# ============================================================================
print("\n[Figure 1] 生成TSCM细胞识别图...")

if HAS_SCANPY:
    h5ad_file = ROOT / "data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad"
    if h5ad_file.exists():
        adata = sc.read_h5ad(h5ad_file)
        
        if "X_umap" in adata.obsm and "cluster_k20" in adata.obs.columns:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            fig.suptitle('Figure 1. Identification of TSCM cells in SLE T-cell populations', 
                        fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
            
            # A: UMAP clusters
            ax1 = axes[0]
            sc.pl.umap(adata, color='cluster_k20', ax=ax1, show=False, 
                      legend_loc='right margin', title='',
                      frameon=False, size=1 if adata.n_obs > 50000 else 10)
            ax1.text(0.05, 0.95, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            ax1.set_title('UMAP visualization of T-cell clusters', fontsize=12, fontweight='bold')
            
            # B: TSCM label
            ax2 = axes[1]
            if "tscm_label" in adata.obs.columns:
                sc.pl.umap(adata, color='tscm_label', ax=ax2, show=False,
                          legend_loc='right margin', title='',
                          frameon=False, size=1 if adata.n_obs > 50000 else 10)
                ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                ax2.set_title('TSCM-high clusters highlighted', fontsize=12, fontweight='bold')
            
            # C: TSCM score
            ax3 = axes[2]
            if "score_tscm" in adata.obs.columns:
                sc.pl.umap(adata, color='score_tscm', ax=ax3, show=False,
                          cmap='Reds', title='',
                          frameon=False, size=1 if adata.n_obs > 50000 else 10)
                ax3.text(0.05, 0.95, 'C', transform=ax3.transAxes, fontsize=16, fontweight='bold', 
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                ax3.set_title('TSCM signature score distribution', fontsize=12, fontweight='bold')
            
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            out_fig = FIG_DIR / "figure1_tscm_identification.tiff"
            buf = io.BytesIO()
            plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
            buf.seek(0)
            img = Image.open(buf)
            img = img.convert('RGB')
            img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
            buf.close()
            print(f"  ✓ 已保存: {out_fig}")
            plt.close()
        else:
            print(f"  ⚠ 缺少必要数据，跳过Figure 1")
    else:
        print(f"  ⚠ 未找到数据文件: {h5ad_file}")
else:
    print(f"  ⚠ scanpy未安装，跳过Figure 1")

# ============================================================================
# Figure 2: TSCM signature definition and functional enrichment
# ============================================================================
print("\n[Figure 2] 生成TSCM signature定义和富集分析图...")

de_file = ROOT / "results/tables/tscm_vs_other_de.full.tsv"
if de_file.exists():
    df_de = pd.read_csv(de_file, sep="\t")
    
    # 筛选显著基因
    sig_up = df_de[(df_de["logfoldchange"] > 0.5) & (df_de["pvals_adj"] < 0.05)].copy()
    sig_down = df_de[(df_de["logfoldchange"] < -0.5) & (df_de["pvals_adj"] < 0.05)].copy()
    ns = df_de[(df_de["pvals_adj"] >= 0.05) | 
               ((df_de["pvals_adj"] < 0.05) & (df_de["logfoldchange"].abs() < 0.5))].copy()
    
    # 创建图表 - 包含火山图和富集分析
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle('Figure 2. TSCM signature definition and functional enrichment', 
                 fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
    
    # A: 火山图
    ax1 = axes[0]
    ax1.scatter(ns["logfoldchange"], -np.log10(ns["pvals_adj"] + 1e-300), 
               c="gray", alpha=0.5, s=15, label="Not significant", edgecolors='none')
    
    if len(sig_down) > 0:
        ax1.scatter(sig_down["logfoldchange"], 
                  -np.log10(sig_down["pvals_adj"] + 1e-300),
                  c="blue", alpha=0.6, s=25, label="Downregulated", edgecolors='black', linewidths=0.5)
    
    if len(sig_up) > 0:
        ax1.scatter(sig_up["logfoldchange"], 
                  -np.log10(sig_up["pvals_adj"] + 1e-300),
                  c="red", alpha=0.6, s=25, label="Upregulated", edgecolors='black', linewidths=0.5)
    
    # 标记关键TSCM marker
    tscm_markers = ["CCR7", "LEF1", "TCF7", "SELL", "IL7R", "CD27", "LRRN3"]
    for marker in tscm_markers:
        if marker in df_de["gene"].values:
            row = df_de[df_de["gene"] == marker].iloc[0]
            ax1.scatter(row["logfoldchange"], 
                      -np.log10(row["pvals_adj"] + 1e-300),
                      c="orange", s=150, marker="*", edgecolors="black", linewidths=1.5, zorder=5)
            ax1.annotate(marker, 
                       (row["logfoldchange"], -np.log10(row["pvals_adj"] + 1e-300)),
                       xytext=(8, 8), textcoords="offset points", fontsize=10, fontweight="bold",
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    
    # 计算x轴范围，使0点居中
    x_max = max(abs(df_de["logfoldchange"].min()), abs(df_de["logfoldchange"].max()))
    x_range = x_max * 1.1
    ax1.set_xlim(-x_range, x_range)
    
    # 添加阈值线
    ax1.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1.5, alpha=0.7, label='p = 0.05')
    ax1.axvline(0, color="black", linestyle="--", linewidth=1.5, alpha=0.7)
    ax1.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    ax1.axvline(-0.5, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    
    ax1.set_xlabel("Log2 Fold Change", fontsize=12, fontweight="bold")
    ax1.set_ylabel("-Log10 Adjusted P-value", fontsize=12, fontweight="bold")
    ax1.text(0.05, 0.95, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax1.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # B: 富集分析（如果有数据）
    ax2 = axes[1]
    enrichment_file = ROOT / "results/tables/tscm_enrichment_go.tsv"
    if enrichment_file.exists():
        df_enrich = pd.read_csv(enrichment_file, sep="\t")
        if len(df_enrich) > 0:
            # 选择top 15
            df_top = df_enrich.head(15).copy()
            df_top = df_top.sort_values('Adjusted P-value', ascending=True)
            
            # 计算-log10(p)
            df_top['-log10(p)'] = -np.log10(df_top['Adjusted P-value'] + 1e-300)
            
            y_pos = np.arange(len(df_top))
            colors_bar = plt.cm.Reds(np.linspace(0.4, 0.9, len(df_top)))
            bars = ax2.barh(y_pos, df_top['-log10(p)'], color=colors_bar)
            
            ax2.set_yticks(y_pos)
            ax2.set_yticklabels(df_top['Term'].str[:50], fontsize=9)
            ax2.set_xlabel("-Log10 Adjusted P-value", fontsize=12, fontweight="bold")
            ax2.set_title("Top enriched GO biological processes", fontsize=12, fontweight="bold")
            ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            ax2.grid(True, axis='x', alpha=0.3, linestyle='--')
        else:
            ax2.text(0.5, 0.5, 'Enrichment data not available', 
                    transform=ax2.transAxes, ha='center', va='center', fontsize=12)
            ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax2.text(0.5, 0.5, 'Enrichment data not available', 
                transform=ax2.transAxes, ha='center', va='center', fontsize=12)
        ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_fig = FIG_DIR / "figure2_tscm_signature_enrichment.tiff"
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
    buf.seek(0)
    img = Image.open(buf)
    img = img.convert('RGB')
    img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
    buf.close()
    print(f"  ✓ 已保存: {out_fig}")
    plt.close()
else:
    print(f"  ⚠ 未找到DE结果文件: {de_file}")

# ============================================================================
# Figure 3: Distribution of TSCM signature scores in validation dataset
# ============================================================================
print("\n[Figure 3] 生成TSCM signature score分布图...")

data_file = OUT_DIR / "GSE88884_tscm_vs_sledai.tsv"
if data_file.exists():
    df = pd.read_csv(data_file, sep="\t", index_col=0)
    df = df[df['tscm_score'].notna()].copy()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle('Figure 3. Distribution of TSCM signature scores in validation dataset', 
                 fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
    
    ax.hist(df['tscm_score'], bins=50, alpha=0.7, color='steelblue', edgecolor='black', linewidth=1)
    ax.axvline(df['tscm_score'].mean(), color='red', linestyle='--', linewidth=2, label=f'Mean = {df["tscm_score"].mean():.3f}')
    ax.axvline(df['tscm_score'].median(), color='green', linestyle='--', linewidth=2, label=f'Median = {df["tscm_score"].median():.3f}')
    
    ax.set_xlabel('TSCM Signature Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title(f'Distribution across {len(df)} samples from GSE88884', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # 添加统计信息
    stats_text = f'Mean = {df["tscm_score"].mean():.3f} ± {df["tscm_score"].std():.3f}\n'
    stats_text += f'Range: {df["tscm_score"].min():.3f} - {df["tscm_score"].max():.3f}'
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_fig = FIG_DIR / "figure3_tscm_score_distribution.tiff"
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
    buf.seek(0)
    img = Image.open(buf)
    img = img.convert('RGB')
    img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
    buf.close()
    print(f"  ✓ 已保存: {out_fig}")
    plt.close()
else:
    print(f"  ⚠ 未找到数据文件: {data_file}")

# ============================================================================
# Figure 4: Association between TSCM signature and SLEDAI disease activity
# 合并A, B, C三个子图
# ============================================================================
print("\n[Figure 4] 生成TSCM signature与SLEDAI关联图（合并A, B, C）...")

if data_file.exists():
    df = pd.read_csv(data_file, sep="\t", index_col=0)
    df = df[df['SLEDAI'].notna() & df['tscm_score'].notna()].copy()
    
    # 计算相关性
    pearson_r, pearson_p = stats.pearsonr(df['tscm_score'], df['SLEDAI'])
    spearman_r, spearman_p = stats.spearmanr(df['tscm_score'], df['SLEDAI'])
    
    # 定义分组
    df['SLEDAI_group'] = pd.cut(
        df['SLEDAI'],
        bins=[0, 4, 8, 12, 100],
        labels=['Low (0-4)', 'Moderate (4-8)', 'High (8-12)', 'Very high (>12)']
    )
    
    # 计算统计
    group_order = ['Low (0-4)', 'Moderate (4-8)', 'High (8-12)', 'Very high (>12)']
    groups_data = [df[df['SLEDAI_group'] == g]['tscm_score'].values for g in group_order if g in df['SLEDAI_group'].values]
    labels_for_plot = [g for g in group_order if g in df['SLEDAI_group'].values]
    
    h_stat, kw_p = kruskal(*groups_data)
    
    # 趋势分析
    group_medians_sledai = []
    group_means_tscm = []
    group_sems = []
    for group in group_order:
        if group in df['SLEDAI_group'].values:
            group_data = df[df['SLEDAI_group'] == group]
            group_medians_sledai.append(group_data['SLEDAI'].median())
            group_means_tscm.append(group_data['tscm_score'].mean())
            group_sems.append(group_data['tscm_score'].std() / np.sqrt(len(group_data)))
    
    all_scores = np.concatenate(groups_data)
    all_ranks = np.concatenate([[i] * len(g) for i, g in enumerate(groups_data)])
    tau, tau_p = kendalltau(all_ranks, all_scores)
    
    # 创建图表 - 三个子图
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Figure 4. Association between TSCM signature and SLEDAI disease activity', 
                 fontsize=14, fontweight='bold', x=0.02, ha='left', y=0.98)
    
    # A: 散点图
    ax1 = axes[0]
    ax1.scatter(df['SLEDAI'], df['tscm_score'], alpha=0.6, s=30, color='steelblue', edgecolors='black', linewidths=0.5)
    
    # 添加回归线
    z = np.polyfit(df['SLEDAI'], df['tscm_score'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(df['SLEDAI'].min(), df['SLEDAI'].max(), 100)
    ax1.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2, label='Regression line')
    
    ax1.set_xlabel('SLEDAI Score', fontsize=12, fontweight='bold')
    ax1.set_ylabel('TSCM Signature Score', fontsize=12, fontweight='bold')
    ax1.text(0.05, 0.95, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax1.text(0.05, 0.85, f'Pearson r = {pearson_r:.3f}, p = {pearson_p:.4f}\nSpearman ρ = {spearman_r:.3f}, p = {spearman_p:.4f}', 
            transform=ax1.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.legend(fontsize=9)
    
    # B: 箱线图
    ax2 = axes[1]
    bp = ax2.boxplot(groups_data, labels=labels_for_plot, patch_artist=True, 
                     widths=0.6, showmeans=True, meanline=True)
    colors = ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3']
    for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax2.set_ylabel('TSCM Signature Score', fontsize=12, fontweight='bold')
    ax2.set_xlabel('SLEDAI Activity Group', fontsize=12, fontweight='bold')
    ax2.text(0.05, 0.95, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.text(0.05, 0.85, f'Kruskal-Wallis: H = {h_stat:.2f}, p = {kw_p:.4f}', 
            transform=ax2.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax2.tick_params(axis='x', rotation=0, labelsize=10)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    # C: 趋势图
    ax3 = axes[2]
    ax3.errorbar(range(len(group_medians_sledai)), group_means_tscm, yerr=group_sems, 
                marker='o', markersize=10, linewidth=2.5, capsize=6, capthick=2.5,
                color='steelblue', markerfacecolor='white', markeredgecolor='steelblue', 
                markeredgewidth=2, elinewidth=2)
    
    ax3.set_xticks(range(len(group_order)))
    ax3.set_xticklabels(group_order, fontsize=10, fontweight='bold')
    ax3.set_ylabel('TSCM Signature Score (Mean ± SEM)', fontsize=12, fontweight='bold')
    ax3.set_xlabel('SLEDAI Activity Group', fontsize=12, fontweight='bold')
    ax3.text(0.05, 0.95, 'C', transform=ax3.transAxes, fontsize=16, fontweight='bold', 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax3.text(0.05, 0.85, f"Kendall's τ = {tau:.4f}, p = {tau_p:.4f}", 
            transform=ax3.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax3.grid(alpha=0.3, linestyle='--', axis='y')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_fig = FIG_DIR / "figure4_tscm_sledai_association.tiff"
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
    buf.seek(0)
    img = Image.open(buf)
    img = img.convert('RGB')
    img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
    buf.close()
    print(f"  ✓ 已保存: {out_fig}")
    plt.close()
else:
    print(f"  ⚠ 未找到数据文件: {data_file}")

print("\n" + "=" * 80)
print("✓ 所有图片修复完成！")
print("=" * 80)
print(f"\n输出目录: {FIG_DIR}")
print("\n生成的图表:")
for fig_file in sorted(FIG_DIR.glob("figure*.tiff")):
    size_mb = fig_file.stat().st_size / (1024 * 1024)
    print(f"  - {fig_file.name}: {size_mb:.2f} MB")


