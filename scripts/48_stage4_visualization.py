#!/usr/bin/env python3
"""
Stage 4H: 生成Stage 4所有分析的可视化图表

输入:
  data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad
  data/qc/GSE137029_sle.tcells.cluster_k20_scores.tsv

输出:
  results/figures/stage4_qc_violin.pdf          # QC指标小提琴图
  results/figures/stage4_umap_clusters.pdf      # UMAP按cluster着色
  results/figures/stage4_umap_tscm.pdf          # UMAP按TSCM label着色
  results/figures/stage4_tscm_scores.pdf        # TSCM打分分布
  results/figures/stage4_cluster_scores_heatmap.pdf  # Cluster打分热图
  results/figures/stage4_marker_genes.pdf       # 关键marker基因表达
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

try:
    import scanpy as sc
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=300, facecolor='white')
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False
    print("[ERROR] scanpy未安装")

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm_labeled.h5ad"
IN_SCORES = ROOT / "data/qc" / "GSE137029_sle.tcells.cluster_k20_scores.tsv"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 输出文件
FIG_QC_VIOLIN = FIG_DIR / "stage4_qc_violin.pdf"
FIG_UMAP_CLUSTERS = FIG_DIR / "stage4_umap_clusters.pdf"
FIG_UMAP_TSCM = FIG_DIR / "stage4_umap_tscm.pdf"
FIG_TSCM_SCORES = FIG_DIR / "stage4_tscm_scores.pdf"
FIG_CLUSTER_SCORES_HEATMAP = FIG_DIR / "stage4_cluster_scores_heatmap.pdf"
FIG_MARKER_GENES = FIG_DIR / "stage4_marker_genes.pdf"
FIG_CLUSTER_COMPOSITION = FIG_DIR / "stage4_cluster_composition.pdf"


def plot_qc_violin(adata, out_file):
    """绘制QC指标小提琴图"""
    print(f"[INFO] 绘制QC指标小提琴图: {out_file}")
    
    if adata.n_obs > 100000:
        # 如果细胞数太多，随机抽样
        print(f"[INFO] 细胞数过多 ({adata.n_obs})，随机抽样100,000个用于可视化")
        np.random.seed(42)
        sample_idx = np.random.choice(adata.n_obs, size=100000, replace=False)
        adata_viz = adata[sample_idx].copy()
    else:
        adata_viz = adata
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # n_genes_by_counts
    sc.pl.violin(adata_viz, keys='n_genes_by_counts', ax=axes[0], show=False, 
                 title='Number of genes', rotation=90)
    
    # total_counts
    sc.pl.violin(adata_viz, keys='total_counts', ax=axes[1], show=False, 
                 title='Total counts', rotation=90)
    
    # pct_counts_mt
    sc.pl.violin(adata_viz, keys='pct_counts_mt', ax=axes[2], show=False, 
                 title='% Mitochondrial', rotation=90)
    
    plt.suptitle("QC Metrics Distribution", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存QC小提琴图: {out_file}")
    plt.close()


def plot_umap_clusters(adata, out_file):
    """绘制UMAP按cluster着色"""
    print(f"[INFO] 绘制UMAP cluster图: {out_file}")
    
    if "X_umap" not in adata.obsm:
        print("[WARN] 未找到UMAP坐标，跳过")
        return
    
    if "cluster_k20" not in adata.obs.columns:
        print("[WARN] 未找到cluster_k20，跳过")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='cluster_k20', ax=ax, show=False, 
               legend_loc='right margin', title='T cell clusters (KMeans k=20)',
               frameon=False, size=1 if adata.n_obs > 50000 else 10)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存UMAP cluster图: {out_file}")
    plt.close()


def plot_umap_tscm(adata, out_file):
    """绘制UMAP按TSCM label着色"""
    print(f"[INFO] 绘制UMAP TSCM图: {out_file}")
    
    if "X_umap" not in adata.obsm:
        print("[WARN] 未找到UMAP坐标，跳过")
        return
    
    if "tscm_label" not in adata.obs.columns:
        print("[WARN] 未找到tscm_label，跳过")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # 左图：TSCM label
    sc.pl.umap(adata, color='tscm_label', ax=axes[0], show=False,
               legend_loc='right margin', title='TSCM Label',
               frameon=False, size=1 if adata.n_obs > 50000 else 10)
    
    # 右图：TSCM score
    if "score_tscm" in adata.obs.columns:
        sc.pl.umap(adata, color='score_tscm', ax=axes[1], show=False,
                   cmap='Reds', title='TSCM Score',
                   frameon=False, size=1 if adata.n_obs > 50000 else 10)
    
    plt.suptitle("TSCM identification in UMAP space", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存UMAP TSCM图: {out_file}")
    plt.close()


def plot_tscm_scores(adata, out_file):
    """绘制TSCM打分分布"""
    print(f"[INFO] 绘制TSCM打分分布图: {out_file}")
    
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if len(score_cols) == 0:
        print("[WARN] 未找到score列，跳过")
        return
    
    n_scores = len(score_cols)
    n_cols = 3
    n_rows = (n_scores + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_scores > 1 else [axes]
    
    for i, score_col in enumerate(score_cols):
        ax = axes[i]
        
        # 按TSCM label分组绘制
        if "tscm_label" in adata.obs.columns:
            for label in ["TSCM_high", "T_other"]:
                mask = adata.obs["tscm_label"] == label
                if mask.sum() > 0:
                    data = adata.obs.loc[mask, score_col]
                    ax.hist(data, bins=50, alpha=0.6, label=label, density=True)
            ax.legend()
        else:
            ax.hist(adata.obs[score_col], bins=50, alpha=0.7, density=True)
        
        ax.set_xlabel(score_col.replace("score_", "").upper(), fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.set_title(score_col.replace("score_", "").upper(), fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.3)
    
    # 隐藏多余的子图
    for i in range(len(score_cols), len(axes)):
        axes[i].axis('off')
    
    plt.suptitle("TSCM and related scores distribution", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存TSCM打分分布图: {out_file}")
    plt.close()


def plot_cluster_scores_heatmap(scores_file, out_file):
    """绘制cluster打分热图"""
    print(f"[INFO] 绘制cluster打分热图: {out_file}")
    
    if not scores_file.exists():
        print(f"[WARN] 未找到打分文件: {scores_file}，跳过")
        return
    
    df = pd.read_csv(scores_file, sep="\t", index_col=0)
    
    # 选择打分列
    score_cols = [c for c in df.columns if c.startswith("score_")]
    if len(score_cols) == 0:
        print("[WARN] 未找到score列，跳过")
        return
    
    # 准备热图数据
    heatmap_data = df[score_cols].T
    
    fig, ax = plt.subplots(figsize=(max(10, len(df)*0.5), max(6, len(score_cols)*0.8)))
    
    sns.heatmap(heatmap_data, annot=False, fmt=".2f", cmap="RdYlBu_r", 
                center=0, vmin=-1, vmax=1, cbar_kws={"label": "Score"},
                ax=ax, linewidths=0.5)
    
    ax.set_title("Cluster scores heatmap", fontsize=14, fontweight="bold")
    ax.set_xlabel("Cluster", fontsize=12)
    ax.set_ylabel("Score type", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存cluster打分热图: {out_file}")
    plt.close()


def plot_marker_genes(adata, out_file):
    """绘制关键marker基因表达"""
    print(f"[INFO] 绘制关键marker基因表达图: {out_file}")
    
    # TSCM关键基因
    key_genes = ["CCR7", "LEF1", "TCF7", "IL7R", "SELL", "CD27", "CD28", 
                 "GZMB", "PRF1", "IFNG", "CD3D", "CD3E"]
    
    # 检查哪些基因在数据中
    available_genes = [g for g in key_genes if g in adata.var_names]
    if len(available_genes) == 0:
        print("[WARN] 未找到任何关键基因，跳过")
        return
    
    print(f"[INFO] 绘制 {len(available_genes)} 个关键基因的表达")
    
    n_genes = len(available_genes)
    n_cols = 4
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    for i, gene in enumerate(available_genes):
        ax = axes[i]
        sc.pl.umap(adata, color=gene, ax=ax, show=False, cmap='Reds', 
                   title=gene, frameon=False, size=1 if adata.n_obs > 50000 else 10)
    
    # 隐藏多余的子图
    for i in range(len(available_genes), len(axes)):
        axes[i].axis('off')
    
    plt.suptitle("Key TSCM marker genes expression in UMAP", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存marker基因表达图: {out_file}")
    plt.close()


def plot_cluster_composition(adata, out_file):
    """绘制cluster组成图"""
    print(f"[INFO] 绘制cluster组成图: {out_file}")
    
    if "cluster_k20" not in adata.obs.columns or "tscm_label" not in adata.obs.columns:
        print("[WARN] 缺少必要列，跳过")
        return
    
    # 计算每个cluster的TSCM比例
    composition = adata.obs.groupby("cluster_k20")["tscm_label"].value_counts(normalize=True).unstack(fill_value=0)
    composition = composition * 100  # 转换为百分比
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # 左图：堆叠条形图
    composition.plot(kind='bar', stacked=True, ax=axes[0], 
                     color=['#FF6B6B', '#4ECDC4'], width=0.8)
    axes[0].set_xlabel("Cluster", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Percentage (%)", fontsize=12, fontweight="bold")
    axes[0].set_title("TSCM composition by cluster", fontsize=13, fontweight="bold")
    axes[0].legend(title="TSCM Label", fontsize=10)
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].grid(True, axis='y', alpha=0.3)
    
    # 右图：TSCM数量
    tscm_counts = adata.obs.groupby("cluster_k20")["tscm_label"].value_counts()
    tscm_counts = tscm_counts.unstack(fill_value=0)
    tscm_counts.plot(kind='bar', ax=axes[1], color=['#FF6B6B', '#4ECDC4'], width=0.8)
    axes[1].set_xlabel("Cluster", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("Number of cells", fontsize=12, fontweight="bold")
    axes[1].set_title("TSCM cell counts by cluster", fontsize=13, fontweight="bold")
    axes[1].legend(title="TSCM Label", fontsize=10)
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存cluster组成图: {out_file}")
    plt.close()


def main():
    print("=" * 60)
    print("Stage 4H: Visualization for all Stage 4 analyses")
    print("=" * 60)
    
    # 1. 加载数据
    print(f"\n[INFO] 读取数据: {IN_H5AD}")
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[INFO] 数据维度: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 2. 绘制所有图表
    print("\n[INFO] 开始生成可视化图表...")
    
    plot_qc_violin(adata, FIG_QC_VIOLIN)
    plot_umap_clusters(adata, FIG_UMAP_CLUSTERS)
    plot_umap_tscm(adata, FIG_UMAP_TSCM)
    plot_tscm_scores(adata, FIG_TSCM_SCORES)
    plot_cluster_scores_heatmap(IN_SCORES, FIG_CLUSTER_SCORES_HEATMAP)
    plot_marker_genes(adata, FIG_MARKER_GENES)
    plot_cluster_composition(adata, FIG_CLUSTER_COMPOSITION)
    
    # 3. 总结
    print("\n" + "=" * 60)
    print("[DONE] Stage 4H (Visualization) 完成")
    print("=" * 60)
    print(f"\n生成的图表文件:")
    print(f"  - QC指标小提琴图: {FIG_QC_VIOLIN}")
    print(f"  - UMAP cluster图: {FIG_UMAP_CLUSTERS}")
    print(f"  - UMAP TSCM图: {FIG_UMAP_TSCM}")
    print(f"  - TSCM打分分布: {FIG_TSCM_SCORES}")
    print(f"  - Cluster打分热图: {FIG_CLUSTER_SCORES_HEATMAP}")
    print(f"  - Marker基因表达: {FIG_MARKER_GENES}")
    print(f"  - Cluster组成图: {FIG_CLUSTER_COMPOSITION}")


if __name__ == "__main__":
    main()

