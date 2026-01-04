#!/usr/bin/env python3
"""
Stage 5D: Trajectory analysis for TSCM - PAGA and pseudotime

输入:
  data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad
    - 需要包含: tscm_label, cluster_k20, PCA, neighbors

输出:
  results/tables/tscm_trajectory_paga.tsv          # PAGA连接性结果
  results/tables/tscm_pseudotime.tsv               # 伪时间分析结果
  results/figures/tscm_paga_graph.pdf              # PAGA图
  results/figures/tscm_trajectory_umap.pdf          # 轨迹UMAP图
  results/figures/tscm_pseudotime_umap.pdf          # 伪时间UMAP图
  results/figures/tscm_trajectory_gene_expression.pdf  # 关键基因沿轨迹表达

说明:
  - 使用PAGA分析TSCM在T细胞分化轨迹中的位置
  - 使用diffusion map和pseudotime分析验证TSCM的干性特征
  - 探索TSCM是否是naive → effector fate转换的源头
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
OUT_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 输出文件
OUT_PAGA = OUT_DIR / "tscm_trajectory_paga.tsv"
OUT_PSEUDOTIME = OUT_DIR / "tscm_pseudotime.tsv"
OUT_TRAJECTORY_SUMMARY = OUT_DIR / "tscm_trajectory_summary.tsv"

# 图表文件
FIG_PAGA = FIG_DIR / "tscm_paga_graph.pdf"
FIG_TRAJECTORY_UMAP = FIG_DIR / "tscm_trajectory_umap.pdf"
FIG_PSEUDOTIME_UMAP = FIG_DIR / "tscm_pseudotime_umap.pdf"
FIG_GENE_EXPRESSION = FIG_DIR / "tscm_trajectory_gene_expression.pdf"


def prepare_data(adata):
    """准备数据用于轨迹分析"""
    print("[INFO] 准备数据用于轨迹分析...")
    
    # 确保有PCA
    if "X_pca" not in adata.obsm:
        print("[INFO] 计算PCA...")
        use_hvg = "highly_variable" in adata.var.columns
        sc.pp.pca(adata, n_comps=50, use_highly_variable=use_hvg, svd_solver="arpack")
    
    # 确保有neighbors
    if "neighbors" not in adata.uns:
        print("[INFO] 计算neighbors...")
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    
    # 确保有UMAP
    if "X_umap" not in adata.obsm:
        print("[INFO] 计算UMAP...")
        sc.tl.umap(adata)
    
    # 确保有cluster_k20
    if "cluster_k20" not in adata.obs.columns:
        raise SystemExit("[ERROR] 未找到cluster_k20，请先运行Stage 4D")
    
    # 确保有tscm_label
    if "tscm_label" not in adata.obs.columns:
        raise SystemExit("[ERROR] 未找到tscm_label，请先运行Stage 4F (脚本46)")
    
    print(f"[INFO] 数据准备完成: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def run_paga(adata):
    """运行PAGA分析"""
    print("\n[INFO] 运行PAGA分析...")
    
    # 使用cluster_k20作为分组
    adata.uns['neighbors']['params']['n_neighbors'] = 15
    adata.uns['neighbors']['params']['n_pcs'] = 50
    
    # 运行PAGA
    sc.tl.paga(adata, groups='cluster_k20', use_rna_velocity=False)
    
    print("[INFO] PAGA分析完成")
    
    # 提取PAGA连接性
    if 'paga' in adata.uns:
        paga_connectivities = adata.uns['paga']['connectivities']
        paga_connectivities_tree = adata.uns['paga']['connectivities_tree']
        
        # 转换为DataFrame
        clusters = adata.obs['cluster_k20'].cat.categories if hasattr(adata.obs['cluster_k20'], 'cat') else sorted(adata.obs['cluster_k20'].unique())
        
        # 保存连接性矩阵
        df_connect = pd.DataFrame(
            paga_connectivities.toarray(),
            index=clusters,
            columns=clusters
        )
        df_connect.to_csv(OUT_PAGA, sep="\t")
        print(f"[INFO] PAGA连接性已保存: {OUT_PAGA}")
        
        return df_connect
    else:
        print("[WARN] PAGA结果未找到")
        return None


def run_diffusion_map(adata):
    """运行diffusion map分析"""
    print("\n[INFO] 运行Diffusion Map分析...")
    
    try:
        sc.tl.diffmap(adata, n_comps=10)
        print("[INFO] Diffusion Map分析完成")
        return True
    except Exception as e:
        print(f"[WARN] Diffusion Map分析失败: {e}")
        return False


def run_pseudotime(adata, root_cluster=None):
    """运行伪时间分析"""
    print("\n[INFO] 运行伪时间分析...")
    
    # 如果没有指定root，使用TSCM cluster作为root
    if root_cluster is None:
        # 找到TSCM_high细胞最多的cluster
        tscm_clusters = adata.obs[adata.obs['tscm_label'] == 'TSCM_high']['cluster_k20'].value_counts()
        if len(tscm_clusters) > 0:
            root_cluster = tscm_clusters.index[0]
            print(f"[INFO] 使用TSCM富集的cluster {root_cluster}作为root")
        else:
            print("[WARN] 未找到TSCM cluster，使用cluster 0作为root")
            root_cluster = '0'
    
    # 使用diffusion map计算伪时间
    if 'X_diffmap' in adata.obsm:
        print("[INFO] 使用Diffusion Map计算伪时间...")
        # 找到root cluster的细胞
        root_mask = adata.obs['cluster_k20'] == str(root_cluster)
        root_cells = np.where(root_mask)[0]
        
        if len(root_cells) > 0:
            # 使用第一个diffusion component作为伪时间
            adata.obs['dpt_pseudotime'] = adata.obsm['X_diffmap'][:, 0]
            # 归一化到0-1
            adata.obs['dpt_pseudotime'] = (adata.obs['dpt_pseudotime'] - adata.obs['dpt_pseudotime'].min()) / (adata.obs['dpt_pseudotime'].max() - adata.obs['dpt_pseudotime'].min())
            print("[INFO] 伪时间计算完成")
            return True
        else:
            print(f"[WARN] 未找到root cluster {root_cluster}的细胞")
            return False
    else:
        print("[WARN] 未找到Diffusion Map，跳过伪时间分析")
        return False


def plot_paga_graph(adata, out_file):
    """绘制PAGA图"""
    print(f"[INFO] 绘制PAGA图: {out_file}")
    
    if 'paga' not in adata.uns:
        print("[WARN] 未找到PAGA结果，跳过PAGA图绘制")
        return
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # 绘制PAGA图
    sc.pl.paga(adata, ax=ax, show=False, node_size_scale=2, edge_width_scale=2)
    
    plt.title("PAGA Graph - T cell differentiation trajectory", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存PAGA图: {out_file}")
    plt.close()


def plot_trajectory_umap(adata, out_file):
    """在UMAP上绘制轨迹"""
    print(f"[INFO] 绘制轨迹UMAP图: {out_file}")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # 左图：按cluster着色
    sc.pl.umap(adata, color='cluster_k20', ax=axes[0], show=False, legend_loc='right margin', title='Clusters')
    
    # 右图：按TSCM label着色
    sc.pl.umap(adata, color='tscm_label', ax=axes[1], show=False, legend_loc='right margin', title='TSCM Label')
    
    plt.suptitle("T cell trajectory in UMAP space", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存轨迹UMAP图: {out_file}")
    plt.close()


def plot_pseudotime_umap(adata, out_file):
    """在UMAP上绘制伪时间"""
    print(f"[INFO] 绘制伪时间UMAP图: {out_file}")
    
    if 'dpt_pseudotime' not in adata.obs.columns:
        print("[WARN] 未找到伪时间数据，跳过伪时间UMAP图")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # 左图：伪时间
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0], show=False, cmap='viridis', title='Pseudotime')
    
    # 右图：TSCM label叠加伪时间
    sc.pl.umap(adata, color='tscm_label', ax=axes[1], show=False, legend_loc='right margin', title='TSCM Label')
    
    plt.suptitle("Pseudotime analysis - TSCM as trajectory origin", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存伪时间UMAP图: {out_file}")
    plt.close()


def plot_gene_expression_along_trajectory(adata, out_file):
    """绘制关键基因沿轨迹的表达"""
    print(f"[INFO] 绘制关键基因沿轨迹表达图: {out_file}")
    
    # TSCM关键基因
    key_genes = ["CCR7", "LEF1", "TCF7", "IL7R", "SELL", "GZMB", "PRF1", "IFNG"]
    
    # 检查哪些基因在数据中
    available_genes = [g for g in key_genes if g in adata.var_names]
    if len(available_genes) == 0:
        print("[WARN] 未找到任何关键基因，跳过基因表达图")
        return
    
    print(f"[INFO] 绘制 {len(available_genes)} 个关键基因的表达")
    
    n_genes = len(available_genes)
    n_cols = 4
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    for i, gene in enumerate(available_genes):
        ax = axes[i]
        sc.pl.umap(adata, color=gene, ax=ax, show=False, cmap='Reds', title=gene, frameon=False)
    
    # 隐藏多余的子图
    for i in range(len(available_genes), len(axes)):
        axes[i].axis('off')
    
    plt.suptitle("Key TSCM marker genes expression along trajectory", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存基因表达图: {out_file}")
    plt.close()


def analyze_trajectory_summary(adata):
    """分析轨迹摘要"""
    print("\n[INFO] 分析轨迹摘要...")
    
    summary = []
    
    # 按cluster分析
    for cluster in sorted(adata.obs['cluster_k20'].unique()):
        mask = adata.obs['cluster_k20'] == cluster
        n_cells = int(mask.sum())
        n_tscm = int((mask & (adata.obs['tscm_label'] == 'TSCM_high')).sum())
        pct_tscm = (n_tscm / n_cells * 100) if n_cells > 0 else 0
        
        # 伪时间（如果有）
        if 'dpt_pseudotime' in adata.obs.columns:
            pseudotime_mean = float(adata.obs.loc[mask, 'dpt_pseudotime'].mean())
            pseudotime_std = float(adata.obs.loc[mask, 'dpt_pseudotime'].std())
        else:
            pseudotime_mean = np.nan
            pseudotime_std = np.nan
        
        summary.append({
            'cluster': cluster,
            'n_cells': n_cells,
            'n_tscm': n_tscm,
            'pct_tscm': pct_tscm,
            'pseudotime_mean': pseudotime_mean,
            'pseudotime_std': pseudotime_std
        })
    
    df_summary = pd.DataFrame(summary)
    df_summary = df_summary.sort_values('pct_tscm', ascending=False)
    df_summary.to_csv(OUT_TRAJECTORY_SUMMARY, sep="\t", index=False)
    print(f"[INFO] 轨迹摘要已保存: {OUT_TRAJECTORY_SUMMARY}")
    
    # 显示top TSCM clusters
    print("\n[INFO] Top TSCM富集的clusters:")
    print(df_summary.head(10).to_string(index=False))
    
    return df_summary


def main():
    print("=" * 60)
    print("Stage 5D: Trajectory Analysis for TSCM")
    print("=" * 60)
    
    # 1. 加载数据
    print(f"\n[INFO] 读取数据: {IN_H5AD}")
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[INFO] 初始维度: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 2. 准备数据
    adata = prepare_data(adata)
    
    # 3. 运行PAGA分析
    df_paga = run_paga(adata)
    
    # 4. 运行Diffusion Map
    has_diffmap = run_diffusion_map(adata)
    
    # 5. 运行伪时间分析
    if has_diffmap:
        run_pseudotime(adata)
    
    # 6. 分析轨迹摘要
    df_summary = analyze_trajectory_summary(adata)
    
    # 7. 保存伪时间结果
    if 'dpt_pseudotime' in adata.obs.columns:
        df_pseudotime = adata.obs[['cluster_k20', 'tscm_label', 'dpt_pseudotime']].copy()
        df_pseudotime.to_csv(OUT_PSEUDOTIME, sep="\t", index=True)
        print(f"[INFO] 伪时间结果已保存: {OUT_PSEUDOTIME}")
    
    # 8. 绘制图表
    print("\n[INFO] 绘制可视化图表...")
    plot_paga_graph(adata, FIG_PAGA)
    plot_trajectory_umap(adata, FIG_TRAJECTORY_UMAP)
    plot_pseudotime_umap(adata, FIG_PSEUDOTIME_UMAP)
    plot_gene_expression_along_trajectory(adata, FIG_GENE_EXPRESSION)
    
    # 9. 保存更新后的数据对象（包含轨迹分析结果）
    out_h5ad = ROOT / "data" / "qc" / "GSE137029_sle.tcells.clustered.tscm_labeled.trajectory.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"[INFO] 已保存带轨迹分析结果的对象: {out_h5ad}")
    
    # 10. 总结
    print("\n" + "=" * 60)
    print("[DONE] Stage 5D (Trajectory Analysis) 完成")
    print("=" * 60)
    print(f"\n输出文件:")
    if df_paga is not None:
        print(f"  - PAGA连接性: {OUT_PAGA}")
    if 'dpt_pseudotime' in adata.obs.columns:
        print(f"  - 伪时间结果: {OUT_PSEUDOTIME}")
    print(f"  - 轨迹摘要: {OUT_TRAJECTORY_SUMMARY}")
    print(f"\n图表文件:")
    print(f"  - PAGA图: {FIG_PAGA}")
    print(f"  - 轨迹UMAP: {FIG_TRAJECTORY_UMAP}")
    print(f"  - 伪时间UMAP: {FIG_PSEUDOTIME_UMAP}")
    print(f"  - 基因表达图: {FIG_GENE_EXPRESSION}")
    
    # 关键发现
    print("\n[INFO] 关键发现:")
    if 'dpt_pseudotime' in adata.obs.columns:
        tscm_pseudotime = adata.obs[adata.obs['tscm_label'] == 'TSCM_high']['dpt_pseudotime'].mean()
        other_pseudotime = adata.obs[adata.obs['tscm_label'] == 'T_other']['dpt_pseudotime'].mean()
        print(f"  - TSCM平均伪时间: {tscm_pseudotime:.3f}")
        print(f"  - Other平均伪时间: {other_pseudotime:.3f}")
        if tscm_pseudotime < other_pseudotime:
            print("  - ✅ TSCM位于轨迹起点，支持干性特征")
        else:
            print("  - ⚠️ TSCM不在轨迹起点，需要进一步分析")


if __name__ == "__main__":
    main()

