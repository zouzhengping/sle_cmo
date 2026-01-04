#!/usr/bin/env python3
"""
Stage 5C: Ligand-receptor interaction analysis for TSCM

输入:
  data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad
    - 需要包含: tscm_label (TSCM_high vs T_other)

输出:
  results/tables/tscm_lr_interactions.tsv          # 配体-受体相互作用结果
  results/tables/tscm_lr_network.tsv              # 相互作用网络
  results/figures/tscm_lr_heatmap.pdf             # 相互作用热图
  results/figures/tscm_lr_network.pdf             # 相互作用网络图
  results/figures/tscm_lr_bubble.pdf              # 相互作用气泡图

说明:
  - 分析TSCM_high和T_other之间的潜在配体-受体相互作用
  - 使用CellChat或简化的配体-受体数据库
  - 即使只有T细胞数据，也可以识别TSCM与其他T细胞亚群的相互作用
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 尝试导入CellChat相关库
try:
    import scanpy as sc
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False
    print("[ERROR] scanpy未安装")

try:
    import cellchat
    HAS_CELLCHAT = True
except ImportError:
    HAS_CELLCHAT = False
    print("[WARN] cellchat未安装，将使用简化的配体-受体分析")
    print("      安装命令: pip install cellchat")

ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = ROOT / "data/qc" / "GSE137029_sle.tcells.clustered.tscm_labeled.h5ad"
OUT_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 输出文件
OUT_LR = OUT_DIR / "tscm_lr_interactions.tsv"
OUT_NETWORK = OUT_DIR / "tscm_lr_network.tsv"
OUT_SIGNALS = OUT_DIR / "tscm_lr_signaling_pathways.tsv"

# 图表文件
FIG_HEATMAP = FIG_DIR / "tscm_lr_heatmap.pdf"
FIG_NETWORK = FIG_DIR / "tscm_lr_network.pdf"
FIG_BUBBLE = FIG_DIR / "tscm_lr_bubble.pdf"
FIG_DOT = FIG_DIR / "tscm_lr_dotplot.pdf"


def load_ligand_receptor_db():
    """加载配体-受体数据库"""
    # 使用CellChat的配体-受体数据库（简化版）
    # 如果CellChat不可用，使用手动整理的常见配体-受体对
    
    # 常见的T细胞相关配体-受体对
    lr_pairs = [
        # T细胞激活相关
        ("CD40LG", "CD40"),  # T细胞-B细胞相互作用
        ("CD28", "CD80"), ("CD28", "CD86"),  # T细胞共刺激
        ("CD27", "CD70"),  # T细胞-T细胞相互作用
        ("ICOS", "ICOSLG"),  # T细胞共刺激
        ("TNF", "TNFRSF1A"), ("TNF", "TNFRSF1B"),  # TNF信号
        ("FASLG", "FAS"),  # 凋亡信号
        
        # 细胞因子-受体
        ("IL2", "IL2RA"), ("IL2", "IL2RB"), ("IL2", "IL2RG"),
        ("IL7", "IL7R"),  # IL-7信号（TSCM关键）
        ("IL15", "IL15RA"), ("IL15", "IL2RB"), ("IL15", "IL2RG"),
        ("IL21", "IL21R"),
        ("IFNG", "IFNGR1"), ("IFNG", "IFNGR2"),
        ("IL4", "IL4R"), ("IL13", "IL4R"),
        ("IL10", "IL10RA"), ("IL10", "IL10RB"),
        ("TGFB1", "TGFBR1"), ("TGFB1", "TGFBR2"),
        
        # 趋化因子-受体
        ("CCL3", "CCR1"), ("CCL3", "CCR5"),
        ("CCL4", "CCR5"),
        ("CCL5", "CCR1"), ("CCL5", "CCR3"), ("CCL5", "CCR5"),
        ("CCL19", "CCR7"), ("CCL21", "CCR7"),  # TSCM关键：淋巴结归巢
        ("CXCL9", "CXCR3"), ("CXCL10", "CXCR3"), ("CXCL11", "CXCR3"),
        ("CXCL12", "CXCR4"),
        ("CCL2", "CCR2"), ("CCL7", "CCR2"),
        
        # 粘附分子
        ("ITGB2", "ICAM1"), ("ITGB2", "ICAM2"),
        ("ITGA4", "VCAM1"),
        ("ITGA9", "VCAM1"),
        ("SELL", "SELL"),  # CD62L (L-selectin)
        
        # WNT信号（TSCM关键）
        ("WNT1", "FZD1"), ("WNT3A", "FZD1"),
        ("WNT5A", "FZD2"), ("WNT5A", "FZD5"),
        
        # Notch信号
        ("DLL1", "NOTCH1"), ("DLL4", "NOTCH1"),
        ("JAG1", "NOTCH1"), ("JAG2", "NOTCH1"),
    ]
    
    return pd.DataFrame(lr_pairs, columns=["ligand", "receptor"])


def calculate_lr_interactions(adata, lr_db):
    """计算配体-受体相互作用强度"""
    print("[INFO] 计算配体-受体相互作用...")
    
    # 获取两组细胞的表达矩阵
    tscm_mask = adata.obs["tscm_label"] == "TSCM_high"
    other_mask = adata.obs["tscm_label"] == "T_other"
    
    # 计算平均表达（使用log1p后的数据）
    if "X" in adata.layers:
        tscm_expr = adata[tscm_mask].layers["X"].mean(axis=0).A1 if hasattr(adata[tscm_mask].layers["X"], 'A') else adata[tscm_mask].layers["X"].mean(axis=0)
        other_expr = adata[other_mask].layers["X"].mean(axis=0).A1 if hasattr(adata[other_mask].layers["X"], 'A') else adata[other_mask].layers["X"].mean(axis=0)
    else:
        tscm_expr = adata[tscm_mask].X.mean(axis=0).A1 if hasattr(adata[tscm_mask].X, 'A') else adata[tscm_mask].X.mean(axis=0)
        other_expr = adata[other_mask].X.mean(axis=0).A1 if hasattr(adata[other_mask].X, 'A') else adata[other_mask].X.mean(axis=0)
    
    # 创建基因表达字典
    gene_expr = {}
    for i, gene in enumerate(adata.var_names):
        gene_expr[gene] = {
            "TSCM_high": tscm_expr[i],
            "T_other": other_expr[i]
        }
    
    # 计算相互作用强度
    interactions = []
    
    for _, row in lr_db.iterrows():
        ligand = row["ligand"]
        receptor = row["receptor"]
        
        # 检查基因是否在数据中
        if ligand not in gene_expr or receptor not in gene_expr:
            continue
        
        # 计算相互作用强度（多种模式）
        # 模式1: TSCM作为配体发送者，T_other作为受体接收者
        score_tscm_to_other = (gene_expr[ligand]["TSCM_high"] * 
                              gene_expr[receptor]["T_other"])
        
        # 模式2: T_other作为配体发送者，TSCM作为受体接收者
        score_other_to_tscm = (gene_expr[ligand]["T_other"] * 
                               gene_expr[receptor]["TSCM_high"])
        
        # 模式3: 双向相互作用
        score_bidirectional = (gene_expr[ligand]["TSCM_high"] * 
                              gene_expr[receptor]["TSCM_high"] +
                              gene_expr[ligand]["T_other"] * 
                              gene_expr[receptor]["T_other"])
        
        # 选择最大相互作用强度
        max_score = max(score_tscm_to_other, score_other_to_tscm, score_bidirectional)
        
        interactions.append({
            "ligand": ligand,
            "receptor": receptor,
            "ligand_TSCM": gene_expr[ligand]["TSCM_high"],
            "ligand_other": gene_expr[ligand]["T_other"],
            "receptor_TSCM": gene_expr[receptor]["TSCM_high"],
            "receptor_other": gene_expr[receptor]["T_other"],
            "score_tscm_to_other": score_tscm_to_other,
            "score_other_to_tscm": score_other_to_tscm,
            "score_bidirectional": score_bidirectional,
            "max_score": max_score,
            "interaction_type": "TSCM_to_other" if max_score == score_tscm_to_other else 
                               "other_to_TSCM" if max_score == score_other_to_tscm else 
                               "bidirectional"
        })
    
    df_interactions = pd.DataFrame(interactions)
    
    # 排序并过滤
    df_interactions = df_interactions.sort_values("max_score", ascending=False)
    
    # 只保留有表达的相互作用（至少一方表达 > 0）
    df_interactions = df_interactions[
        (df_interactions["ligand_TSCM"] > 0) | 
        (df_interactions["ligand_other"] > 0) |
        (df_interactions["receptor_TSCM"] > 0) | 
        (df_interactions["receptor_other"] > 0)
    ]
    
    print(f"[INFO] 识别到 {len(df_interactions)} 个配体-受体相互作用")
    
    return df_interactions


def plot_lr_heatmap(df_interactions, out_file, top_n=50):
    """绘制配体-受体相互作用热图"""
    print(f"[INFO] 绘制相互作用热图: {out_file}")
    
    # 选择top N相互作用
    df_top = df_interactions.head(top_n).copy()
    
    # 准备热图数据
    heatmap_data = df_top[["ligand_TSCM", "ligand_other", 
                           "receptor_TSCM", "receptor_other"]].values
    
    fig, ax = plt.subplots(figsize=(12, max(8, len(df_top) * 0.3)))
    
    # 绘制热图
    sns.heatmap(heatmap_data, 
                yticklabels=[f"{row['ligand']} → {row['receptor']}" 
                            for _, row in df_top.iterrows()],
                xticklabels=["Ligand\nTSCM", "Ligand\nOther", 
                            "Receptor\nTSCM", "Receptor\nOther"],
                cmap="YlOrRd", 
                annot=False,
                fmt=".2f",
                cbar_kws={"label": "Expression Level"})
    
    ax.set_title("Top Ligand-Receptor Interactions\n(TSCM_high vs T_other)", 
                fontsize=14, fontweight="bold")
    ax.set_ylabel("Ligand-Receptor Pair", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存热图: {out_file}")
    plt.close()


def plot_lr_bubble(df_interactions, out_file, top_n=30):
    """绘制配体-受体相互作用气泡图"""
    print(f"[INFO] 绘制相互作用气泡图: {out_file}")
    
    # 选择top N相互作用
    df_top = df_interactions.head(top_n).copy()
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(df_top) * 0.25)))
    
    # 准备数据
    y_pos = np.arange(len(df_top))
    x_pos = df_top["max_score"].values
    sizes = df_top["max_score"].values * 10  # 调整气泡大小
    
    # 根据相互作用类型着色
    colors = []
    for itype in df_top["interaction_type"]:
        if itype == "TSCM_to_other":
            colors.append("red")
        elif itype == "other_to_TSCM":
            colors.append("blue")
        else:
            colors.append("purple")
    
    # 绘制气泡
    scatter = ax.scatter(x_pos, y_pos, s=sizes, c=colors, alpha=0.6, edgecolors="black")
    
    # 设置标签
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{row['ligand']} → {row['receptor']}" 
                        for _, row in df_top.iterrows()], fontsize=9)
    ax.set_xlabel("Interaction Score", fontsize=12, fontweight="bold")
    ax.set_title("Top Ligand-Receptor Interactions\n(TSCM_high vs T_other)", 
                fontsize=14, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label='TSCM → Other'),
        Patch(facecolor='blue', label='Other → TSCM'),
        Patch(facecolor='purple', label='Bidirectional')
    ]
    ax.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存气泡图: {out_file}")
    plt.close()


def plot_lr_network(df_interactions, out_file, top_n=20):
    """绘制配体-受体相互作用网络图（简化版）"""
    print(f"[INFO] 绘制相互作用网络图: {out_file}")
    
    # 选择top N相互作用
    df_top = df_interactions.head(top_n).copy()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # 简化的网络可视化：显示配体和受体的连接
    # 这里使用简单的文本布局，更复杂的网络图可以使用networkx
    
    # 提取所有配体和受体
    all_ligands = set(df_top["ligand"].unique())
    all_receptors = set(df_top["receptor"].unique())
    
    # 创建连接列表
    connections = []
    for _, row in df_top.iterrows():
        connections.append({
            "ligand": row["ligand"],
            "receptor": row["receptor"],
            "score": row["max_score"],
            "type": row["interaction_type"]
        })
    
    # 简单的可视化：显示top相互作用
    y_pos = 0
    for i, conn in enumerate(connections[:top_n]):
        ax.text(0.1, 1 - i*0.04, f"{conn['ligand']} → {conn['receptor']}", 
               fontsize=9, va="center")
        ax.text(0.6, 1 - i*0.04, f"Score: {conn['score']:.2f}", 
               fontsize=8, va="center", style="italic")
        ax.text(0.8, 1 - i*0.04, f"({conn['type']})", 
               fontsize=7, va="center", alpha=0.7)
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Top Ligand-Receptor Interactions\n(TSCM_high vs T_other)", 
                fontsize=14, fontweight="bold", pad=20)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存网络图: {out_file}")
    plt.close()


def identify_signaling_pathways(df_interactions):
    """识别信号通路"""
    print("[INFO] 识别信号通路...")
    
    # 定义信号通路和相关的配体-受体对
    pathways = {
        "T cell activation": ["CD28", "CD80", "CD86", "CD27", "CD70", "ICOS", "ICOSLG"],
        "IL-7 signaling": ["IL7", "IL7R"],
        "IL-15 signaling": ["IL15", "IL15RA", "IL2RB", "IL2RG"],
        "IL-2 signaling": ["IL2", "IL2RA", "IL2RB", "IL2RG"],
        "TNF signaling": ["TNF", "TNFRSF1A", "TNFRSF1B"],
        "FAS signaling": ["FASLG", "FAS"],
        "CCR7 homing": ["CCL19", "CCL21", "CCR7"],
        "CXCR3 signaling": ["CXCL9", "CXCL10", "CXCL11", "CXCR3"],
        "WNT signaling": ["WNT1", "WNT3A", "WNT5A", "FZD1", "FZD2", "FZD5"],
        "Notch signaling": ["DLL1", "DLL4", "JAG1", "JAG2", "NOTCH1"],
        "IFN signaling": ["IFNG", "IFNGR1", "IFNGR2"],
    }
    
    pathway_scores = []
    
    for pathway, genes in pathways.items():
        # 查找包含这些基因的相互作用
        pathway_lr = df_interactions[
            (df_interactions["ligand"].isin(genes)) | 
            (df_interactions["receptor"].isin(genes))
        ]
        
        if len(pathway_lr) > 0:
            total_score = pathway_lr["max_score"].sum()
            avg_score = pathway_lr["max_score"].mean()
            n_interactions = len(pathway_lr)
            
            pathway_scores.append({
                "pathway": pathway,
                "n_interactions": n_interactions,
                "total_score": total_score,
                "avg_score": avg_score,
                "top_interaction": pathway_lr.iloc[0]["ligand"] + " → " + pathway_lr.iloc[0]["receptor"],
                "genes_involved": ", ".join(set(pathway_lr["ligand"].tolist() + pathway_lr["receptor"].tolist()))
            })
    
    df_pathways = pd.DataFrame(pathway_scores)
    df_pathways = df_pathways.sort_values("total_score", ascending=False)
    
    print(f"[INFO] 识别到 {len(df_pathways)} 个信号通路")
    
    return df_pathways


def main():
    print("=" * 60)
    print("Stage 5C: Ligand-Receptor Interaction Analysis")
    print("=" * 60)
    
    # 1. 加载数据
    print(f"[INFO] 读取数据: {IN_H5AD}")
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[INFO] 细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    
    if "tscm_label" not in adata.obs.columns:
        raise SystemExit("[ERROR] obs中未发现tscm_label，请先运行Stage 4F (脚本46)")
    
    # 检查tscm_label分布
    print("\n[INFO] TSCM label分布:")
    print(adata.obs["tscm_label"].value_counts())
    
    # 2. 加载配体-受体数据库
    lr_db = load_ligand_receptor_db()
    print(f"\n[INFO] 配体-受体数据库: {len(lr_db)} 对")
    
    # 3. 计算相互作用
    df_interactions = calculate_lr_interactions(adata, lr_db)
    
    # 保存结果
    df_interactions.to_csv(OUT_LR, sep="\t", index=False)
    print(f"[INFO] 相互作用结果已保存: {OUT_LR}")
    
    # 4. 识别信号通路
    df_pathways = identify_signaling_pathways(df_interactions)
    df_pathways.to_csv(OUT_SIGNALS, sep="\t", index=False)
    print(f"[INFO] 信号通路结果已保存: {OUT_SIGNALS}")
    
    # 显示top信号通路
    print("\n[INFO] Top 10 信号通路:")
    print(df_pathways.head(10).to_string(index=False))
    
    # 5. 绘制图表
    print("\n[INFO] 绘制可视化图表...")
    plot_lr_heatmap(df_interactions, FIG_HEATMAP, top_n=50)
    plot_lr_bubble(df_interactions, FIG_BUBBLE, top_n=30)
    plot_lr_network(df_interactions, FIG_NETWORK, top_n=20)
    
    # 6. 总结
    print("\n" + "=" * 60)
    print("[DONE] Stage 5C (Ligand-Receptor Analysis) 完成")
    print("=" * 60)
    print(f"\n输出文件:")
    print(f"  - 相互作用结果: {OUT_LR}")
    print(f"  - 信号通路结果: {OUT_SIGNALS}")
    print(f"  - 热图: {FIG_HEATMAP}")
    print(f"  - 气泡图: {FIG_BUBBLE}")
    print(f"  - 网络图: {FIG_NETWORK}")
    print(f"\n识别到 {len(df_interactions)} 个配体-受体相互作用")
    print(f"识别到 {len(df_pathways)} 个信号通路")


if __name__ == "__main__":
    main()

