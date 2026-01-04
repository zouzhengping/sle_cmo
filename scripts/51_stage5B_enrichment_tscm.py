#!/usr/bin/env python3
"""
Stage 5B: Pathway enrichment analysis for TSCM signature

输入:
  results/tables/tscm_vs_other_de.full.tsv
    - 包含: gene, logfoldchange, pvals, pvals_adj, scores

输出:
  results/tables/tscm_enrichment_go.tsv          # GO富集结果
  results/tables/tscm_enrichment_kegg.tsv        # KEGG富集结果
  results/tables/tscm_enrichment_reactome.tsv    # Reactome富集结果
  results/figures/tscm_volcano_plot.pdf          # 火山图
  results/figures/tscm_enrichment_barplot.pdf     # 富集条形图
  results/figures/tscm_enrichment_dotplot.pdf     # 富集气泡图
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 尝试导入富集分析库
try:
    import gseapy as gp
    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False
    print("[WARN] gseapy未安装，将使用简化的富集分析。")
    print("      安装命令: pip install gseapy")

try:
    from goatools import GOEnrichmentStudy
    from goatools.obo_parser import GODag
    HAS_GOATOOLS = True
except ImportError:
    HAS_GOATOOLS = False
    print("[WARN] goatools未安装，GO富集分析将使用gseapy。")

ROOT = Path(__file__).resolve().parents[1]
IN_DE = ROOT / "results" / "tables" / "tscm_vs_other_de.full.tsv"
OUT_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 输出文件
OUT_GO = OUT_DIR / "tscm_enrichment_go.tsv"
OUT_KEGG = OUT_DIR / "tscm_enrichment_kegg.tsv"
OUT_REACTOME = OUT_DIR / "tscm_enrichment_reactome.tsv"
OUT_SIG_GENES = OUT_DIR / "tscm_significant_genes.tsv"

# 图表文件
FIG_VOLCANO = FIG_DIR / "tscm_volcano_plot.pdf"
FIG_BARPLOT = FIG_DIR / "tscm_enrichment_barplot.pdf"
FIG_DOTPLOT = FIG_DIR / "tscm_enrichment_dotplot.pdf"


def load_de_results(de_file):
    """加载差异表达结果"""
    print(f"[INFO] 读取差异表达结果: {de_file}")
    df = pd.read_csv(de_file, sep="\t")
    print(f"[INFO] 总基因数: {len(df)}")
    
    # 过滤显著基因
    sig_up = df[(df["pvals_adj"] < 0.05) & (df["logfoldchange"] > 0)]
    sig_down = df[(df["pvals_adj"] < 0.05) & (df["logfoldchange"] < 0)]
    
    print(f"[INFO] 显著上调基因: {len(sig_up)}")
    print(f"[INFO] 显著下调基因: {len(sig_down)}")
    
    return df, sig_up, sig_down


def plot_volcano(df, sig_up, sig_down, out_file):
    """绘制火山图"""
    print(f"[INFO] 绘制火山图: {out_file}")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 非显著基因
    ns = df[(df["pvals_adj"] >= 0.05) | 
            ((df["pvals_adj"] < 0.05) & (df["logfoldchange"].abs() < 0.5))]
    ax.scatter(ns["logfoldchange"], -np.log10(ns["pvals_adj"] + 1e-300), 
               c="gray", alpha=0.5, s=10, label="Not significant")
    
    # 显著下调
    if len(sig_down) > 0:
        ax.scatter(sig_down["logfoldchange"], 
                  -np.log10(sig_down["pvals_adj"] + 1e-300),
                  c="blue", alpha=0.6, s=20, label="Downregulated")
    
    # 显著上调
    if len(sig_up) > 0:
        ax.scatter(sig_up["logfoldchange"], 
                  -np.log10(sig_up["pvals_adj"] + 1e-300),
                  c="red", alpha=0.6, s=20, label="Upregulated")
    
    # 标记关键TSCM marker
    tscm_markers = ["CCR7", "LEF1", "TCF7", "SELL", "IL7R", "CD27", "LRRN3"]
    for marker in tscm_markers:
        if marker in df["gene"].values:
            row = df[df["gene"] == marker].iloc[0]
            ax.scatter(row["logfoldchange"], 
                      -np.log10(row["pvals_adj"] + 1e-300),
                      c="orange", s=100, marker="*", edgecolors="black", linewidths=1)
            ax.annotate(marker, 
                       (row["logfoldchange"], -np.log10(row["pvals_adj"] + 1e-300)),
                       xytext=(5, 5), textcoords="offset points", fontsize=9, fontweight="bold")
    
    # 添加阈值线
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.3)
    ax.axvline(-0.5, color="black", linestyle="--", linewidth=1, alpha=0.3)
    
    ax.set_xlabel("Log2 Fold Change", fontsize=12, fontweight="bold")
    ax.set_ylabel("-Log10 Adjusted P-value", fontsize=12, fontweight="bold")
    ax.set_title("TSCM vs Other T cells - Volcano Plot", fontsize=14, fontweight="bold")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存火山图: {out_file}")
    plt.close()


def enrichment_with_gseapy(gene_list, gene_sets, out_file, name="Enrichment"):
    """使用gseapy进行富集分析"""
    if not HAS_GSEAPY:
        print(f"[WARN] gseapy未安装，跳过{name}富集分析")
        return None
    
    try:
        print(f"[INFO] 进行{name}富集分析 (gseapy)...")
        print(f"      输入基因数: {len(gene_list)}")
        
        # 使用enrichr进行富集分析
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism='Human',
            outdir=None,
            cutoff=0.05
        )
        
        if enr is not None and hasattr(enr, 'results') and len(enr.results) > 0:
            results = enr.results
            results.to_csv(out_file, sep="\t", index=False)
            print(f"[INFO] {name}富集结果已保存: {out_file}")
            print(f"      富集通路数: {len(results)}")
            return results
        else:
            print(f"[WARN] {name}富集分析未找到显著结果")
            return None
            
    except Exception as e:
        print(f"[ERROR] {name}富集分析失败: {e}")
        return None


def simple_enrichment_ora(gene_list, background_genes, gene_sets_dict, out_file, name="Enrichment"):
    """简化的ORA富集分析（当gseapy不可用时）"""
    print(f"[INFO] 进行简化的{name} ORA分析...")
    print(f"      输入基因数: {len(gene_list)}")
    print(f"      背景基因数: {len(background_genes)}")
    
    results = []
    
    for term, genes_in_term in gene_sets_dict.items():
        # 计算重叠
        overlap = set(gene_list) & set(genes_in_term)
        if len(overlap) == 0:
            continue
        
        # 超几何检验
        # M: 背景中属于该term的基因数
        # N: 背景总基因数
        # n: 输入基因数
        # k: 输入基因中属于该term的基因数
        M = len(set(background_genes) & set(genes_in_term))
        N = len(background_genes)
        n = len(gene_list)
        k = len(overlap)
        
        if M == 0 or k == 0:
            continue
        
        # 超几何分布p值
        pvalue = stats.hypergeom.sf(k-1, N, M, n)
        
        # 计算富集倍数
        expected = (M / N) * n
        enrichment_ratio = k / expected if expected > 0 else 0
        
        results.append({
            "Term": term,
            "Overlap": f"{len(overlap)}/{len(genes_in_term)}",
            "P-value": pvalue,
            "Adjusted P-value": pvalue,  # 简化版本，不进行多重检验校正
            "Odds Ratio": enrichment_ratio,
            "Combined Score": -np.log10(pvalue + 1e-300) * enrichment_ratio,
            "Genes": ",".join(sorted(overlap))
        })
    
    if results:
        df_results = pd.DataFrame(results)
        df_results = df_results.sort_values("P-value")
        df_results.to_csv(out_file, sep="\t", index=False)
        print(f"[INFO] {name}富集结果已保存: {out_file}")
        print(f"      富集通路数: {len(df_results)}")
        return df_results
    else:
        print(f"[WARN] {name}富集分析未找到显著结果")
        return None


def plot_enrichment_barplot(enrichment_results, out_file, top_n=20, title="Enrichment Analysis"):
    """绘制富集分析条形图"""
    if enrichment_results is None or len(enrichment_results) == 0:
        print(f"[WARN] 无富集结果，跳过条形图绘制")
        return
    
    print(f"[INFO] 绘制富集条形图: {out_file}")
    
    # 选择top N
    df = enrichment_results.head(top_n).copy()
    
    # 获取p值列名
    pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
    if pval_col not in df.columns:
        pval_col = df.columns[df.columns.str.contains("P", case=False)][0]
    
    # 获取term名称
    term_col = "Term" if "Term" in df.columns else df.columns[0]
    
    # 排序
    df = df.sort_values(pval_col, ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.3)))
    
    # 计算-log10(p)
    df["-log10(p)"] = -np.log10(df[pval_col] + 1e-300)
    
    # 绘制条形图
    colors = plt.cm.Reds(np.linspace(0.4, 0.9, len(df)))
    bars = ax.barh(range(len(df)), df["-log10(p)"], color=colors)
    
    # 设置y轴标签
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df[term_col].str[:60], fontsize=9)  # 截断长名称
    
    ax.set_xlabel("-Log10 P-value", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存条形图: {out_file}")
    plt.close()


def plot_enrichment_dotplot(enrichment_results, out_file, top_n=20, title="Enrichment Analysis"):
    """绘制富集分析气泡图"""
    if enrichment_results is None or len(enrichment_results) == 0:
        print(f"[WARN] 无富集结果，跳过气泡图绘制")
        return
    
    print(f"[INFO] 绘制富集气泡图: {out_file}")
    
    # 选择top N
    df = enrichment_results.head(top_n).copy()
    
    # 获取列名
    pval_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
    if pval_col not in df.columns:
        pval_col = df.columns[df.columns.str.contains("P", case=False)][0]
    
    term_col = "Term" if "Term" in df.columns else df.columns[0]
    
    # 计算基因数（从Overlap列提取）
    if "Overlap" in df.columns:
        df["n_genes"] = df["Overlap"].str.split("/").str[0].astype(int)
    elif "Genes" in df.columns:
        df["n_genes"] = df["Genes"].str.count(",") + 1
    else:
        df["n_genes"] = 1
    
    # 计算-log10(p)
    df["-log10(p)"] = -np.log10(df[pval_col] + 1e-300)
    
    # 排序
    df = df.sort_values("-log10(p)", ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.3)))
    
    # 绘制气泡图
    scatter = ax.scatter(df["-log10(p)"], range(len(df)), 
                        s=df["n_genes"] * 20, 
                        c=df["-log10(p)"], 
                        cmap="Reds", 
                        alpha=0.7, 
                        edgecolors="black", 
                        linewidths=1)
    
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df[term_col].str[:60], fontsize=9)
    ax.set_xlabel("-Log10 P-value", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    
    # 添加颜色条
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("-Log10 P-value", fontsize=10)
    
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"[INFO] 已保存气泡图: {out_file}")
    plt.close()


def main():
    print("=" * 60)
    print("Stage 5B: TSCM Pathway Enrichment Analysis")
    print("=" * 60)
    
    # 1. 加载DE结果
    df_de, sig_up, sig_down = load_de_results(IN_DE)
    
    # 保存显著基因列表
    sig_genes = pd.concat([sig_up, sig_down])
    sig_genes.to_csv(OUT_SIG_GENES, sep="\t", index=False)
    print(f"[INFO] 显著基因列表已保存: {OUT_SIG_GENES}")
    
    # 2. 绘制火山图
    plot_volcano(df_de, sig_up, sig_down, FIG_VOLCANO)
    
    # 3. 准备基因列表（仅显著上调基因，用于富集分析）
    sig_up_genes = sig_up["gene"].tolist()
    background_genes = df_de["gene"].tolist()
    
    print(f"\n[INFO] 准备进行富集分析...")
    print(f"      显著上调基因数: {len(sig_up_genes)}")
    print(f"      背景基因数: {len(background_genes)}")
    
    # 4. 富集分析
    enrichment_results = {}
    
    if HAS_GSEAPY:
        # 使用gseapy进行富集分析
        print("\n[INFO] 使用gseapy进行富集分析...")
        
        # GO富集
        go_results = enrichment_with_gseapy(
            sig_up_genes, 
            ['GO_Biological_Process_2021', 'GO_Molecular_Function_2021', 'GO_Cellular_Component_2021'],
            OUT_GO,
            name="GO"
        )
        enrichment_results['GO'] = go_results
        
        # KEGG富集
        kegg_results = enrichment_with_gseapy(
            sig_up_genes,
            ['KEGG_2021_Human'],
            OUT_KEGG,
            name="KEGG"
        )
        enrichment_results['KEGG'] = kegg_results
        
        # Reactome富集
        reactome_results = enrichment_with_gseapy(
            sig_up_genes,
            ['Reactome_2022'],
            OUT_REACTOME,
            name="Reactome"
        )
        enrichment_results['Reactome'] = reactome_results
        
    else:
        print("\n[WARN] gseapy未安装，跳过自动富集分析")
        print("      建议安装: pip install gseapy")
        print("      或手动使用在线工具（如Enrichr）进行富集分析")
    
    # 5. 绘制富集结果图表
    print("\n[INFO] 绘制富集分析图表...")
    
    for db_name, results in enrichment_results.items():
        if results is not None and hasattr(results, '__len__') and len(results) > 0:
            plot_enrichment_barplot(
                results, 
                FIG_DIR / f"tscm_enrichment_{db_name.lower()}_barplot.pdf",
                top_n=20,
                title=f"TSCM Signature - {db_name} Enrichment"
            )
            plot_enrichment_dotplot(
                results,
                FIG_DIR / f"tscm_enrichment_{db_name.lower()}_dotplot.pdf",
                top_n=20,
                title=f"TSCM Signature - {db_name} Enrichment"
            )
    
    # 6. 总结
    print("\n" + "=" * 60)
    print("[DONE] Stage 5B (Pathway Enrichment) 完成")
    print("=" * 60)
    print(f"\n输出文件:")
    print(f"  - 显著基因列表: {OUT_SIG_GENES}")
    print(f"  - 火山图: {FIG_VOLCANO}")
    if enrichment_results.get('GO'):
        print(f"  - GO富集结果: {OUT_GO}")
    if enrichment_results.get('KEGG'):
        print(f"  - KEGG富集结果: {OUT_KEGG}")
    if enrichment_results.get('Reactome'):
        print(f"  - Reactome富集结果: {OUT_REACTOME}")
    print(f"\n图表文件保存在: {FIG_DIR}")


if __name__ == "__main__":
    main()

