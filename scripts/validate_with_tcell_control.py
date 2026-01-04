#!/usr/bin/env python3
"""
使用T细胞作为对照的TSCM Signature验证脚本

改进点：
1. 使用T细胞marker基因来"加权"或"过滤"bulk表达数据
2. 比较SLE和Healthy中T细胞相关的TSCM signature
3. 更接近原始分析（TSCM_high vs T_other都是在T细胞内部）

输入:
  - GSE211700或GSE61635的bulk表达数据
  - TSCM signature（来自GSE137029 TSCM_high vs T_other）

输出:
  - T细胞加权的TSCM signature score
  - SLE vs Healthy的T细胞层面比较
  - 区分比例变化vs细胞内在变化的分析
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 输入文件
SIGNATURE_DE = ROOT / "results/tables/tscm_vs_other_de.full.tsv"

# T细胞marker基因（用于识别和加权T细胞信号）
TCELL_MARKERS = [
    "CD3D", "CD3E", "CD3G",  # CD3复合体
    "CD3Z", "CD247",         # CD3 zeta链
    "TRAC", "TRBC1", "TRBC2",  # T细胞受体
    "CD2", "CD5", "CD7",     # T细胞表面标记
]

# 输出目录
OUT_DIR = ROOT / "results/validation"
FIG_DIR = ROOT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

def load_tscm_signature():
    """加载TSCM signature"""
    print("[步骤1] 加载TSCM Signature...")
    
    if not SIGNATURE_DE.exists():
        raise FileNotFoundError(f"未找到TSCM signature文件: {SIGNATURE_DE}")
    
    df_de = pd.read_csv(SIGNATURE_DE, sep="\t")
    print(f"  加载了 {len(df_de)} 个基因的DE结果")
    
    # 提取显著上调的基因作为TSCM signature
    sig_up = df_de[
        (df_de["logfoldchange"] > 0.5) & 
        (df_de["pvals_adj"] < 0.05)
    ].copy()
    
    tscm_signature = sig_up.nlargest(100, "logfoldchange")["gene"].tolist()
    print(f"  定义了 {len(tscm_signature)} 个TSCM signature基因")
    
    return tscm_signature, df_de

def calculate_tcell_score(expr_matrix, tcell_markers):
    """计算T细胞score（用于加权）"""
    # 找到可用的T细胞marker
    available_markers = [g for g in tcell_markers if g.upper() in expr_matrix.columns.str.upper().values]
    
    if len(available_markers) == 0:
        print(f"  ⚠ 未找到T细胞marker基因，使用所有基因")
        return pd.Series(1.0, index=expr_matrix.index)
    
    # 提取T细胞marker的表达
    marker_expr = []
    for marker in available_markers:
        # 找到匹配的列（不区分大小写）
        matched_cols = expr_matrix.columns[expr_matrix.columns.str.upper() == marker.upper()]
        if len(matched_cols) > 0:
            marker_expr.append(expr_matrix[matched_cols[0]])
    
    if len(marker_expr) == 0:
        return pd.Series(1.0, index=expr_matrix.index)
    
    # 计算平均表达作为T细胞score
    tcell_expr = pd.concat(marker_expr, axis=1).mean(axis=1)
    
    # 标准化到0-1范围（用于加权）
    tcell_score = (tcell_expr - tcell_expr.min()) / (tcell_expr.max() - tcell_expr.min() + 1e-10)
    
    print(f"  使用了 {len(available_markers)}/{len(tcell_markers)} 个T细胞marker")
    print(f"  T细胞score范围: {tcell_score.min():.3f} - {tcell_score.max():.3f}")
    
    return tcell_score

def calculate_weighted_tscm_score(expr_matrix, tscm_signature, tcell_score, method='weighted_mean'):
    """
    计算T细胞加权的TSCM signature score
    
    Parameters:
    -----------
    method : str
        'weighted_mean': 使用T细胞score加权TSCM基因表达
        'tcell_filtered': 只使用T细胞score高的样本
        'both': 同时计算两种方法
    """
    print(f"\n[步骤3] 计算T细胞加权的TSCM Signature Score (方法: {method})...")
    
    # 处理基因名匹配
    expr_cols_upper = expr_matrix.columns.str.upper()
    tscm_signature_upper = [g.upper() for g in tscm_signature]
    
    matched_genes = [g for g in tscm_signature_upper if g in expr_cols_upper.values]
    print(f"  匹配到 {len(matched_genes)}/{len(tscm_signature)} 个signature基因")
    
    if len(matched_genes) == 0:
        raise ValueError("没有可用的signature基因")
    
    # 提取signature基因的表达
    signature_expr = []
    for sig_gene in matched_genes:
        matched_cols = expr_matrix.columns[expr_matrix.columns.str.upper() == sig_gene]
        if len(matched_cols) > 0:
            signature_expr.append(expr_matrix[matched_cols[0]])
    
    if len(signature_expr) == 0:
        raise ValueError("无法提取signature基因表达")
    
    df_sig_expr = pd.concat(signature_expr, axis=1)
    df_sig_expr.columns = matched_genes[:len(df_sig_expr.columns)]
    
    if method == 'weighted_mean':
        # 方法1: 使用T细胞score加权每个基因的表达
        weighted_expr = df_sig_expr.multiply(tcell_score, axis=0)
        tscm_score = weighted_expr.mean(axis=1)
        
    elif method == 'tcell_filtered':
        # 方法2: 只使用T细胞score高的样本（top 50%）
        threshold = tcell_score.quantile(0.5)
        high_tcell_mask = tcell_score >= threshold
        tscm_score = df_sig_expr[high_tcell_mask].mean(axis=1)
        # 对于低T细胞的样本，设为NaN
        tscm_score = tscm_score.reindex(expr_matrix.index)
        
    elif method == 'both':
        # 同时计算两种方法
        weighted_expr = df_sig_expr.multiply(tcell_score, axis=0)
        score_weighted = weighted_expr.mean(axis=1)
        
        threshold = tcell_score.quantile(0.5)
        high_tcell_mask = tcell_score >= threshold
        score_filtered = df_sig_expr[high_tcell_mask].mean(axis=1)
        score_filtered = score_filtered.reindex(expr_matrix.index)
        
        return {
            'weighted': score_weighted,
            'filtered': score_filtered,
            'tcell_score': tcell_score
        }
    
    return {
        method: tscm_score,
        'tcell_score': tcell_score
    }

def validate_dataset(dataset_id, expr_file, metadata_file):
    """验证单个数据集"""
    print("=" * 80)
    print(f"{dataset_id} TSCM Signature 验证（T细胞对照）")
    print("=" * 80)
    
    # 加载数据
    print(f"\n[步骤2] 加载{dataset_id}数据...")
    df_expr = pd.read_csv(expr_file, index_col=0, sep="\t")
    print(f"  表达矩阵: {df_expr.shape[0]} 个样本 x {df_expr.shape[1]} 个基因")
    
    df_meta = pd.read_csv(metadata_file, sep="\t")
    print(f"  样本metadata: {len(df_meta)} 个样本")
    
    # 匹配样本和metadata
    condition_map = {}
    for _, row in df_meta.iterrows():
        sample_id = str(row.get("sample_id", "")).strip('"')
        geo_acc = str(row.get("geo_accession", "")).strip('"')
        condition = str(row.get("condition", "Unknown")).strip('"')
        
        for idx in df_expr.index:
            idx_str = str(idx).strip('"')
            if sample_id and sample_id in idx_str:
                condition_map[idx] = condition
            elif geo_acc and geo_acc in idx_str:
                condition_map[idx] = condition
    
    # 如果匹配失败，从样本名推断
    if len(condition_map) == 0 or all(v == "Unknown" for v in condition_map.values()):
        print("  ⚠ 无法通过metadata匹配，从样本名推断condition...")
        for idx in df_expr.index:
            idx_str = str(idx).upper().strip('"')
            if "CTRL" in idx_str or "CONTROL" in idx_str or "HEALTHY" in idx_str:
                condition_map[idx] = "Healthy"
            elif "SLE" in idx_str or "LUPUS" in idx_str:
                condition_map[idx] = "SLE"
            elif "LN" in idx_str:
                condition_map[idx] = "SLE"  # LN也是SLE的一种
            else:
                condition_map[idx] = "Unknown"
    
    df_meta_expr = pd.DataFrame({
        'sample_id': df_expr.index,
        'condition': [condition_map.get(idx, "Unknown") for idx in df_expr.index]
    }, index=df_expr.index)
    
    print(f"  Condition分布: {df_meta_expr['condition'].value_counts().to_dict()}")
    
    # 加载TSCM signature
    tscm_signature, df_de = load_tscm_signature()
    
    # 计算T细胞score
    print(f"\n[步骤2.5] 计算T细胞score...")
    tcell_score = calculate_tcell_score(df_expr, TCELL_MARKERS)
    df_meta_expr['tcell_score'] = tcell_score.values
    
    # 计算T细胞加权的TSCM score
    results = calculate_weighted_tscm_score(
        df_expr, tscm_signature, tcell_score, method='both'
    )
    
    df_meta_expr['tscm_score_weighted'] = results['weighted'].values
    df_meta_expr['tscm_score_filtered'] = results['filtered'].values
    
    # 统计分析
    print(f"\n[步骤4] 统计分析...")
    
    conditions = df_meta_expr['condition'].unique()
    stats_results = {}
    
    for method in ['weighted', 'filtered']:
        score_col = f'tscm_score_{method}'
        print(f"\n  方法: {method}")
        
        for cond in conditions:
            if cond == "Unknown":
                continue
            
            scores = df_meta_expr[df_meta_expr['condition'] == cond][score_col].dropna()
            if len(scores) > 0:
                stats_results[f"{cond}_{method}"] = {
                    'n': len(scores),
                    'mean': scores.mean(),
                    'std': scores.std(),
                    'median': scores.median()
                }
                print(f"    {cond}: n={len(scores)}, mean={scores.mean():.3f}±{scores.std():.3f}")
        
        # 统计检验
        if 'SLE' in conditions and 'Healthy' in conditions:
            scores_sle = df_meta_expr[df_meta_expr['condition'] == 'SLE'][score_col].dropna()
            scores_healthy = df_meta_expr[df_meta_expr['condition'] == 'Healthy'][score_col].dropna()
            
            if len(scores_sle) > 0 and len(scores_healthy) > 0:
                t_stat, t_pval = stats.ttest_ind(scores_sle, scores_healthy)
                u_stat, u_pval = stats.mannwhitneyu(scores_sle, scores_healthy, alternative='two-sided')
                fc = scores_sle.mean() / scores_healthy.mean() if scores_healthy.mean() > 0 else np.nan
                
                stats_results[f'comparison_{method}'] = {
                    't_test_pvalue': t_pval,
                    'mannwhitney_pvalue': u_pval,
                    'fold_change': fc
                }
                
                print(f"    SLE vs Healthy: p={t_pval:.4f}, FC={fc:.3f}")
    
    # 可视化
    print(f"\n[步骤5] 生成可视化...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'{dataset_id} TSCM Signature Validation (T-cell Control)', fontsize=16, fontweight='bold')
    
    # 1. T细胞score分布
    ax1 = axes[0, 0]
    for cond in conditions:
        if cond != "Unknown":
            scores = df_meta_expr[df_meta_expr['condition'] == cond]['tcell_score']
            ax1.hist(scores, alpha=0.6, label=cond, bins=20)
    ax1.set_xlabel('T-cell Score', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('T-cell Score Distribution', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. 加权TSCM score箱线图
    ax2 = axes[0, 1]
    plot_data = []
    plot_labels = []
    for cond in conditions:
        if cond != "Unknown":
            scores = df_meta_expr[df_meta_expr['condition'] == cond]['tscm_score_weighted'].dropna()
            if len(scores) > 0:
                plot_data.append(scores.values)
                plot_labels.append(cond)
    
    if plot_data:
        bp = ax2.boxplot(plot_data, labels=plot_labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax2.set_ylabel('TSCM Score (T-cell Weighted)', fontsize=12)
        ax2.set_title('TSCM Score by Condition (Weighted)', fontsize=13, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
    
    # 3. 过滤后TSCM score箱线图
    ax3 = axes[1, 0]
    plot_data = []
    plot_labels = []
    for cond in conditions:
        if cond != "Unknown":
            scores = df_meta_expr[df_meta_expr['condition'] == cond]['tscm_score_filtered'].dropna()
            if len(scores) > 0:
                plot_data.append(scores.values)
                plot_labels.append(cond)
    
    if plot_data:
        bp = ax3.boxplot(plot_data, labels=plot_labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightgreen')
        ax3.set_ylabel('TSCM Score (T-cell Filtered)', fontsize=12)
        ax3.set_title('TSCM Score by Condition (Filtered)', fontsize=13, fontweight='bold')
        ax3.grid(axis='y', alpha=0.3)
    
    # 4. T细胞score vs TSCM score散点图
    ax4 = axes[1, 1]
    colors = {'SLE': 'red', 'Healthy': 'blue'}
    for cond in conditions:
        if cond != "Unknown":
            mask = df_meta_expr['condition'] == cond
            ax4.scatter(
                df_meta_expr[mask]['tcell_score'],
                df_meta_expr[mask]['tscm_score_weighted'],
                label=cond, alpha=0.6, s=50, c=colors.get(cond, 'gray')
            )
    ax4.set_xlabel('T-cell Score', fontsize=12)
    ax4.set_ylabel('TSCM Score (Weighted)', fontsize=12)
    ax4.set_title('T-cell Score vs TSCM Score', fontsize=13, fontweight='bold')
    ax4.legend()
    ax4.grid(alpha=0.3)
    
    plt.tight_layout()
    out_fig = FIG_DIR / f"{dataset_id}_tscm_validation_tcell_control.pdf"
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    print(f"  ✓ 已保存图片: {out_fig}")
    
    # 保存结果
    out_scores = OUT_DIR / f"{dataset_id}_tscm_scores_tcell_control.tsv"
    df_meta_expr.to_csv(out_scores, sep="\t", index=True)
    print(f"  ✓ 已保存scores: {out_scores}")
    
    # 生成报告
    print(f"\n[步骤6] 生成验证报告...")
    report_content = f"""# {dataset_id} TSCM Signature 验证报告（T细胞对照）

## 执行信息
- **执行时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
- **数据集**: {dataset_id}
- **验证方法**: TSCM Signature Score with T-cell Control
- **改进**: 使用T细胞marker基因加权，更接近原始分析（TSCM_high vs T_other都在T细胞内部）

---

## 一、方法改进

### 原始验证方法的问题
- 比较SLE vs Healthy的整个PBMC
- 包含B细胞、NK细胞、单核细胞等其他细胞类型
- 无法区分TSCM比例变化vs TSCM细胞内在变化

### 新方法
1. **T细胞score计算**: 使用T细胞marker基因（CD3D, CD3E, CD3G等）计算每个样本的T细胞score
2. **加权方法**: 使用T细胞score加权TSCM signature基因的表达
3. **过滤方法**: 只使用T细胞score高的样本（top 50%）进行分析

---

## 二、T细胞Score分析

### T细胞Marker基因
{', '.join(TCELL_MARKERS)}

### T细胞Score统计
"""
    
    for cond in conditions:
        if cond != "Unknown":
            scores = df_meta_expr[df_meta_expr['condition'] == cond]['tcell_score']
            report_content += f"""
#### {cond}
- 样本数: {len(scores)}
- 均值: {scores.mean():.3f} ± {scores.std():.3f}
- 范围: {scores.min():.3f} - {scores.max():.3f}
"""
    
    report_content += f"""
---

## 三、TSCM Signature Score结果

### 方法1: T细胞加权（Weighted）
使用T细胞score加权每个TSCM signature基因的表达

"""
    
    if 'SLE_weighted' in stats_results and 'Healthy_weighted' in stats_results:
        sle_stats = stats_results['SLE_weighted']
        healthy_stats = stats_results['Healthy_weighted']
        comp = stats_results.get('comparison_weighted', {})
        
        report_content += f"""
#### SLE
- 样本数: {sle_stats['n']}
- 均值: {sle_stats['mean']:.3f} ± {sle_stats['std']:.3f}

#### Healthy
- 样本数: {healthy_stats['n']}
- 均值: {healthy_stats['mean']:.3f} ± {healthy_stats['std']:.3f}

#### 比较
"""
        
        fc_val = comp.get('fold_change', np.nan)
        pval = comp.get('t_test_pvalue', np.nan)
        
        if isinstance(fc_val, (int, float)) and not np.isnan(fc_val):
            report_content += f"- Fold Change: {fc_val:.3f}\n"
        else:
            report_content += "- Fold Change: N/A\n"
        
        if isinstance(pval, (int, float)) and not np.isnan(pval):
            report_content += f"- t-test p值: {pval:.4f}\n"
        else:
            report_content += "- t-test p值: N/A\n"
    
    report_content += f"""
### 方法2: T细胞过滤（Filtered）
只使用T细胞score高的样本（top 50%）

"""
    
    if 'SLE_filtered' in stats_results and 'Healthy_filtered' in stats_results:
        sle_stats = stats_results['SLE_filtered']
        healthy_stats = stats_results['Healthy_filtered']
        comp = stats_results.get('comparison_filtered', {})
        
        report_content += f"""
#### SLE
- 样本数: {sle_stats['n']}
- 均值: {sle_stats['mean']:.3f} ± {sle_stats['std']:.3f}

#### Healthy
- 样本数: {healthy_stats['n']}
- 均值: {healthy_stats['mean']:.3f} ± {healthy_stats['std']:.3f}

#### 比较
"""
        
        fc_val = comp.get('fold_change', np.nan)
        pval = comp.get('t_test_pvalue', np.nan)
        
        if isinstance(fc_val, (int, float)) and not np.isnan(fc_val):
            report_content += f"- Fold Change: {fc_val:.3f}\n"
        else:
            report_content += "- Fold Change: N/A\n"
        
        if isinstance(pval, (int, float)) and not np.isnan(pval):
            report_content += f"- t-test p值: {pval:.4f}\n"
        else:
            report_content += "- t-test p值: N/A\n"
    
    report_content += f"""
---

## 四、结论

通过使用T细胞作为对照，我们可以：
1. **更准确地比较T细胞层面的TSCM signature**
2. **减少其他细胞类型的干扰**
3. **更接近原始分析**（TSCM_high vs T_other都在T细胞内部）

### 关键发现
"""
    
    if 'comparison_weighted' in stats_results:
        comp = stats_results['comparison_weighted']
        if comp.get('fold_change', 0) > 1 and comp.get('t_test_pvalue', 1) < 0.05:
            report_content += "- ✅ **与原始发现一致**: SLE患者的T细胞中TSCM signature显著高于Healthy对照\n"
        elif comp.get('fold_change', 0) < 1 and comp.get('t_test_pvalue', 1) < 0.05:
            report_content += "- ⚠️ **结果不一致**: SLE患者的T细胞中TSCM signature显著低于Healthy对照\n"
        else:
            report_content += "- ℹ️ **无显著差异**: SLE和Healthy之间TSCM signature无显著差异\n"
    
    report_content += f"""
---

## 五、输出文件

- **Score数据**: `{out_scores}`
- **可视化图片**: `{out_fig}`

---

*报告生成时间: {pd.Timestamp.now()}*
"""
    
    out_report = OUT_DIR / f"{dataset_id}_tscm_validation_tcell_control_report.md"
    with open(out_report, 'w', encoding='utf-8') as f:
        f.write(report_content)
    
    print(f"  ✓ 已生成验证报告: {out_report}")
    
    print("\n" + "=" * 80)
    print("验证完成！")
    print("=" * 80)
    
    return df_meta_expr, stats_results

def main():
    """主函数：验证GSE211700和GSE61635"""
    
    datasets = [
        {
            'id': 'GSE211700',
            'expr_file': ROOT / 'data/raw/GSE211700/GSE211700.txt_expression_matrix.tsv',
            'metadata_file': ROOT / 'data/raw/GSE211700/GSE211700.txt_metadata.tsv'
        },
        {
            'id': 'GSE61635',
            'expr_file': ROOT / 'data/raw/GSE61635/GSE61635_expression_matrix_gene_symbol.tsv',
            'metadata_file': ROOT / 'data/raw/GSE61635/GSE61635.txt_metadata.tsv'
        }
    ]
    
    all_results = {}
    
    for dataset in datasets:
        if dataset['expr_file'].exists() and dataset['metadata_file'].exists():
            try:
                df_results, stats = validate_dataset(
                    dataset['id'],
                    dataset['expr_file'],
                    dataset['metadata_file']
                )
                all_results[dataset['id']] = {'data': df_results, 'stats': stats}
            except Exception as e:
                print(f"\n✗ {dataset['id']}验证失败: {e}")
                import traceback
                traceback.print_exc()
        else:
            print(f"\n⚠ {dataset['id']}数据文件不存在，跳过")
    
    # 生成综合报告
    if len(all_results) > 0:
        print("\n" + "=" * 80)
        print("综合验证结果")
        print("=" * 80)
        
        for dataset_id, results in all_results.items():
            print(f"\n{dataset_id}:")
            stats = results['stats']
            if 'comparison_weighted' in stats:
                comp = stats['comparison_weighted']
                print(f"  T细胞加权方法: FC={comp.get('fold_change', 'N/A'):.3f}, p={comp.get('t_test_pvalue', 'N/A'):.4f}")

if __name__ == "__main__":
    main()

