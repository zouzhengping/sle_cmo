#!/usr/bin/env python3
"""
TSCM细胞比例Deconvolution分析

使用单细胞参考数据（GSE137029）作为signature矩阵，
对bulk数据（GSE211700, GSE61635）进行deconvolution，
估计TSCM细胞在PBMC中的比例。

方法：
1. 从GSE137029单细胞数据构建TSCM和其他细胞类型的signature矩阵
2. 使用CIBERSORTx或其他deconvolution方法估计比例
3. 比较SLE vs Healthy中TSCM细胞比例的变化
4. 区分比例变化vs细胞内在变化
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
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
SC_REF_H5AD = ROOT / "data/qc/GSE137029_sle.tcells.clustered.tscm_labeled.h5ad"

# 输出目录
OUT_DIR = ROOT / "results/deconvolution"
FIG_DIR = ROOT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TSCM细胞比例Deconvolution分析")
print("=" * 80)

# ============================================================================
# 步骤1: 从单细胞数据构建signature矩阵
# ============================================================================
print("\n[步骤1] 从单细胞数据构建signature矩阵...")

if not SC_REF_H5AD.exists():
    raise FileNotFoundError(f"未找到单细胞参考数据: {SC_REF_H5AD}")

print(f"  读取单细胞参考数据: {SC_REF_H5AD}")
adata_sc = sc.read_h5ad(SC_REF_H5AD)
print(f"  单细胞数据: {adata_sc.n_obs} 个细胞 x {adata_sc.n_vars} 个基因")

if "tscm_label" not in adata_sc.obs.columns:
    raise ValueError("单细胞数据中未找到tscm_label，请先运行Stage 4F")

# 检查TSCM_high和T_other的细胞数
tscm_counts = adata_sc.obs["tscm_label"].value_counts()
print(f"  细胞类型分布:")
for label, count in tscm_counts.items():
    print(f"    {label}: {count} 个细胞")

# 构建signature矩阵
# 方法1: 使用平均表达作为signature
print("\n[步骤1.1] 计算各细胞类型的平均表达...")

# 获取表达矩阵（使用log1p后的数据或raw counts）
if adata_sc.raw is not None:
    print("  使用raw counts计算signature")
    X = adata_sc.raw.X
    gene_names = adata_sc.raw.var_names
else:
    print("  使用当前X计算signature")
    X = adata_sc.X
    gene_names = adata_sc.var_names

# 转换为dense矩阵（如果必要）
if hasattr(X, 'toarray'):
    X = X.toarray()

# 按细胞类型分组计算平均表达
signature_dict = {}
for label in adata_sc.obs["tscm_label"].unique():
    mask = adata_sc.obs["tscm_label"] == label
    if mask.sum() > 0:
        # 计算平均表达
        mean_expr = np.mean(X[mask], axis=0)
        signature_dict[label] = mean_expr
        print(f"    {label}: {mask.sum()} 个细胞, 平均表达范围: {mean_expr.min():.3f} - {mean_expr.max():.3f}")

# 创建signature矩阵（基因 x 细胞类型）
df_signature = pd.DataFrame(
    signature_dict,
    index=gene_names
)
print(f"  Signature矩阵: {df_signature.shape[0]} 个基因 x {df_signature.shape[1]} 个细胞类型")

# 保存signature矩阵
out_sig_file = OUT_DIR / "tscm_signature_matrix.tsv"
df_signature.to_csv(out_sig_file, sep="\t", index=True)
print(f"  ✓ 已保存signature矩阵: {out_sig_file}")

# ============================================================================
# 步骤2: 简化deconvolution方法（基于相关性）
# ============================================================================
def simple_deconvolution(bulk_expr, signature_matrix, method='correlation'):
    """
    简化的deconvolution方法
    
    方法：
    1. correlation: 计算bulk表达与signature的相关性
    2. regression: 使用线性回归估计比例
    """
    print(f"\n[步骤2] 执行deconvolution (方法: {method})...")
    
    # 匹配基因
    common_genes = bulk_expr.columns.intersection(signature_matrix.index)
    print(f"  共同基因数: {len(common_genes)}/{len(signature_matrix.index)}")
    
    if len(common_genes) < 100:
        raise ValueError("共同基因数太少，无法进行deconvolution")
    
    # 提取共同基因的表达
    bulk_common = bulk_expr[common_genes].T
    sig_common = signature_matrix.loc[common_genes]
    
    if method == 'correlation':
        # 方法1: 使用相关性估计比例
        # 计算每个样本与每个细胞类型signature的相关性
        proportions = []
        for sample in bulk_common.columns:
            sample_expr = bulk_common[sample].values
            corrs = {}
            for cell_type in sig_common.columns:
                sig_expr = sig_common[cell_type].values
                # 计算Pearson相关系数
                corr, _ = stats.pearsonr(sample_expr, sig_expr)
                corrs[cell_type] = corr
            
            # 将相关性转换为比例（归一化）
            corr_values = np.array(list(corrs.values()))
            # 处理负相关（设为0）
            corr_values = np.maximum(corr_values, 0)
            # 归一化
            if corr_values.sum() > 0:
                proportions.append(corr_values / corr_values.sum())
            else:
                proportions.append(np.ones(len(corrs)) / len(corrs))
        
        df_proportions = pd.DataFrame(
            proportions,
            index=bulk_common.columns,
            columns=sig_common.columns
        )
        
    elif method == 'regression':
        # 方法2: 使用非负最小二乘回归
        from scipy.optimize import nnls
        
        proportions = []
        for sample in bulk_common.columns:
            sample_expr = bulk_common[sample].values
            # 使用NNLS估计比例
            coeffs, _ = nnls(sig_common.values, sample_expr)
            # 归一化
            if coeffs.sum() > 0:
                coeffs = coeffs / coeffs.sum()
            proportions.append(coeffs)
        
        df_proportions = pd.DataFrame(
            proportions,
            index=bulk_common.columns,
            columns=sig_common.columns
        )
    
    print(f"  ✓ 估计了 {len(df_proportions)} 个样本的细胞类型比例")
    print(f"  比例范围:")
    for col in df_proportions.columns:
        print(f"    {col}: {df_proportions[col].min():.3f} - {df_proportions[col].max():.3f} (均值: {df_proportions[col].mean():.3f})")
    
    return df_proportions

# ============================================================================
# 步骤3: 对验证数据集进行deconvolution
# ============================================================================
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
    print("\n" + "=" * 80)
    print(f"处理数据集: {dataset['id']}")
    print("=" * 80)
    
    if not dataset['expr_file'].exists():
        print(f"  ⚠ 表达文件不存在，跳过: {dataset['expr_file']}")
        continue
    
    # 加载bulk表达数据
    print(f"\n[步骤3] 加载{dataset['id']}数据...")
    df_bulk = pd.read_csv(dataset['expr_file'], index_col=0, sep="\t")
    print(f"  Bulk表达矩阵: {df_bulk.shape[0]} 个样本 x {df_bulk.shape[1]} 个基因")
    
    # 加载metadata
    if dataset['metadata_file'].exists():
        df_meta = pd.read_csv(dataset['metadata_file'], sep="\t")
        # 匹配condition
        condition_map = {}
        for _, row in df_meta.iterrows():
            sample_id = str(row.get("sample_id", "")).strip('"')
            geo_acc = str(row.get("geo_accession", "")).strip('"')
            condition = str(row.get("condition", "Unknown")).strip('"')
            
            for idx in df_bulk.index:
                idx_str = str(idx).strip('"')
                if sample_id and sample_id in idx_str:
                    condition_map[idx] = condition
                elif geo_acc and geo_acc in idx_str:
                    condition_map[idx] = condition
        
        # 如果匹配失败，从样本名推断
        if len(condition_map) == 0 or all(v == "Unknown" for v in condition_map.values()):
            for idx in df_bulk.index:
                idx_str = str(idx).upper().strip('"')
                if "CTRL" in idx_str or "CONTROL" in idx_str or "HEALTHY" in idx_str:
                    condition_map[idx] = "Healthy"
                elif "SLE" in idx_str or "LUPUS" in idx_str:
                    condition_map[idx] = "SLE"
                elif "LN" in idx_str:
                    condition_map[idx] = "SLE"
                else:
                    condition_map[idx] = "Unknown"
    else:
        condition_map = {idx: "Unknown" for idx in df_bulk.index}
    
    # 执行deconvolution
    try:
        df_proportions = simple_deconvolution(df_bulk, df_signature, method='correlation')
        
        # 添加condition信息
        df_proportions['condition'] = [condition_map.get(idx, "Unknown") for idx in df_proportions.index]
        
        # 保存结果
        out_prop_file = OUT_DIR / f"{dataset['id']}_cell_proportions.tsv"
        df_proportions.to_csv(out_prop_file, sep="\t", index=True)
        print(f"  ✓ 已保存细胞比例: {out_prop_file}")
        
        # 统计分析
        print(f"\n[步骤4] 统计分析{dataset['id']}...")
        
        if 'TSCM_high' in df_proportions.columns:
            tscm_prop = df_proportions['TSCM_high']
            
            for cond in df_proportions['condition'].unique():
                if cond != "Unknown":
                    props = tscm_prop[df_proportions['condition'] == cond]
                    print(f"  {cond}: TSCM比例 = {props.mean():.3f} ± {props.std():.3f} (n={len(props)})")
            
            # 统计检验
            if 'SLE' in df_proportions['condition'].values and 'Healthy' in df_proportions['condition'].values:
                props_sle = tscm_prop[df_proportions['condition'] == 'SLE']
                props_healthy = tscm_prop[df_proportions['condition'] == 'Healthy']
                
                if len(props_sle) > 0 and len(props_healthy) > 0:
                    t_stat, t_pval = stats.ttest_ind(props_sle, props_healthy)
                    u_stat, u_pval = stats.mannwhitneyu(props_sle, props_healthy, alternative='two-sided')
                    fc = props_sle.mean() / props_healthy.mean() if props_healthy.mean() > 0 else np.nan
                    
                    print(f"\n  SLE vs Healthy TSCM比例比较:")
                    print(f"    SLE: {props_sle.mean():.3f} ± {props_sle.std():.3f}")
                    print(f"    Healthy: {props_healthy.mean():.3f} ± {props_healthy.std():.3f}")
                    print(f"    Fold Change: {fc:.3f}")
                    print(f"    t-test p值: {t_pval:.4f}")
                    print(f"    Mann-Whitney U p值: {u_pval:.4f}")
        
        # 可视化
        print(f"\n[步骤5] 生成可视化{dataset['id']}...")
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'{dataset["id"]} TSCM Cell Proportion Deconvolution', fontsize=16, fontweight='bold')
        
        # 1. TSCM比例箱线图
        ax1 = axes[0, 0]
        if 'TSCM_high' in df_proportions.columns:
            plot_data = []
            plot_labels = []
            for cond in df_proportions['condition'].unique():
                if cond != "Unknown":
                    props = df_proportions[df_proportions['condition'] == cond]['TSCM_high']
                    if len(props) > 0:
                        plot_data.append(props.values)
                        plot_labels.append(cond)
            
            if plot_data:
                bp = ax1.boxplot(plot_data, labels=plot_labels, patch_artist=True)
                for patch in bp['boxes']:
                    patch.set_facecolor('lightblue')
                ax1.set_ylabel('TSCM Proportion', fontsize=12)
                ax1.set_title('TSCM Cell Proportion by Condition', fontsize=13, fontweight='bold')
                ax1.grid(axis='y', alpha=0.3)
        
        # 2. 所有细胞类型比例堆叠图
        ax2 = axes[0, 1]
        if len(df_proportions.columns) > 1:
            prop_cols = [col for col in df_proportions.columns if col != 'condition']
            df_plot = df_proportions.sort_values('condition')[prop_cols]
            
            # 按condition分组
            conditions = df_proportions['condition'].unique()
            condition_order = [c for c in ['Healthy', 'SLE'] if c in conditions] + [c for c in conditions if c not in ['Healthy', 'SLE']]
            
            plot_data = []
            for cond in condition_order:
                if cond != "Unknown":
                    mask = df_proportions['condition'] == cond
                    props_mean = df_proportions[mask][prop_cols].mean()
                    plot_data.append(props_mean.values)
            
            if plot_data:
                ax2.bar(range(len(condition_order)), [d[0] if len(d) > 0 else 0 for d in plot_data], 
                       label=prop_cols[0] if len(prop_cols) > 0 else 'TSCM_high')
                bottom = [d[0] if len(d) > 0 else 0 for d in plot_data]
                for i in range(1, len(prop_cols)):
                    ax2.bar(range(len(condition_order)), 
                           [d[i] if len(d) > i else 0 for d in plot_data],
                           bottom=bottom, label=prop_cols[i])
                    bottom = [bottom[j] + (plot_data[j][i] if len(plot_data[j]) > i else 0) 
                             for j in range(len(condition_order))]
                
                ax2.set_xticks(range(len(condition_order)))
                ax2.set_xticklabels(condition_order)
                ax2.set_ylabel('Mean Cell Proportion', fontsize=12)
                ax2.set_title('Cell Type Proportions by Condition', fontsize=13, fontweight='bold')
                ax2.legend()
                ax2.grid(axis='y', alpha=0.3)
        
        # 3. TSCM比例分布
        ax3 = axes[1, 0]
        if 'TSCM_high' in df_proportions.columns:
            for cond in df_proportions['condition'].unique():
                if cond != "Unknown":
                    props = df_proportions[df_proportions['condition'] == cond]['TSCM_high']
                    ax3.hist(props, alpha=0.6, label=cond, bins=20)
            ax3.set_xlabel('TSCM Proportion', fontsize=12)
            ax3.set_ylabel('Frequency', fontsize=12)
            ax3.set_title('TSCM Proportion Distribution', fontsize=13, fontweight='bold')
            ax3.legend()
            ax3.grid(axis='y', alpha=0.3)
        
        # 4. TSCM比例 vs TSCM signature score（如果可用）
        ax4 = axes[1, 1]
        # 尝试加载TSCM signature score
        score_file = ROOT / f"results/validation/{dataset['id']}_tscm_scores_tcell_control.tsv"
        if score_file.exists():
            df_scores = pd.read_csv(score_file, index_col=0, sep="\t")
            if 'tscm_score_weighted' in df_scores.columns:
                # 合并数据
                df_merge = df_proportions.merge(
                    df_scores[['tscm_score_weighted']],
                    left_index=True, right_index=True, how='inner'
                )
                
                colors = {'SLE': 'red', 'Healthy': 'blue'}
                for cond in df_merge['condition'].unique():
                    if cond != "Unknown":
                        mask = df_merge['condition'] == cond
                        ax4.scatter(
                            df_merge[mask]['TSCM_high'],
                            df_merge[mask]['tscm_score_weighted'],
                            label=cond, alpha=0.6, s=50, c=colors.get(cond, 'gray')
                        )
                ax4.set_xlabel('TSCM Proportion', fontsize=12)
                ax4.set_ylabel('TSCM Signature Score (Weighted)', fontsize=12)
                ax4.set_title('TSCM Proportion vs Signature Score', fontsize=13, fontweight='bold')
                ax4.legend()
                ax4.grid(alpha=0.3)
        
        plt.tight_layout()
        out_fig = FIG_DIR / f"{dataset['id']}_deconvolution_tscm_proportion.pdf"
        plt.savefig(out_fig, dpi=300, bbox_inches='tight')
        print(f"  ✓ 已保存图片: {out_fig}")
        
        all_results[dataset['id']] = {
            'proportions': df_proportions,
            'stats': {
                'sle_mean': props_sle.mean() if 'props_sle' in locals() else np.nan,
                'healthy_mean': props_healthy.mean() if 'props_healthy' in locals() else np.nan,
                'fold_change': fc if 'fc' in locals() else np.nan,
                'pvalue': t_pval if 't_pval' in locals() else np.nan
            }
        }
        
    except Exception as e:
        print(f"  ✗ Deconvolution失败: {e}")
        import traceback
        traceback.print_exc()

# ============================================================================
# 步骤6: 生成综合报告
# ============================================================================
print("\n" + "=" * 80)
print("生成综合报告")
print("=" * 80)

report_content = f"""# TSCM细胞比例Deconvolution分析报告

## 执行时间
{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 方法

### Signature矩阵构建
- **来源**: GSE137029单细胞数据
- **细胞类型**: TSCM_high, T_other
- **方法**: 计算各细胞类型的平均表达作为signature

### Deconvolution方法
- **方法**: 基于相关性的简化deconvolution
- **原理**: 计算bulk表达与signature的相关性，归一化后作为比例估计

---

## 结果

"""

for dataset_id, results in all_results.items():
    stats = results['stats']
    report_content += f"""
### {dataset_id}

#### TSCM细胞比例估计
- **SLE**: {stats['sle_mean']:.3f}
- **Healthy**: {stats['healthy_mean']:.3f}
- **Fold Change**: {stats['fold_change']:.3f}
- **p值**: {stats['pvalue']:.4f}

#### 解释
"""
    
    if stats['fold_change'] < 1 and stats['pvalue'] < 0.05:
        report_content += f"- ✅ **SLE患者中TSCM细胞比例显著降低** (p = {stats['pvalue']:.4f})\n"
        report_content += "- 这解释了为什么bulk数据中SLE的TSCM signature score更低\n"
        report_content += "- 支持比例变化是主要因素的假设\n"
    elif stats['fold_change'] > 1 and stats['pvalue'] < 0.05:
        report_content += f"- ✅ **SLE患者中TSCM细胞比例显著升高** (p = {stats['pvalue']:.4f})\n"
    else:
        report_content += f"- ⚠️ **SLE和Healthy之间TSCM细胞比例无显著差异** (p = {stats['pvalue']:.4f})\n"

report_content += f"""
---

## 结论

### 主要发现
1. **Deconvolution方法可以估计TSCM细胞比例**
2. **比例变化vs内在变化的区分**:
   - 如果比例显著降低 → 支持比例变化是主要因素
   - 如果比例无显著差异 → 可能主要是内在变化

### 与之前结果的一致性
- T细胞过滤方法显示差异显著减小（FC接近1）
- Deconvolution方法可以进一步验证比例变化的假设

---

## 输出文件

"""

for dataset_id in all_results.keys():
    report_content += f"""
### {dataset_id}
- 细胞比例: `results/deconvolution/{dataset_id}_cell_proportions.tsv`
- 可视化: `results/figures/{dataset_id}_deconvolution_tscm_proportion.pdf`
"""

report_content += f"""
---

*报告生成时间: {pd.Timestamp.now()}*
"""

out_report = OUT_DIR / "deconvolution_tscm_proportion_report.md"
with open(out_report, 'w', encoding='utf-8') as f:
    f.write(report_content)

print(f"  ✓ 已生成综合报告: {out_report}")

print("\n" + "=" * 80)
print("Deconvolution分析完成！")
print("=" * 80)

