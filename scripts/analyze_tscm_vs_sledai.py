#!/usr/bin/env python3
"""
TSCM Signature与SLEDAI疾病活动指数相关性分析

分析TSCM signature score与SLEDAI的关系，探索：
1. TSCM signature与SLEDAI的相关性
2. 按SLEDAI分组比较TSCM signature
3. TSCM signature与疾病复发/Flare的关系
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

# 输出目录
OUT_DIR = ROOT / "results/validation/sledai"
FIG_DIR = ROOT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TSCM Signature与SLEDAI相关性分析")
print("=" * 80)

# ============================================================================
# 步骤1: 加载TSCM Signature
# ============================================================================
print("\n[步骤1] 加载TSCM Signature...")

if not SIGNATURE_DE.exists():
    raise FileNotFoundError(f"未找到TSCM signature文件: {SIGNATURE_DE}")

df_de = pd.read_csv(SIGNATURE_DE, sep="\t")
sig_up = df_de[
    (df_de["logfoldchange"] > 0.5) & 
    (df_de["pvals_adj"] < 0.05)
].copy()

tscm_signature = sig_up.nlargest(100, "logfoldchange")["gene"].tolist()
print(f"  定义了 {len(tscm_signature)} 个TSCM signature基因")

# ============================================================================
# 步骤2: 检查包含SLEDAI的数据集
# ============================================================================
print("\n[步骤2] 检查包含SLEDAI的数据集...")

sledai_datasets = []

# GSE49454
gse49454_meta = ROOT / "data/metadata/sledai/GSE49454_full_metadata.tsv"
if gse49454_meta.exists():
    df_meta = pd.read_csv(gse49454_meta, sep="\t", index_col=0)
    if "sledai:ch1" in df_meta.columns:
        sledai_datasets.append({
            'id': 'GSE49454',
            'metadata_file': gse49454_meta,
            'sledai_col': 'sledai:ch1',
            'expr_file': None  # 需要下载
        })
        print(f"  ✓ GSE49454: 找到SLEDAI列 (sledai:ch1)")
        print(f"    样本数: {len(df_meta)}")
        # 提取SLEDAI数值范围
        sledai_str = df_meta['sledai:ch1'].dropna().astype(str)
        import re
        sledai_values = []
        for val in sledai_str:
            numbers = re.findall(r'\d+\.?\d*', val)
            if numbers:
                sledai_values.append(float(numbers[0]))
        if sledai_values:
            print(f"    SLEDAI范围: {min(sledai_values):.1f} - {max(sledai_values):.1f}")

# GSE88884
gse88884_meta = ROOT / "data/metadata/sledai/GSE88884_full_metadata.tsv"
if gse88884_meta.exists():
    df_meta = pd.read_csv(gse88884_meta, sep="\t", index_col=0)
    if "sledai_at_baseline:ch1" in df_meta.columns:
        sledai_datasets.append({
            'id': 'GSE88884',
            'metadata_file': gse88884_meta,
            'sledai_col': 'sledai_at_baseline:ch1',
            'expr_file': None  # 需要下载
        })
        print(f"  ✓ GSE88884: 找到SLEDAI列 (sledai_at_baseline:ch1)")
        print(f"    样本数: {len(df_meta)}")

# GSE72747
gse72747_meta = ROOT / "data/metadata/sledai/GSE72747_full_metadata.tsv"
if gse72747_meta.exists():
    df_meta = pd.read_csv(gse72747_meta, sep="\t", index_col=0)
    if "sledai:ch1" in df_meta.columns:
        sledai_datasets.append({
            'id': 'GSE72747',
            'metadata_file': gse72747_meta,
            'sledai_col': 'sledai:ch1',
            'expr_file': None  # 需要下载
        })
        print(f"  ✓ GSE72747: 找到SLEDAI列 (sledai:ch1)")
        print(f"    样本数: {len(df_meta)}")

if len(sledai_datasets) == 0:
    print("  ✗ 未找到包含SLEDAI的数据集")
    sys.exit(0)

# ============================================================================
# 步骤3: 分析每个数据集
# ============================================================================
def extract_sledai_value(sledai_str):
    """从SLEDAI字符串中提取数值"""
    if pd.isna(sledai_str):
        return np.nan
    sledai_str = str(sledai_str)
    # 尝试提取数字
    import re
    numbers = re.findall(r'\d+\.?\d*', sledai_str)
    if numbers:
        return float(numbers[0])
    return np.nan

def analyze_tscm_vs_sledai(dataset_info):
    """分析TSCM signature与SLEDAI的关系"""
    dataset_id = dataset_info['id']
    print(f"\n" + "=" * 80)
    print(f"分析 {dataset_id}")
    print("=" * 80)
    
    # 加载metadata
    df_meta = pd.read_csv(dataset_info['metadata_file'], sep="\t", index_col=0)
    
    # 提取SLEDAI值
    sledai_values = df_meta[dataset_info['sledai_col']].apply(extract_sledai_value)
    df_meta['SLEDAI'] = sledai_values
    
    # 过滤有效SLEDAI值
    valid_mask = df_meta['SLEDAI'].notna()
    df_valid = df_meta[valid_mask].copy()
    
    print(f"\n[步骤3.1] SLEDAI数据概览")
    print(f"  总样本数: {len(df_meta)}")
    print(f"  有效SLEDAI样本数: {len(df_valid)}")
    print(f"  SLEDAI范围: {df_valid['SLEDAI'].min():.1f} - {df_valid['SLEDAI'].max():.1f}")
    print(f"  SLEDAI均值: {df_valid['SLEDAI'].mean():.2f} ± {df_valid['SLEDAI'].std():.2f}")
    
    # 检查是否有表达数据
    # 这里先检查已下载的数据，如果没有则提示需要下载
    expr_files = [
        ROOT / f"data/raw/{dataset_id}/{dataset_id}_expression_matrix_gene_symbol.tsv",
        ROOT / f"data/raw/{dataset_id}/{dataset_id}.txt_expression_matrix.tsv",
        ROOT / f"data/raw/{dataset_id}/{dataset_id}_series_matrix.txt.gz",
    ]
    
    df_expr = None
    for expr_file in expr_files:
        if expr_file.exists():
            if expr_file.suffix == '.gz':
                # 需要解析GEO matrix
                print(f"  发现GEO matrix文件，需要解析...")
                # 这里可以调用之前的解析函数
            else:
                try:
                    df_expr = pd.read_csv(expr_file, index_col=0, sep="\t")
                    print(f"  ✓ 加载表达数据: {df_expr.shape}")
                    break
                except:
                    continue
    
    if df_expr is None:
        print(f"  ⚠ 未找到表达数据，需要先下载")
        print(f"  建议: 运行下载脚本获取{dataset_id}的表达数据")
        return None
    
    # 计算TSCM signature score
    print(f"\n[步骤3.2] 计算TSCM Signature Score...")
    
    # 匹配基因（表达矩阵的行名是基因，列名是样本）
    # 处理gene symbol格式（可能包含 "NM_xxx // GENE // description" 或 "--- // --- // GENE" 格式）
    def extract_gene_symbol(gene_str):
        """从gene symbol字符串中提取实际的基因名"""
        if pd.isna(gene_str):
            return None
        gene_str = str(gene_str).strip()
        # 如果是 "//" 分隔的格式
        if ' // ' in gene_str:
            parts = [p.strip() for p in gene_str.split(' // ')]
            # 优先取第二个部分（通常是gene symbol）
            if len(parts) >= 2 and parts[1] and parts[1] != '---' and parts[1] != 'NULL':
                gene = parts[1]
                # 如果第二个部分是标准基因名（字母数字组合，不含特殊字符）
                if gene and not gene.startswith('LOC') and not gene.startswith('OTTHUM'):
                    return gene.upper()
            # 如果第二个部分不合适，尝试其他部分
            for part in parts:
                part = part.strip()
                if part and part != '---' and part != 'NULL' and not part.startswith('NM_') and not part.startswith('NR_') and not part.startswith('ENST') and not part.startswith('LOC') and not part.startswith('OTTHUM'):
                    # 检查是否是标准基因名格式（字母开头，可能包含数字和连字符）
                    if part and part[0].isalpha() and len(part) <= 20:
                        return part.upper()
        return gene_str.upper()
    
    # 提取表达矩阵中的基因名
    expr_genes = [extract_gene_symbol(g) for g in df_expr.index]
    expr_genes_dict = {g: idx for idx, g in enumerate(expr_genes) if g is not None}
    
    # 匹配signature基因
    tscm_signature_upper = [g.upper() for g in tscm_signature]
    matched_genes = []
    matched_indices = []
    for sig_gene in tscm_signature_upper:
        if sig_gene in expr_genes_dict:
            matched_genes.append(sig_gene)
            matched_indices.append(expr_genes_dict[sig_gene])
    
    if len(matched_genes) == 0:
        print(f"  ✗ 无法匹配signature基因")
        print(f"    尝试匹配的基因数: {len(tscm_signature_upper)}")
        print(f"    表达矩阵中的基因数: {len(expr_genes_dict)}")
        print(f"    示例表达矩阵基因: {list(expr_genes_dict.keys())[:5]}")
        return None
    
    # 提取signature基因表达
    df_sig_expr = df_expr.iloc[matched_indices]
    tscm_scores = df_sig_expr.mean(axis=0)  # 按样本（列）求平均
    
    print(f"  使用了 {len(matched_genes)}/{len(tscm_signature)} 个signature基因")
    print(f"  TSCM Score范围: {tscm_scores.min():.3f} - {tscm_scores.max():.3f}")
    
    # 合并数据（tscm_scores的index是样本ID）
    df_analysis = pd.DataFrame({
        'sample_id': tscm_scores.index,
        'tscm_score': tscm_scores.values
    })
    df_analysis.set_index('sample_id', inplace=True)
    
    # 匹配SLEDAI
    # GSE88884的样本ID格式可能是CEL文件名，需要从metadata中匹配
    df_analysis['SLEDAI'] = np.nan
    
    # 尝试多种匹配方式
    for sample_id in df_analysis.index:
        sample_id_str = str(sample_id)
        # 移除可能的文件扩展名
        if sample_id_str.endswith('.CEL'):
            sample_id_str = sample_id_str[:-4]
        elif sample_id_str.endswith('.CEL.gz'):
            sample_id_str = sample_id_str[:-7]
        
        # 在metadata中查找匹配
        for meta_idx in df_valid.index:
            meta_idx_str = str(meta_idx)
            # 检查是否包含样本ID
            if sample_id_str in meta_idx_str or meta_idx_str in sample_id_str:
                df_analysis.loc[sample_id, 'SLEDAI'] = df_valid.loc[meta_idx, 'SLEDAI']
                break
        
        # 如果还没匹配到，尝试从supplementary_file列匹配
        if pd.isna(df_analysis.loc[sample_id, 'SLEDAI']):
            if 'supplementary_file' in df_valid.columns:
                for meta_idx in df_valid.index:
                    supp_file = str(df_valid.loc[meta_idx, 'supplementary_file'])
                    if sample_id_str in supp_file:
                        df_analysis.loc[sample_id, 'SLEDAI'] = df_valid.loc[meta_idx, 'SLEDAI']
                        break
    
    df_analysis = df_analysis[df_analysis['SLEDAI'].notna()].copy()
    
    if len(df_analysis) == 0:
        print(f"  ✗ 无法匹配SLEDAI数据")
        return None
    
    print(f"  匹配了 {len(df_analysis)} 个样本的SLEDAI数据")
    
    # 相关性分析
    print(f"\n[步骤3.3] TSCM Signature与SLEDAI相关性分析...")
    
    # Pearson相关系数
    pearson_r, pearson_p = stats.pearsonr(df_analysis['tscm_score'], df_analysis['SLEDAI'])
    # Spearman相关系数
    spearman_r, spearman_p = stats.spearmanr(df_analysis['tscm_score'], df_analysis['SLEDAI'])
    
    print(f"  Pearson相关系数: r = {pearson_r:.3f}, p = {pearson_p:.4f}")
    print(f"  Spearman相关系数: ρ = {spearman_r:.3f}, p = {spearman_p:.4f}")
    
    # 按SLEDAI分组
    print(f"\n[步骤3.4] 按SLEDAI分组分析...")
    
    # 定义分组
    df_analysis['SLEDAI_group'] = pd.cut(
        df_analysis['SLEDAI'],
        bins=[0, 4, 8, 12, 100],
        labels=['低活动度 (0-4)', '中活动度 (4-8)', '高活动度 (8-12)', '极高活动度 (>12)']
    )
    
    print(f"  SLEDAI分组分布:")
    for group, count in df_analysis['SLEDAI_group'].value_counts().items():
        if pd.notna(group):
            scores = df_analysis[df_analysis['SLEDAI_group'] == group]['tscm_score']
            print(f"    {group}: n={count}, TSCM score = {scores.mean():.3f} ± {scores.std():.3f}")
    
    # 统计检验（低活动度 vs 高活动度）
    low_activity = df_analysis[df_analysis['SLEDAI'] <= 4]['tscm_score']
    high_activity = df_analysis[df_analysis['SLEDAI'] > 8]['tscm_score']
    
    if len(low_activity) > 0 and len(high_activity) > 0:
        t_stat, t_pval = stats.ttest_ind(low_activity, high_activity)
        print(f"\n  低活动度 (SLEDAI ≤ 4) vs 高活动度 (SLEDAI > 8):")
        print(f"    低活动度: {low_activity.mean():.3f} ± {low_activity.std():.3f} (n={len(low_activity)})")
        print(f"    高活动度: {high_activity.mean():.3f} ± {high_activity.std():.3f} (n={len(high_activity)})")
        print(f"    t-test: p = {t_pval:.4f}")
    
    # 可视化
    print(f"\n[步骤3.5] 生成可视化...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'{dataset_id} TSCM Signature vs SLEDAI', fontsize=16, fontweight='bold')
    
    # 1. 散点图：TSCM score vs SLEDAI
    ax1 = axes[0, 0]
    ax1.scatter(df_analysis['SLEDAI'], df_analysis['tscm_score'], alpha=0.6, s=50)
    
    # 添加回归线
    z = np.polyfit(df_analysis['SLEDAI'], df_analysis['tscm_score'], 1)
    p = np.poly1d(z)
    ax1.plot(df_analysis['SLEDAI'], p(df_analysis['SLEDAI']), "r--", alpha=0.8, linewidth=2)
    
    ax1.set_xlabel('SLEDAI Score', fontsize=12)
    ax1.set_ylabel('TSCM Signature Score', fontsize=12)
    ax1.set_title(f'TSCM Score vs SLEDAI (r={pearson_r:.3f}, p={pearson_p:.4f})', fontsize=13, fontweight='bold')
    ax1.grid(alpha=0.3)
    
    # 2. 按SLEDAI分组箱线图
    ax2 = axes[0, 1]
    plot_data = []
    plot_labels = []
    for group in df_analysis['SLEDAI_group'].cat.categories:
        if pd.notna(group):
            scores = df_analysis[df_analysis['SLEDAI_group'] == group]['tscm_score']
            if len(scores) > 0:
                plot_data.append(scores.values)
                plot_labels.append(str(group))
    
    if plot_data:
        bp = ax2.boxplot(plot_data, labels=plot_labels, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax2.set_ylabel('TSCM Signature Score', fontsize=12)
        ax2.set_title('TSCM Score by SLEDAI Group', fontsize=13, fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        ax2.grid(axis='y', alpha=0.3)
    
    # 3. SLEDAI分布
    ax3 = axes[1, 0]
    ax3.hist(df_analysis['SLEDAI'], bins=20, alpha=0.6, edgecolor='black')
    ax3.set_xlabel('SLEDAI Score', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('SLEDAI Distribution', fontsize=13, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # 4. TSCM score分布
    ax4 = axes[1, 1]
    ax4.hist(df_analysis['tscm_score'], bins=20, alpha=0.6, edgecolor='black')
    ax4.set_xlabel('TSCM Signature Score', fontsize=12)
    ax4.set_ylabel('Frequency', fontsize=12)
    ax4.set_title('TSCM Score Distribution', fontsize=13, fontweight='bold')
    ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    out_fig = FIG_DIR / f"{dataset_id}_tscm_vs_sledai.pdf"
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    print(f"  ✓ 已保存图片: {out_fig}")
    
    # 保存结果
    out_file = OUT_DIR / f"{dataset_id}_tscm_vs_sledai.tsv"
    df_analysis.to_csv(out_file, sep="\t", index=True)
    print(f"  ✓ 已保存结果: {out_file}")
    
    # 生成报告
    report_content = f"""# {dataset_id} TSCM Signature与SLEDAI相关性分析报告

## 执行时间
{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 数据集信息
- **数据集**: {dataset_id}
- **有效样本数**: {len(df_analysis)}
- **SLEDAI范围**: {df_analysis['SLEDAI'].min():.1f} - {df_analysis['SLEDAI'].max():.1f}

---

## 一、相关性分析

### Pearson相关系数
- **r**: {pearson_r:.3f}
- **p值**: {pearson_p:.4f}
- **解释**: {"显著正相关" if pearson_p < 0.05 and pearson_r > 0 else "显著负相关" if pearson_p < 0.05 and pearson_r < 0 else "无显著相关"}

### Spearman相关系数（非参数）
- **ρ**: {spearman_r:.3f}
- **p值**: {spearman_p:.4f}
- **解释**: {"显著正相关" if spearman_p < 0.05 and spearman_r > 0 else "显著负相关" if spearman_p < 0.05 and spearman_r < 0 else "无显著相关"}

---

## 二、按SLEDAI分组分析

### SLEDAI分组
- **低活动度 (0-4)**: SLEDAI ≤ 4
- **中活动度 (4-8)**: 4 < SLEDAI ≤ 8
- **高活动度 (8-12)**: 8 < SLEDAI ≤ 12
- **极高活动度 (>12)**: SLEDAI > 12

### 各组TSCM Signature Score

"""
    
    for group in df_analysis['SLEDAI_group'].cat.categories:
        if pd.notna(group):
            scores = df_analysis[df_analysis['SLEDAI_group'] == group]['tscm_score']
            if len(scores) > 0:
                report_content += f"""
#### {group}
- 样本数: {len(scores)}
- TSCM Score: {scores.mean():.3f} ± {scores.std():.3f}
- 范围: {scores.min():.3f} - {scores.max():.3f}
"""
    
    if len(low_activity) > 0 and len(high_activity) > 0:
        report_content += f"""
### 低活动度 vs 高活动度比较

- **低活动度 (SLEDAI ≤ 4)**: {low_activity.mean():.3f} ± {low_activity.std():.3f} (n={len(low_activity)})
- **高活动度 (SLEDAI > 8)**: {high_activity.mean():.3f} ± {high_activity.std():.3f} (n={len(high_activity)})
- **t-test p值**: {t_pval:.4f}
- **解释**: {"显著差异" if t_pval < 0.05 else "无显著差异"}
"""
    
    report_content += f"""
---

## 三、结论

### 主要发现
"""
    
    if pearson_p < 0.05:
        if pearson_r > 0:
            report_content += "- ✅ **TSCM signature与SLEDAI显著正相关**: TSCM signature score越高，疾病活动度越高\n"
            report_content += "- **生物学意义**: TSCM细胞可能参与疾病活动过程\n"
        else:
            report_content += "- ✅ **TSCM signature与SLEDAI显著负相关**: TSCM signature score越高，疾病活动度越低\n"
            report_content += "- **生物学意义**: TSCM细胞可能具有保护作用，或在高活动度时被抑制\n"
    else:
        report_content += "- ⚠️ **TSCM signature与SLEDAI无显著相关**: TSCM signature score与疾病活动度无关\n"
        report_content += "- **生物学意义**: TSCM signature可能主要反映疾病状态而非活动度\n"
    
    report_content += f"""
---

## 四、输出文件

- **分析结果**: `{out_file}`
- **可视化图片**: `{out_fig}`

---

*报告生成时间: {pd.Timestamp.now()}*
"""
    
    out_report = OUT_DIR / f"{dataset_id}_tscm_vs_sledai_report.md"
    with open(out_report, 'w', encoding='utf-8') as f:
        f.write(report_content)
    
    print(f"  ✓ 已生成报告: {out_report}")
    
    return {
        'dataset_id': dataset_id,
        'n_samples': len(df_analysis),
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p
    }

# 分析所有数据集
all_results = []

for dataset_info in sledai_datasets:
    try:
        result = analyze_tscm_vs_sledai(dataset_info)
        if result:
            all_results.append(result)
    except Exception as e:
        print(f"\n✗ {dataset_info['id']}分析失败: {e}")
        import traceback
        traceback.print_exc()

# 生成综合报告
if len(all_results) > 0:
    print("\n" + "=" * 80)
    print("综合结果")
    print("=" * 80)
    
    for result in all_results:
        print(f"\n{result['dataset_id']}:")
        print(f"  样本数: {result['n_samples']}")
        print(f"  Pearson r = {result['pearson_r']:.3f}, p = {result['pearson_p']:.4f}")
        print(f"  Spearman ρ = {result['spearman_r']:.3f}, p = {result['spearman_p']:.4f}")

print("\n" + "=" * 80)
print("分析完成！")
print("=" * 80)

