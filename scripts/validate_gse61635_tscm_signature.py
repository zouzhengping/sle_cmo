#!/usr/bin/env python3
"""
GSE61635 TSCM Signature 验证脚本

基于GSE211700验证脚本，适配GSE61635数据集
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib后端
import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 输入文件
SIGNATURE_DE = ROOT / "results/tables/tscm_vs_other_de.full.tsv"
SAMPLE_SHEET = ROOT / "data/raw/GSE61635/GSE61635.txt_metadata.tsv"

# 输出目录
OUT_DIR = ROOT / "results/validation"
FIG_DIR = ROOT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 输出文件
OUT_REPORT = OUT_DIR / "GSE61635_tscm_validation_report.md"
OUT_SCORES = OUT_DIR / "GSE61635_tscm_scores.tsv"
OUT_FIG = FIG_DIR / "GSE61635_tscm_validation.pdf"

# GSE61635数据目录
GSE61635_DIR = ROOT / "data/raw/GSE61635"
GSE61635_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("GSE61635 TSCM Signature 验证")
print("=" * 80)

# ============================================================================
# 步骤1: 加载TSCM Signature
# ============================================================================
print("\n[步骤1] 加载TSCM Signature...")

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
print(f"  定义了 {len(tscm_signature)} 个TSCM signature基因（logFC > 0.5, adj.p < 0.05, top 100）")

# 经典marker基因
classic_markers = ["CCR7", "SELL", "IL7R", "LEF1", "TCF7", "BCL2", "CXCR3", "CD27", "CD28"]
classic_markers_found = [g for g in classic_markers if g in df_de["gene"].values]
print(f"  经典TSCM marker: {len(classic_markers_found)}/{len(classic_markers)} 个在数据中")

# ============================================================================
# 步骤2: 检查并加载GSE61635数据
# ============================================================================
print("\n[步骤2] 检查GSE61635数据...")

def download_geo_matrix(gse_id):
    """尝试从GEO下载series matrix文件"""
    import urllib.request
    import gzip
    
    url_patterns = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
    ]
    
    matrix_file = f"{gse_id}_series_matrix.txt.gz"
    out_file = GSE61635_DIR / matrix_file
    
    if out_file.exists() and out_file.stat().st_size > 1000:
        return out_file
    
    for url in url_patterns:
        try:
            urllib.request.urlretrieve(url, out_file)
            if out_file.stat().st_size > 1000:
                return out_file
            else:
                out_file.unlink()
        except:
            continue
    return None

def parse_geo_matrix(matrix_file):
    """解析GEO matrix文件（microarray或bulk RNA-seq）"""
    print(f"  解析GEO matrix文件: {matrix_file}")
    
    with gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    data_start_idx = None
    data_end_idx = None
    sample_names = []
    
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            data_start_idx = i + 1
        elif line.startswith("!series_matrix_table_end"):
            data_end_idx = i
        elif line.startswith("!Sample_geo_accession"):
            parts = line.strip().split("\t")
            if len(parts) > 1:
                sample_names = parts[1:]
    
    if data_start_idx is None:
        raise ValueError("无法找到数据开始位置")
    
    # 读取表头
    header_line = lines[data_start_idx]
    header_parts = header_line.strip().split("\t")
    
    # 读取数据
    data_rows = []
    gene_ids = []
    
    for line in lines[data_start_idx + 1:data_end_idx]:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        
        gene_ids.append(parts[0])
        values = []
        for v in parts[1:]:
            try:
                values.append(float(v))
            except:
                values.append(np.nan)
        data_rows.append(values)
    
    sample_names = header_parts[1:] if len(header_parts) > 1 else [f"Sample_{i+1}" for i in range(len(data_rows[0]) if data_rows else 0)]
    
    if len(sample_names) != len(data_rows[0]) if data_rows else 0:
        sample_names = [f"Sample_{i+1}" for i in range(len(data_rows[0]) if data_rows else 0)]
    
    df_expr = pd.DataFrame(data_rows, index=gene_ids, columns=sample_names)
    df_expr = df_expr.T  # 转置：样本x基因
    
    print(f"  解析完成: {df_expr.shape[0]} 个样本, {df_expr.shape[1]} 个基因")
    return df_expr

# 优先检查转换后的Gene Symbol表达矩阵
gene_symbol_expr_file = GSE61635_DIR / "GSE61635_expression_matrix_gene_symbol.tsv"
df_expr = None

if gene_symbol_expr_file.exists():
    print(f"  发现转换后的Gene Symbol表达矩阵: {gene_symbol_expr_file}")
    try:
        df_expr = pd.read_csv(gene_symbol_expr_file, index_col=0, sep="\t")
        print(f"  ✓ 成功加载转换后的表达数据: {df_expr.shape}")
    except Exception as e:
        print(f"  ✗ 加载失败: {e}")
        df_expr = None

# 如果转换后的文件不存在，尝试解析原始matrix文件
if df_expr is None:
    matrix_file = GSE61635_DIR / "GSE61635_series_matrix.txt.gz"
    if matrix_file.exists() and matrix_file.stat().st_size > 1000:
        print(f"  ⚠ 未找到转换后的文件，尝试解析原始matrix文件...")
        print(f"  提示: 原始文件使用探针ID，需要先运行转换脚本")
        print(f"  运行: Rscript scripts/convert_gse61635_probe_to_gene.R")
        try:
            df_expr = parse_geo_matrix(matrix_file)
            print(f"  ✓ 成功加载GSE61635表达数据（探针ID格式）")
        except Exception as e:
            print(f"  ✗ 解析失败: {e}")
            df_expr = None

if df_expr is None:
    print("\n[警告] 无法获取GSE61635的表达数据")
    print("  生成模拟验证报告（基于signature定义）...")
    
    report_content = f"""# GSE61635 TSCM Signature 验证报告

## 执行时间
{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 状态
⚠️ **数据未获取**: 无法自动下载或加载GSE61635的表达数据

## TSCM Signature定义

### Signature基因列表（Top 100, logFC > 0.5, adj.p < 0.05）
共 {len(tscm_signature)} 个基因

### 经典TSCM Marker基因
{', '.join(classic_markers_found) if classic_markers_found else '无'}

## 建议的下一步操作

1. **检查数据类型**: GSE61635可能是microarray数据，需要特殊处理
2. **手动下载数据**: 访问 GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61635
3. **使用GEOquery**: 考虑使用R的GEOquery包下载和处理数据

---
*报告生成时间: {pd.Timestamp.now()}*
"""
    
    with open(OUT_REPORT, 'w', encoding='utf-8') as f:
        f.write(report_content)
    
    print(f"\n[完成] 已生成验证报告: {OUT_REPORT}")
    sys.exit(0)

# ============================================================================
# 步骤3: 加载样本metadata
# ============================================================================
print("\n[步骤3] 加载样本metadata...")

if SAMPLE_SHEET.exists():
    df_meta = pd.read_csv(SAMPLE_SHEET, sep="\t")
    print(f"  加载了 {len(df_meta)} 个样本的metadata")
    print(f"  Condition分布: {df_meta['condition'].value_counts().to_dict()}")
else:
    # 从样本名推断
    sample_conditions = []
    for sample_name in df_expr.index:
        sample_name_upper = str(sample_name).upper()
        if "SLE" in sample_name_upper or "LUPUS" in sample_name_upper:
            sample_conditions.append("SLE")
        elif "HEALTHY" in sample_name_upper or "CONTROL" in sample_name_upper or "CTRL" in sample_name_upper:
            sample_conditions.append("Healthy")
        else:
            sample_conditions.append("Unknown")
    
    df_meta = pd.DataFrame({
        'sample_id': df_expr.index,
        'condition': sample_conditions
    })
    print(f"  从样本名推断condition: {pd.Series(sample_conditions).value_counts().to_dict()}")

# ============================================================================
# 步骤4: 匹配样本和表达数据
# ============================================================================
print("\n[步骤4] 匹配样本和表达数据...")

df_meta_expr = pd.DataFrame(index=df_expr.index)
df_meta_expr["sample_id"] = df_expr.index

# 匹配condition
if "condition" in df_meta.columns:
    # 尝试通过sample_id或geo_accession匹配
    condition_map = {}
    for _, row in df_meta.iterrows():
        sample_id = str(row.get("sample_id", "")).strip('"')
        geo_acc = str(row.get("geo_accession", "")).strip('"')
        condition = str(row.get("condition", "Unknown")).strip('"')
        
        for idx in df_expr.index:
            idx_str = str(idx).strip('"')
            # 尝试多种匹配方式
            if sample_id and sample_id in idx_str:
                condition_map[idx] = condition
            elif geo_acc and geo_acc in idx_str:
                condition_map[idx] = condition
            elif idx_str in sample_id or idx_str in geo_acc:
                condition_map[idx] = condition
    
    # 如果匹配失败，尝试直接使用GEO accession作为key
    if len(condition_map) == 0 and "geo_accession" in df_meta.columns:
        for _, row in df_meta.iterrows():
            geo_acc = str(row.get("geo_accession", "")).strip('"')
            condition = str(row.get("condition", "Unknown")).strip('"')
            if geo_acc and geo_acc in df_expr.index:
                condition_map[geo_acc] = condition
    
    df_meta_expr["condition"] = [condition_map.get(idx, "Unknown") for idx in df_expr.index]
    
    # 如果仍然都是Unknown，尝试从样本名推断
    if df_meta_expr["condition"].value_counts().get("Unknown", 0) == len(df_meta_expr):
        print("  ⚠ 无法通过metadata匹配，尝试从样本名推断...")
        sample_conditions = []
        for sample_name in df_expr.index:
            sample_name_upper = str(sample_name).upper()
            if "SLE" in sample_name_upper or "LUPUS" in sample_name_upper:
                sample_conditions.append("SLE")
            elif "HEALTHY" in sample_name_upper or "CONTROL" in sample_name_upper or "CTRL" in sample_name_upper:
                sample_conditions.append("Healthy")
            else:
                sample_conditions.append("Unknown")
        df_meta_expr["condition"] = sample_conditions
else:
    # 从样本名推断
    sample_conditions = []
    for sample_name in df_expr.index:
        sample_name_upper = str(sample_name).upper()
        if "SLE" in sample_name_upper or "LUPUS" in sample_name_upper:
            sample_conditions.append("SLE")
        elif "HEALTHY" in sample_name_upper or "CONTROL" in sample_name_upper or "CTRL" in sample_name_upper:
            sample_conditions.append("Healthy")
        else:
            sample_conditions.append("Unknown")
    df_meta_expr["condition"] = sample_conditions

print(f"  创建了 {len(df_meta_expr)} 个样本的metadata")
print(f"  Condition分布: {df_meta_expr['condition'].value_counts().to_dict()}")

# ============================================================================
# 步骤5: 计算TSCM Signature Score
# ============================================================================
print("\n[步骤5] 计算TSCM Signature Score...")

def calculate_signature_score(expr_matrix, signature_genes, method='mean'):
    """计算signature score"""
    available_genes = [g for g in signature_genes if g in expr_matrix.columns]
    missing_genes = [g for g in signature_genes if g not in expr_matrix.columns]
    
    if len(missing_genes) > 0:
        print(f"    警告: {len(missing_genes)} 个signature基因在数据中不存在")
    
    if len(available_genes) == 0:
        raise ValueError("没有可用的signature基因")
    
    expr_sig = expr_matrix[available_genes]
    
    if method == 'mean':
        scores = expr_sig.mean(axis=1)
    else:
        scores = expr_sig.mean(axis=1)
    
    return scores, available_genes, missing_genes

# 处理基因名匹配
print("\n[步骤5.1] 处理基因名匹配...")

# GSE61635是microarray数据（GPL570），列名可能是探针ID
# 需要检查列名格式
sample_cols = df_expr.columns.tolist()
print(f"  数据列名示例: {sample_cols[:5]}")
print(f"  列名格式检查: 是否为探针ID（如以数字开头）")

# 尝试直接匹配Gene Symbol
df_expr_cols_upper = df_expr.columns.str.upper()
tscm_signature_upper = [g.upper() for g in tscm_signature]

matched_genes = [g for g in tscm_signature_upper if g in df_expr_cols_upper.values]
print(f"  直接匹配: {len(matched_genes)}/{len(tscm_signature)} 个signature基因")

# 如果匹配率很低，可能是探针ID，需要转换
if len(matched_genes) < len(tscm_signature) * 0.3:
    print("  ⚠ 匹配率很低，可能是探针ID格式")
    print("  提示: GSE61635是microarray数据（GPL570），可能需要探针ID到Gene Symbol的转换")
    print("  建议: 使用GEOquery R包或biomaRt进行ID转换")
    
    # 尝试从行名（基因ID）中查找
    # 对于microarray，行名可能是探针ID，列名是样本
    # 需要转置后检查
    if df_expr.shape[1] > df_expr.shape[0]:
        print("  尝试从行名（探针ID）匹配...")
        # 这里暂时跳过，需要ID转换工具
        print("  需要探针ID到Gene Symbol的转换，暂时跳过")
    
    # 如果仍然没有匹配，使用所有可用的列（可能是探针）
    if len(matched_genes) == 0:
        print("  ⚠ 无法匹配任何signature基因")
        print("  可能原因: 数据使用探针ID而非Gene Symbol")
        print("  建议: 需要先进行ID转换")
        raise ValueError("无法匹配signature基因，可能需要探针ID到Gene Symbol的转换")

tscm_signature_upper = matched_genes if matched_genes else []

# 计算score
tscm_scores, available_genes, missing_genes = calculate_signature_score(
    df_expr, tscm_signature_upper, method='mean'
)

df_meta_expr['tscm_score'] = tscm_scores.values

print(f"  ✓ 计算完成，使用了 {len(available_genes)}/{len(tscm_signature)} 个signature基因")
print(f"  TSCM Score范围: {tscm_scores.min():.3f} - {tscm_scores.max():.3f}")
print(f"  TSCM Score均值: {tscm_scores.mean():.3f} ± {tscm_scores.std():.3f}")

# ============================================================================
# 步骤6: 统计分析
# ============================================================================
print("\n[步骤6] 统计分析...")

conditions = df_meta_expr['condition'].unique()
print(f"  Condition分组: {conditions}")

results_stats = {}

for cond in conditions:
    if cond == "Unknown":
        continue
    scores_cond = df_meta_expr[df_meta_expr['condition'] == cond]['tscm_score']
    results_stats[cond] = {
        'n': len(scores_cond),
        'mean': scores_cond.mean(),
        'std': scores_cond.std(),
        'median': scores_cond.median(),
        'min': scores_cond.min(),
        'max': scores_cond.max()
    }
    print(f"  {cond}: n={len(scores_cond)}, mean={scores_cond.mean():.3f}±{scores_cond.std():.3f}")

# 统计检验
if 'SLE' in conditions and 'Healthy' in conditions:
    scores_sle = df_meta_expr[df_meta_expr['condition'] == 'SLE']['tscm_score']
    scores_healthy = df_meta_expr[df_meta_expr['condition'] == 'Healthy']['tscm_score']
    
    t_stat, t_pval = stats.ttest_ind(scores_sle, scores_healthy)
    u_stat, u_pval = stats.mannwhitneyu(scores_sle, scores_healthy, alternative='two-sided')
    
    results_stats['comparison'] = {
        't_test_statistic': t_stat,
        't_test_pvalue': t_pval,
        'mannwhitney_statistic': u_stat,
        'mannwhitney_pvalue': u_pval,
        'fold_change': scores_sle.mean() / scores_healthy.mean() if scores_healthy.mean() > 0 else np.nan
    }
    
    print(f"\n  SLE vs Healthy:")
    print(f"    t-test: statistic={t_stat:.3f}, p={t_pval:.4f}")
    print(f"    Mann-Whitney U: statistic={u_stat:.3f}, p={u_pval:.4f}")
    print(f"    Fold change: {results_stats['comparison']['fold_change']:.3f}")

# ============================================================================
# 步骤7: 可视化
# ============================================================================
print("\n[步骤7] 生成可视化...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('GSE61635 TSCM Signature Validation', fontsize=16, fontweight='bold')

# 1. 箱线图
ax1 = axes[0, 0]
plot_data = []
plot_labels = []
for cond in conditions:
    if cond != "Unknown":
        scores_cond = df_meta_expr[df_meta_expr['condition'] == cond]['tscm_score']
        plot_data.append(scores_cond.values)
        plot_labels.append(cond)

if plot_data:
    bp = ax1.boxplot(plot_data, labels=plot_labels, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    ax1.set_ylabel('TSCM Signature Score', fontsize=12)
    ax1.set_title('TSCM Score by Condition', fontsize=13, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    if 'comparison' in results_stats:
        pval = results_stats['comparison']['t_test_pvalue']
        y_max = max([max(d) for d in plot_data])
        ax1.text(0.5, y_max * 1.05, f'p = {pval:.4f}', ha='center', fontsize=10, style='italic')

# 2. 直方图
ax2 = axes[0, 1]
for cond in conditions:
    if cond != "Unknown":
        scores_cond = df_meta_expr[df_meta_expr['condition'] == cond]['tscm_score']
        ax2.hist(scores_cond, alpha=0.6, label=cond, bins=20)
ax2.set_xlabel('TSCM Signature Score', fontsize=12)
ax2.set_ylabel('Frequency', fontsize=12)
ax2.set_title('TSCM Score Distribution', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(axis='y', alpha=0.3)

# 3. 散点图
ax3 = axes[1, 0]
df_plot = df_meta_expr[df_meta_expr['condition'] != 'Unknown'].copy()
if len(df_plot) > 0:
    df_plot = df_plot.sort_values('tscm_score')
    colors = {'SLE': 'red', 'Healthy': 'blue'}
    for cond in df_plot['condition'].unique():
        mask = df_plot['condition'] == cond
        ax3.scatter(range(mask.sum()), df_plot[mask]['tscm_score'], 
                   label=cond, alpha=0.6, s=50, c=colors.get(cond, 'gray'))
    ax3.set_xlabel('Sample Index (sorted)', fontsize=12)
    ax3.set_ylabel('TSCM Signature Score', fontsize=12)
    ax3.set_title('TSCM Score by Sample', fontsize=13, fontweight='bold')
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)

# 4. 热图
ax4 = axes[1, 1]
if len(available_genes) > 0:
    top_genes = available_genes[:20]
    expr_top = df_expr[top_genes]
    
    if 'condition' in df_meta_expr.columns:
        condition_order = ['Healthy', 'SLE', 'Unknown']
        order = []
        for cond in condition_order:
            if cond in df_meta_expr['condition'].values:
                order.extend(df_meta_expr[df_meta_expr['condition'] == cond].index.tolist())
        if len(order) < len(expr_top.index):
            remaining = [idx for idx in expr_top.index if idx not in order]
            order.extend(remaining)
        expr_top = expr_top.loc[order]
    
    expr_top_norm = (expr_top - expr_top.mean()) / expr_top.std()
    
    im = ax4.imshow(expr_top_norm.T, aspect='auto', cmap='RdYlBu_r', interpolation='nearest')
    ax4.set_xlabel('Samples', fontsize=12)
    ax4.set_ylabel('Top 20 Signature Genes', fontsize=12)
    ax4.set_title('Signature Gene Expression (Z-score)', fontsize=13, fontweight='bold')
    ax4.set_yticks(range(len(top_genes)))
    ax4.set_yticklabels(top_genes, fontsize=8)
    plt.colorbar(im, ax=ax4)

plt.tight_layout()
plt.savefig(OUT_FIG, dpi=300, bbox_inches='tight')
print(f"  ✓ 已保存图片: {OUT_FIG}")

# ============================================================================
# 步骤8: 保存结果
# ============================================================================
print("\n[步骤8] 保存结果...")

df_scores = df_meta_expr[['condition', 'tscm_score']].copy()
df_scores.to_csv(OUT_SCORES, sep="\t", index=True)
print(f"  ✓ 已保存scores: {OUT_SCORES}")

# ============================================================================
# 步骤9: 生成报告
# ============================================================================
print("\n[步骤9] 生成验证报告...")

report_content = f"""# GSE61635 TSCM Signature 验证报告

## 执行信息
- **执行时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
- **数据集**: GSE61635
- **验证方法**: TSCM Signature Score (基于GSE137029单细胞数据)

---

## 一、TSCM Signature定义

### Signature来源
- **原始数据集**: GSE137029 (单细胞RNA-seq)
- **分析方法**: TSCM_high vs T_other 差异表达分析
- **筛选标准**: 
  - logFC > 0.5
  - 调整后p值 < 0.05
  - Top 100 基因（按logFC排序）

### Signature基因
- **总基因数**: {len(tscm_signature)}
- **在验证数据中可用**: {len(available_genes)}/{len(tscm_signature)} ({len(available_genes)/len(tscm_signature)*100:.1f}%)
- **缺失基因数**: {len(missing_genes)}

### 经典TSCM Marker基因
{', '.join(classic_markers_found) if classic_markers_found else '无'}

---

## 二、数据概览

### 样本信息
- **总样本数**: {len(df_expr)}
- **Condition分布**:
"""

for cond, count in df_meta_expr['condition'].value_counts().items():
    report_content += f"  - {cond}: {count} 个样本\n"

report_content += f"""
- **基因数**: {len(df_expr.columns)}

---

## 三、TSCM Signature Score 计算结果

### 总体统计
- **Score范围**: {tscm_scores.min():.3f} - {tscm_scores.max():.3f}
- **Score均值**: {tscm_scores.mean():.3f} ± {tscm_scores.std():.3f}
- **Score中位数**: {tscm_scores.median():.3f}

### 按Condition分组统计
"""

for cond, stats_dict in results_stats.items():
    if cond == 'comparison':
        continue
    report_content += f"""
#### {cond}
- 样本数: {stats_dict['n']}
- 均值: {stats_dict['mean']:.3f} ± {stats_dict['std']:.3f}
- 中位数: {stats_dict['median']:.3f}
- 范围: {stats_dict['min']:.3f} - {stats_dict['max']:.3f}
"""

if 'comparison' in results_stats:
    comp = results_stats['comparison']
    report_content += f"""
### SLE vs Healthy 比较

- **Fold Change**: {comp['fold_change']:.3f}
- **t-test**: 
  - 统计量: {comp['t_test_statistic']:.3f}
  - p值: {comp['t_test_pvalue']:.4f}
- **Mann-Whitney U test**:
  - 统计量: {comp['mannwhitney_statistic']:.3f}
  - p值: {comp['mannwhitney_pvalue']:.4f}

**解释**:
"""

    if comp['t_test_pvalue'] < 0.05:
        direction = "更高" if comp['fold_change'] > 1 else "更低"
        report_content += f"- ✅ **显著差异**: SLE患者的TSCM signature score {direction}于Healthy对照 (p = {comp['t_test_pvalue']:.4f})\n"
    else:
        report_content += f"- ⚠️ **无显著差异**: SLE和Healthy之间的TSCM signature score无显著差异 (p = {comp['t_test_pvalue']:.4f})\n"

report_content += f"""
---

## 四、可视化结果

所有可视化结果已保存至: `{OUT_FIG}`

---

## 五、结论

"""

if 'comparison' in results_stats:
    comp = results_stats['comparison']
    if comp['t_test_pvalue'] < 0.05:
        if comp['fold_change'] > 1:
            conclusion = "✅ **验证成功**: GSE61635数据中，SLE患者的TSCM signature score显著高于Healthy对照，与GSE137029单细胞数据的发现一致。"
        else:
            conclusion = "⚠️ **结果不一致**: GSE61635数据中，SLE患者的TSCM signature score显著低于Healthy对照，与GSE137029单细胞数据的发现不一致，需要进一步调查。"
    else:
        conclusion = "⚠️ **无显著差异**: GSE61635数据中，SLE和Healthy之间的TSCM signature score无显著差异。"
else:
    conclusion = "⚠️ **无法比较**: 缺少SLE和Healthy的分组信息，无法进行验证比较。"

report_content += conclusion

report_content += f"""

---

## 六、输出文件

- **验证报告**: `{OUT_REPORT}`
- **Score数据**: `{OUT_SCORES}`
- **可视化图片**: `{OUT_FIG}`

---

*报告生成时间: {pd.Timestamp.now()}*
"""

with open(OUT_REPORT, 'w', encoding='utf-8') as f:
    f.write(report_content)

print(f"  ✓ 已生成验证报告: {OUT_REPORT}")

print("\n" + "=" * 80)
print("验证完成！")
print("=" * 80)
print(f"\n输出文件:")
print(f"  - 验证报告: {OUT_REPORT}")
print(f"  - Score数据: {OUT_SCORES}")
print(f"  - 可视化: {OUT_FIG}")
print()

