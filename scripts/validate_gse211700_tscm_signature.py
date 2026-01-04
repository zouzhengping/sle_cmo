#!/usr/bin/env python3
"""
GSE211700 TSCM Signature 验证脚本

自动完成以下步骤：
1. 检查并下载GSE211700数据（从GEO或SRA）
2. 加载TSCM signature（从GSE137029的DE结果）
3. 计算TSCM signature score
4. 比较SLE vs Healthy的差异
5. 生成验证报告

输出:
  results/validation/GSE211700_tscm_validation_report.md
  results/validation/GSE211700_tscm_scores.tsv
  results/figures/GSE211700_tscm_validation.pdf
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

# 设置matplotlib后端（避免GUI依赖）
import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 输入文件
SIGNATURE_DE = ROOT / "results/tables/tscm_vs_other_de.full.tsv"
SAMPLE_SHEET = ROOT / "data/metadata/sample_sheet.sorted.tsv"

# 输出目录
OUT_DIR = ROOT / "results/validation"
FIG_DIR = ROOT / "results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 输出文件
OUT_REPORT = OUT_DIR / "GSE211700_tscm_validation_report.md"
OUT_SCORES = OUT_DIR / "GSE211700_tscm_scores.tsv"
OUT_FIG = FIG_DIR / "GSE211700_tscm_validation.pdf"

# GSE211700数据目录
GSE211700_DIR = ROOT / "data/raw/GSE211700"
GSE211700_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("GSE211700 TSCM Signature 验证")
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
# 使用logFC > 0.5 和 pvals_adj < 0.05 作为阈值
sig_up = df_de[
    (df_de["logfoldchange"] > 0.5) & 
    (df_de["pvals_adj"] < 0.05)
].copy()

# 按logFC排序，取top 100作为核心signature
tscm_signature = sig_up.nlargest(100, "logfoldchange")["gene"].tolist()
print(f"  定义了 {len(tscm_signature)} 个TSCM signature基因（logFC > 0.5, adj.p < 0.05, top 100）")

# 也提取经典marker基因
classic_markers = ["CCR7", "SELL", "IL7R", "LEF1", "TCF7", "BCL2", "CXCR3", "CD27", "CD28"]
classic_markers_found = [g for g in classic_markers if g in df_de["gene"].values]
print(f"  经典TSCM marker: {len(classic_markers_found)}/{len(classic_markers)} 个在数据中")


# ============================================================================
# 步骤2: 检查并下载GSE211700数据
# ============================================================================
print("\n[步骤2] 检查GSE211700数据...")

def download_geo_matrix(gse_id):
    """尝试从GEO下载处理好的表达矩阵"""
    import urllib.request
    import gzip
    import shutil
    
    # 尝试多种URL格式
    url_patterns = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:7]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
    ]
    
    matrix_file = f"{gse_id}_series_matrix.txt.gz"
    out_file = GSE211700_DIR / matrix_file
    
    if out_file.exists():
        print(f"  发现已存在的matrix文件: {out_file}")
        return out_file
    
    for url in url_patterns:
        try:
            print(f"  尝试从GEO下载: {url}")
            urllib.request.urlretrieve(url, out_file)
            if out_file.stat().st_size > 0:
                print(f"  下载成功: {out_file}")
                return out_file
            else:
                out_file.unlink()  # 删除空文件
        except Exception as e:
            print(f"  尝试失败: {e}")
            continue
    
    print(f"  所有GEO下载尝试均失败")
    return None

def parse_geo_matrix(matrix_file):
    """解析GEO matrix文件"""
    import gzip
    
    print(f"  解析GEO matrix文件: {matrix_file}")
    
    # 读取文件
    with gzip.open(matrix_file, 'rt') as f:
        lines = f.readlines()
    
    # 查找数据开始位置
    data_start = None
    sample_names = []
    gene_ids = []
    
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            data_start = i + 1
        elif line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
            # 提取样本信息
            parts = line.strip().split("\t")
            if len(parts) > 1:
                sample_names = parts[1:]
        elif line.startswith("ID_REF"):
            # 这是表头行
            if not sample_names:
                sample_names = line.strip().split("\t")[1:]
    
    if data_start is None:
        raise ValueError("无法找到数据开始位置")
    
    # 读取表达数据
    data_rows = []
    for line in lines[data_start:]:
        if line.startswith("!series_matrix_table_end"):
            break
        parts = line.strip().split("\t")
        if len(parts) > 1:
            gene_ids.append(parts[0])
            # 转换为float，处理缺失值
            values = []
            for v in parts[1:]:
                try:
                    values.append(float(v))
                except:
                    values.append(np.nan)
            data_rows.append(values)
    
    # 创建DataFrame
    df_expr = pd.DataFrame(data_rows, index=gene_ids, columns=sample_names[:len(data_rows[0])])
    df_expr = df_expr.T  # 转置：样本x基因
    
    print(f"  解析完成: {df_expr.shape[0]} 个样本, {df_expr.shape[1]} 个基因")
    return df_expr

# 优先检查FPKM表达矩阵（包含实际表达数据）
fpkm_expr_file = GSE211700_DIR / "GSE211700.txt_expression_matrix.tsv"
if fpkm_expr_file.exists():
    print(f"  发现FPKM表达矩阵: {fpkm_expr_file}")
    try:
        df_expr = pd.read_csv(fpkm_expr_file, index_col=0, sep="\t")
        print(f"  ✓ 成功加载FPKM表达数据: {df_expr.shape}")
    except Exception as e:
        print(f"  ✗ 加载FPKM文件失败: {e}")
        df_expr = None
else:
    print(f"  ⚠ 未找到FPKM表达矩阵，尝试下载...")
    df_expr = None

# 如果FPKM文件不存在，尝试从GEO下载
if df_expr is None:
    matrix_file = download_geo_matrix("GSE211700")
    
    if matrix_file and matrix_file.exists():
        try:
            df_expr = parse_geo_matrix(matrix_file)
            print(f"  ✓ 成功加载GSE211700表达数据")
        except Exception as e:
            print(f"  ✗ 解析matrix文件失败: {e}")
            df_expr = None
    else:
        print(f"  ⚠ 无法从GEO下载matrix文件，尝试其他方法...")
        df_expr = None

# 如果无法从GEO获取，尝试检查是否有其他格式的数据
if df_expr is None:
    # 优先检查FPKM表达矩阵
    fpkm_expr_file = GSE211700_DIR / "GSE211700.txt_expression_matrix.tsv"
    if fpkm_expr_file.exists():
        print(f"  发现FPKM表达矩阵: {fpkm_expr_file}")
        try:
            df_expr = pd.read_csv(fpkm_expr_file, index_col=0, sep="\t")
            print(f"  ✓ 成功加载FPKM表达数据: {df_expr.shape}")
        except Exception as e:
            print(f"  ✗ 加载失败: {e}")
            df_expr = None
    
    # 如果FPKM文件不存在，检查其他格式
    if df_expr is None:
        count_files = [f for f in (list(GSE211700_DIR.glob("*.tsv")) + list(GSE211700_DIR.glob("*.csv"))) 
                       if "expression_matrix" in f.name or "count" in f.name.lower()]
        if count_files:
            print(f"  发现已处理的表达文件: {count_files[0]}")
            try:
                df_expr = pd.read_csv(count_files[0], index_col=0, sep="\t")
                print(f"  ✓ 成功加载表达数据: {df_expr.shape}")
            except Exception as e:
                print(f"  ✗ 加载失败: {e}")
                df_expr = None

if df_expr is None:
    print("\n[警告] 无法获取GSE211700的表达数据")
    print("  可能的原因:")
    print("  1. 数据需要从SRA下载FASTQ并定量（需要较长时间）")
    print("  2. GEO matrix文件格式不标准")
    print("  3. 需要手动下载并放置到 data/raw/GSE211700/ 目录")
    print("\n  生成模拟验证报告（基于signature定义）...")
    
    # 生成一个基于signature的报告
    report_content = f"""# GSE211700 TSCM Signature 验证报告

## 执行时间
{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 状态
⚠️ **数据未获取**: 无法自动下载或加载GSE211700的表达数据

## TSCM Signature定义

### Signature基因列表（Top 100, logFC > 0.5, adj.p < 0.05）
共 {len(tscm_signature)} 个基因：

{', '.join(tscm_signature[:20])}...
（完整列表见 results/tables/tscm_vs_other_de.full.tsv）

### 经典TSCM Marker基因
{', '.join(classic_markers_found) if classic_markers_found else '无'}

## 建议的下一步操作

1. **手动下载数据**:
   - 访问 GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211700
   - 下载处理好的表达矩阵（如果有）
   - 或下载原始FASTQ文件并进行定量

2. **数据格式要求**:
   - 表达矩阵：TSV或CSV格式
   - 行：基因（Gene Symbol或Ensembl ID）
   - 列：样本ID
   - 放置位置：`data/raw/GSE211700/`

3. **重新运行验证**:
   ```bash
   python scripts/validate_gse211700_tscm_signature.py
   ```

## Signature统计

- 总signature基因数: {len(tscm_signature)}
- 经典marker基因数: {len(classic_markers_found)}
- 来源: GSE137029 TSCM_high vs T_other DE分析

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

# 尝试从sample_sheet获取metadata
df_meta = pd.read_csv(SAMPLE_SHEET, sep="\t")
gse211700_meta = df_meta[df_meta["dataset_id"] == "GSE211700"].copy()

# 如果condition列为空，尝试从样本名称推断
if gse211700_meta["condition"].isna().all() or (gse211700_meta["condition"] == "").all():
    print("  ⚠ condition信息缺失，尝试从样本名称推断...")
    # 从表达数据的样本名推断condition
    if df_expr is not None:
        sample_conditions = []
        for sample_name in df_expr.index:
            sample_name_upper = str(sample_name).upper()
            if "CTRL" in sample_name_upper or "CONTROL" in sample_name_upper or "HEALTHY" in sample_name_upper:
                sample_conditions.append("Healthy")
            elif "SLE" in sample_name_upper:
                sample_conditions.append("SLE")
            elif "LN" in sample_name_upper:
                sample_conditions.append("SLE")  # LN也是SLE的一种
            else:
                sample_conditions.append("Unknown")
        
        # 创建样本名到condition的映射
        sample_to_condition = dict(zip(df_expr.index, sample_conditions))
        print(f"  从样本名推断condition: {pd.Series(sample_conditions).value_counts().to_dict()}")
    else:
        gse211700_meta["condition"] = "Unknown"

print(f"  加载了 {len(gse211700_meta)} 个样本的metadata")
print(f"  Condition分布: {gse211700_meta['condition'].value_counts().to_dict()}")


# ============================================================================
# 步骤4: 匹配样本和表达数据
# ============================================================================
print("\n[步骤4] 匹配样本和表达数据...")

# 创建metadata DataFrame（直接从样本名推断condition）
df_meta_expr = pd.DataFrame(index=df_expr.index)
df_meta_expr["sample_id"] = df_expr.index

# 从样本名推断condition
sample_conditions = []
for sample_name in df_expr.index:
    sample_name_upper = str(sample_name).upper()
    if "CTRL" in sample_name_upper or "CONTROL" in sample_name_upper or "HEALTHY" in sample_name_upper:
        sample_conditions.append("Healthy")
    elif "SLE" in sample_name_upper:
        sample_conditions.append("SLE")
    elif "LN" in sample_name_upper:
        sample_conditions.append("SLE")  # LN也是SLE的一种
    else:
        sample_conditions.append("Unknown")

df_meta_expr["condition"] = sample_conditions

print(f"  创建了 {len(df_meta_expr)} 个样本的metadata")
print(f"  Condition分布: {pd.Series(sample_conditions).value_counts().to_dict()}")


# ============================================================================
# 步骤5: 计算TSCM Signature Score
# ============================================================================
print("\n[步骤5] 计算TSCM Signature Score...")

def calculate_signature_score(expr_matrix, signature_genes, method='mean'):
    """
    计算signature score
    
    Parameters:
    -----------
    expr_matrix : pd.DataFrame
        表达矩阵（样本x基因）
    signature_genes : list
        signature基因列表
    method : str
        'mean': 平均表达
        'ssgsea': 类似SSGSEA的方法（简化版）
    """
    # 找到数据中存在的signature基因
    available_genes = [g for g in signature_genes if g in expr_matrix.columns]
    missing_genes = [g for g in signature_genes if g not in expr_matrix.columns]
    
    if len(missing_genes) > 0:
        print(f"    警告: {len(missing_genes)} 个signature基因在数据中不存在")
        print(f"    缺失基因示例: {missing_genes[:5]}")
    
    if len(available_genes) == 0:
        raise ValueError("没有可用的signature基因")
    
    # 提取signature基因的表达
    expr_sig = expr_matrix[available_genes]
    
    if method == 'mean':
        # 简单平均
        scores = expr_sig.mean(axis=1)
    elif method == 'ssgsea':
        # 简化的SSGSEA方法
        # 1. 按样本排序
        expr_ranked = expr_sig.rank(axis=1, method='average', pct=True)
        # 2. 计算加权得分
        scores = expr_ranked[available_genes].mean(axis=1)
    else:
        raise ValueError(f"未知的方法: {method}")
    
    return scores, available_genes, missing_genes

# 处理基因名匹配
print("\n[步骤5.1] 处理基因名匹配...")

# FPKM文件已经使用Official_Symbol作为列名，直接匹配即可
df_expr_cols_upper = df_expr.columns.str.upper()
tscm_signature_upper = [g.upper() for g in tscm_signature]

# 直接匹配
matched_genes = [g for g in tscm_signature_upper if g in df_expr_cols_upper.values]

print(f"  匹配到 {len(matched_genes)}/{len(tscm_signature)} 个signature基因 ({len(matched_genes)/len(tscm_signature)*100:.1f}%)")

if len(matched_genes) < len(tscm_signature) * 0.5:
    missing = [g for g in tscm_signature_upper if g not in df_expr_cols_upper.values]
    print(f"  缺失的signature基因示例: {missing[:10]}")

# 使用匹配到的基因
tscm_signature_upper = matched_genes

# 计算score
tscm_scores, available_genes, missing_genes = calculate_signature_score(
    df_expr, tscm_signature_upper, method='mean'
)

# 添加到metadata
df_meta_expr['tscm_score'] = tscm_scores.values

print(f"  ✓ 计算完成，使用了 {len(available_genes)}/{len(tscm_signature)} 个signature基因")
print(f"  TSCM Score范围: {tscm_scores.min():.3f} - {tscm_scores.max():.3f}")
print(f"  TSCM Score均值: {tscm_scores.mean():.3f} ± {tscm_scores.std():.3f}")


# ============================================================================
# 步骤6: 统计分析
# ============================================================================
print("\n[步骤6] 统计分析...")

# 按condition分组
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

# 如果有SLE和Healthy，进行统计检验
if 'SLE' in conditions and 'Healthy' in conditions:
    scores_sle = df_meta_expr[df_meta_expr['condition'] == 'SLE']['tscm_score']
    scores_healthy = df_meta_expr[df_meta_expr['condition'] == 'Healthy']['tscm_score']
    
    # t检验
    t_stat, t_pval = stats.ttest_ind(scores_sle, scores_healthy)
    
    # Mann-Whitney U检验（非参数）
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
fig.suptitle('GSE211700 TSCM Signature Validation', fontsize=16, fontweight='bold')

# 1. 箱线图：按condition分组
ax1 = axes[0, 0]
if len(conditions) > 1:
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
        
        # 添加统计检验结果
        if 'comparison' in results_stats:
            pval = results_stats['comparison']['t_test_pvalue']
            y_max = max([max(d) for d in plot_data])
            ax1.text(0.5, y_max * 1.05, f'p = {pval:.4f}', 
                    ha='center', fontsize=10, style='italic')

# 2. 直方图：score分布
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

# 3. 散点图：样本排序
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

# 4. Signature基因表达热图（top 20）
ax4 = axes[1, 1]
if len(available_genes) > 0:
    top_genes = available_genes[:20]  # top 20
    expr_top = df_expr[top_genes]
    
    # 按condition排序
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
    
    # 标准化（Z-score）
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

# 保存scores
df_scores = df_meta_expr[['condition', 'tscm_score']].copy()
df_scores.to_csv(OUT_SCORES, sep="\t", index=True)
print(f"  ✓ 已保存scores: {OUT_SCORES}")


# ============================================================================
# 步骤9: 生成报告
# ============================================================================
print("\n[步骤9] 生成验证报告...")

report_content = f"""# GSE211700 TSCM Signature 验证报告

## 执行信息
- **执行时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
- **数据集**: GSE211700 (Bulk RNA-seq PBMC)
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
{', '.join(classic_markers_found) if classic_markers_found else '无可用数据'}

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

包含以下图表:
1. **箱线图**: 按condition分组的TSCM score分布
2. **直方图**: TSCM score的总体分布
3. **散点图**: 各样本的TSCM score排序
4. **热图**: Top 20 signature基因的表达模式

---

## 五、结论

"""

# 生成结论
if 'comparison' in results_stats:
    comp = results_stats['comparison']
    if comp['t_test_pvalue'] < 0.05:
        if comp['fold_change'] > 1:
            conclusion = "✅ **验证成功**: GSE211700数据中，SLE患者的TSCM signature score显著高于Healthy对照，与GSE137029单细胞数据的发现一致。"
        else:
            conclusion = "⚠️ **结果不一致**: GSE211700数据中，SLE患者的TSCM signature score显著低于Healthy对照，与GSE137029单细胞数据的发现不一致，需要进一步调查。"
    else:
        conclusion = "⚠️ **无显著差异**: GSE211700数据中，SLE和Healthy之间的TSCM signature score无显著差异，可能是由于bulk RNA-seq数据的细胞类型混合效应。"
else:
    conclusion = "⚠️ **无法比较**: 缺少SLE和Healthy的分组信息，无法进行验证比较。"

report_content += conclusion

report_content += f"""

### 局限性
1. **数据来源差异**: GSE211700是bulk RNA-seq数据，而signature来自单细胞数据，可能存在细胞类型混合效应
2. **样本量**: 当前分析的样本量可能不足以检测到细微差异
3. **批次效应**: 不同数据集的批次效应可能影响结果

### 建议
1. 如果可能，使用单细胞数据或分选的T细胞数据进行验证
2. 考虑使用deconvolution方法估计TSCM细胞比例
3. 结合其他独立数据集进行meta分析

---

## 六、输出文件

- **验证报告**: `{OUT_REPORT}`
- **Score数据**: `{OUT_SCORES}`
- **可视化图片**: `{OUT_FIG}`

---

*报告生成时间: {pd.Timestamp.now()}*
*脚本: scripts/validate_gse211700_tscm_signature.py*
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

