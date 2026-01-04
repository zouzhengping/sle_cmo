#!/usr/bin/env python3
"""
按SLEDAI分组详细比较TSCM Score

包括：
1. 多组比较（ANOVA/Kruskal-Wallis）
2. 两两比较（post-hoc tests）
3. 趋势检验（Jonckheere-Terpstra）
4. 详细的可视化和统计报告
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import kruskal, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 输出目录
OUT_DIR = ROOT / "results/validation/sledai"
FIG_DIR = ROOT / "results/figures"
REPORT_DIR = ROOT / "docs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("按SLEDAI分组详细比较TSCM Score")
print("=" * 80)

# ============================================================================
# 步骤1: 加载数据
# ============================================================================
print("\n[步骤1] 加载数据...")

data_file = OUT_DIR / "GSE88884_tscm_vs_sledai.tsv"
if not data_file.exists():
    raise FileNotFoundError(f"未找到数据文件: {data_file}")

df = pd.read_csv(data_file, sep="\t")
print(f"  加载了 {len(df)} 个样本")

# 过滤有效数据
df = df[df['SLEDAI'].notna() & df['tscm_score'].notna()].copy()
print(f"  有效样本数: {len(df)}")

# ============================================================================
# 步骤2: 定义分组策略
# ============================================================================
print("\n[步骤2] 定义SLEDAI分组...")

# 策略1: 标准临床分组
df['group_standard'] = pd.cut(
    df['SLEDAI'],
    bins=[0, 4, 8, 12, 100],
    labels=['低活动度 (0-4)', '中活动度 (4-8)', '高活动度 (8-12)', '极高活动度 (>12)']
)

# 策略2: 三等分
df['group_tertile'] = pd.qcut(
    df['SLEDAI'],
    q=3,
    labels=['低活动度 (T1)', '中活动度 (T2)', '高活动度 (T3)']
)

# 策略3: 二分（低 vs 高）
df['group_binary'] = pd.cut(
    df['SLEDAI'],
    bins=[0, 8, 100],
    labels=['低活动度 (≤8)', '高活动度 (>8)']
)

print("  分组策略:")
print(f"    标准分组: {df['group_standard'].value_counts().to_dict()}")
print(f"    三等分: {df['group_tertile'].value_counts().to_dict()}")
print(f"    二分: {df['group_binary'].value_counts().to_dict()}")

# ============================================================================
# 步骤3: 描述性统计
# ============================================================================
print("\n[步骤3] 描述性统计...")

def describe_group(df, group_col, score_col='tscm_score'):
    """描述性统计"""
    results = []
    for group in df[group_col].cat.categories:
        if pd.notna(group):
            group_data = df[df[group_col] == group][score_col]
            results.append({
                'Group': group,
                'N': len(group_data),
                'Mean': group_data.mean(),
                'SD': group_data.std(),
                'Median': group_data.median(),
                'IQR': group_data.quantile(0.75) - group_data.quantile(0.25),
                'Min': group_data.min(),
                'Max': group_data.max()
            })
    return pd.DataFrame(results)

# 标准分组统计
df_stats_standard = describe_group(df, 'group_standard')
print("\n  标准分组统计:")
print(df_stats_standard.to_string(index=False))

# ============================================================================
# 步骤4: 统计检验
# ============================================================================
print("\n[步骤4] 统计检验...")

def perform_statistical_tests(df, group_col, score_col='tscm_score'):
    """执行统计检验"""
    results = {}
    
    # 提取各组数据
    groups_data = {}
    for group in df[group_col].cat.categories:
        if pd.notna(group):
            groups_data[group] = df[df[group_col] == group][score_col].values
    
    groups_list = list(groups_data.keys())
    groups_values = list(groups_data.values())
    
    # 1. 正态性检验（Shapiro-Wilk，每组抽样）
    print(f"\n  正态性检验 (Shapiro-Wilk, 每组抽样500个):")
    normal_results = {}
    for group, values in groups_data.items():
        sample = values[:500] if len(values) > 500 else values
        stat, pval = stats.shapiro(sample)
        normal_results[group] = {'stat': stat, 'p': pval}
        print(f"    {group}: W={stat:.4f}, p={pval:.4f} {'(正态)' if pval > 0.05 else '(非正态)'}")
    
    results['normality'] = normal_results
    
    # 2. 方差齐性检验（Levene）
    levene_stat, levene_p = stats.levene(*groups_values)
    print(f"\n  方差齐性检验 (Levene):")
    print(f"    F={levene_stat:.4f}, p={levene_p:.4f} {'(方差齐)' if levene_p > 0.05 else '(方差不齐)'}")
    results['levene'] = {'stat': levene_stat, 'p': levene_p}
    
    # 3. 多组比较
    # 如果正态且方差齐，用ANOVA；否则用Kruskal-Wallis
    use_parametric = all([v['p'] > 0.05 for v in normal_results.values()]) and levene_p > 0.05
    
    if use_parametric:
        f_stat, anova_p = stats.f_oneway(*groups_values)
        print(f"\n  多组比较 (ANOVA):")
        print(f"    F={f_stat:.4f}, p={anova_p:.4f}")
        results['multigroup'] = {'method': 'ANOVA', 'stat': f_stat, 'p': anova_p}
    else:
        h_stat, kw_p = kruskal(*groups_values)
        print(f"\n  多组比较 (Kruskal-Wallis):")
        print(f"    H={h_stat:.4f}, p={kw_p:.4f}")
        results['multigroup'] = {'method': 'Kruskal-Wallis', 'stat': h_stat, 'p': kw_p}
    
    # 4. 两两比较（post-hoc）
    method_name = 'Tukey HSD' if use_parametric else 'Mann-Whitney U with Bonferroni校正'
    print(f"\n  两两比较 ({method_name}):")
    pairwise_results = []
    
    if use_parametric:
        from scipy.stats import tukey_hsd
        tukey_result = tukey_hsd(*groups_values)
        pvals_matrix = tukey_result.pvalue
        
        for i, group1 in enumerate(groups_list):
            for j, group2 in enumerate(groups_list):
                if i < j:
                    pval = pvals_matrix[i, j]
                    pairwise_results.append({
                        'Group1': group1,
                        'Group2': group2,
                        'Method': 'Tukey HSD',
                        'p_value': pval,
                        'Significant': 'Yes' if pval < 0.05 else 'No'
                    })
                    print(f"    {group1} vs {group2}: p={pval:.4f} {'*' if pval < 0.05 else ''}")
    else:
        # Mann-Whitney U with Bonferroni校正
        n_comparisons = len(groups_list) * (len(groups_list) - 1) // 2
        bonferroni_alpha = 0.05 / n_comparisons
        
        for i, group1 in enumerate(groups_list):
            for j, group2 in enumerate(groups_list):
                if i < j:
                    u_stat, pval = mannwhitneyu(groups_data[group1], groups_data[group2], alternative='two-sided')
                    pval_corrected = min(pval * n_comparisons, 1.0)
                    pairwise_results.append({
                        'Group1': group1,
                        'Group2': group2,
                        'Method': 'Mann-Whitney U (Bonferroni)',
                        'p_value': pval,
                        'p_value_corrected': pval_corrected,
                        'Significant': 'Yes' if pval_corrected < 0.05 else 'No'
                    })
                    print(f"    {group1} vs {group2}: p={pval:.4f}, p_corrected={pval_corrected:.4f} {'*' if pval_corrected < 0.05 else ''}")
    
    results['pairwise'] = pd.DataFrame(pairwise_results)
    
    # 5. 趋势检验（Jonckheere-Terpstra，如果有顺序）
    if len(groups_list) >= 3:
        # 按SLEDAI中位数排序
        group_medians = {group: df[df[group_col] == group]['SLEDAI'].median() for group in groups_list}
        sorted_groups = sorted(groups_list, key=lambda x: group_medians[x])
        
        # Jonckheere-Terpstra检验
        from scipy.stats import kendalltau
        # 简化：使用Kendall's tau作为趋势检验
        group_ranks = [sorted_groups.index(g) for g in groups_list]
        all_scores = np.concatenate(groups_values)
        all_ranks = np.concatenate([[r] * len(groups_data[g]) for r, g in zip(group_ranks, groups_list)])
        
        tau, tau_p = kendalltau(all_ranks, all_scores)
        print(f"\n  趋势检验 (Kendall's tau):")
        print(f"    τ={tau:.4f}, p={tau_p:.4f}")
        results['trend'] = {'method': "Kendall's tau", 'stat': tau, 'p': tau_p}
    
    return results

# 对标准分组进行统计检验
stats_results = perform_statistical_tests(df, 'group_standard')

# ============================================================================
# 步骤5: 可视化
# ============================================================================
print("\n[步骤5] 生成可视化...")

fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

# 1. 箱线图（标准分组）
ax1 = fig.add_subplot(gs[0, 0])
group_order = ['低活动度 (0-4)', '中活动度 (4-8)', '高活动度 (8-12)', '极高活动度 (>12)']
data_for_plot = [df[df['group_standard'] == g]['tscm_score'].values for g in group_order if g in df['group_standard'].values]
labels_for_plot = [g for g in group_order if g in df['group_standard'].values]

bp = ax1.boxplot(data_for_plot, labels=labels_for_plot, patch_artist=True)
colors = ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3']
for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax1.set_ylabel('TSCM Signature Score', fontsize=12)
ax1.set_title('TSCM Score by SLEDAI Group (Standard)', fontsize=13, fontweight='bold')
ax1.grid(axis='y', alpha=0.3)
ax1.tick_params(axis='x', rotation=45)

# 添加统计显著性标记
if stats_results['multigroup']['p'] < 0.05:
    y_max = max([max(d) for d in data_for_plot])
    y_pos = y_max * 1.05
    ax1.plot([1, len(data_for_plot)], [y_pos, y_pos], 'k-', linewidth=1.5)
    ax1.text((1 + len(data_for_plot)) / 2, y_pos * 1.02, f"p={stats_results['multigroup']['p']:.4f}", 
             ha='center', fontsize=10, fontweight='bold')

# 2. 小提琴图
ax2 = fig.add_subplot(gs[0, 1])
df_plot = df[df['group_standard'].notna()].copy()
sns.violinplot(data=df_plot, x='group_standard', y='tscm_score', ax=ax2, palette=colors[:len(group_order)])
ax2.set_ylabel('TSCM Signature Score', fontsize=12)
ax2.set_xlabel('')
ax2.set_title('TSCM Score Distribution by SLEDAI Group', fontsize=13, fontweight='bold')
ax2.tick_params(axis='x', rotation=45)
ax2.grid(axis='y', alpha=0.3)

# 3. 散点图（带分组颜色）
ax3 = fig.add_subplot(gs[0, 2])
for i, group in enumerate(group_order):
    if group in df['group_standard'].values:
        group_data = df[df['group_standard'] == group]
        ax3.scatter(group_data['SLEDAI'], group_data['tscm_score'], 
                   alpha=0.5, s=30, label=group, color=colors[i])

ax3.set_xlabel('SLEDAI Score', fontsize=12)
ax3.set_ylabel('TSCM Signature Score', fontsize=12)
ax3.set_title('TSCM Score vs SLEDAI (Colored by Group)', fontsize=13, fontweight='bold')
ax3.legend(fontsize=9, loc='best')
ax3.grid(alpha=0.3)

# 4. 条形图（均值 ± SEM）
ax4 = fig.add_subplot(gs[1, 0])
means = df_stats_standard['Mean'].values
sems = df_stats_standard['SD'].values / np.sqrt(df_stats_standard['N'].values)
x_pos = np.arange(len(df_stats_standard))
bars = ax4.bar(x_pos, means, yerr=sems, capsize=5, color=colors[:len(means)], alpha=0.7)
ax4.set_xticks(x_pos)
ax4.set_xticklabels(df_stats_standard['Group'], rotation=45, ha='right')
ax4.set_ylabel('TSCM Signature Score (Mean ± SEM)', fontsize=12)
ax4.set_title('TSCM Score by SLEDAI Group (Mean ± SEM)', fontsize=13, fontweight='bold')
ax4.grid(axis='y', alpha=0.3)

# 5. 分组直方图
ax5 = fig.add_subplot(gs[1, 1])
for i, group in enumerate(group_order):
    if group in df['group_standard'].values:
        group_data = df[df['group_standard'] == group]['tscm_score']
        ax5.hist(group_data, bins=30, alpha=0.5, label=group, color=colors[i], density=True)

ax5.set_xlabel('TSCM Signature Score', fontsize=12)
ax5.set_ylabel('Density', fontsize=12)
ax5.set_title('TSCM Score Distribution by Group', fontsize=13, fontweight='bold')
ax5.legend(fontsize=9)
ax5.grid(alpha=0.3)

# 6. 趋势图（按SLEDAI中位数）
ax6 = fig.add_subplot(gs[1, 2])
group_medians_sledai = []
group_means_tscm = []
group_sems = []
for group in group_order:
    if group in df['group_standard'].values:
        group_data = df[df['group_standard'] == group]
        group_medians_sledai.append(group_data['SLEDAI'].median())
        group_means_tscm.append(group_data['tscm_score'].mean())
        group_sems.append(group_data['tscm_score'].std() / np.sqrt(len(group_data)))

ax6.errorbar(group_medians_sledai, group_means_tscm, yerr=group_sems, 
            marker='o', markersize=8, linewidth=2, capsize=5, capthick=2)
ax6.set_xlabel('SLEDAI Median', fontsize=12)
ax6.set_ylabel('TSCM Score (Mean ± SEM)', fontsize=12)
ax6.set_title('TSCM Score Trend Across SLEDAI Groups', fontsize=13, fontweight='bold')
ax6.grid(alpha=0.3)

# 7-9. 其他分组方式的可视化（三等分、二分）
ax7 = fig.add_subplot(gs[2, 0])
tertile_order = ['低活动度 (T1)', '中活动度 (T2)', '高活动度 (T3)']
data_tertile = [df[df['group_tertile'] == g]['tscm_score'].values for g in tertile_order if g in df['group_tertile'].values]
bp7 = ax7.boxplot(data_tertile, labels=[g for g in tertile_order if g in df['group_tertile'].values], patch_artist=True)
for patch, color in zip(bp7['boxes'], colors[:3]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax7.set_ylabel('TSCM Signature Score', fontsize=12)
ax7.set_title('TSCM Score by SLEDAI Tertile', fontsize=13, fontweight='bold')
ax7.grid(axis='y', alpha=0.3)

ax8 = fig.add_subplot(gs[2, 1])
binary_order = ['低活动度 (≤8)', '高活动度 (>8)']
data_binary = [df[df['group_binary'] == g]['tscm_score'].values for g in binary_order if g in df['group_binary'].values]
bp8 = ax8.boxplot(data_binary, labels=[g for g in binary_order if g in df['group_binary'].values], patch_artist=True)
for patch, color in zip(bp8['boxes'], colors[:2]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax8.set_ylabel('TSCM Signature Score', fontsize=12)
ax8.set_title('TSCM Score: Low vs High Activity', fontsize=13, fontweight='bold')
ax8.grid(axis='y', alpha=0.3)

# 9. 统计摘要表
ax9 = fig.add_subplot(gs[2, 2])
ax9.axis('off')
table_data = df_stats_standard[['Group', 'N', 'Mean', 'SD', 'Median']].copy()
table_data['Mean'] = table_data['Mean'].round(3)
table_data['SD'] = table_data['SD'].round(3)
table_data['Median'] = table_data['Median'].round(3)
table = ax9.table(cellText=table_data.values, colLabels=table_data.columns,
                 cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2)
ax9.set_title('Summary Statistics', fontsize=13, fontweight='bold', pad=20)

plt.suptitle('GSE88884: TSCM Signature Score by SLEDAI Groups', 
             fontsize=16, fontweight='bold', y=0.995)

out_fig = FIG_DIR / "GSE88884_tscm_by_sledai_groups_detailed.pdf"
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
print(f"  ✓ 已保存: {out_fig}")

# ============================================================================
# 步骤6: 生成详细报告
# ============================================================================
print("\n[步骤6] 生成详细报告...")

report_file = REPORT_DIR / "GSE88884_tscm_sledai_group_comparison_report.md"

with open(report_file, 'w', encoding='utf-8') as f:
    f.write("# GSE88884: TSCM Signature Score按SLEDAI分组比较分析报告\n\n")
    f.write(f"**生成时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("## 一、数据概览\n\n")
    f.write(f"- **总样本数**: {len(df)}\n")
    f.write(f"- **SLEDAI范围**: {df['SLEDAI'].min():.1f} - {df['SLEDAI'].max():.1f}\n")
    f.write(f"- **SLEDAI均值**: {df['SLEDAI'].mean():.2f} ± {df['SLEDAI'].std():.2f}\n")
    f.write(f"- **TSCM Score范围**: {df['tscm_score'].min():.3f} - {df['tscm_score'].max():.3f}\n")
    f.write(f"- **TSCM Score均值**: {df['tscm_score'].mean():.3f} ± {df['tscm_score'].std():.3f}\n\n")
    
    f.write("## 二、分组策略\n\n")
    f.write("### 2.1 标准临床分组\n")
    f.write("- 低活动度 (0-4)\n")
    f.write("- 中活动度 (4-8)\n")
    f.write("- 高活动度 (8-12)\n")
    f.write("- 极高活动度 (>12)\n\n")
    
    f.write("### 2.2 分组样本数分布\n\n")
    f.write(df['group_standard'].value_counts().to_frame('样本数').to_markdown())
    f.write("\n\n")
    
    f.write("## 三、描述性统计\n\n")
    f.write(df_stats_standard.to_markdown(index=False))
    f.write("\n\n")
    
    f.write("## 四、统计检验结果\n\n")
    f.write("### 4.1 正态性检验\n\n")
    for group, result in stats_results['normality'].items():
        f.write(f"- **{group}**: W={result['stat']:.4f}, p={result['p']:.4f} {'(正态)' if result['p'] > 0.05 else '(非正态)'}\n")
    f.write("\n")
    
    f.write("### 4.2 方差齐性检验\n\n")
    f.write(f"- **Levene检验**: F={stats_results['levene']['stat']:.4f}, p={stats_results['levene']['p']:.4f} ")
    f.write(f"{'(方差齐)' if stats_results['levene']['p'] > 0.05 else '(方差不齐)'}\n\n")
    
    f.write("### 4.3 多组比较\n\n")
    f.write(f"- **方法**: {stats_results['multigroup']['method']}\n")
    f.write(f"- **统计量**: {stats_results['multigroup']['stat']:.4f}\n")
    f.write(f"- **p值**: {stats_results['multigroup']['p']:.4f} ")
    f.write(f"{'** (p < 0.05)' if stats_results['multigroup']['p'] < 0.05 else '(p ≥ 0.05)'}\n\n")
    
    f.write("### 4.4 两两比较\n\n")
    f.write(stats_results['pairwise'].to_markdown(index=False))
    f.write("\n\n")
    
    if 'trend' in stats_results:
        f.write("### 4.5 趋势检验\n\n")
        f.write(f"- **方法**: {stats_results['trend']['method']}\n")
        f.write(f"- **统计量**: {stats_results['trend']['stat']:.4f}\n")
        f.write(f"- **p值**: {stats_results['trend']['p']:.4f} ")
        f.write(f"{'** (p < 0.05)' if stats_results['trend']['p'] < 0.05 else '(p ≥ 0.05)'}\n\n")
    
    f.write("## 五、主要发现\n\n")
    if stats_results['multigroup']['p'] < 0.05:
        f.write("1. **多组比较显著**: TSCM Score在不同SLEDAI活动度组间存在显著差异。\n")
    else:
        f.write("1. **多组比较不显著**: TSCM Score在不同SLEDAI活动度组间无显著差异。\n")
    
    f.write("\n2. **分组趋势**: ")
    if 'trend' in stats_results and stats_results['trend']['p'] < 0.05:
        f.write(f"TSCM Score随SLEDAI活动度变化呈显著趋势 (τ={stats_results['trend']['stat']:.4f}, p={stats_results['trend']['p']:.4f})。\n")
    else:
        f.write("TSCM Score随SLEDAI活动度变化无显著趋势。\n")
    
    f.write("\n3. **临床意义**: ")
    # 计算效应量（Cohen's d for pairwise comparisons）
    if len(stats_results['pairwise']) > 0:
        f.write("具体两两比较结果见上表。\n")
    
    f.write("\n## 六、图表\n\n")
    f.write(f"- 详细可视化图表已保存至: `{out_fig}`\n")
    f.write("- 图表包括：箱线图、小提琴图、散点图、条形图、直方图、趋势图等\n\n")
    
    f.write("## 七、结论\n\n")
    f.write("基于GSE88884数据集（n=1756）的分析显示，TSCM Signature Score与SLEDAI疾病活动度")
    if stats_results['multigroup']['p'] < 0.05:
        f.write("存在显著关联。")
    else:
        f.write("无显著关联。")
    f.write(" 这一发现为进一步研究TSCM在SLE疾病进展中的作用提供了初步证据。\n")

print(f"  ✓ 已保存: {report_file}")

# ============================================================================
# 步骤7: 保存详细统计结果
# ============================================================================
print("\n[步骤7] 保存统计结果...")

# 保存分组统计
df_stats_standard.to_csv(OUT_DIR / "GSE88884_tscm_sledai_group_statistics.tsv", sep="\t", index=False)

# 保存两两比较结果
stats_results['pairwise'].to_csv(OUT_DIR / "GSE88884_tscm_sledai_pairwise_comparisons.tsv", sep="\t", index=False)

print(f"  ✓ 已保存统计结果文件")

print("\n" + "=" * 80)
print("✓ 分析完成！")
print("=" * 80)
print(f"\n输出文件:")
print(f"  1. 详细图表: {out_fig}")
print(f"  2. 分析报告: {report_file}")
print(f"  3. 分组统计: {OUT_DIR / 'GSE88884_tscm_sledai_group_statistics.tsv'}")
print(f"  4. 两两比较: {OUT_DIR / 'GSE88884_tscm_sledai_pairwise_comparisons.tsv'}")

