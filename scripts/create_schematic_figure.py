#!/usr/bin/env python3
"""
创建研究流程图（Schematic Figure）
SLE PBMC (单细胞) → T细胞 → TSCM cluster → pseudo-bulk → DE → 富集 → 构建TSCM signature → 外部队列验证
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, ConnectionPatch
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
from PIL import Image
import io

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results/figures/plos_one_format_optimized"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 设置字体
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# 定义颜色（仅用于边框和文本，填充透明）
color_pbmc = '#E8F4F8'
color_tcell = '#B8E0D2'
color_tscm = '#F4A261'
color_analysis = '#E76F51'
color_validation = '#2A9D8F'
color_text = '#264653'

# 方框尺寸（缩小）
box_width = 1.2
box_height = 1.2

# 步骤1: SLE PBMC (单细胞)
box1 = FancyBboxPatch((0.5, 7), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box1)
ax.text(1.1, 7.6, 'SLE PBMC\n(scRNA-seq)', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(1.1, 7.2, 'GSE137029', ha='center', va='center', 
        fontsize=8, style='italic', color=color_text)

# 箭头1（调整位置以适应缩小后的方框）
arrow1 = FancyArrowPatch((1.7, 7.6), (2.5, 7.6), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow1)

# 步骤2: T细胞分离
box2 = FancyBboxPatch((2.5, 7), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box2)
ax.text(3.1, 7.6, 'T-cell\nIsolation', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)

# 箭头2
arrow2 = FancyArrowPatch((3.7, 7.6), (4.5, 7.6), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow2)

# 步骤3: TSCM cluster识别
box3 = FancyBboxPatch((4.5, 7), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box3)
ax.text(5.1, 7.6, 'TSCM\nCluster\nIdentification', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(5.1, 7.2, 'Clusters 5, 13, 16, 17', ha='center', va='center', 
        fontsize=7, color=color_text)

# 箭头3 (向下)
arrow3 = FancyArrowPatch((5.1, 7), (5.1, 5.9), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow3)

# 步骤4: Pseudo-bulk
box4 = FancyBboxPatch((4.2, 5.4), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box4)
ax.text(4.8, 6.0, 'Pseudo-bulk\nAggregation', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(4.8, 5.6, 'TSCM vs Others', ha='center', va='center', 
        fontsize=8, color=color_text)

# 箭头4
arrow4 = FancyArrowPatch((4.2, 6.0), (3.4, 6.0), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow4)

# 步骤5: 差异表达分析
box5 = FancyBboxPatch((2.2, 5.4), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box5)
ax.text(2.8, 6.0, 'Differential\nExpression', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(2.8, 5.6, 'Wilcoxon test', ha='center', va='center', 
        fontsize=8, color=color_text)

# 箭头5
arrow5 = FancyArrowPatch((2.2, 6.0), (1.4, 6.0), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow5)

# 步骤6: 富集分析
box6 = FancyBboxPatch((0.2, 5.4), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box6)
ax.text(0.8, 6.0, 'Pathway\nEnrichment', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(0.8, 5.6, 'GO/KEGG/Reactome', ha='center', va='center', 
        fontsize=8, color=color_text)

# 箭头6 (向下)
arrow6 = FancyArrowPatch((0.8, 5.4), (0.8, 4.2), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow6)

# 步骤7: 构建TSCM signature
box7 = FancyBboxPatch((0.2, 3.7), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box7)
ax.text(0.8, 4.3, 'TSCM\nSignature\nDefinition', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(0.8, 3.9, 'Core: 49 genes\nFull: 386 genes', ha='center', va='center', 
        fontsize=7, color=color_text)

# 箭头7
arrow7 = FancyArrowPatch((1.4, 4.0), (2.2, 4.0), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow7)

# 步骤8: 外部队列验证
box8 = FancyBboxPatch((2.2, 3.7), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box8)
ax.text(2.8, 4.3, 'External\nValidation', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(2.8, 3.9, 'GSE88884\n(n=1,756)', ha='center', va='center', 
        fontsize=8, color=color_text)

# 箭头8
arrow8 = FancyArrowPatch((3.4, 4.0), (4.2, 4.0), 
                         arrowstyle='->', lw=2, color=color_text)
ax.add_patch(arrow8)

# 步骤9: SLEDAI关联分析
box9 = FancyBboxPatch((4.2, 3.7), box_width, box_height, 
                      boxstyle="round,pad=0.1", 
                      facecolor='none',  # 透明填充
                      edgecolor=color_text, 
                      linewidth=2)
ax.add_patch(box9)
ax.text(4.8, 4.3, 'SLEDAI\nAssociation', ha='center', va='center', 
        fontsize=10, fontweight='bold', color=color_text)
ax.text(4.8, 3.9, 'Correlation\nAnalysis', ha='center', va='center', 
        fontsize=8, color=color_text)

# 添加标题
ax.text(5, 9.5, 'Figure S1. Study Workflow and Analysis Pipeline', 
        ha='center', va='center', fontsize=14, fontweight='bold', color=color_text)

# 添加说明文字
ax.text(7.5, 8.5, 'Key Steps:', ha='left', va='top', 
        fontsize=11, fontweight='bold', color=color_text)
ax.text(7.5, 8, '1. Single-cell RNA-seq data from SLE PBMC', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 7.5, '2. T-cell isolation and clustering', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 7, '3. TSCM-enriched cluster identification', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 6.5, '4. Pseudo-bulk aggregation for DE analysis', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 6, '5. Pathway enrichment analysis', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 5.5, '6. TSCM signature definition (core vs full)', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 5, '7. Validation in independent bulk dataset', ha='left', va='top', 
        fontsize=9, color=color_text)
ax.text(7.5, 4.5, '8. Association with SLEDAI disease activity', ha='left', va='top', 
        fontsize=9, color=color_text)

# 添加图例（所有填充透明，仅用边框颜色区分）
legend_elements = [
    mpatches.Patch(facecolor='none', edgecolor=color_text, linewidth=2, label='Data Source / Cell Isolation'),
    mpatches.Patch(facecolor='none', edgecolor=color_text, linewidth=2, label='TSCM Identification'),
    mpatches.Patch(facecolor='none', edgecolor=color_text, linewidth=2, label='Analysis'),
    mpatches.Patch(facecolor='none', edgecolor=color_text, linewidth=2, label='Validation')
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9)

plt.tight_layout()

# 保存图片
out_fig = FIG_DIR / "figure_s1_workflow_schematic.tiff"
buf = io.BytesIO()
plt.savefig(buf, dpi=300, bbox_inches='tight', format='png')
buf.seek(0)
img = Image.open(buf)
img = img.convert('RGB')
img.save(out_fig, 'TIFF', compression='tiff_lzw', dpi=(300, 300))
buf.close()
print(f"已保存流程图: {out_fig}")

plt.close()


