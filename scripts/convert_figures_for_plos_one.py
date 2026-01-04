#!/usr/bin/env python3
"""
Convert figures to PLOS ONE format (TIFF, 300 DPI, RGB)

PLOS ONE requirements:
- Format: TIFF (preferred) or EPS
- Resolution: Color ≥300 DPI, Grayscale ≥500 DPI, Line art ≥1000 DPI
- Color mode: RGB
- File size: <20 MB
- Fonts: Arial, Helvetica, or Times New Roman, minimum 8pt
"""

import sys
import os
from pathlib import Path
from PIL import Image
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results/figures"
OUT_DIR = ROOT / "results/figures/plos_one_format"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("Converting figures to PLOS ONE format")
print("=" * 80)

# Required figures for paper
REQUIRED_FIGURES = {
    'Figure1': [
        'stage4_umap_clusters.pdf',  # UMAP with clusters
        'stage4_tscm_scores.pdf',     # TSCM scores
        'stage4_marker_genes.pdf'     # Marker genes
    ],
    'Figure2': [
        'tscm_volcano_plot.pdf',      # Volcano plot
        'tscm_enrichment_go_barplot.pdf',  # GO enrichment
        'tscm_enrichment_kegg_barplot.pdf'  # KEGG enrichment
    ],
    'Figure3': [
        # Distribution plot (may need to generate)
    ],
    'Figure4': [
        'GSE88884_tscm_vs_sledai.pdf',  # Scatter + box plots
        'GSE88884_tscm_by_sledai_groups_detailed.pdf'  # Detailed group comparison
    ]
}

def convert_pdf_to_tiff(pdf_path, output_path, dpi=300):
    """Convert PDF to TIFF using ImageMagick or pdftoppm"""
    try:
        # Try ImageMagick first
        cmd = [
            'convert',
            '-density', str(dpi),
            '-quality', '100',
            '-colorspace', 'RGB',
            pdf_path,
            output_path
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    
    try:
        # Try pdftoppm + convert
        temp_png = output_path.with_suffix('.temp.png')
        cmd = [
            'pdftoppm',
            '-png',
            '-r', str(dpi),
            pdf_path,
            str(temp_png).replace('.temp.png', '')
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            # Find the generated PNG
            png_files = list(temp_png.parent.glob(temp_png.stem + '*.png'))
            if png_files:
                img = Image.open(png_files[0])
                if img.mode != 'RGB':
                    img = img.convert('RGB')
                img.save(output_path, 'TIFF', compression='tiff_lzw', dpi=(dpi, dpi))
                # Clean up
                for png in png_files:
                    png.unlink()
                return True
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    
    return False

def check_figure_requirements(tiff_path):
    """Check if figure meets PLOS ONE requirements"""
    try:
        img = Image.open(tiff_path)
        width, height = img.size
        dpi_x, dpi_y = img.info.get('dpi', (0, 0))
        
        # Calculate actual DPI if not in metadata
        if dpi_x == 0:
            # Assume A4 size (8.27 x 11.69 inches) for estimation
            dpi_x = width / 8.27
            dpi_y = height / 11.69
        
        issues = []
        if dpi_x < 300 or dpi_y < 300:
            issues.append(f"DPI too low: {dpi_x:.1f} x {dpi_y:.1f} (required: ≥300)")
        
        if img.mode != 'RGB':
            issues.append(f"Color mode: {img.mode} (required: RGB)")
        
        file_size_mb = tiff_path.stat().st_size / (1024 * 1024)
        if file_size_mb > 20:
            issues.append(f"File size too large: {file_size_mb:.1f} MB (required: <20 MB)")
        
        return {
            'valid': len(issues) == 0,
            'dpi': (dpi_x, dpi_y),
            'mode': img.mode,
            'size_mb': file_size_mb,
            'dimensions': (width, height),
            'issues': issues
        }
    except Exception as e:
        return {'valid': False, 'error': str(e)}

# Convert figures
converted_files = []
failed_files = []

print("\n[步骤1] 转换必需图表...")

for fig_name, pdf_files in REQUIRED_FIGURES.items():
    if not pdf_files:
        continue
    
    print(f"\n处理 {fig_name}:")
    for pdf_file in pdf_files:
        pdf_path = FIG_DIR / pdf_file
        if not pdf_path.exists():
            print(f"  ⚠ {pdf_file}: 文件不存在")
            failed_files.append(pdf_file)
            continue
        
        # Output filename
        tiff_name = pdf_file.replace('.pdf', '.tiff')
        tiff_path = OUT_DIR / tiff_name
        
        print(f"  转换: {pdf_file} -> {tiff_name}")
        
        if convert_pdf_to_tiff(pdf_path, tiff_path, dpi=300):
            # Check requirements
            check_result = check_figure_requirements(tiff_path)
            if check_result.get('valid'):
                print(f"    ✓ 转换成功")
                print(f"      DPI: {check_result['dpi'][0]:.1f} x {check_result['dpi'][1]:.1f}")
                print(f"      大小: {check_result['size_mb']:.2f} MB")
                converted_files.append((pdf_file, tiff_path))
            else:
                print(f"    ⚠ 转换成功但可能不符合要求:")
                for issue in check_result.get('issues', []):
                    print(f"      - {issue}")
                converted_files.append((pdf_file, tiff_path))
        else:
            print(f"    ✗ 转换失败")
            failed_files.append(pdf_file)

print("\n" + "=" * 80)
print("转换完成")
print("=" * 80)
print(f"\n成功转换: {len(converted_files)} 个文件")
print(f"失败: {len(failed_files)} 个文件")

if converted_files:
    print("\n转换后的文件:")
    for pdf_file, tiff_path in converted_files:
        size_mb = tiff_path.stat().st_size / (1024 * 1024)
        print(f"  {tiff_path.name} ({size_mb:.2f} MB)")

if failed_files:
    print("\n失败的文件:")
    for f in failed_files:
        print(f"  {f}")

print(f"\n输出目录: {OUT_DIR}")
print("\n注意: 请检查转换后的TIFF文件是否符合PLOS ONE要求")
print("如果DPI不足，可能需要使用更高分辨率的源文件或重新生成图表")

