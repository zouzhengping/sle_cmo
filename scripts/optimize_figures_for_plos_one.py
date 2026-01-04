#!/usr/bin/env python3
"""
优化图表文件大小，确保符合PLOS ONE要求（<20 MB, 300 DPI, RGB）

PLOS ONE要求:
- 格式: TIFF
- 分辨率: ≥300 DPI (彩色)
- 颜色模式: RGB
- 文件大小: <20 MB
- 压缩: LZW
"""

import sys
from pathlib import Path
from PIL import Image
import subprocess

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "results/figures/plos_one_format"
OUTPUT_DIR = ROOT / "results/figures/plos_one_format_optimized"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_SIZE_MB = 20
TARGET_DPI = 300

print("=" * 80)
print("优化图表文件大小 - PLOS ONE格式")
print("=" * 80)
print(f"输入目录: {INPUT_DIR}")
print(f"输出目录: {OUTPUT_DIR}")
print(f"目标大小: <{TARGET_SIZE_MB} MB")
print(f"目标DPI: {TARGET_DPI}")
print()

def optimize_tiff(input_path, output_path):
    """优化TIFF文件"""
    try:
        img = Image.open(input_path)
        original_size_mb = input_path.stat().st_size / (1024 * 1024)
        
        # 转换为RGB（如果还不是）
        if img.mode != 'RGB':
            print(f"  转换颜色模式: {img.mode} -> RGB")
            img = img.convert('RGB')
        
        # 获取当前DPI信息
        dpi_x, dpi_y = img.info.get('dpi', (72, 72))
        
        # 保存时使用LZW压缩，设置DPI为300
        img.save(output_path, 'TIFF', compression='tiff_lzw', dpi=(TARGET_DPI, TARGET_DPI))
        
        new_size_mb = output_path.stat().st_size / (1024 * 1024)
        
        # 检查结果
        status = "✓" if new_size_mb < TARGET_SIZE_MB else "⚠"
        print(f"{status} {input_path.name}")
        print(f"   原始: {original_size_mb:.1f} MB (DPI: {dpi_x:.0f})")
        print(f"   优化: {new_size_mb:.1f} MB (DPI: {TARGET_DPI})")
        
        if new_size_mb > TARGET_SIZE_MB:
            print(f"   ⚠️  仍然超过{TARGET_SIZE_MB}MB，可能需要进一步优化（降低分辨率或裁剪）")
        elif new_size_mb < original_size_mb * 0.5:
            print(f"   ✓ 文件大小减少了 {(1 - new_size_mb/original_size_mb)*100:.1f}%")
        
        return new_size_mb < TARGET_SIZE_MB
        
    except Exception as e:
        print(f"✗ {input_path.name}: 错误 - {e}")
        return False

# 处理所有TIFF文件
tiff_files = list(INPUT_DIR.glob("*.tiff"))
if not tiff_files:
    print(f"未找到TIFF文件在 {INPUT_DIR}")
    sys.exit(1)

print(f"找到 {len(tiff_files)} 个TIFF文件\n")

success_count = 0
fail_count = 0

for tiff_file in sorted(tiff_files):
    output_file = OUTPUT_DIR / tiff_file.name
    if optimize_tiff(tiff_file, output_file):
        success_count += 1
    else:
        fail_count += 1
    print()

print("=" * 80)
print("优化完成")
print("=" * 80)
print(f"成功: {success_count}/{len(tiff_files)}")
print(f"失败: {fail_count}/{len(tiff_files)}")
print(f"\n优化后的文件保存在: {OUTPUT_DIR}")
print("\n注意:")
print("  - 如果文件仍然过大，可能需要:")
print("    1. 降低DPI到300（最低要求）")
print("    2. 裁剪不必要的空白区域")
print("    3. 将大图拆分为多个子图")

