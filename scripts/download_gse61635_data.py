#!/usr/bin/env python3
"""
GSE61635 数据下载和处理脚本

支持单细胞RNA-seq数据的下载和处理
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import gzip
import urllib.request
import time
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/raw/GSE61635"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = ROOT / "logs/validation/gse61635_download.log"
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

def log(message):
    """记录日志"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}\n"
    print(message)
    with open(LOG_FILE, 'a') as f:
        f.write(log_msg)

def download_geo_matrix(gse_id):
    """下载GEO series matrix文件"""
    log("=" * 60)
    log("下载GEO series matrix文件")
    log("=" * 60)
    
    url_patterns = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
    ]
    
    matrix_file = f"{gse_id}_series_matrix.txt.gz"
    out_file = OUT_DIR / matrix_file
    
    if out_file.exists() and out_file.stat().st_size > 1000:
        log(f"✓ 文件已存在: {out_file} ({out_file.stat().st_size / 1024 / 1024:.2f} MB)")
        return out_file
    
    for url in url_patterns:
        try:
            log(f"尝试下载: {url}")
            cmd = [
                "wget",
                "--tries=5",
                "--timeout=60",
                "--limit-rate=2m",
                "--continue",
                "-O", str(out_file),
                url
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode == 0 and out_file.exists() and out_file.stat().st_size > 1000:
                log(f"✓ 下载成功: {out_file.stat().st_size / 1024 / 1024:.2f} MB")
                return out_file
            else:
                if out_file.exists():
                    out_file.unlink()
                log(f"✗ 下载失败: {result.stderr}")
        except Exception as e:
            log(f"✗ 下载异常: {e}")
    
    return None

def check_supplementary_files(gse_id):
    """检查并下载补充文件（单细胞数据通常在suppl目录）"""
    log("=" * 60)
    log("检查补充文件（单细胞数据）")
    log("=" * 60)
    
    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/{gse_id}/suppl/"
    
    # 常见的单细胞数据文件格式
    possible_files = [
        f"{gse_id}_*.h5",
        f"{gse_id}_*.h5ad",
        f"{gse_id}_*.mtx.gz",
        f"{gse_id}_*.tsv.gz",
        f"{gse_id}_*.csv.gz",
    ]
    
    log(f"补充文件URL: {suppl_url}")
    log("提示: 单细胞数据可能需要手动下载或使用其他工具（如GEOquery）")
    
    return None

def parse_geo_matrix_for_metadata(matrix_file):
    """从series matrix文件中提取metadata"""
    log("=" * 60)
    log("解析GEO matrix文件提取metadata")
    log("=" * 60)
    
    if not matrix_file or not matrix_file.exists():
        return None
    
    try:
        with gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        sample_info = {}
        sample_titles = []
        sample_conditions = []
        sample_geo_accessions = []
        
        for line in lines:
            if line.startswith("!Sample_title"):
                parts = line.strip().split("\t")
                if len(parts) > 1:
                    sample_titles = parts[1:]
            elif line.startswith("!Sample_geo_accession"):
                parts = line.strip().split("\t")
                if len(parts) > 1:
                    sample_geo_accessions = parts[1:]
            elif "!Sample_characteristics" in line or "diagnosis" in line.lower() or "condition" in line.lower():
                parts = line.strip().split("\t")
                if len(parts) > 1:
                    for part in parts[1:]:
                        part_upper = part.upper()
                        if "SLE" in part_upper or "LUPUS" in part_upper:
                            sample_conditions.append("SLE")
                        elif "HEALTHY" in part_upper or "CONTROL" in part_upper or "CTRL" in part_upper:
                            sample_conditions.append("Healthy")
                        else:
                            sample_conditions.append("")
        
        # 创建metadata DataFrame
        if sample_titles:
            n_samples = len(sample_titles)
            df_meta = pd.DataFrame({
                'sample_id': sample_titles[:n_samples],
                'geo_accession': sample_geo_accessions[:n_samples] if sample_geo_accessions else [''] * n_samples,
                'condition': sample_conditions[:n_samples] if sample_conditions else ['Unknown'] * n_samples,
            })
            
            log(f"✓ 提取了 {len(df_meta)} 个样本的metadata")
            log(f"  Condition分布: {df_meta['condition'].value_counts().to_dict()}")
            
            # 保存metadata
            out_meta_file = OUT_DIR / f"{Path(matrix_file).stem.replace('_series_matrix', '')}_metadata.tsv"
            df_meta.to_csv(out_meta_file, sep="\t", index=False)
            log(f"✓ 已保存metadata: {out_meta_file}")
            
            return df_meta
        
    except Exception as e:
        log(f"✗ 解析失败: {e}")
        import traceback
        log(traceback.format_exc())
    
    return None

def main():
    log("=" * 60)
    log("GSE61635 数据下载和处理")
    log("=" * 60)
    
    gse_id = "GSE61635"
    
    # 下载series matrix文件（用于获取metadata）
    matrix_file = download_geo_matrix(gse_id)
    
    # 提取metadata
    if matrix_file:
        df_meta = parse_geo_matrix_for_metadata(matrix_file)
    
    # 检查补充文件
    check_supplementary_files(gse_id)
    
    log("=" * 60)
    log("下载完成")
    log("=" * 60)
    log("")
    log("注意: 单细胞数据可能需要:")
    log("1. 使用GEOquery R包下载")
    log("2. 或手动从GEO网站下载")
    log("3. 或使用其他工具（如SRA toolkit）")
    log("")
    log("下一步: 检查 data/raw/GSE61635/ 目录中的文件")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

