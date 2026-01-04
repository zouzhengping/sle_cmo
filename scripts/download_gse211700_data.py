#!/usr/bin/env python3
"""
GSE211700 数据下载和处理脚本

支持多种下载方式：
1. 从GEO下载series matrix文件
2. 使用GEOparse解析
3. 从SRA下载FASTQ（如果需要）

输出处理好的表达矩阵到 data/raw/GSE211700/
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
OUT_DIR = ROOT / "data/raw/GSE211700"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = ROOT / "logs/validation/gse211700_download.log"
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

def log(message):
    """记录日志"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}\n"
    print(message)
    with open(LOG_FILE, 'a') as f:
        f.write(log_msg)

def download_geo_matrix_method1(gse_id):
    """方法1: 直接从GEO FTP下载"""
    log("=" * 60)
    log("方法1: 从GEO FTP下载series matrix文件")
    log("=" * 60)
    
    # 尝试多种URL格式
    url_patterns = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
    ]
    
    out_file = OUT_DIR / f"{gse_id}_series_matrix.txt.gz"
    
    if out_file.exists() and out_file.stat().st_size > 1000:
        log(f"✓ 文件已存在: {out_file} ({out_file.stat().st_size / 1024 / 1024:.2f} MB)")
        return out_file
    
    for url in url_patterns:
        try:
            log(f"尝试下载: {url}")
            urllib.request.urlretrieve(url, out_file)
            if out_file.stat().st_size > 1000:
                log(f"✓ 下载成功: {out_file} ({out_file.stat().st_size / 1024 / 1024:.2f} MB)")
                return out_file
            else:
                out_file.unlink()
                log(f"✗ 下载的文件太小，可能不完整")
        except Exception as e:
            log(f"✗ 下载失败: {e}")
            continue
    
    return None

def download_geo_matrix_method2(gse_id):
    """方法2: 使用wget下载（更稳定）"""
    log("=" * 60)
    log("方法2: 使用wget下载")
    log("=" * 60)
    
    out_file = OUT_DIR / f"{gse_id}_series_matrix.txt.gz"
    
    if out_file.exists() and out_file.stat().st_size > 1000:
        log(f"✓ 文件已存在")
        return out_file
    
    urls = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
    ]
    
    for url in urls:
        try:
            log(f"使用wget下载: {url}")
            cmd = [
                "wget",
                "--tries=3",
                "--timeout=30",
                "--limit-rate=2m",
                "-O", str(out_file),
                url
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0 and out_file.exists() and out_file.stat().st_size > 1000:
                log(f"✓ wget下载成功: {out_file.stat().st_size / 1024 / 1024:.2f} MB")
                return out_file
            else:
                if out_file.exists():
                    out_file.unlink()
                log(f"✗ wget下载失败: {result.stderr}")
        except Exception as e:
            log(f"✗ wget下载异常: {e}")
            continue
    
    return None

def parse_geo_matrix_correct(matrix_file):
    """正确解析GEO matrix文件"""
    log("=" * 60)
    log("解析GEO matrix文件")
    log("=" * 60)
    
    log(f"读取文件: {matrix_file}")
    
    # 读取文件内容
    with gzip.open(matrix_file, 'rt', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    log(f"文件总行数: {len(lines)}")
    
    # 查找关键位置
    data_start_idx = None
    data_end_idx = None
    sample_titles = []
    sample_geo_accessions = []
    sample_conditions = []
    
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            data_start_idx = i + 1
        elif line.startswith("!series_matrix_table_end"):
            data_end_idx = i
        elif line.startswith("!Sample_title"):
            parts = line.strip().split("\t")
            if len(parts) > 1:
                sample_titles = parts[1:]
                log(f"找到 {len(sample_titles)} 个样本标题")
        elif line.startswith("!Sample_geo_accession"):
            parts = line.strip().split("\t")
            if len(parts) > 1:
                sample_geo_accessions = parts[1:]
                log(f"找到 {len(sample_geo_accessions)} 个GEO accession")
        elif line.startswith("!Sample_characteristics_ch1"):
            # 尝试提取condition信息
            parts = line.strip().split("\t")
            if len(parts) > 1:
                for part in parts[1:]:
                    if "disease" in part.lower() or "condition" in part.lower() or "sle" in part.lower():
                        sample_conditions.append(part)
                    elif "control" in part.lower() or "healthy" in part.lower() or "normal" in part.lower():
                        sample_conditions.append("Healthy")
                    else:
                        sample_conditions.append("")
    
    if data_start_idx is None:
        raise ValueError("无法找到数据开始位置")
    
    log(f"数据开始行: {data_start_idx}, 结束行: {data_end_idx}")
    
    # 读取表头行（ID_REF）
    header_line = lines[data_start_idx]
    header_parts = header_line.strip().split("\t")
    log(f"表头列数: {len(header_parts)}")
    
    # 读取数据
    data_rows = []
    gene_ids = []
    
    for line in lines[data_start_idx + 1:data_end_idx]:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        
        gene_id = parts[0]
        gene_ids.append(gene_id)
        
        # 转换数值
        values = []
        for v in parts[1:]:
            try:
                values.append(float(v))
            except:
                values.append(np.nan)
        data_rows.append(values)
    
    log(f"解析了 {len(gene_ids)} 个基因, {len(data_rows[0]) if data_rows else 0} 个样本")
    
    # 创建DataFrame
    sample_names = header_parts[1:] if len(header_parts) > 1 else [f"Sample_{i+1}" for i in range(len(data_rows[0]) if data_rows else 0)]
    
    # 确保样本名称数量匹配
    if len(sample_names) != len(data_rows[0]) if data_rows else 0:
        log(f"警告: 样本名称数量({len(sample_names)})与数据列数({len(data_rows[0]) if data_rows else 0})不匹配，使用默认名称")
        sample_names = [f"Sample_{i+1}" for i in range(len(data_rows[0]) if data_rows else 0)]
    
    df_expr = pd.DataFrame(data_rows, index=gene_ids, columns=sample_names)
    df_expr = df_expr.T  # 转置：样本x基因
    
    log(f"✓ 表达矩阵创建成功: {df_expr.shape[0]} 个样本 x {df_expr.shape[1]} 个基因")
    
    # 保存处理好的表达矩阵
    out_expr_file = OUT_DIR / f"{Path(matrix_file).stem.replace('_series_matrix', '')}_expression_matrix.tsv"
    df_expr.to_csv(out_expr_file, sep="\t", index=True)
    log(f"✓ 已保存表达矩阵: {out_expr_file}")
    
    # 保存样本信息
    if sample_titles or sample_geo_accessions:
        df_meta = pd.DataFrame({
            'sample_id': sample_names,
            'geo_accession': sample_geo_accessions[:len(sample_names)] if sample_geo_accessions else [''] * len(sample_names),
            'title': sample_titles[:len(sample_names)] if sample_titles else [''] * len(sample_names),
            'condition': sample_conditions[:len(sample_names)] if sample_conditions else ['Unknown'] * len(sample_names),
        })
        
        out_meta_file = OUT_DIR / f"{Path(matrix_file).stem.replace('_series_matrix', '')}_metadata.tsv"
        df_meta.to_csv(out_meta_file, sep="\t", index=False)
        log(f"✓ 已保存样本信息: {out_meta_file}")
    
    return df_expr, out_expr_file

def download_fpkm_file(gse_id):
    """下载FPKM表达数据文件"""
    log("=" * 60)
    log("下载FPKM表达数据文件")
    log("=" * 60)
    
    # FPKM文件URL（从series matrix中提取）
    fpkm_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/{gse_id}/suppl/{gse_id}_Transcript_FPKM.txt.gz"
    out_file = OUT_DIR / f"{gse_id}_Transcript_FPKM.txt.gz"
    
    if out_file.exists() and out_file.stat().st_size > 1000:
        log(f"✓ FPKM文件已存在: {out_file} ({out_file.stat().st_size / 1024 / 1024:.2f} MB)")
        return out_file
    
    # 尝试下载
    try:
        log(f"下载FPKM文件: {fpkm_url}")
        cmd = [
            "wget",
            "--tries=5",
            "--timeout=60",
            "--limit-rate=2m",
            "--continue",
            "-O", str(out_file),
            fpkm_url
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0 and out_file.exists() and out_file.stat().st_size > 1000:
            log(f"✓ FPKM文件下载成功: {out_file.stat().st_size / 1024 / 1024:.2f} MB")
            return out_file
        else:
            log(f"✗ FPKM文件下载失败: {result.stderr}")
            if out_file.exists():
                out_file.unlink()
    except Exception as e:
        log(f"✗ FPKM文件下载异常: {e}")
    
    return None

def parse_fpkm_file(fpkm_file):
    """解析FPKM文件（处理包含metadata行的特殊格式）"""
    log("=" * 60)
    log("解析FPKM文件")
    log("=" * 60)
    
    log(f"读取文件: {fpkm_file}")
    
    try:
        # 直接读取，第一行是表头
        df_raw = pd.read_csv(fpkm_file, sep="\t", compression='gzip')
        log(f"原始数据形状: {df_raw.shape}")
        log(f"列名: {list(df_raw.columns[:5])}...")
        
        # 检查列名
        if "Official_Symbol" not in df_raw.columns:
            raise ValueError("未找到Official_Symbol列")
        
        # 提取基因符号和表达数据
        # 样本列从"Gene_type"之后开始（或从第5列开始，如果前4列是ID信息）
        sample_cols = [col for col in df_raw.columns if col not in 
                      ["Transcript_id", "Gene_id", "Official_Symbol", "Gene_type"]]
        
        log(f"找到 {len(sample_cols)} 个样本列")
        
        # 过滤掉Official_Symbol为"--"或空的基因
        valid_mask = (df_raw["Official_Symbol"] != "--") & (df_raw["Official_Symbol"].notna())
        df_valid = df_raw[valid_mask].copy()
        log(f"有效基因数（Official_Symbol不为'--'）: {len(df_valid)}/{len(df_raw)}")
        
        # 提取基因符号和表达数据
        gene_names = df_valid["Official_Symbol"].values
        expr_data = df_valid[sample_cols].values
        
        # 创建表达矩阵（基因x样本）
        df_expr_genes = pd.DataFrame(expr_data, index=gene_names, columns=sample_cols)
        
        # 如果有重复的基因名，取平均值
        if df_expr_genes.index.duplicated().any():
            n_duplicates = df_expr_genes.index.duplicated().sum()
            log(f"发现 {n_duplicates} 个重复基因名，进行聚合（取平均值）...")
            df_expr_genes = df_expr_genes.groupby(df_expr_genes.index).mean()
            log(f"聚合后: {df_expr_genes.shape[0]} 个唯一基因")
        
        # 转置：样本x基因
        df_expr = df_expr_genes.T
        
        log(f"✓ 处理完成: {df_expr.shape[0]} 个样本 x {df_expr.shape[1]} 个基因")
        
        # 保存
        out_expr_file = OUT_DIR / f"{Path(fpkm_file).stem.replace('_Transcript_FPKM', '')}_expression_matrix.tsv"
        df_expr.to_csv(out_expr_file, sep="\t", index=True)
        log(f"✓ 已保存表达矩阵: {out_expr_file}")
        
        return df_expr, out_expr_file
            
    except Exception as e:
        log(f"✗ 解析FPKM文件失败: {e}")
        import traceback
        log(traceback.format_exc())
        return None, None

def main():
    log("=" * 60)
    log("GSE211700 数据下载和处理")
    log("=" * 60)
    
    gse_id = "GSE211700"
    
    # 优先尝试下载FPKM文件（包含实际表达数据）
    fpkm_file = download_fpkm_file(gse_id)
    
    if fpkm_file and fpkm_file.exists():
        try:
            df_expr, expr_file = parse_fpkm_file(fpkm_file)
            if df_expr is not None:
                log("=" * 60)
                log("✓ 数据处理完成")
                log("=" * 60)
                log(f"表达矩阵: {expr_file}")
                log(f"矩阵维度: {df_expr.shape[0]} 样本 x {df_expr.shape[1]} 个基因/转录本")
                return 0
        except Exception as e:
            log(f"✗ 处理FPKM文件失败: {e}")
    
    # 如果FPKM文件失败，尝试series matrix
    log("")
    log("尝试下载series matrix文件...")
    matrix_file = None
    
    # 方法1: 直接下载
    matrix_file = download_geo_matrix_method1(gse_id)
    
    # 方法2: wget下载
    if not matrix_file or not matrix_file.exists() or matrix_file.stat().st_size < 1000:
        matrix_file = download_geo_matrix_method2(gse_id)
    
    if not matrix_file or not matrix_file.exists() or matrix_file.stat().st_size < 1000:
        log("=" * 60)
        log("所有下载方法均失败")
        log("=" * 60)
        log("")
        log("请尝试手动下载:")
        log("1. 访问: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211700")
        log("2. 下载 'Series Matrix File(s)' 或 'Supplementary Files' 中的 FPKM 文件")
        log(f"3. 放置到: {OUT_DIR}/")
        log("4. 重新运行此脚本")
        return 1
    
    # 解析文件
    try:
        df_expr, expr_file = parse_geo_matrix_correct(matrix_file)
        if df_expr.shape[1] > 0:  # 确保有基因数据
            log("=" * 60)
            log("✓ 数据处理完成")
            log("=" * 60)
            log(f"表达矩阵: {expr_file}")
            log(f"矩阵维度: {df_expr.shape[0]} 样本 x {df_expr.shape[1]} 基因")
            return 0
        else:
            log("✗ Series matrix文件不包含表达数据，请下载FPKM文件")
            return 1
    except Exception as e:
        log(f"✗ 解析失败: {e}")
        import traceback
        log(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())

