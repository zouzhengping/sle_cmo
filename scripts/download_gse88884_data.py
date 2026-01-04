#!/usr/bin/env python3
"""
GSE88884 数据下载和处理脚本

GSE88884是Affymetrix HTA2.0 microarray数据，包含SLEDAI信息
需要：
1. 下载series matrix文件
2. 转换为gene symbol表达矩阵
3. 提取SLEDAI metadata
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
import re

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/raw/GSE88884"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = ROOT / "logs/validation/gse88884_download.log"
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
    
    out_file = OUT_DIR / f"{gse_id}_series_matrix.txt.gz"
    
    if out_file.exists() and out_file.stat().st_size > 1000:
        log(f"✓ 文件已存在: {out_file} ({out_file.stat().st_size / 1024 / 1024:.2f} MB)")
        return out_file
    
    urls = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE888nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_series_matrix.txt.gz",
    ]
    
    for url in urls:
        try:
            log(f"尝试下载: {url}")
            cmd = [
                "wget",
                "--tries=3",
                "--timeout=60",
                "--limit-rate=2m",
                "-O", str(out_file),
                url
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode == 0 and out_file.exists() and out_file.stat().st_size > 1000:
                log(f"✓ 下载成功: {out_file} ({out_file.stat().st_size / 1024 / 1024:.2f} MB)")
                return out_file
            else:
                if out_file.exists():
                    out_file.unlink()
                log(f"✗ 下载失败或文件不完整")
        except Exception as e:
            log(f"✗ 下载失败: {e}")
            continue
    
    return None

def parse_geo_matrix(matrix_file):
    """解析GEO matrix文件，提取表达数据和metadata"""
    log("=" * 60)
    log("解析GEO matrix文件")
    log("=" * 60)
    
    log(f"读取文件: {matrix_file}")
    
    # 读取文件
    with gzip.open(matrix_file, 'rt') as f:
        lines = f.readlines()
    
    # 找到数据开始行
    data_start_idx = None
    metadata = {}
    
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("!Series_"):
            # 提取metadata
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.replace("!Series_", "").strip()
                value = value.strip().strip('"')
                metadata[key] = value
        elif line.startswith("!Sample_"):
            # 提取sample metadata
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.replace("!Sample_", "").strip()
                value = value.strip().strip('"')
                if key not in metadata:
                    metadata[key] = []
                metadata[key].append(value)
        elif line.startswith("ID_REF") or (line.startswith("\"ID_REF\"") and data_start_idx is None):
            data_start_idx = i
            break
    
    if data_start_idx is None:
        raise ValueError("无法找到数据开始行")
    
    log(f"数据开始行: {data_start_idx}")
    
    # 读取表达数据
    log("读取表达数据...")
    df_expr = pd.read_csv(
        matrix_file,
        sep="\t",
        skiprows=data_start_idx,
        compression='gzip',
        index_col=0
    )
    
    # 清理列名（去除引号）
    df_expr.columns = df_expr.columns.str.strip().str.strip('"')
    
    log(f"表达矩阵维度: {df_expr.shape}")
    log(f"样本数: {len(df_expr.columns)}")
    
    return df_expr, metadata

def convert_probe_to_gene_r(df_expr, gse_id):
    """使用R的GEOquery包将probe ID转换为gene symbol"""
    log("=" * 60)
    log("使用R转换probe ID到gene symbol")
    log("=" * 60)
    
    # 保存probe表达矩阵（临时）
    temp_expr_file = OUT_DIR / f"{gse_id}_probe_expr.tsv"
    df_expr.to_csv(temp_expr_file, sep="\t")
    log(f"已保存probe表达矩阵: {temp_expr_file}")
    
    # R脚本
    r_script = f"""
library(GEOquery)
library(Biobase)
library(dplyr)
library(data.table)

gse_id <- "{gse_id}"
gse_file <- file.path("{OUT_DIR}", paste0(gse_id, "_series_matrix.txt.gz"))
temp_expr_file <- "{temp_expr_file}"
out_expr_file <- file.path("{OUT_DIR}", paste0(gse_id, "_expression_matrix_gene_symbol.tsv"))
out_meta_file <- file.path("{OUT_DIR}", paste0(gse_id, "_metadata.tsv"))

cat("\\n==========================================\\n")
cat("转换", gse_id, "probe ID到gene symbol\\n")
cat("==========================================\\n\\n")

tryCatch({{
    # 加载GEO数据
    cat("加载GEO数据...\\n")
    gset <- getGEO(filename = gse_file, GSEMatrix = TRUE, AnnotGPL = TRUE)
    gset <- gset[[1]]
    
    # 获取表达矩阵
    expr_matrix <- exprs(gset)
    cat("表达矩阵维度:", dim(expr_matrix), "\\n")
    
    # 获取feature data（probe annotation）
    feature_data <- fData(gset)
    cat("Feature data列数:", ncol(feature_data), "\\n")
    cat("Feature data列名:", paste(head(colnames(feature_data), 10), collapse = ", "), "...\\n")
    
    # 查找Gene Symbol列
    gene_symbol_cols <- grep("Gene Symbol|gene symbol|Symbol|SYMBOL", colnames(feature_data), ignore.case = TRUE, value = TRUE)
    if (length(gene_symbol_cols) == 0) {{
        # 尝试其他可能的列名
        gene_symbol_cols <- grep("Gene|gene", colnames(feature_data), ignore.case = TRUE, value = TRUE)
    }}
    
    if (length(gene_symbol_cols) > 0) {{
        cat("找到Gene Symbol列:", paste(gene_symbol_cols, collapse = ", "), "\\n")
        gene_symbol_col <- gene_symbol_cols[1]
    }} else {{
        cat("未找到Gene Symbol列，使用ID\\n")
        gene_symbol_col <- "ID"
    }}
    
    # 提取Gene Symbol
    feature_data$gene_symbol <- ifelse(
        !is.na(feature_data[[gene_symbol_col]]) & feature_data[[gene_symbol_col]] != "",
        as.character(feature_data[[gene_symbol_col]]),
        as.character(feature_data$ID)
    )
    
    # 处理多个gene symbol（用分号分隔的情况）
    feature_data$gene_symbol <- sapply(feature_data$gene_symbol, function(x) {{
        if (is.na(x) || x == "") return(NA)
        # 取第一个gene symbol
        strsplit(as.character(x), split = "///|;|,")[[1]][1]
    }})
    
    # 过滤掉没有gene symbol的probe
    valid_probes <- !is.na(feature_data$gene_symbol) & feature_data$gene_symbol != "" & feature_data$gene_symbol != "---"
    cat("有效probe数:", sum(valid_probes), "/", length(valid_probes), "\\n")
    
    expr_matrix_filtered <- expr_matrix[valid_probes, ]
    feature_data_filtered <- feature_data[valid_probes, ]
    
    # 转换为data.table进行高效聚合
    dt_expr <- as.data.table(expr_matrix_filtered, keep.rownames = "probe_id")
    dt_expr$gene_symbol <- feature_data_filtered$gene_symbol
    
    # 按gene symbol聚合（取平均值）
    cat("按gene symbol聚合...\\n")
    gene_expr_aggregated <- dt_expr %>%
        group_by(gene_symbol) %>%
        summarise(across(starts_with("GSM"), mean, na.rm = TRUE)) %>%
        ungroup()
    
    # 转换为data.frame并设置行名
    gene_expr_final <- as.data.frame(gene_expr_aggregated)
    rownames(gene_expr_final) <- gene_expr_final$gene_symbol
    gene_expr_final$gene_symbol <- NULL
    
    # 转置为samples x genes
    gene_expr_final <- t(gene_expr_final)
    
    cat("最终表达矩阵维度:", dim(gene_expr_final), "\\n")
    cat("基因数:", ncol(gene_expr_final), "\\n")
    cat("样本数:", nrow(gene_expr_final), "\\n")
    
    # 保存表达矩阵
    write.table(gene_expr_final, file = out_expr_file, sep = "\\t", quote = FALSE, row.names = TRUE)
    cat("\\n已保存gene symbol表达矩阵:", out_expr_file, "\\n")
    
    # 保存metadata
    pdata <- pData(gset)
    write.table(pdata, file = out_meta_file, sep = "\\t", quote = FALSE, row.names = TRUE)
    cat("已保存metadata:", out_meta_file, "\\n")
    
    cat("\\n✓ 转换完成\\n")
    
}}, error = function(e) {{
    cat("✗ 转换失败:", e$message, "\\n")
    traceback()
    quit(status = 1)
}})
"""
    
    # 保存R脚本
    r_script_file = OUT_DIR / "convert_probe_to_gene.R"
    with open(r_script_file, 'w') as f:
        f.write(r_script)
    
    log(f"已创建R脚本: {r_script_file}")
    
    # 运行R脚本
    log("运行R脚本...")
    result = subprocess.run(
        ["Rscript", str(r_script_file)],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        timeout=600
    )
    
    if result.stdout:
        log(result.stdout)
    if result.stderr:
        log(f"R脚本stderr: {result.stderr}")
    
    if result.returncode != 0:
        raise RuntimeError(f"R脚本执行失败: {result.returncode}")
    
    # 检查输出文件
    out_expr_file = OUT_DIR / f"{gse_id}_expression_matrix_gene_symbol.tsv"
    if out_expr_file.exists():
        log(f"✓ 已生成gene symbol表达矩阵: {out_expr_file}")
        return out_expr_file
    else:
        raise FileNotFoundError(f"未找到输出文件: {out_expr_file}")

def main():
    """主函数"""
    gse_id = "GSE88884"
    
    log("=" * 80)
    log(f"GSE88884 数据下载和处理")
    log("=" * 80)
    log(f"输出目录: {OUT_DIR}")
    log(f"日志文件: {LOG_FILE}")
    
    try:
        # 步骤1: 下载series matrix
        matrix_file = download_geo_matrix(gse_id)
        if matrix_file is None:
            raise FileNotFoundError("无法下载series matrix文件")
        
        # 步骤2: 解析matrix文件
        df_expr, metadata = parse_geo_matrix(matrix_file)
        
        # 步骤3: 转换为gene symbol（使用R）
        log("\n" + "=" * 60)
        log("开始转换probe ID到gene symbol")
        log("=" * 60)
        
        out_expr_file = convert_probe_to_gene_r(df_expr, gse_id)
        
        # 验证输出
        log("\n" + "=" * 60)
        log("验证输出文件")
        log("=" * 60)
        
        df_final = pd.read_csv(out_expr_file, sep="\t", index_col=0, nrows=5)
        log(f"✓ 最终表达矩阵维度: {df_final.shape}")
        log(f"✓ 样本数: {len(df_final)}")
        log(f"✓ 基因数: {len(df_final.columns)}")
        
        log("\n" + "=" * 80)
        log("✓ GSE88884数据下载和处理完成！")
        log("=" * 80)
        log(f"表达矩阵: {out_expr_file}")
        log(f"Metadata: {OUT_DIR / f'{gse_id}_metadata.tsv'}")
        
    except Exception as e:
        log(f"\n✗ 错误: {e}")
        import traceback
        log(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()

