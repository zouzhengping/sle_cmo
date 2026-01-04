#!/usr/bin/env python3
"""
GSE88884 数据预处理脚本
将预处理文件转换为gene symbol表达矩阵
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import gzip
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data/raw/GSE88884"
OUT_DIR = DATA_DIR
LOG_FILE = ROOT / "logs/validation/gse88884_processing.log"
PID_FILE = ROOT / "logs/validation/gse88884_processing.pid"

def log(message):
    """记录日志"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}\n"
    print(message, flush=True)
    with open(LOG_FILE, 'a', buffering=1) as f:  # 行缓冲
        f.write(log_msg)
        f.flush()

def check_if_running():
    """检查是否已在运行"""
    if PID_FILE.exists():
        old_pid = int(PID_FILE.read_text().strip())
        try:
            os.kill(old_pid, 0)  # 检查进程是否存在
            log(f"处理进程已在运行 (PID: {old_pid})")
            return True
        except OSError:
            PID_FILE.unlink()
    return False

def save_pid():
    """保存PID"""
    PID_FILE.parent.mkdir(parents=True, exist_ok=True)
    PID_FILE.write_text(str(os.getpid()))

def load_preprocessed_data():
    """加载预处理文件"""
    log("=" * 60)
    log("加载预处理文件")
    log("=" * 60)
    
    preprocessed_file = DATA_DIR / "GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz"
    
    if not preprocessed_file.exists():
        raise FileNotFoundError(f"预处理文件不存在: {preprocessed_file}")
    
    log(f"读取文件: {preprocessed_file}")
    log(f"文件大小: {preprocessed_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    # 读取文件（第一列是transcript cluster ID）
    log("读取数据...")
    df = pd.read_csv(
        preprocessed_file,
        sep="\t",
        index_col=0,
        compression='gzip',
        low_memory=False
    )
    
    log(f"数据维度: {df.shape}")
    log(f"行数（transcript clusters）: {len(df)}")
    log(f"列数（样本）: {len(df.columns)}")
    log(f"前5个transcript cluster ID: {list(df.index[:5])}")
    
    return df

def convert_transcript_to_gene_symbol(df):
    """将transcript cluster ID转换为gene symbol"""
    log("=" * 60)
    log("转换transcript cluster ID到gene symbol")
    log("=" * 60)
    
    # 检查是否已有注释文件
    annot_file = OUT_DIR / "GPL17586_annotation.tsv"
    
    if not annot_file.exists():
        # 使用R脚本获取注释
        log("使用R脚本获取GPL17586注释信息...")
        r_script = f'''
library(GEOquery)
library(Biobase)

gpl_id <- "GPL17586"
gpl_file <- file.path("{DATA_DIR}", paste0(gpl_id, ".soft.gz"))
out_file <- "{annot_file}"

cat("下载/加载GPL17586...\\n")
if (!file.exists(gpl_file)) {{
    gpl <- getGEO(gpl_id, destdir = "{DATA_DIR}")
}} else {{
    gpl <- getGEO(filename = gpl_file)
}}

annot <- Table(gpl)
cat("注释表维度:", dim(annot), "\\n")

# 查找Gene Symbol列
gene_symbol_cols <- grep("Gene Symbol|gene symbol|Symbol|SYMBOL|gene_assignment", 
                          colnames(annot), ignore.case = TRUE, value = TRUE)
if (length(gene_symbol_cols) > 0) {{
    gene_symbol_col <- gene_symbol_cols[1]
    cat("使用列:", gene_symbol_col, "\\n")
}} else {{
    gene_symbol_col <- "ID"
    cat("使用ID列\\n")
}}

# 提取ID和Gene Symbol
annot_subset <- annot[, c("ID", gene_symbol_col)]
colnames(annot_subset) <- c("transcript_cluster_id", "gene_symbol")

# 处理gene symbol（取第一个）
annot_subset$gene_symbol <- sapply(annot_subset$gene_symbol, function(x) {{
    if (is.na(x) || x == "") return(NA)
    strsplit(as.character(x), split = "///|;|,")[[1]][1]
}})

# 过滤有效注释
annot_subset <- annot_subset[!is.na(annot_subset$gene_symbol) & 
                              annot_subset$gene_symbol != "" & 
                              annot_subset$gene_symbol != "---", ]

cat("有效注释数:", nrow(annot_subset), "\\n")

write.table(annot_subset, file = out_file, sep = "\\t", quote = FALSE, row.names = FALSE)
cat("已保存:", out_file, "\\n")
'''
        
        r_script_file = OUT_DIR / "get_gpl17586_annotation.R"
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        log(f"运行R脚本...")
        result = subprocess.run(
            ["Rscript", str(r_script_file)],
            capture_output=True,
            text=True,
            cwd=str(ROOT),
            timeout=600
        )
        
        if result.stdout:
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    log(f"R: {line}")
        if result.stderr and result.stderr.strip():
            log(f"R stderr: {result.stderr}")
        
        if result.returncode != 0:
            log("⚠ R脚本执行失败，将尝试其他方法")
    
    # 读取注释文件
    if annot_file.exists():
        log(f"读取注释文件: {annot_file}")
        df_annot = pd.read_csv(annot_file, sep="\t")
        log(f"注释表维度: {df_annot.shape}")
        log(f"注释列: {list(df_annot.columns)}")
        
        # 提取transcript cluster ID
        # 预处理文件中的ID格式是 "TC01000774.hg.1:1"，注释文件中是 "TC01000001.hg.1"
        # 需要提取到 ".hg.1" 之前的部分来匹配
        log("提取transcript cluster ID...")
        # 从 "TC01000774.hg.1:1" 提取 "TC01000774.hg.1"
        df['transcript_cluster_id'] = df.index.str.split(':').str[0]
        log(f"示例ID: {df['transcript_cluster_id'].head().tolist()}")
        
        # 合并注释
        log("合并注释信息...")
        df_merged = df.merge(
            df_annot,
            left_on='transcript_cluster_id',
            right_on='transcript_cluster_id',
            how='left'
        )
        
        # 过滤有gene symbol的行
        df_with_symbol = df_merged[df_merged['gene_symbol'].notna() & 
                                   (df_merged['gene_symbol'] != '')]
        log(f"有gene symbol的行数: {len(df_with_symbol)}/{len(df)}")
        
        if len(df_with_symbol) == 0:
            log("⚠ 没有匹配到gene symbol，使用transcript cluster ID作为行名")
            return df
        
        # 移除辅助列
        sample_cols = [col for col in df.columns if col != 'transcript_cluster_id']
        df_final = df_with_symbol.set_index('gene_symbol')[sample_cols]
        
        # 聚合重复的gene symbol（取平均值）
        log("聚合重复的gene symbol...")
        df_aggregated = df_final.groupby(df_final.index).mean()
        
        log(f"最终表达矩阵维度: {df_aggregated.shape}")
        log(f"基因数: {len(df_aggregated)}")
        log(f"样本数: {len(df_aggregated.columns)}")
        
        return df_aggregated
    else:
        log("⚠ 未找到注释文件，使用transcript cluster ID作为行名")
        # 如果无法获取注释，至少保存原始数据
        return df

def save_expression_matrix(df_expr):
    """保存表达矩阵"""
    log("=" * 60)
    log("保存表达矩阵")
    log("=" * 60)
    
    out_file = OUT_DIR / "GSE88884_expression_matrix_gene_symbol.tsv"
    
    log(f"保存到: {out_file}")
    df_expr.to_csv(out_file, sep="\t")
    
    log(f"✓ 已保存: {out_file}")
    log(f"  文件大小: {out_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    return out_file

def main():
    """主函数"""
    log("=" * 80)
    log("GSE88884 数据预处理")
    log("=" * 80)
    log(f"工作目录: {ROOT}")
    log(f"数据目录: {DATA_DIR}")
    log(f"日志文件: {LOG_FILE}")
    log(f"PID: {os.getpid()}")
    
    save_pid()
    
    try:
        # 步骤1: 加载预处理数据
        df = load_preprocessed_data()
        
        # 步骤2: 转换为gene symbol
        df_expr = convert_transcript_to_gene_symbol(df)
        
        # 步骤3: 保存表达矩阵
        out_file = save_expression_matrix(df_expr)
        
        log("\n" + "=" * 80)
        log("✓ GSE88884数据预处理完成！")
        log("=" * 80)
        log(f"输出文件: {out_file}")
        log(f"表达矩阵维度: {df_expr.shape}")
        
    except Exception as e:
        log(f"\n✗ 错误: {e}")
        import traceback
        log(traceback.format_exc())
        sys.exit(1)
    finally:
        if PID_FILE.exists():
            PID_FILE.unlink()

if __name__ == "__main__":
    main()

