#!/bin/bash
# GSE88884 后台下载和处理脚本

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs/validation"
PID_FILE="$LOG_DIR/gse88884_download.pid"
LOG_FILE="$LOG_DIR/gse88884_download_$(date +%Y%m%d_%H%M%S).log"

# 创建日志目录
mkdir -p "$LOG_DIR"

# 检查是否已经在运行
if [ -f "$PID_FILE" ]; then
    OLD_PID=$(cat "$PID_FILE")
    if ps -p "$OLD_PID" > /dev/null 2>&1; then
        echo "下载进程已在运行 (PID: $OLD_PID)"
        echo "日志文件: $(ls -t $LOG_DIR/gse88884_download_*.log 2>/dev/null | head -1)"
        exit 0
    else
        rm -f "$PID_FILE"
    fi
fi

echo "=========================================="
echo "GSE88884 后台下载和处理"
echo "=========================================="
echo "日志文件: $LOG_FILE"
echo "PID文件: $PID_FILE"
echo ""

# 后台运行下载脚本
cd "$PROJECT_ROOT"

nohup bash -c '
    echo "PID: $$" > "$PID_FILE"
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] 开始下载GSE88884" >> "$LOG_FILE"
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] PID: $$" >> "$LOG_FILE"
    
    # 检查预处理文件是否已下载
    PREPROCESSED_FILE="$PROJECT_ROOT/data/raw/GSE88884/GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz"
    
    if [ -f "$PREPROCESSED_FILE" ] && [ $(stat -c%s "$PREPROCESSED_FILE") -gt 1000000 ]; then
        echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ✓ 预处理文件已存在: $(du -h "$PREPROCESSED_FILE" | cut -f1)" >> "$LOG_FILE"
    else
        echo "[$(date +%Y-%m-%d\ %H:%M:%S)] 开始下载预处理文件..." >> "$LOG_FILE"
        wget -q --show-progress --tries=3 --timeout=60 --limit-rate=2m \
            -O "$PREPROCESSED_FILE" \
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE88nnn/GSE88884/suppl/GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz" \
            >> "$LOG_FILE" 2>&1
        
        if [ $? -eq 0 ] && [ -f "$PREPROCESSED_FILE" ]; then
            echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ✓ 下载完成: $(du -h "$PREPROCESSED_FILE" | cut -f1)" >> "$LOG_FILE"
        else
            echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ✗ 下载失败" >> "$LOG_FILE"
            rm -f "$PID_FILE"
            exit 1
        fi
    fi
    
    # 处理数据（如果需要）
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] 开始处理数据..." >> "$LOG_FILE"
    
    # 检查是否需要转换为gene symbol格式
    GENE_SYMBOL_FILE="$PROJECT_ROOT/data/raw/GSE88884/GSE88884_expression_matrix_gene_symbol.tsv"
    
    if [ ! -f "$GENE_SYMBOL_FILE" ]; then
        echo "[$(date +%Y-%m-%d\ %H:%M:%S)] 需要转换数据格式..." >> "$LOG_FILE"
        # 这里可以添加数据转换逻辑
        # 由于预处理文件可能已经是gene symbol格式，先检查一下
        zcat "$PREPROCESSED_FILE" | head -5 >> "$LOG_FILE" 2>&1
    else
        echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ✓ 表达矩阵已存在" >> "$LOG_FILE"
    fi
    
    echo "[$(date +%Y-%m-%d\ %H:%M:%S)] ✓ 处理完成" >> "$LOG_FILE"
    rm -f "$PID_FILE"
' > /dev/null 2>&1 &

NEW_PID=$!
echo "$NEW_PID" > "$PID_FILE"

echo "✓ 后台下载已启动"
echo "  PID: $NEW_PID"
echo "  日志: $LOG_FILE"
echo ""
echo "查看进度: tail -f $LOG_FILE"
echo "检查状态: bash $SCRIPT_DIR/check_gse88884_download.sh"

