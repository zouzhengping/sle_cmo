#!/bin/bash
# 检查GSE88884处理状态

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs/validation"
PID_FILE="$LOG_DIR/gse88884_processing.pid"

echo "=========================================="
echo "GSE88884 处理状态检查"
echo "=========================================="
echo ""

# 检查PID文件
if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if ps -p "$PID" > /dev/null 2>&1; then
        echo "✓ 处理进程正在运行"
        echo "  PID: $PID"
        echo "  运行时间: $(ps -o etime= -p $PID | tr -d ' ')"
    else
        echo "✗ 进程已结束 (PID: $PID)"
        rm -f "$PID_FILE"
    fi
else
    echo "✗ 未找到运行中的处理进程"
fi

echo ""

# 检查文件
PREPROCESSED_FILE="$PROJECT_ROOT/data/raw/GSE88884/GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz"
EXPR_MATRIX_FILE="$PROJECT_ROOT/data/raw/GSE88884/GSE88884_expression_matrix_gene_symbol.tsv"

echo "文件状态:"
if [ -f "$PREPROCESSED_FILE" ]; then
    SIZE=$(du -h "$PREPROCESSED_FILE" | cut -f1)
    echo "  ✓ 预处理文件: $SIZE"
else
    echo "  ✗ 预处理文件: 不存在"
fi

if [ -f "$EXPR_MATRIX_FILE" ]; then
    SIZE=$(du -h "$EXPR_MATRIX_FILE" | cut -f1)
    LINES=$(wc -l < "$EXPR_MATRIX_FILE" 2>/dev/null || echo "?")
    echo "  ✓ 表达矩阵: $SIZE (约 $LINES 行)"
else
    echo "  ✗ 表达矩阵: 不存在"
fi

echo ""

# 显示最新日志
LATEST_LOG=$(ls -t "$LOG_DIR"/gse88884_processing_*.log 2>/dev/null | head -1)
if [ -n "$LATEST_LOG" ]; then
    echo "最新日志 (最后15行):"
    echo "---"
    tail -15 "$LATEST_LOG"
    echo ""
    echo "查看完整日志: tail -f $LATEST_LOG"
fi

