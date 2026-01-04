#!/bin/bash
# 检查TSCM与SLEDAI分析状态

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs/validation"
PID_FILE="$LOG_DIR/analyze_tscm_sledai.pid"

echo "=========================================="
echo "TSCM与SLEDAI分析状态检查"
echo "=========================================="
echo ""

# 检查PID文件
if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if ps -p "$PID" > /dev/null 2>&1; then
        echo "✓ 分析进程正在运行"
        echo "  PID: $PID"
        echo "  运行时间: $(ps -o etime= -p $PID | tr -d ' ')"
    else
        echo "✗ 进程已结束 (PID: $PID)"
        rm -f "$PID_FILE"
    fi
else
    echo "✗ 未找到运行中的分析进程"
fi

echo ""

# 检查输出文件
OUT_DIR="$PROJECT_ROOT/results/validation/sledai"
FIG_DIR="$PROJECT_ROOT/results/figures"

echo "输出文件状态:"
if [ -d "$OUT_DIR" ]; then
    RESULT_FILES=$(find "$OUT_DIR" -name "*.tsv" -o -name "*.txt" 2>/dev/null | wc -l)
    echo "  ✓ 结果目录: $OUT_DIR ($RESULT_FILES 个文件)"
    ls -lh "$OUT_DIR"/*.tsv 2>/dev/null | head -5
else
    echo "  ✗ 结果目录: 不存在"
fi

if [ -d "$FIG_DIR" ]; then
    FIG_FILES=$(find "$FIG_DIR" -name "*sledai*" -o -name "*SLEDAI*" 2>/dev/null | wc -l)
    echo "  ✓ 图表目录: $FIG_DIR ($FIG_FILES 个相关文件)"
    ls -lh "$FIG_DIR"/*sledai* "$FIG_DIR"/*SLEDAI* 2>/dev/null | head -5
fi

echo ""

# 显示最新日志
LATEST_LOG=$(ls -t "$LOG_DIR"/analyze_tscm_sledai_*.log 2>/dev/null | head -1)
if [ -n "$LATEST_LOG" ]; then
    echo "最新日志 (最后20行):"
    echo "---"
    tail -20 "$LATEST_LOG"
    echo ""
    echo "查看完整日志: tail -f $LATEST_LOG"
fi

