#!/bin/bash
# TSCM与SLEDAI相关性分析后台运行脚本

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs/validation"
PID_FILE="$LOG_DIR/analyze_tscm_sledai.pid"
LOG_FILE="$LOG_DIR/analyze_tscm_sledai_$(date +%Y%m%d_%H%M%S).log"

# 创建日志目录
mkdir -p "$LOG_DIR"

# 检查是否已经在运行
if [ -f "$PID_FILE" ]; then
    OLD_PID=$(cat "$PID_FILE")
    if ps -p "$OLD_PID" > /dev/null 2>&1; then
        echo "分析进程已在运行 (PID: $OLD_PID)"
        echo "日志文件: $(ls -t $LOG_DIR/analyze_tscm_sledai_*.log 2>/dev/null | head -1)"
        exit 0
    else
        rm -f "$PID_FILE"
    fi
fi

echo "=========================================="
echo "TSCM与SLEDAI相关性分析（后台运行）"
echo "=========================================="
echo "日志文件: $LOG_FILE"
echo "PID文件: $PID_FILE"
echo ""

# 后台运行分析脚本
cd "$PROJECT_ROOT"

nohup python3 "$SCRIPT_DIR/analyze_tscm_vs_sledai.py" >> "$LOG_FILE" 2>&1 &

NEW_PID=$!
echo "$NEW_PID" > "$PID_FILE"

echo "✓ 后台分析已启动"
echo "  PID: $NEW_PID"
echo "  日志: $LOG_FILE"
echo ""
echo "查看进度: tail -f $LOG_FILE"
echo "检查状态: bash $SCRIPT_DIR/check_analyze_tscm_sledai.sh"

