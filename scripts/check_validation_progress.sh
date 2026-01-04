#!/bin/bash
# 查看GSE211700验证进度

PROJECT_ROOT="/home/ug2217/projects/sle_tscm_bulk"
LOG_DIR="${PROJECT_ROOT}/logs/validation"
PID_FILE="${PROJECT_ROOT}/tmp/gse211700_validation.pid"
REPORT="${PROJECT_ROOT}/results/validation/GSE211700_tscm_validation_report.md"
SCORES="${PROJECT_ROOT}/results/validation/GSE211700_tscm_scores.tsv"
FIG="${PROJECT_ROOT}/results/figures/GSE211700_tscm_validation.pdf"

echo "=========================================="
echo "GSE211700 TSCM Signature 验证进度"
echo "=========================================="
echo ""

# 1. 检查进程状态
echo "【1】进程状态"
if [[ -f "$PID_FILE" ]]; then
    PID=$(cat "$PID_FILE")
    if ps -p "$PID" > /dev/null 2>&1; then
        echo "  ✓ 验证进程正在运行"
        echo "    PID: $PID"
        echo "    进程信息:"
        ps -p "$PID" -o pid,user,%cpu,%mem,etime,cmd | tail -1 | sed 's/^/      /'
        
        # 检查是否有Python进程
        PYTHON_PID=$(pgrep -f "validate_gse211700_tscm_signature.py" | head -1)
        if [[ -n "$PYTHON_PID" ]]; then
            echo "    Python进程: $PYTHON_PID"
            echo "    CPU使用率: $(ps -p $PYTHON_PID -o %cpu --no-headers | xargs)%"
            echo "    内存使用: $(ps -p $PYTHON_PID -o %mem --no-headers | xargs)%"
        fi
    else
        echo "  ⚠ PID文件存在但进程未运行（可能已完成或异常退出）"
        rm -f "$PID_FILE"
    fi
else
    echo "  ℹ 未发现运行中的验证进程"
fi

echo ""

# 2. 检查输出文件
echo "【2】输出文件状态"
files_created=0
files_total=3

if [[ -f "$REPORT" ]]; then
    echo "  ✓ 验证报告已生成"
    echo "    位置: $REPORT"
    echo "    大小: $(du -h "$REPORT" | cut -f1)"
    echo "    修改时间: $(stat -c %y "$REPORT" | cut -d'.' -f1)"
    files_created=$((files_created + 1))
else
    echo "  ⏳ 验证报告未生成"
fi

if [[ -f "$SCORES" ]]; then
    echo "  ✓ Score数据已生成"
    echo "    位置: $SCORES"
    echo "    样本数: $(wc -l < "$SCORES" | xargs)"
    files_created=$((files_created + 1))
else
    echo "  ⏳ Score数据未生成"
fi

if [[ -f "$FIG" ]]; then
    echo "  ✓ 可视化图片已生成"
    echo "    位置: $FIG"
    echo "    大小: $(du -h "$FIG" | cut -f1)"
    files_created=$((files_created + 1))
else
    echo "  ⏳ 可视化图片未生成"
fi

echo ""
echo "  完成度: $files_created/$files_total 文件已生成"

echo ""

# 3. 查看最新日志
echo "【3】最新日志"
if [[ -d "$LOG_DIR" ]]; then
    LATEST_LOG=$(ls -t "$LOG_DIR"/gse211700_validation_*.log 2>/dev/null | head -1)
    if [[ -n "$LATEST_LOG" && -f "$LATEST_LOG" ]]; then
        echo "  日志文件: $LATEST_LOG"
        echo "  文件大小: $(du -h "$LATEST_LOG" | cut -f1)"
        echo "  最后修改: $(stat -c %y "$LATEST_LOG" | cut -d'.' -f1)"
        echo ""
        echo "  --- 最后20行日志 ---"
        tail -20 "$LATEST_LOG" | sed 's/^/    /'
        
        # 检查是否有错误
        if grep -qi "error\|失败\|失败\|exception" "$LATEST_LOG"; then
            echo ""
            echo "  ⚠ 发现错误信息:"
            grep -i "error\|失败\|失败\|exception" "$LATEST_LOG" | tail -5 | sed 's/^/    /'
        fi
        
        # 检查是否完成
        if grep -qi "验证完成\|DONE\|完成" "$LATEST_LOG"; then
            echo ""
            echo "  ✓ 验证已完成"
        fi
    else
        echo "  ℹ 暂无日志文件"
    fi
else
    echo "  ℹ 日志目录不存在"
fi

echo ""

# 4. 检查数据文件
echo "【4】数据文件状态"
DATA_DIR="${PROJECT_ROOT}/data/raw/GSE211700"
if [[ -d "$DATA_DIR" ]]; then
    file_count=$(find "$DATA_DIR" -type f 2>/dev/null | wc -l)
    if [[ $file_count -gt 0 ]]; then
        echo "  ✓ 发现 $file_count 个数据文件"
        echo "    目录: $DATA_DIR"
        echo "    文件列表:"
        ls -lh "$DATA_DIR" | tail -5 | sed 's/^/      /'
    else
        echo "  ⏳ 数据目录存在但无文件"
        echo "    提示: 需要下载GSE211700数据到 $DATA_DIR"
    fi
else
    echo "  ⏳ 数据目录不存在"
    echo "    提示: 需要创建目录并下载数据"
fi

echo ""

# 5. 总体状态
echo "【5】总体状态"
if [[ -f "$REPORT" ]]; then
    # 检查报告内容判断状态
    if grep -q "数据未获取\|无法获取" "$REPORT"; then
        echo "  ⚠ 等待数据: 验证脚本已运行，但需要数据文件才能完成完整验证"
        echo "     下一步: 下载GSE211700数据到 data/raw/GSE211700/"
    elif grep -q "验证成功\|验证完成" "$REPORT"; then
        echo "  ✓ 验证完成: 已成功完成TSCM signature验证"
    else
        echo "  ℹ 验证进行中或状态未知"
    fi
else
    echo "  ⏳ 验证尚未开始或刚刚启动"
fi

echo ""
echo "=========================================="
echo "快速命令"
echo "=========================================="
echo "  查看完整报告: cat $REPORT"
echo "  查看完整日志: tail -f $LATEST_LOG"
echo "  重新运行验证: bash ${PROJECT_ROOT}/scripts/run_validate_gse211700.sh"
echo "=========================================="

