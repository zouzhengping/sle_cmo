#!/bin/bash
# 快速查看GSE211700验证结果

PROJECT_ROOT="/home/ug2217/projects/sle_tscm_bulk"
REPORT="${PROJECT_ROOT}/results/validation/GSE211700_tscm_validation_report.md"
SCORES="${PROJECT_ROOT}/results/validation/GSE211700_tscm_scores.tsv"
FIG="${PROJECT_ROOT}/results/figures/GSE211700_tscm_validation.pdf"
LOG_DIR="${PROJECT_ROOT}/logs/validation"

echo "=========================================="
echo "GSE211700 TSCM Signature 验证结果查看"
echo "=========================================="
echo ""

# 检查报告文件
if [[ -f "$REPORT" ]]; then
    echo "✓ 验证报告已生成"
    echo "  位置: $REPORT"
    echo ""
    echo "--- 报告摘要 ---"
    head -50 "$REPORT"
    echo ""
    echo "... (完整报告请查看: $REPORT)"
else
    echo "✗ 验证报告尚未生成"
    echo "  可能原因:"
    echo "    1. 验证仍在进行中"
    echo "    2. 验证失败（请查看日志）"
fi

echo ""
echo "--- 其他文件 ---"

# 检查scores文件
if [[ -f "$SCORES" ]]; then
    echo "✓ Score数据: $SCORES"
    echo "  样本数: $(wc -l < "$SCORES" | xargs)"
else
    echo "✗ Score数据未生成"
fi

# 检查图片
if [[ -f "$FIG" ]]; then
    echo "✓ 可视化图片: $FIG"
    echo "  文件大小: $(du -h "$FIG" | cut -f1)"
else
    echo "✗ 可视化图片未生成"
fi

echo ""
echo "--- 最新日志 ---"
if [[ -d "$LOG_DIR" ]]; then
    LATEST_LOG=$(ls -t "$LOG_DIR"/gse211700_validation_*.log 2>/dev/null | head -1)
    if [[ -n "$LATEST_LOG" ]]; then
        echo "最新日志: $LATEST_LOG"
        echo ""
        tail -20 "$LATEST_LOG"
    else
        echo "暂无日志文件"
    fi
else
    echo "日志目录不存在"
fi

echo ""
echo "=========================================="
echo "提示: 使用以下命令查看完整报告"
echo "  cat $REPORT"
echo "  或"
echo "  less $REPORT"
echo "=========================================="

