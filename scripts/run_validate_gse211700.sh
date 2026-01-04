#!/bin/bash
# GSE211700验证脚本 - 后台运行，自动错误处理

set -euo pipefail

PROJECT_ROOT="/home/ug2217/projects/sle_tscm_bulk"
SCRIPT="${PROJECT_ROOT}/scripts/validate_gse211700_tscm_signature.py"
LOG_DIR="${PROJECT_ROOT}/logs/validation"
LOG_FILE="${LOG_DIR}/gse211700_validation_$(date +%Y%m%d_%H%M%S).log"
PID_FILE="${PROJECT_ROOT}/tmp/gse211700_validation.pid"
MAX_RETRIES=3

# 创建必要目录
mkdir -p "$LOG_DIR" "$(dirname "$PID_FILE")"

# 检查是否已在运行
if [[ -f "$PID_FILE" ]]; then
    EXISTING_PID=$(cat "$PID_FILE")
    if ps -p "$EXISTING_PID" > /dev/null 2>&1; then
        echo "[INFO] 验证任务已在运行，PID: $EXISTING_PID"
        echo "[INFO] 日志文件: $LOG_FILE"
        exit 0
    else
        echo "[WARN] 发现旧的PID文件但进程不存在，删除旧PID文件..."
        rm -f "$PID_FILE"
    fi
fi

# 记录PID
echo "$$" > "$PID_FILE"

# 清理函数
cleanup() {
    rm -f "$PID_FILE"
    echo "[INFO] 清理完成，PID文件已删除" >> "$LOG_FILE"
}
trap cleanup EXIT

# 日志函数
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log "=========================================="
log "GSE211700 TSCM Signature 验证开始"
log "=========================================="
log "脚本: $SCRIPT"
log "日志文件: $LOG_FILE"
log "PID: $$"
log ""

# 运行验证脚本（带重试机制）
retry_count=0
success=false

while [[ $retry_count -lt $MAX_RETRIES ]]; do
    log "尝试 $((retry_count + 1))/$MAX_RETRIES"
    
    # 运行Python脚本
    if python3 "$SCRIPT" >> "$LOG_FILE" 2>&1; then
        success=true
        log "✓ 验证成功完成"
        break
    else
        retry_count=$((retry_count + 1))
        if [[ $retry_count -lt $MAX_RETRIES ]]; then
            log "✗ 验证失败，等待10秒后重试..."
            sleep 10
        else
            log "✗ 验证失败，已达到最大重试次数"
        fi
    fi
done

log ""
log "=========================================="
if $success; then
    log "验证完成 - 成功"
    log "查看报告: ${PROJECT_ROOT}/results/validation/GSE211700_tscm_validation_report.md"
else
    log "验证完成 - 失败（已重试${MAX_RETRIES}次）"
    log "请查看日志文件: $LOG_FILE"
fi
log "=========================================="

exit $([ "$success" = true ] && echo 0 || echo 1)

