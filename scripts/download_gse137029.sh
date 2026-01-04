#!/bin/bash
set -euo pipefail

########################################
# GSE137029 sle.mtx 稳健下载脚本（修正版）
########################################

# 配置参数
GSE="GSE137029"
# ✅ 改成 HTTPS，避免 FTP 被屏蔽
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/${GSE}/suppl"
FILE="${GSE}_sle.mtx.gz"

PROJECT_ROOT="/home/ug2217/projects/sle_tscm_bulk"
OUT_DIR="${PROJECT_ROOT}/data/raw/GSE137029/sle"
LOG_DIR="${PROJECT_ROOT}/logs/download"
MANIFEST="${PROJECT_ROOT}/data/metadata/download_manifest.tsv"
PID_FILE="${PROJECT_ROOT}/tmp/${GSE}_mtx_download.pid"
NOTIFY_FILE="${PROJECT_ROOT}/tmp/download_complete_notification.txt"

# 资源限制配置
MAX_RETRIES=1000      # 最大重试次数
LIMIT_RATE="3m"       # 下载速度限制（3MB/s）
TIMEOUT=1200          # 连接/整体超时（秒）
READ_TIMEOUT=300      # 读取超时（秒）

# 全局重试计数（⚠️ 不要在函数里再 local）
retry_count=0

# 创建必要目录
mkdir -p "$OUT_DIR" "$LOG_DIR" "$(dirname "$PID_FILE")" "$(dirname "$NOTIFY_FILE")"

# 下载URL和文件路径
URL="${BASE_URL}/${FILE}"
OUT_FILE="${OUT_DIR}/${FILE}"
LOG_FILE="${LOG_DIR}/${FILE}_robust.log"
TEMP_FILE="${OUT_FILE}.tmp"

# 检查是否已经在运行
if [[ -f "$PID_FILE" ]]; then
    EXISTING_PID=$(cat "$PID_FILE")
    if ps -p "$EXISTING_PID" > /dev/null 2>&1; then
        echo "[ERROR] 下载任务已经在运行，PID: $EXISTING_PID"
        echo "[INFO] 如需重启，请先停止该进程或删除 PID 文件: $PID_FILE"
        exit 1
    else
        echo "[WARN] 发现旧的 PID 文件但进程不存在，删除旧 PID 文件..."
        rm -f "$PID_FILE"
    fi
fi

# 记录当前 PID
echo "$$" > "$PID_FILE"

# 退出清理
cleanup() {
    rm -f "$PID_FILE"
    echo "[INFO] 清理完成，PID 文件已删除"
}
trap cleanup EXIT

# 日志初始化
: > "$LOG_FILE"

log() {
    local level="$1"
    local message="$2"
    local ts
    ts="$(date +'%Y-%m-%d %H:%M:%S')"
    echo "[$ts] [$level] $message" | tee -a "$LOG_FILE"
}

log "INFO" "启动稳健的 MTX 文件下载: $FILE"
log "INFO" "下载 URL: $URL"
log "INFO" "输出文件: $OUT_FILE"
log "INFO" "资源限制: 下载速度=$LIMIT_RATE, 最大重试=$MAX_RETRIES"
log "INFO" "PID: $$, 日志文件: $LOG_FILE"

# 若已有完整文件，直接退出
if [[ -f "$OUT_FILE" ]]; then
    if gzip -t "$OUT_FILE" 2>/dev/null; then
        log "INFO" "检测到已存在且完整的文件，无需重新下载。"
        exit 0
    else
        log "WARN" "检测到损坏文件，删除后重新下载: $OUT_FILE"
        rm -f "$OUT_FILE"
    fi
fi

# 稳健下载函数
robust_download() {
    local success=false

    while [[ $retry_count -lt $MAX_RETRIES && "$success" == false ]]; do
        retry_count=$((retry_count + 1))
        log "INFO" "下载尝试 #$retry_count / $MAX_RETRIES"

        wget \
          -t 3 \
          -T "$TIMEOUT" \
          --read-timeout="$READ_TIMEOUT" \
          --limit-rate="$LIMIT_RATE" \
          --no-dns-cache \
          --max-redirect=5 \
          --retry-connrefused \
          --continue \
          --user-agent="Robust MTX Downloader/1.0" \
          -O "$TEMP_FILE" "$URL" \
          >> "$LOG_FILE" 2>&1

        if [[ $? -eq 0 ]]; then
            log "INFO" "下载尝试 #$retry_count 成功"
            success=true
        else
            log "WARN" "下载尝试 #$retry_count 失败，30 秒后重试..."
            sleep 30
        fi
    done

    $success && return 0 || return 1
}

# 执行下载
if ! robust_download; then
    log "ERROR" "超过最大重试次数 ($MAX_RETRIES)，下载仍未成功，退出。"
    exit 1
fi

# 验证文件完整性
log "INFO" "开始验证文件完整性: $TEMP_FILE"
if ! gzip -t "$TEMP_FILE"; then
    log "ERROR" "gzip 完整性验证失败: $TEMP_FILE"
    exit 1
fi

# 重命名临时文件
mv "$TEMP_FILE" "$OUT_FILE"
log "INFO" "文件重命名完成: $OUT_FILE"

# 计算文件信息
SHA256=$(sha256sum "$OUT_FILE" | awk '{print $1}')
SIZE=$(stat -c %s "$OUT_FILE")
HUMAN_SIZE=$(du -h "$OUT_FILE" | awk '{print $1}')

log "INFO" "文件信息:"
log "INFO" "  SHA256: $SHA256"
log "INFO" "  大小: $SIZE bytes ($HUMAN_SIZE)"

# 初始化 manifest（如不存在）
if [[ ! -f "$MANIFEST" ]]; then
    echo -e "file\turl\tsha256\tsize_bytes\ttimestamp\tstatus" > "$MANIFEST"
fi

# 更新 manifest
log "INFO" "更新下载清单: $MANIFEST"
grep -v -P "^${FILE}\t" "$MANIFEST" > "${MANIFEST}.tmp" || true
echo -e "${FILE}\t${URL}\t${SHA256}\t${SIZE}\t$(date +'%Y-%m-%d %H:%M:%S')\tCOMPLETED" >> "${MANIFEST}.tmp"
mv "${MANIFEST}.tmp" "$MANIFEST"

# 写入通知文件
cat > "$NOTIFY_FILE" << EOF
========================================
MTX 文件下载完成通知
========================================
文件: $FILE
状态: 下载成功，文件完整
时间: $(date +'%Y-%m-%d %H:%M:%S')
路径: $OUT_FILE
大小: $HUMAN_SIZE ($SIZE bytes)
SHA256: $SHA256
下载尝试次数: $retry_count
========================================
EOF

log "INFO" "========================================"
log "INFO" "MTX 文件下载完成！"
log "INFO" "文件: $FILE"
log "INFO" "状态: 下载成功，文件完整"
log "INFO" "时间: $(date +'%Y-%m-%d %H:%M:%S')"
log "INFO" "路径: $OUT_FILE"
log "INFO" "大小: $HUMAN_SIZE ($SIZE bytes)"
log "INFO" "SHA256: $SHA256"
log "INFO" "下载尝试次数: $retry_count"
log "INFO" "通知文件: $NOTIFY_FILE"
log "INFO" "========================================"

echo ""
echo "========================================"
echo "MTX 文件下载完成！"
echo "文件: $FILE"
echo "状态: 下载成功，文件完整"
echo "时间: $(date +'%Y-%m-%d %H:%M:%S')"
echo "路径: $OUT_FILE"
echo "大小: $HUMAN_SIZE ($SIZE bytes)"
echo "SHA256: $SHA256"
echo "下载尝试次数: $retry_count"
echo "通知文件: $NOTIFY_FILE"
echo "========================================"
echo ""

