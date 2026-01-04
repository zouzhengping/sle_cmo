#!/bin/bash
# GitHub 上传脚本
# 使用方法: bash upload_commands.sh YOUR_USERNAME REPO_NAME

if [ $# -lt 2 ]; then
    echo "使用方法: bash upload_commands.sh YOUR_USERNAME REPO_NAME"
    echo "例如: bash upload_commands.sh Raghu150999 sle_tscm_analysis"
    exit 1
fi

USERNAME=$1
REPO_NAME=$2

echo "=== 初始化 Git 仓库 ==="
git init

echo "=== 添加所有文件 ==="
git add .

echo "=== 创建提交 ==="
git commit -m "Initial commit: TSCM signature analysis code for SLE manuscript"

echo "=== 设置主分支 ==="
git branch -M main

echo "=== 添加远程仓库 ==="
git remote add origin https://github.com/${USERNAME}/${REPO_NAME}.git

echo "=== 推送到 GitHub ==="
echo "注意: 如果要求输入密码，请使用 GitHub Personal Access Token"
git push -u origin main

echo "=== 完成！==="
echo "你的代码已上传到: https://github.com/${USERNAME}/${REPO_NAME}"
