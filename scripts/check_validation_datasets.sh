#!/bin/bash
# 检查SLE验证数据集的可用性

echo "=========================================="
echo "SLE TSCM Signature 验证数据集检查"
echo "=========================================="
echo ""

# 候选数据集列表
datasets=(
    "GSE61635:单细胞RNA-seq:SLE和Healthy对照"
    "GSE88884:单细胞RNA-seq:SLE患者PBMC"
    "GSE49454:Bulk RNA-seq:SLE和Healthy对照"
    "GSE50635:Bulk RNA-seq:SLE患者PBMC"
    "GSE81622:Bulk RNA-seq:SLE和Healthy对照"
)

echo "候选数据集列表:"
echo ""

for dataset in "${datasets[@]}"; do
    IFS=':' read -r gse_id tech_type description <<< "$dataset"
    echo "【$gse_id】"
    echo "  技术: $tech_type"
    echo "  描述: $description"
    echo "  GEO链接: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$gse_id"
    echo ""
done

echo "=========================================="
echo "推荐验证顺序"
echo "=========================================="
echo ""
echo "1. 第一优先级: GSE61635（单细胞，如果确认有Healthy对照）"
echo "2. 第二优先级: GSE88884（单细胞，如果确认有Healthy对照）"
echo "3. 第三优先级: GSE49454（Bulk RNA-seq，备选）"
echo ""
echo "=========================================="
echo "详细评估报告"
echo "=========================================="
echo "请查看: docs/sle_validation_datasets_comprehensive_search.md"
echo ""

