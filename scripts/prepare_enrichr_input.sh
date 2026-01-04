#!/bin/bash
# 准备Enrichr输入的基因列表文件

cd "$(dirname "$0")/.."

DE_FILE="results/tables/tscm_vs_other_de.full.tsv"
OUT_UP="results/tables/tscm_upregulated_genes.txt"
OUT_DOWN="results/tables/tscm_downregulated_genes.txt"
OUT_ALL_SIG="results/tables/tscm_all_significant_genes.txt"

echo "[INFO] 从DE结果中提取基因列表..."

# 提取显著上调基因（p_adj < 0.05, logFC > 0）
awk -F'\t' 'NR>1 && $4<0.05 && $2>0 {print $1}' "$DE_FILE" > "$OUT_UP"
echo "[INFO] 显著上调基因数: $(wc -l < "$OUT_UP")"
echo "      已保存到: $OUT_UP"

# 提取显著下调基因（p_adj < 0.05, logFC < 0）
awk -F'\t' 'NR>1 && $4<0.05 && $2<0 {print $1}' "$DE_FILE" > "$OUT_DOWN"
echo "[INFO] 显著下调基因数: $(wc -l < "$OUT_DOWN")"
echo "      已保存到: $OUT_DOWN"

# 提取所有显著基因
awk -F'\t' 'NR>1 && $4<0.05 {print $1}' "$DE_FILE" > "$OUT_ALL_SIG"
echo "[INFO] 所有显著基因数: $(wc -l < "$OUT_ALL_SIG")"
echo "      已保存到: $OUT_ALL_SIG"

echo ""
echo "[INFO] 前10个上调基因:"
head -10 "$OUT_UP"

echo ""
echo "[DONE] 基因列表已准备完成！"
echo ""
echo "使用Enrichr进行富集分析:"
echo "1. 访问 https://maayanlab.cloud/Enrichr/"
echo "2. 上传文件: $OUT_UP (用于上调基因富集)"
echo "3. 选择数据库: GO, KEGG, Reactome"
echo "4. 下载结果并保存到 results/tables/"

