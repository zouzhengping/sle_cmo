#!/bin/bash
# 将所有paper相关的md文件转换为docx格式

cd /home/ug2217/projects/sle_tscm_bulk

OUT_DIR="docs/word_format"
mkdir -p "$OUT_DIR"

# 需要转换的文件列表
FILES=(
    "docs/paper_title_page_and_statements.md"
    "docs/paper_abstract.md"
    "docs/paper_introduction.md"
    "docs/paper_methods_results.md"
    "docs/paper_discussion.md"
    "docs/paper_references.md"
)

echo "开始转换Markdown文件为Word格式..."

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        basename=$(basename "$file" .md)
        output="$OUT_DIR/${basename}.docx"
        echo "转换: $file -> $output"
        pandoc "$file" -o "$output" --reference-doc=/usr/share/pandoc/data/reference.docx 2>/dev/null || \
        pandoc "$file" -o "$output" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo "  ✓ 成功"
        else
            echo "  ✗ 失败"
        fi
    else
        echo "  ⚠ 文件不存在: $file"
    fi
done

echo ""
echo "转换完成！文件保存在: $OUT_DIR"



