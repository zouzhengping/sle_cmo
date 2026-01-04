#!/usr/bin/env python3
"""
检查GSE137029和其他数据集的SLEDAI metadata

从GEO下载完整的样本metadata，检查是否有SLEDAI、疾病活动度等信息
"""

import sys
from pathlib import Path
import pandas as pd
import subprocess
import json

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data/metadata/sledai"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("检查SLEDAI metadata")
print("=" * 80)

# 使用R的GEOquery包检查metadata
r_script = f"""
library(GEOquery)

# 检查GSE137029
cat("\\n==========================================\\n")
cat("检查 GSE137029\\n")
cat("==========================================\\n\\n")

result <- tryCatch({{
    gset <- getGEO('GSE137029', GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = '{OUT_DIR}')
    gset <- gset[[1]]
    
    pdata <- pData(gset)
    cat("样本数:", nrow(pdata), "\\n")
    cat("Metadata列数:", ncol(pdata), "\\n\\n")
    
    cat("所有列名:\\n")
    print(colnames(pdata))
    
    cat("\\n\\n查找SLEDAI相关列:\\n")
    sledai_cols <- grep("SLEDAI|sledai|DAI|dai|activity|Activity|flare|Flare", colnames(pdata), ignore.case = TRUE, value = TRUE)
    if (length(sledai_cols) > 0) {{
        cat("找到SLEDAI相关列:\\n")
        print(sledai_cols)
        cat("\\n这些列的内容示例:\\n")
        for (col in sledai_cols) {{
            cat("\\n", col, ":\\n")
            print(table(pdata[[col]], useNA = "always"))
        }}
    }} else {{
        cat("未找到SLEDAI相关列\\n")
    }}
    
    # 保存完整metadata
    write.table(pdata, file = '{OUT_DIR}/GSE137029_full_metadata.tsv', sep = "\\t", quote = FALSE, row.names = TRUE)
    cat("\\n\\n已保存完整metadata到: {OUT_DIR}/GSE137029_full_metadata.tsv\\n")
    
    return(TRUE)
}}, error = function(e) {{
    cat("错误:", e$message, "\\n")
    return(FALSE)
}})

# 检查其他数据集
datasets <- c("GSE88884", "GSE49454", "GSE81622", "GSE72747")

for (gse_id in datasets) {{
    cat("\\n\\n==========================================\\n")
    cat("检查", gse_id, "\\n")
    cat("==========================================\\n\\n")
    
    result <- tryCatch({{
        gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = '{OUT_DIR}')
        gset <- gset[[1]]
        
        pdata <- pData(gset)
        cat("样本数:", nrow(pdata), "\\n")
        
        sledai_cols <- grep("SLEDAI|sledai|DAI|dai|activity|Activity|flare|Flare", colnames(pdata), ignore.case = TRUE, value = TRUE)
        if (length(sledai_cols) > 0) {{
            cat("✓ 找到SLEDAI相关列:", paste(sledai_cols, collapse = ", "), "\\n")
            write.table(pdata, file = paste0('{OUT_DIR}/', gse_id, '_full_metadata.tsv'), sep = "\\t", quote = FALSE, row.names = TRUE)
        }} else {{
            cat("✗ 未找到SLEDAI相关列\\n")
        }}
        return(TRUE)
    }}, error = function(e) {{
        cat("✗ 无法获取数据:", e$message, "\\n")
        return(FALSE)
    }})
}}
"""

# 保存R脚本
r_script_file = OUT_DIR / "check_sledai_metadata.R"
with open(r_script_file, 'w') as f:
    f.write(r_script)

print(f"\n[步骤1] 创建R检查脚本: {r_script_file}")
print(f"[步骤2] 运行R脚本检查metadata...\n")

# 运行R脚本
result = subprocess.run(
    ["Rscript", str(r_script_file)],
    capture_output=True,
    text=True,
    cwd=str(ROOT)
)

print(result.stdout)
if result.stderr:
    print("错误信息:", result.stderr)

# 检查生成的文件
print("\n" + "=" * 80)
print("检查生成的文件")
print("=" * 80)

metadata_files = list(OUT_DIR.glob("*_full_metadata.tsv"))
if metadata_files:
    print(f"\n找到 {len(metadata_files)} 个metadata文件:")
    for f in metadata_files:
        print(f"  - {f.name}")
        try:
            df = pd.read_csv(f, sep="\t", index_col=0, nrows=5)
            print(f"    列数: {len(df.columns)}")
            # 查找SLEDAI相关列
            sledai_cols = [col for col in df.columns if any(keyword in col.upper() for keyword in ['SLEDAI', 'DAI', 'ACTIVITY', 'FLARE'])]
            if sledai_cols:
                print(f"    ✓ 找到SLEDAI相关列: {', '.join(sledai_cols)}")
            else:
                print(f"    ✗ 未找到SLEDAI相关列")
        except Exception as e:
            print(f"    ✗ 读取失败: {e}")
else:
    print("\n未找到metadata文件")

print("\n" + "=" * 80)
print("检查完成")
print("=" * 80)

