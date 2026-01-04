#!/usr/bin/env Rscript
# GSE61635 探针ID到Gene Symbol转换脚本

library(GEOquery)
library(dplyr)
library(data.table)

# 设置路径
project_root <- "/home/ug2217/projects/sle_tscm_bulk"
data_dir <- file.path(project_root, "data/raw/GSE61635")
out_dir <- data_dir

# 输出文件
matrix_file <- file.path(data_dir, "GSE61635_series_matrix.txt.gz")
out_expr_file <- file.path(out_dir, "GSE61635_expression_matrix_gene_symbol.tsv")
out_probe_map_file <- file.path(out_dir, "GSE61635_probe_to_gene_map.tsv")

cat("============================================================\n")
cat("GSE61635 探针ID到Gene Symbol转换\n")
cat("============================================================\n\n")

# 方法1: 使用GEOquery下载并获取注释
cat("[方法1] 使用GEOquery下载GSE61635...\n")
tryCatch({
    gset <- getGEO('GSE61635', GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = data_dir)
    gset <- gset[[1]]
    
    # 提取表达矩阵
    expr <- exprs(gset)
    cat(sprintf("  表达矩阵维度: %d 探针 x %d 样本\n", nrow(expr), ncol(expr)))
    
    # 获取探针注释
    fdata <- fData(gset)
    cat(sprintf("  注释信息列: %s\n", paste(colnames(fdata), collapse=", ")))
    
    # 查找Gene Symbol列
    gene_symbol_col <- NULL
    possible_cols <- c("Gene Symbol", "Gene.symbol", "SYMBOL", "Symbol", "gene_symbol")
    for (col in possible_cols) {
        if (col %in% colnames(fdata)) {
            gene_symbol_col <- col
            break
        }
    }
    
    if (is.null(gene_symbol_col)) {
        # 尝试查找包含"symbol"的列
        symbol_cols <- grep("symbol|Symbol", colnames(fdata), ignore.case = TRUE, value = TRUE)
        if (length(symbol_cols) > 0) {
            gene_symbol_col <- symbol_cols[1]
        }
    }
    
    if (is.null(gene_symbol_col)) {
        cat("  ⚠ 未找到Gene Symbol列，尝试使用ID列\n")
        gene_symbol_col <- colnames(fdata)[1]
    }
    
    cat(sprintf("  使用列作为Gene Symbol: %s\n", gene_symbol_col))
    
    # 创建探针到基因的映射
    probe_to_gene <- data.frame(
        probe_id = rownames(expr),
        gene_symbol = fdata[[gene_symbol_col]],
        stringsAsFactors = FALSE
    )
    
    # 清理Gene Symbol（去除空值、多个基因用分号分隔的情况）
    probe_to_gene$gene_symbol <- gsub(" /// ", ";", probe_to_gene$gene_symbol)
    probe_to_gene$gene_symbol <- gsub(" // ", ";", probe_to_gene$gene_symbol)
    probe_to_gene <- probe_to_gene[probe_to_gene$gene_symbol != "" & !is.na(probe_to_gene$gene_symbol), ]
    
    cat(sprintf("  有效探针数: %d/%d\n", nrow(probe_to_gene), nrow(expr)))
    
    # 对于多个基因对应一个探针的情况，取第一个基因
    probe_to_gene$gene_symbol <- sapply(strsplit(probe_to_gene$gene_symbol, ";"), function(x) x[1])
    
    # 聚合：多个探针对应一个基因时，取平均值
    expr_df <- as.data.frame(expr)
    expr_df$probe_id <- rownames(expr_df)
    
    # 合并探针和基因信息
    expr_with_gene <- merge(expr_df, probe_to_gene, by = "probe_id", all.x = TRUE)
    
    # 只保留有Gene Symbol的探针
    expr_with_gene <- expr_with_gene[!is.na(expr_with_gene$gene_symbol) & expr_with_gene$gene_symbol != "", ]
    
    # 按基因聚合（取平均值）
    sample_cols <- setdiff(colnames(expr_with_gene), c("probe_id", "gene_symbol"))
    expr_gene <- expr_with_gene %>%
        group_by(gene_symbol) %>%
        summarise(across(all_of(sample_cols), mean, na.rm = TRUE)) %>%
        as.data.frame()
    
    rownames(expr_gene) <- expr_gene$gene_symbol
    expr_gene$gene_symbol <- NULL
    
    # 转置：样本x基因
    expr_gene_t <- t(expr_gene)
    
    cat(sprintf("  转换后: %d 个样本 x %d 个基因\n", nrow(expr_gene_t), ncol(expr_gene_t)))
    
    # 保存表达矩阵
    write.table(expr_gene_t, out_expr_file, sep = "\t", quote = FALSE, col.names = NA)
    cat(sprintf("  ✓ 已保存表达矩阵: %s\n", out_expr_file))
    
    # 保存探针映射
    write.table(probe_to_gene, out_probe_map_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  ✓ 已保存探针映射: %s\n", out_probe_map_file))
    
    cat("\n============================================================\n")
    cat("✓ 转换完成\n")
    cat("============================================================\n")
    
}, error = function(e) {
    cat(sprintf("  ✗ 错误: %s\n", e$message))
    cat("  尝试方法2: 直接解析series matrix文件...\n")
    
    # 方法2: 如果GEOquery失败，尝试下载平台注释文件
    cat("  需要手动下载GPL570平台注释文件\n")
    cat("  URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570\n")
})

