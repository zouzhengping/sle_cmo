#!/usr/bin/env Rscript
# GSE88884 数据下载和处理（使用R GEOquery）

library(GEOquery)
library(Biobase)
library(dplyr)
library(data.table)

# 设置参数
gse_id <- "GSE88884"
out_dir <- file.path(getwd(), "data/raw/GSE88884")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n==========================================\n")
cat("下载和处理 GSE88884\n")
cat("==========================================\n\n")

cat("输出目录:", out_dir, "\n\n")

# 下载GEO数据
cat("[步骤1] 下载GEO数据...\n")
tryCatch({
    gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = out_dir)
    gset <- gset[[1]]
    
    cat("✓ 下载成功\n")
    cat("样本数:", nrow(pData(gset)), "\n")
    cat("Probe数:", nrow(exprs(gset)), "\n")
    
    # 保存原始series matrix（如果还没有）
    matrix_file <- file.path(out_dir, paste0(gse_id, "_series_matrix.txt.gz"))
    if (!file.exists(matrix_file)) {
        # GEOquery已经下载了，但可能在不同的位置
        cat("注意: series matrix可能已下载到GEOquery缓存目录\n")
    }
    
    # 获取表达矩阵
    expr_matrix <- exprs(gset)
    cat("表达矩阵维度:", dim(expr_matrix), "\n")
    
    # 获取feature data
    feature_data <- fData(gset)
    cat("Feature data列数:", ncol(feature_data), "\n")
    cat("Feature data列名:", paste(head(colnames(feature_data), 10), collapse = ", "), "...\n")
    
    # 查找Gene Symbol列
    gene_symbol_cols <- grep("Gene Symbol|gene symbol|Symbol|SYMBOL", colnames(feature_data), ignore.case = TRUE, value = TRUE)
    if (length(gene_symbol_cols) == 0) {
        # 尝试其他可能的列名
        gene_symbol_cols <- grep("Gene|gene", colnames(feature_data), ignore.case = TRUE, value = TRUE)
    }
    
    if (length(gene_symbol_cols) > 0) {
        cat("找到Gene Symbol列:", paste(gene_symbol_cols, collapse = ", "), "\n")
        gene_symbol_col <- gene_symbol_cols[1]
    } else {
        cat("未找到Gene Symbol列，使用ID\n")
        gene_symbol_col <- "ID"
    }
    
    # 提取Gene Symbol
    cat("\n[步骤2] 提取Gene Symbol...\n")
    feature_data$gene_symbol <- ifelse(
        !is.na(feature_data[[gene_symbol_col]]) & feature_data[[gene_symbol_col]] != "",
        as.character(feature_data[[gene_symbol_col]]),
        as.character(feature_data$ID)
    )
    
    # 处理多个gene symbol（用分号分隔的情况）
    feature_data$gene_symbol <- sapply(feature_data$gene_symbol, function(x) {
        if (is.na(x) || x == "") return(NA)
        # 取第一个gene symbol
        strsplit(as.character(x), split = "///|;|,")[[1]][1]
    })
    
    # 过滤掉没有gene symbol的probe
    valid_probes <- !is.na(feature_data$gene_symbol) & 
                    feature_data$gene_symbol != "" & 
                    feature_data$gene_symbol != "---" &
                    !grepl("^\\d+$", feature_data$gene_symbol)  # 排除纯数字
    
    cat("有效probe数:", sum(valid_probes), "/", length(valid_probes), "\n")
    
    expr_matrix_filtered <- expr_matrix[valid_probes, ]
    feature_data_filtered <- feature_data[valid_probes, ]
    
    # 转换为data.table进行高效聚合
    cat("\n[步骤3] 按gene symbol聚合...\n")
    dt_expr <- as.data.table(expr_matrix_filtered, keep.rownames = "probe_id")
    dt_expr$gene_symbol <- feature_data_filtered$gene_symbol
    
    # 按gene symbol聚合（取平均值）
    gene_expr_aggregated <- dt_expr %>%
        group_by(gene_symbol) %>%
        summarise(across(starts_with("GSM"), mean, na.rm = TRUE)) %>%
        ungroup()
    
    # 转换为data.frame并设置行名
    gene_expr_final <- as.data.frame(gene_expr_aggregated)
    rownames(gene_expr_final) <- gene_expr_final$gene_symbol
    gene_expr_final$gene_symbol <- NULL
    
    # 转置为samples x genes
    gene_expr_final <- t(gene_expr_final)
    
    cat("最终表达矩阵维度:", dim(gene_expr_final), "\n")
    cat("基因数:", ncol(gene_expr_final), "\n")
    cat("样本数:", nrow(gene_expr_final), "\n")
    
    # 保存表达矩阵
    out_expr_file <- file.path(out_dir, paste0(gse_id, "_expression_matrix_gene_symbol.tsv"))
    write.table(gene_expr_final, file = out_expr_file, sep = "\t", quote = FALSE, row.names = TRUE)
    cat("\n✓ 已保存gene symbol表达矩阵:", out_expr_file, "\n")
    
    # 保存metadata
    pdata <- pData(gset)
    out_meta_file <- file.path(out_dir, paste0(gse_id, "_metadata.tsv"))
    write.table(pdata, file = out_meta_file, sep = "\t", quote = FALSE, row.names = TRUE)
    cat("✓ 已保存metadata:", out_meta_file, "\n")
    
    # 检查SLEDAI列
    sledai_cols <- grep("SLEDAI|sledai", colnames(pdata), ignore.case = TRUE, value = TRUE)
    if (length(sledai_cols) > 0) {
        cat("\n✓ 找到SLEDAI列:", paste(sledai_cols, collapse = ", "), "\n")
        # 显示SLEDAI统计
        for (col in sledai_cols) {
            sledai_vals <- pdata[[col]]
            sledai_numeric <- as.numeric(gsub("[^0-9.]", "", as.character(sledai_vals)))
            sledai_numeric <- sledai_numeric[!is.na(sledai_numeric)]
            if (length(sledai_numeric) > 0) {
                cat("  ", col, ": 范围", min(sledai_numeric), "-", max(sledai_numeric), 
                    ", 均值", round(mean(sledai_numeric), 2), "\n")
            }
        }
    }
    
    cat("\n==========================================\n")
    cat("✓ GSE88884数据下载和处理完成！\n")
    cat("==========================================\n")
    
}, error = function(e) {
    cat("\n✗ 错误:", e$message, "\n")
    traceback()
    quit(status = 1)
})

