#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — clusterProfiler 官方 GO/KEGG 富集（CONTRACT §10.3）
#
# 工具：clusterProfiler 4.14.6 + org.Hs.eg.db / org.Mm.eg.db（官方实现）
# 输入：bundle（基因表达 + obs 分群），params：
#   source    go|kegg|both（默认 both）
#   species   human|mouse（默认 human）
#   gene_list 逗号/空格/换行分隔的基因符号（留空 = 用各簇 top markers）
#   top_n     每簇取 top 基因（默认 50）
#   cluster_col 分群列（默认自动：clusters/cluster/cell_type/seurat_clusters）
#   pvalue_cutoff / qvalue_cutoff（默认 0.05 / 0.2）
# 输出：result.json（GO/KEGG 富集表 + 条形图/泡图）、obs_extra 无
# 退出码：0 成功；3 缺前置（无 clusterProfiler / 无参照基因组 / 基因无法映射）
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

# ---------------------------------------------------------------------------
# 向量化 Welch t 检验：从稀疏矩阵算各簇 one-vs-rest top markers
# 返回 data.frame(cluster, gene, score, logFC, pval, padj, pct_in, pct_out)
# ---------------------------------------------------------------------------
sc_welch_markers <- function(mat, labels, levels, top_n = 50L) {
  mat <- methods::as(mat, "CsparseMatrix")
  n <- nrow(mat)
  gene_names <- colnames(mat)
  total_sum <- as.numeric(Matrix::colSums(mat))
  total_sq <- as.numeric(Matrix::colSums(mat * mat))
  total_nnz <- as.numeric(Matrix::colSums(mat != 0))
  # log 数据的 expm1 列和（用于 logFC）
  expm1_matrix <- mat
  expm1_matrix@x <- expm1(expm1_matrix@x)
  total_exp <- as.numeric(Matrix::colSums(expm1_matrix))

  out <- list()
  for (lev in levels) {
    in_group <- labels == lev & !is.na(labels)
    n_g <- sum(in_group)
    n_r <- n - n_g
    if (n_g < 2 || n_r < 2) next
    sub <- mat[in_group, , drop = FALSE]
    sum_g <- as.numeric(Matrix::colSums(sub))
    sq_g <- as.numeric(Matrix::colSums(sub * sub))
    nnz_g <- as.numeric(Matrix::colSums(sub != 0))
    sum_r <- total_sum - sum_g
    sq_r <- pmax(total_sq - sq_g, 0)
    nnz_r <- total_nnz - nnz_g

    mean_g <- sum_g / n_g
    mean_r <- sum_r / n_r
    var_g <- pmax(sq_g - sum_g^2 / n_g, 0) / max(n_g - 1, 1)
    var_r <- pmax(sq_r - sum_r^2 / n_r, 0) / max(n_r - 1, 1)
    se <- sqrt(var_g / n_g + var_r / n_r)
    t_stat <- ifelse(se > 0, (mean_g - mean_r) / pmax(se, 1e-30), 0)
    t_stat[!is.finite(t_stat)] <- 0
    df_num <- (var_g / n_g + var_r / n_r)^2
    df_den <- (var_g / n_g)^2 / max(n_g - 1, 1) + (var_r / n_r)^2 / max(n_r - 1, 1)
    df <- ifelse(df_den > 0, df_num / pmax(df_den, 1e-30), 1)
    df[!is.finite(df)] <- 1
    pval <- 2 * stats::pt(-abs(t_stat), df = pmax(df, 1))
    pval[!is.finite(pval)] <- 1
    padj <- stats::p.adjust(pval, method = "BH")

    exp_sub <- as.numeric(Matrix::colSums(expm1_matrix[in_group, , drop = FALSE]))
    exp_r <- pmax(total_exp - exp_sub, 0)
    mean_g_raw <- exp_sub / n_g
    mean_r_raw <- exp_r / n_r
    logfc <- log2((mean_g_raw + 1e-9) / (mean_r_raw + 1e-9))
    logfc[!is.finite(logfc)] <- 0

    ord <- order(-t_stat, seq_along(t_stat))
    take <- ord[seq_len(min(as.integer(top_n), length(ord)))]
    out[[length(out) + 1L]] <- data.frame(
      cluster = lev,
      gene = gene_names[take],
      score = sc_round(t_stat[take], 4),
      logFC = sc_round(logfc[take], 4),
      pval = sc_round(pval[take], 6),
      padj = sc_round(padj[take], 6),
      pct_in = sc_round(nnz_g[take] / n_g, 4),
      pct_out = sc_round(nnz_r[take] / n_r, 4),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

sc_run("clusterprofiler_ora", "clusterProfiler 官方 GO/KEGG 富集",
       "clusterProfiler 4.14.6（官方 enrichGO/enrichKEGG + org.Hs.eg.db）",
       function(ctx) {
         params <- ctx$params
         sc_require("clusterProfiler")
         species <- tolower(as.character(params$species %||% "human"))
         if (!species %in% c("human", "mouse")) species <- "human"
         orgdb_name <- if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
         sc_require(orgdb_name)
         orgdb <- get(orgdb_name, envir = asNamespace(orgdb_name))
         kegg_org <- if (species == "human") "hsa" else "mmu"
         source_kind <- tolower(as.character(params$source %||% "both"))
         if (!source_kind %in% c("go", "kegg", "both")) source_kind <- "both"
         top_n <- as.integer(params$top_n %||% 50)
         if (is.na(top_n) || top_n < 5) top_n <- 5L
         p_cut <- as.numeric(params$pvalue_cutoff %||% 0.05)
         q_cut <- as.numeric(params$qvalue_cutoff %||% 0.2)
         if (!is.finite(p_cut)) p_cut <- 0.05
         if (!is.finite(q_cut)) q_cut <- 0.2

         b <- sc_load(ctx, need_counts = FALSE, max_cells = 15000, max_genes = 5000)

         # ---- 1. 基因列表
         raw_list <- trimws(as.character(params$gene_list %||% ""))
         marker_df <- NULL
         gene_source <- NULL
         if (nzchar(raw_list)) {
           genes_in <- unique(toupper(unlist(strsplit(raw_list, "[,;[:space:]]+"))))
           genes_in <- genes_in[nzchar(genes_in)]
           gene_source <- sprintf("params.gene_list（%d 个基因）", length(genes_in))
         } else {
           cluster_col <- sc_pick_column(
             b$obs,
             c(as.character(params$cluster_col %||% ""), "clusters", "cluster", "seurat_clusters", "cell_type"),
             what = "分群列（或提供 params.gene_list）"
           )
           groups <- sc_groups(b$obs, cluster_col)
           small <- sc_drop_small_groups(ctx, groups$labels, min_cells = 3)
           if (sum(small$keep) < 20) {
             sc_prerequisite("可用细胞太少，无法计算标记基因；请提供 params.gene_list")
           }
           keep <- small$keep
           marker_df <- sc_welch_markers(b$data[keep, , drop = FALSE], small$labels[keep],
                                         small$levels, top_n = top_n)
           if (is.null(marker_df) || nrow(marker_df) == 0) {
             sc_prerequisite("没有算出任何标记基因；请提供 params.gene_list")
           }
           marker_df <- marker_df[marker_df$score > 0, , drop = FALSE]
           genes_in <- unique(marker_df$gene)
           gene_source <- sprintf("各簇 top %d 标记基因（one-vs-rest Welch t，共 %d 个唯一基因，%d 个簇）",
                                  top_n, length(genes_in), length(unique(marker_df$cluster)))
           ctx$notes <- c(ctx$notes, sprintf("标记基因由数据现算（cluster 列：%s），未使用 cache", cluster_col))
         }
         genes_in <- unique(genes_in[nzchar(genes_in)])
         if (length(genes_in) == 0) sc_prerequisite("基因列表为空")

         # ---- 2. SYMBOL -> ENTREZID
         gene_map <- tryCatch(
           suppressMessages(AnnotationDbi::select(
             orgdb, keys = genes_in, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL"
           )),
           error = function(e) NULL
         )
         gene_map <- gene_map[!is.na(gene_map$ENTREZID) & nzchar(gene_map$ENTREZID), , drop = FALSE]
         gene_map <- gene_map[!duplicated(gene_map$SYMBOL), , drop = FALSE]
         if (nrow(gene_map) == 0) {
           sc_prerequisite(sprintf("基因符号无法映射到 %s 的 Entrez ID（检查基因名是否 %s 全转录组符号）",
                                   orgdb_name, if (species == "human") "人类" else "小鼠"))
         }
         entrez <- unique(as.character(gene_map$ENTREZID))
         ctx$notes <- c(ctx$notes, sprintf("SYMBOL→ENTREZID 映射成功 %d/%d 个基因（%s）",
                                           length(entrez), length(genes_in), orgdb_name))

         res <- sc_result(ctx, metrics = structure(list(), names = character(0)))
         res$title <- "clusterProfiler 官方 GO/KEGG 富集"
         res$method <- sprintf("clusterProfiler %s（官方 enrichGO/enrichKEGG + %s）",
                               sc_version("clusterProfiler"), orgdb_name)
         bullets <- c(sprintf("输入基因：%s", gene_source),
                      sprintf("映射到 Entrez ID：%d 个", length(entrez)))
         metrics <- list(n_input_genes = length(entrez))

         # ---- 3. GO BP
         go_df <- NULL
         if (source_kind %in% c("go", "both")) {
           go_res <- tryCatch(
             suppressWarnings(clusterProfiler::enrichGO(
               gene = entrez, OrgDb = orgdb, keyType = "ENTREZID", ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = p_cut, qvalueCutoff = q_cut,
               readable = TRUE, minGSSize = 5, maxGSSize = 800
             )),
             error = function(e) {
               ctx$warnings <<- c(ctx$warnings, sprintf("enrichGO 失败：%s", conditionMessage(e)))
               NULL
             }
           )
           if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
             go_df <- as.data.frame(go_res)
             go_df$ONTOLOGY <- "BP"
             # 追加 MF/CC 前若干条（信号通路分析通常先看 BP，但官方 GO 三项都可取）
             for (ont in c("MF", "CC")) {
               extra <- tryCatch(
                 suppressWarnings(clusterProfiler::enrichGO(
                   gene = entrez, OrgDb = orgdb, keyType = "ENTREZID", ont = ont,
                   pAdjustMethod = "BH", pvalueCutoff = p_cut, qvalueCutoff = q_cut,
                   readable = TRUE, minGSSize = 5, maxGSSize = 800
                 )),
                 error = function(e) NULL
               )
               if (!is.null(extra) && nrow(as.data.frame(extra)) > 0) {
                 df <- as.data.frame(extra)
                 df$ONTOLOGY <- ont
                 go_df <- rbind(go_df, df)
               }
             }
             go_df <- go_df[order(go_df$p.adjust, go_df$pvalue), , drop = FALSE]
             go_df <- go_df[!duplicated(go_df$ID), , drop = FALSE]
           }
           if (is.null(go_df)) {
             bullets <- c(bullets, sprintf("GO 富集：无显著条目（p<%.2f, q<%.2f）", p_cut, q_cut))
             metrics$go_significant <- 0L
           } else {
             metrics$go_significant <- sum(go_df$p.adjust < p_cut, na.rm = TRUE)
             bullets <- c(bullets, sprintf("GO（BP/MF/CC）显著条目：%d 个（p.adjust < %.2f）",
                                           metrics$go_significant, p_cut))
           }
         }

         # ---- 4. KEGG（官方 enrichKEGG，需要联网下载 KEGG 数据）
         kegg_df <- NULL
         if (source_kind %in% c("kegg", "both")) {
           kegg_res <- tryCatch(
             suppressWarnings(clusterProfiler::enrichKEGG(
               gene = entrez, organism = kegg_org, keyType = "kegg",
               pAdjustMethod = "BH", pvalueCutoff = p_cut, qvalueCutoff = q_cut,
               minGSSize = 5, maxGSSize = 800
             )),
             error = function(e) {
               ctx$warnings <<- c(ctx$warnings, sprintf("enrichKEGG 失败（可能是网络/KEGG 服务不可达）：%s",
                                                        conditionMessage(e)))
               NULL
             }
           )
           if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
             kegg_df <- as.data.frame(kegg_res)
             # KEGG 结果描述形如 "hsa04110: Cell cycle"，去掉编号前缀便于阅读
             kegg_df$Description <- sub("^[a-z]{3,4}[0-9]+:\\s*", "", kegg_df$Description)
             kegg_df <- kegg_df[order(kegg_df$p.adjust, kegg_df$pvalue), , drop = FALSE]
             metrics$kegg_significant <- sum(kegg_df$p.adjust < p_cut, na.rm = TRUE)
             bullets <- c(bullets, sprintf("KEGG 显著通路：%d 条（p.adjust < %.2f）",
                                           metrics$kegg_significant, p_cut))
           } else {
             metrics$kegg_significant <- NA_integer_
             bullets <- c(bullets, "KEGG 富集不可用（KEGG 在线数据下载失败或网络受限，已跳过）")
             ctx$warnings <- c(ctx$warnings, "KEGG 富集未获得结果：clusterProfiler 需要联网访问 KEGG REST 服务")
           }
         }

         # ---- 5. 表
         if (!is.null(go_df)) {
           res$tables <- c(res$tables, list(sc_table(
             "go_enrichment", "GO 富集（BP/MF/CC）", go_df,
             sc_cols(c("ONTOLOGY", "Ontology", "string"), c("ID", "GO ID", "string"),
                     c("Description", "Term", "string"), c("GeneRatio", "Gene ratio", "string"),
                     c("BgRatio", "Background ratio", "string"), c("pvalue", "P value", "number"),
                     c("p.adjust", "P adjust", "number"), c("qvalue", "Q value", "number"),
                     c("Count", "Count", "number"), c("geneID", "Genes", "string")),
             max_rows = 50,
             note = sprintf("clusterProfiler enrichGO（BH 校正，minGSSize=5）")
           )))
         }
         if (!is.null(kegg_df)) {
           res$tables <- c(res$tables, list(sc_table(
             "kegg_enrichment", "KEGG 通路富集", kegg_df,
             sc_cols(c("ID", "KEGG ID", "string"), c("Description", "Pathway", "string"),
                     c("GeneRatio", "Gene ratio", "string"), c("BgRatio", "Background ratio", "string"),
                     c("pvalue", "P value", "number"), c("p.adjust", "P adjust", "number"),
                     c("Count", "Count", "number"), c("geneID", "Genes", "string")),
             max_rows = 50, note = "clusterProfiler enrichKEGG（BH 校正）"
           )))
         }
         if (!is.null(marker_df)) {
           res$tables <- c(res$tables, list(sc_table(
             "cluster_markers", "各簇 top 标记基因", marker_df,
             sc_cols(c("cluster", "Cluster", "string"), c("gene", "Gene", "string"),
                     c("score", "t statistic", "number"), c("logFC", "logFC (log data)", "number"),
                     c("padj", "P adjust", "number"), c("pct_in", "Pct in", "number"),
                     c("pct_out", "Pct out", "number")),
             max_rows = 300,
             note = "one-vs-rest Welch t 检验 + BH（log 标准化矩阵上的 logFC）"
           )))
         }

         # ---- 6. 图（英文白底）
         plot_df <- go_df
         plot_title <- "GO enrichment (BP/MF/CC)"
         if (source_kind == "kegg") {
           plot_df <- kegg_df
           plot_title <- "KEGG pathway enrichment"
         }
         if (is.null(plot_df) && !is.null(kegg_df)) {
           plot_df <- kegg_df
           plot_title <- "KEGG pathway enrichment"
         }
         if (is.null(plot_df) && !is.null(go_df)) {
           plot_df <- go_df
           plot_title <- "GO enrichment (BP/MF/CC)"
         }
         if (!is.null(plot_df) && nrow(plot_df) > 0) {
           top <- head(plot_df[order(plot_df$p.adjust), , drop = FALSE], 12)
           top <- top[order(top$p.adjust, decreasing = TRUE), , drop = FALSE]
           label <- substr(top$Description, 1, 58)
           bar_path <- sc_png(ctx, "clusterprofiler_ora_bar", width = 1500, height = 900, draw = function() {
             graphics::par(bg = "white", mar = c(4.5, 19, 3, 1.2), las = 1, cex = 0.85)
             vals <- -log10(pmax(top$p.adjust, 1e-300))
             cols <- sc_palette(nrow(top))
             graphics::barplot(vals, horiz = TRUE, names.arg = label, col = cols,
                               border = NA, xlab = "-log10(adjusted P)", main = plot_title,
                               cex.names = 0.78)
             graphics::abline(v = -log10(p_cut), lty = 2, col = "grey40")
             graphics::mtext(sprintf("dashed line: p.adjust = %.2f", p_cut), side = 3,
                             line = 0.2, cex = 0.7, adj = 1)
           })
           res$figures <- c(res$figures, list(sc_figure("ora_bar", "富集条形图", bar_path)))

           # 泡图（Count vs -log10 p.adjust）
           dot_path <- sc_png(ctx, "clusterprofiler_ora_dot", width = 1500, height = 900, draw = function() {
             graphics::par(bg = "white", mar = c(4.5, 19, 3, 1.2), las = 1, cex = 0.85)
             vals <- -log10(pmax(top$p.adjust, 1e-300))
             counts <- pmax(top$Count, 1)
             cols <- grDevices::colorRampPalette(c("#4E79A7", "#F28E2B", "#E15759"))(100)
             idx <- cut(vals, breaks = seq(min(vals), max(vals), length.out = 101),
                        include.lowest = TRUE, labels = FALSE)
             sizes <- 0.8 + 3.2 * (counts - min(counts)) / max(diff(range(counts)), 1)
             graphics::plot(vals, seq_len(nrow(top)), pch = 21, bg = cols[idx],
                            col = "grey30", cex = sizes, xlab = "-log10(adjusted P)",
                            ylab = "", yaxt = "n", main = paste0(plot_title, " (dot plot)"),
                            xlim = range(vals) + c(-0.2, 0.6) * diff(range(vals) + 0.5))
             graphics::axis(2, at = seq_len(nrow(top)), labels = label, las = 2, cex.axis = 0.78, tick = FALSE)
             graphics::abline(v = -log10(p_cut), lty = 2, col = "grey40")
             graphics::legend("bottomright", legend = c(min(counts), max(counts)),
                              pch = 21, pt.cex = c(0.8, 4.0), col = "grey30", bg = "grey80",
                              title = "Count", bty = "n", cex = 0.7)
           })
           res$figures <- c(res$figures, list(sc_figure("ora_dot", "富集泡图", dot_path)))
         }

         res$summary$bullets <- bullets
         res$summary$metrics <- metrics
         res$notes <- c(res$notes,
                        sprintf("富集分析使用官方 clusterProfiler（enrichGO ont=BP/MF/CC、enrichKEGG organism=%s）", kegg_org),
                        sprintf("物种=%s；显著性阈值 p<%.3f, q<%.2f（BH 校正）", species, p_cut, q_cut))
         res
       })
