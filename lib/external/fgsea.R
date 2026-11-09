#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell — fgsea 官方预排序 GSEA（CONTRACT §10.3）
#
# 工具：fgsea 1.32.4（官方 fgseaMultilevel）
# 输入：bundle + params：
#   ranking   pca1 | dispersion | custom（默认 pca1）
#   gmt_file  用户 GMT（留空 = 内置精选基因集）
#   gene_stats custom 模式：JSON 对象 {gene: score}
#   gene_list  custom 模式：逗号/空格分隔基因列表（按给定顺序排成降序排名）
#   n_perm_simple 置换次数（默认 10000）
# 输出：GSEA 表（NES/padj/leading edge）+ 富集曲线图 + 排序分布图
# 说明：pca1 = 每个基因与第一主成分（bundle 的 pca.csv 或现场 PCA）的相关性；
#       dispersion = 基因在 log-normalized 数据上的方差。
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

# 与 PC1 的相关性（稀疏矩阵友好：genes × cells）
sc_gene_pc1_cor <- function(expr_cells_genes, pc1) {
  m <- methods::as(Matrix::t(expr_cells_genes), "CsparseMatrix")  # genes × cells
  n <- ncol(m)
  pc <- as.numeric(pc1)
  pc <- pc - mean(pc)
  sd_pc <- stats::sd(pc)
  if (!is.finite(sd_pc) || sd_pc == 0) return(NULL)
  if (inherits(m, "sparseMatrix")) {
    proj <- as.numeric(m %*% pc) / n
    mean_g <- as.numeric(Matrix::rowSums(m)) / n
    mean_g2 <- as.numeric(Matrix::rowSums(m * m)) / n
    sd_g <- sqrt(pmax(mean_g2 - mean_g^2, 0))
  } else {
    mm <- as.matrix(m)
    proj <- rowMeans(mm * pc)
    mean_g <- rowMeans(mm)
    sd_g <- apply(mm, 1, stats::sd)
  }
  cor_vals <- proj / (pmax(sd_g, 1e-12) * sd_pc)
  cor_vals[!is.finite(cor_vals)] <- 0
  stats::setNames(cor_vals, rownames(m))
}

sc_gene_custom_stats <- function(ctx, params) {
  raw <- params$gene_stats
  if (is.list(raw) && length(raw) > 0 && !is.null(names(raw))) {
    values <- suppressWarnings(as.numeric(unlist(raw)))
    out <- stats::setNames(values, names(raw))
    out <- out[is.finite(out)]
    if (length(out) > 0) return(out)
  }
  gene_list <- as.character(params$gene_list %||% "")
  if (nzchar(gene_list)) {
    genes <- unlist(strsplit(gene_list, "[,;[:space:]]+"))
    genes <- genes[nzchar(genes)]
    genes <- unique(genes)
    if (length(genes) > 1) {
      # 列表顺序即排名：第 1 个基因分数最高
      return(stats::setNames(seq(length(genes), 1), genes))
    }
  }
  NULL
}

sc_run("fgsea", "fgsea 官方预排序 GSEA", "fgsea 1.32.4（官方 fgseaMultilevel）",
       function(ctx) {
         params <- ctx$params
         sc_require("fgsea")
         ranking <- tolower(as.character(params$ranking %||% "pca1"))
         if (!ranking %in% c("pca1", "dispersion", "custom")) ranking <- "pca1"
         n_perm <- as.integer(params$n_perm_simple %||% 10000)
         if (is.na(n_perm) || n_perm < 1000) n_perm <- 10000L

         b <- sc_load(ctx, max_cells = 15000, max_genes = 5000)

         # ---- 1. 基因排序分数
         stats_vec <- NULL
         ranking_label <- NULL
         if (ranking == "pca1") {
           pc1 <- NULL
           if (!is.null(b$pca) && ncol(b$pca) >= 1) {
             pc1 <- b$pca[, 1]
             ranking_label <- "PC1（bundle 的 pca.csv）"
           } else {
             # 现场 PCA：取高变基因，prcomp 第一主成分
             var <- sc_gene_variance(b$data)
             top <- order(-var, seq_along(var))[seq_len(min(2000L, length(var)))]
             sub <- b$data[, sort(top), drop = FALSE]
             dense <- as.matrix(sub)
             mu <- colMeans(dense)
             sdv <- apply(dense, 2, stats::sd)
             sdv[sdv == 0] <- 1
             pca <- stats::prcomp((dense - rep(mu, each = nrow(dense))) / rep(sdv, each = nrow(dense)),
                                  center = FALSE, scale. = FALSE, rank. = 1)
             pc1 <- pca$x[, 1]
             ranking_label <- "PC1（现场 PCA，2,000 个高变基因）"
             ctx$notes <- c(ctx$notes, "bundle 无 pca.csv，已现场计算第一主成分作为排序依据")
           }
           stats_vec <- sc_gene_pc1_cor(b$data, pc1)
           if (is.null(stats_vec)) sc_prerequisite("无法计算 PC1 相关性（主成分方差为 0）")
         } else if (ranking == "dispersion") {
           var <- sc_gene_variance(b$data)
           stats_vec <- stats::setNames(var, b$genes)
           ranking_label <- "基因方差（log-normalized 数据）"
         } else {
           stats_vec <- sc_gene_custom_stats(ctx, params)
           if (is.null(stats_vec) || length(stats_vec) < 20) {
             sc_prerequisite("custom 排序需要 params.gene_stats（JSON 对象）或 params.gene_list（基因列表）")
           }
           stats_vec <- stats_vec[is.finite(stats_vec)]
           # 只保留数据里存在的基因
           stats_vec <- stats_vec[names(stats_vec) %in% b$genes]
           if (length(stats_vec) < 20) sc_prerequisite("custom 排序基因与数据基因重叠太少")
           ranking_label <- "自定义基因分数"
         }
         stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(names(stats_vec))]
         stats_vec <- stats_vec[!duplicated(names(stats_vec))]
         stats_vec <- sort(stats_vec, decreasing = TRUE)
         if (length(stats_vec) < 50) sc_prerequisite(sprintf("可用于排序的基因太少（%d 个）", length(stats_vec)))

         # ---- 2. 基因集
         gs <- sc_gene_sets(ctx, b$genes, max_sets = 200L)
         pathways <- gs$sets
         pathways <- lapply(pathways, function(genes) intersect(genes, names(stats_vec)))
         pathways <- pathways[vapply(pathways, length, integer(1)) >= 5]
         if (length(pathways) == 0) sc_prerequisite("没有与排序基因重叠 >= 5 的基因集")

         # ---- 3. 官方 fgseaMultilevel
         set.seed(ctx$seed)
         fg <- tryCatch(
           suppressWarnings(fgsea::fgseaMultilevel(
             pathways = pathways, stats = stats_vec, minSize = 5, maxSize = 500,
             nPermSimple = n_perm, eps = 0,
             nproc = 1, BPPARAM = NULL
           )),
           error = function(e) {
             # 退化到经典 fgsea（官方实现，仍属于 fgsea 包）
             ctx$warnings <<- c(ctx$warnings, sprintf("fgseaMultilevel 失败，退化为 fgsea 经典置换：%s",
                                                      conditionMessage(e)))
             tryCatch(suppressWarnings(fgsea::fgsea(pathways = pathways, stats = stats_vec,
                                                    minSize = 5, maxSize = 500,
                                                    nperm = max(n_perm, 1000), nproc = 1)),
                      error = function(e2) sc_prerequisite(sprintf("fgsea 运行失败：%s", conditionMessage(e2))))
           }
         )
         fg <- as.data.frame(fg)
         if (nrow(fg) == 0) sc_prerequisite("fgsea 没有返回任何结果（基因集过小或排序值与基因集无重叠）")
         fg <- fg[order(fg$padj, fg$pval, -abs(fg$NES)), , drop = FALSE]
         fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(utils::head(as.character(x), 10), collapse = ","), character(1))
         fg$pathway <- as.character(fg$pathway)
         n_sig <- sum(fg$padj < 0.05, na.rm = TRUE)
         metrics <- list(n_genes = length(stats_vec), n_sets = nrow(fg),
                         n_significant = n_sig,
                         n_pos = sum(fg$NES > 0 & fg$padj < 0.05, na.rm = TRUE),
                         n_neg = sum(fg$NES < 0 & fg$padj < 0.05, na.rm = TRUE))

         # ---- 4. 富集曲线图（前 4 条）
         top <- utils::head(fg, 4)
         curve_path <- sc_png(ctx, "fgsea_curves", width = 1450, height = 1100, draw = function() {
           graphics::par(bg = "white", mfrow = c(2, 2), mar = c(4.2, 4.4, 2.8, 1.0), las = 1, cex = 0.8)
           sorted_names <- names(stats_vec)
           for (i in seq_len(nrow(top))) {
             pathway <- top$pathway[i]
             genes <- pathways[[pathway]]
             hits <- which(sorted_names %in% genes)
             es <- rep(0, length(sorted_names))
             if (length(hits) > 0) {
               weight <- abs(stats_vec[hits])
               weight <- weight / sum(weight)
               sign <- sign(stats_vec[hits])
               contributes <- rep(0, length(sorted_names))
               if (sum(sign >= 0) > 0) contributes[hits[sign >= 0]] <- weight[sign >= 0] / max(sum(sign >= 0), 1)
               if (sum(sign < 0) > 0) contributes[hits[sign < 0]] <- -weight[sign < 0] / max(sum(sign < 0), 1)
               running <- cumsum(contributes)
               es <- running
             }
             x <- seq_along(es) / length(es)
             plot(x, es, type = "l", lwd = 2, col = "#2166AC",
                  xlab = "Ranked genes (high -> low statistic)", ylab = "Running enrichment score",
                  main = sprintf("%s\nNES=%.2f  padj=%.3g", substr(pathway, 1, 46), top$NES[i], top$padj[i]),
                  cex.main = 0.72)
             abline(h = 0, col = "grey60", lty = 2)
             rug(x[hits] / length(es), col = "#B2182B", ticksize = 0.03)
           }
         })

         # ---- 5. 排序分布图
         rank_path <- sc_png(ctx, "fgsea_ranking", width = 1400, height = 700, draw = function() {
           graphics::par(bg = "white", mar = c(4.2, 4.6, 3.0, 1.2), las = 1, cex = 0.85)
           plot(seq_along(stats_vec), stats_vec, pch = 16, cex = 0.35,
                col = ifelse(stats_vec > 0, "#B2182B", "#2166AC"),
                xlab = "Rank (descending)", ylab = "Ranking statistic",
                main = sprintf("Ranking metric: %s", ranking_label))
           abline(h = 0, col = "grey50", lty = 2)
           abline(v = which(stats_vec > 0)[sum(stats_vec > 0)], col = "grey50", lty = 3)
         })

         # ---- 6. 表
         table_df <- data.frame(
           pathway = fg$pathway,
           size = as.integer(fg$size),
           ES = sc_round(fg$ES, 4),
           NES = sc_round(fg$NES, 4),
           pval = sc_round(fg$pval, 8),
           padj = sc_round(fg$padj, 8),
           leadingEdge = fg$leadingEdge,
           stringsAsFactors = FALSE, check.names = FALSE
         )
         top_path <- fg$pathway[1]
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("排序依据：%s（%d 个基因）", ranking_label, length(stats_vec)),
                            sprintf("基因集来源：%s（%d 条，每条 >= 5 个基因与排序交集）", gs$source, nrow(fg)),
                            sprintf("显著富集（padj<0.05）：%d 条（上调 %d / 下调 %d）",
                                    n_sig, metrics$n_pos, metrics$n_neg),
                            sprintf("最强富集：%s（NES=%.2f, padj=%.2g）",
                                    substr(top_path, 1, 40), fg$NES[1], fg$padj[1])
                          ),
                          metrics = metrics)
         res$title <- "fgsea 官方预排序 GSEA"
         res$method <- sprintf("fgsea %s（官方 fgseaMultilevel，nPermSimple=%d）",
                               sc_version("fgsea"), n_perm)
         res$tables <- list(sc_table(
           "gsea_results", "GSEA 富集结果", table_df,
           sc_cols(c("pathway", "Gene set", "string"), c("size", "Size", "number"),
                   c("ES", "ES", "number"), c("NES", "NES", "number"),
                   c("pval", "P value", "number"), c("padj", "P adjust", "number"),
                   c("leadingEdge", "Leading edge (top 10)", "string")),
           max_rows = 80, note = "fgseaMultilevel（BH 校正，seed 固定）"
         ))
         res$figures <- list(
           sc_figure("gsea_curves", "GSEA 富集曲线（top 4）", curve_path),
           sc_figure("ranking_distribution", "排序统计量分布", rank_path)
         )
         res$notes <- c(res$notes, gs$source,
                        "pca1 排序 = 基因表达与 PC1 的相关性；dispersion = 基因方差")
         res
       })
