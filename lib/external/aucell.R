#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — AUCell 官方基因集活性打分（CONTRACT §10.3）
#
# 工具：AUCell 1.28.0（官方 AUCell_buildRankings + AUCell_calcAUC）
# 输入：bundle（data.mtx 为 log-normalized）+ params：
#   gmt_file  用户 GMT（留空 = 内置精选基因集）
#   n_sets    使用前 N 条基因集（默认 12）
#   cluster_col 热图分群列（默认自动 clusters/cell_type）
# 输出：每细胞 AUC 写 obs_extra.csv（列名 auc_<基因集>）+ 簇均值热图 + 基因集汇总表
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

sc_run("aucell", "AUCell 官方基因集打分", "AUCell 1.28.0（官方 AUCell_buildRankings + AUCell_calcAUC）",
       function(ctx) {
         params <- ctx$params
         sc_require("AUCell")
         n_sets <- as.integer(params$n_sets %||% 12)
         if (is.na(n_sets) || n_sets < 1) n_sets <- 12L

         b <- sc_load(ctx, max_cells = 15000, max_genes = 5000)
         gs <- sc_gene_sets(ctx, b$genes, max_sets = n_sets)
         set_names <- names(gs$sets)
         if (length(set_names) > 40) {
           gs$sets <- gs$sets[seq_len(40)]
           set_names <- names(gs$sets)
         }

         # ---- AUCell 官方流程（AUCell 要求 genes × cells）
         expr <- b$data
         if (inherits(expr, "sparseMatrix")) {
           expr_mat <- methods::as(Matrix::t(expr), "CsparseMatrix")
         } else {
           expr_mat <- t(as.matrix(expr))
         }
         rownames(expr_mat) <- b$genes
         colnames(expr_mat) <- b$cells
         rankings <- tryCatch(
           suppressWarnings(AUCell::AUCell_buildRankings(expr_mat, plotStats = FALSE,
                                                        verbose = FALSE, nCores = 1)),
           error = function(e) sc_prerequisite(sprintf("AUCell_buildRankings 失败：%s", conditionMessage(e)))
         )
         auc_max_rank <- max(1, min(ceiling(0.05 * nrow(rankings)), ncol(rankings)))
         auc <- tryCatch(
           suppressWarnings(AUCell::AUCell_calcAUC(gs$sets, rankings, nCores = 1,
                                                   aucMaxRank = auc_max_rank, verbose = FALSE)),
           error = function(e) sc_prerequisite(sprintf("AUCell_calcAUC 失败：%s", conditionMessage(e)))
         )
         auc_mat <- t(as.matrix(AUCell::getAUC(auc)))  # cells × sets
         colnames(auc_mat) <- set_names

         # ---- obs_extra
         safe_names <- make.names(set_names, unique = TRUE)
         obs_cols <- list()
         for (i in seq_along(set_names)) {
           obs_cols[[paste0("auc_", safe_names[i])]] <- as.numeric(auc_mat[, i])
         }
         sc_write_obs(ctx, obs_cols, b$cells)

         # ---- 表：基因集汇总
         summary_df <- data.frame(
           gene_set = set_names,
           n_genes = as.integer(vapply(gs$sets, length, integer(1))),
           mean_AUC = sc_round(colMeans(auc_mat), 4),
           max_AUC = sc_round(apply(auc_mat, 2, max), 4),
           stringsAsFactors = FALSE,
           check.names = FALSE
         )
         cluster_col <- sc_pick_column(
           b$obs,
           c(as.character(params$cluster_col %||% ""), "clusters", "cluster", "seurat_clusters", "cell_type"),
           what = "分群列", required = FALSE
         )
         cluster_levels <- NULL
         cluster_labels <- NULL
         if (!is.null(cluster_col)) {
           grp <- sc_groups(b$obs, cluster_col)
           cluster_labels <- grp$labels
           cluster_levels <- grp$levels
           keep <- !is.na(cluster_labels)
           if (sum(keep) > 10 && length(cluster_levels) > 1) {
             agg <- rowsum(auc_mat[keep, , drop = FALSE], group = cluster_labels[keep], reorder = TRUE)
             n_per <- as.numeric(table(cluster_labels[keep])[rownames(agg)])
             agg <- agg / n_per
             summary_df$top_cluster <- cluster_levels[apply(agg, 2, which.max)]
             # 热图：按平均 AUC 取前 12 条基因集
             ord <- order(-colMeans(agg), seq_along(set_names))
             top <- ord[seq_len(min(12, length(ord)))]
             mat <- t(agg[, top, drop = FALSE])
             rownames(mat) <- set_names[top]
             colnames(mat) <- rownames(agg)
             rng <- range(mat, na.rm = TRUE)
             if (diff(rng) == 0) mat <- mat + 1e-9
             heat_path <- sc_png(ctx, "aucell_heatmap", width = 1500, height = 850, draw = function() {
               sc_heatmap(mat, title = "AUCell gene-set activity (cluster mean)",
                          value_label = "AUC", row_title = "Gene set",
                          col_title = cluster_col, cex = 0.72)
             })
             ctx$figures_heat <- heat_path
           }
         }

         # ---- 图：每细胞 AUC 分布（top 6 基因集）
         top_sets <- order(-colMeans(auc_mat))[seq_len(min(6, ncol(auc_mat)))]
         dist_path <- sc_png(ctx, "aucell_distribution", width = 1400, height = 800, draw = function() {
           graphics::par(bg = "white", mar = c(4.2, 7.5, 3, 1.2), las = 1, cex = 0.85)
           vals <- lapply(top_sets, function(i) auc_mat[, i])
           names(vals) <- set_names[top_sets]
           graphics::boxplot(vals, horizontal = TRUE, col = sc_palette(length(vals)),
                             border = "grey30", xlab = "AUCell AUC per cell",
                             main = "AUCell activity distribution (top gene sets)",
                             cex.axis = 0.75)
         })

         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("基因集来源：%s", gs$source),
                            sprintf("参与打分基因集：%d 条（每条与数据重叠 >= 5 个基因）", length(set_names)),
                            sprintf("每细胞 AUC 已写入 obs_extra.csv（%d 个 obs 列，可嵌入图着色）", length(obs_cols)),
                            if (!is.null(cluster_col)) sprintf("簇均值热图：以 %s 分群", cluster_col) else "未找到分群列，热图仅给出分布图"
                          ),
                          metrics = list(
                            n_cells = nrow(auc_mat), n_sets = length(set_names),
                            auc_max_rank = auc_max_rank,
                            mean_auc_top_set = round(max(colMeans(auc_mat)), 4)
                          ))
         res$title <- "AUCell 官方基因集打分"
         res$method <- sprintf("AUCell %s（官方 AUCell_buildRankings + AUCell_calcAUC）", sc_version("AUCell"))
         res$figures <- list(sc_figure("aucell_distribution", "AUCell 活性分布", dist_path))
         if (!is.null(ctx$figures_heat)) {
           res$figures <- c(res$figures, list(sc_figure("aucell_heatmap", "AUCell 基因集活性热图", ctx$figures_heat)))
         }
         res$tables <- list(sc_table(
           "geneset_auc", "基因集 AUC 汇总", summary_df,
           sc_cols(c("gene_set", "Gene set", "string"), c("n_genes", "Genes", "number"),
                   c("mean_AUC", "Mean AUC", "number"), c("max_AUC", "Max AUC", "number"),
                   c("top_cluster", "Top cluster", "string")),
           max_rows = 60, note = sprintf("AUCell（aucMaxRank=%d，normAUC）", auc_max_rank)
         ))
         res$notes <- c(res$notes, gs$source,
                        "obs_extra 列名格式：auc_<基因集名（英文句点替换）>")
         res
       })
