#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell — GSVA 官方通路活性（CONTRACT §10.3）
#
# 工具：GSVA 2.0.7（官方 gsva(gsvaParam/ssgseaParam)）
# 输入：bundle（data.mtx = log-normalized）+ params：
#   method     gsva | ssgsea（默认 gsva）
#   gmt_file   用户 GMT（留空 = 内置精选基因集）
#   n_sets     使用前 N 条基因集（默认 12）
#   cluster_col 热图分群列（默认自动）
# 输出：每细胞通路分数 obs_extra.csv（列名 gsva_<通路>）+ 簇均值热图 + 通路汇总表
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

sc_run("gsva", "GSVA 官方通路活性", "GSVA 2.0.7（官方 gsva/ssgsea）",
       function(ctx) {
         params <- ctx$params
         sc_require("GSVA")
         method <- tolower(as.character(params$method %||% "gsva"))
         if (!method %in% c("gsva", "ssgsea")) method <- "gsva"
         n_sets <- as.integer(params$n_sets %||% 12)
         if (is.na(n_sets) || n_sets < 1) n_sets <- 12L

         b <- sc_load(ctx, max_cells = 15000, max_genes = 5000)
         gs <- sc_gene_sets(ctx, b$genes, max_sets = n_sets)
         expr <- b$data  # cells × genes（log-normalized）
         gene_sets <- lapply(gs$sets, function(genes) intersect(genes, b$genes))
         gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 5]
         if (length(gene_sets) == 0) sc_prerequisite("没有可用基因集（过滤后为空）")

         # GSVA 需要 genes × samples 的矩阵
         expr_gs <- methods::as(Matrix::t(expr), "CsparseMatrix")  # genes × cells
         param <- if (method == "gsva") {
           tryCatch(GSVA::gsvaParam(expr_gs, gene_sets, kcdf = "Gaussian", sparse = TRUE),
                    error = function(e) {
                      GSVA::gsvaParam(as.matrix(expr_gs), gene_sets, kcdf = "Gaussian", sparse = FALSE)
                    })
         } else {
           tryCatch(GSVA::ssgseaParam(expr_gs, gene_sets, normalize = TRUE),
                    error = function(e) {
                      GSVA::ssgseaParam(as.matrix(expr_gs), gene_sets, normalize = TRUE)
                    })
         }
         scores <- tryCatch(
           suppressWarnings(GSVA::gsva(param, verbose = FALSE)),
           error = function(e) sc_prerequisite(sprintf("GSVA 计算失败：%s", conditionMessage(e)))
         )
         scores <- as.matrix(scores)
         set_names <- rownames(scores)
         cell_names <- colnames(scores)
         # 对齐 bundle 细胞顺序
         idx <- match(b$cells, cell_names)
         if (anyNA(idx)) sc_prerequisite("GSVA 结果细胞名与 barcodes.tsv 不一致")
         scores <- scores[, idx, drop = FALSE]
         score_mat <- t(scores)  # cells × sets

         safe_names <- make.names(set_names, unique = TRUE)
         obs_cols <- list()
         for (i in seq_along(set_names)) obs_cols[[paste0("gsva_", safe_names[i])]] <- as.numeric(score_mat[, i])
         sc_write_obs(ctx, obs_cols, b$cells)

         # ---- 分群均值
         cluster_col <- sc_pick_column(
           b$obs,
           c(as.character(params$cluster_col %||% ""), "clusters", "cluster", "seurat_clusters", "cell_type"),
           what = "分群列", required = FALSE
         )
         summary_df <- data.frame(
           pathway = set_names,
           n_genes = as.integer(vapply(gene_sets[set_names], length, integer(1))),
           mean_score = sc_round(colMeans(score_mat), 4),
           sd_score = sc_round(apply(score_mat, 2, stats::sd), 4),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         heat_path <- NULL
         if (!is.null(cluster_col)) {
           grp <- sc_groups(b$obs, cluster_col)
           keep <- !is.na(grp$labels)
           if (sum(keep) > 10 && length(grp$levels) > 1) {
             agg <- rowsum(score_mat[keep, , drop = FALSE], group = grp$labels[keep], reorder = TRUE)
             n_per <- as.numeric(table(grp$labels[keep])[rownames(agg)])
             agg <- agg / n_per
             summary_df$top_cluster <- grp$levels[apply(agg, 2, which.max)]
             # 热图行 = 按分群间标准差排序的前 20 条通路，列 = 簇
             ord <- order(-apply(agg, 2, stats::sd), seq_along(set_names))
             top <- ord[seq_len(min(20, length(ord)))]
             mat <- t(agg[, top, drop = FALSE])
             rownames(mat) <- set_names[top]
             colnames(mat) <- rownames(agg)
             heat_path <- sc_png(ctx, "gsva_cluster_heatmap", width = 1550, height = 900, draw = function() {
               sc_heatmap(mat, title = sprintf("%s pathway activity (cluster mean)", toupper(method)),
                          value_label = "score", row_title = "Pathway", col_title = cluster_col,
                          cex = 0.7)
             })
           }
         }
         if (is.null(heat_path)) {
           # 无分群：抽样 120 个细胞画热图
           ord <- order(rowMeans(score_mat))
           take <- ord[round(seq(1, length(ord), length.out = min(120, length(ord))))]
           mat <- t(score_mat[take, , drop = FALSE])
           rownames(mat) <- set_names
           colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
           heat_path <- sc_png(ctx, "gsva_cell_heatmap", width = 1500, height = 800, draw = function() {
             sc_heatmap(mat, title = sprintf("%s pathway activity (sampled cells)", toupper(method)),
                        value_label = "score", row_title = "Pathway", col_title = "Cells (sorted)",
                        cex = 0.7)
           })
         }

         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("GSVA 方法：%s（官方 GSVA %s）", method, sc_version("GSVA")),
                            sprintf("基因集来源：%s（%d 条通路）", gs$source, length(set_names)),
                            sprintf("每细胞通路分数已写入 obs_extra.csv（%d 列）", length(obs_cols)),
                            if (!is.null(cluster_col)) sprintf("热图为 %s 各簇平均活性", cluster_col) else "热图为抽样细胞活性"
                          ),
                          metrics = list(
                            n_cells = nrow(score_mat), n_sets = length(set_names),
                            method = method,
                            score_mad = round(stats::mad(as.numeric(score_mat)), 4)
                          ))
         res$title <- "GSVA 官方通路活性"
         res$method <- sprintf("GSVA %s（官方 gsva/%s，表达矩阵 log-normalized）", sc_version("GSVA"), method)
         res$figures <- list(sc_figure("gsva_heatmap", "GSVA 通路活性热图", heat_path))
         res$tables <- list(sc_table(
           "pathway_scores", "通路活性汇总", summary_df,
           sc_cols(c("pathway", "Pathway", "string"), c("n_genes", "Genes", "number"),
                   c("mean_score", "Mean score", "number"), c("sd_score", "SD", "number"),
                   c("top_cluster", "Top cluster", "string")),
           max_rows = 60, note = sprintf("GSVA %s；obs_extra 列名 = gsva_<通路>", method)
         ))
         res$notes <- c(res$notes, gs$source)
         res
       })
