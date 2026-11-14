#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — monocle3 官方轨迹（CONTRACT §10.3）
#
# 工具：monocle3 1.3.7（官方 new_cell_data_set → preprocess_cds → reduce_dimension(UMAP)
#       → cluster_cells → learn_graph → order_cells）
# 输入：bundle（counts.mtx 优先，缺失用 data.mtx）+ params：
#   resolution   cluster_cells 分辨率（默认 0.8）
#   n_top_genes  preprocess_cds 高变基因数（默认 2000）
#   root_cell    可选；根细胞名或根节点名（默认自动：最大连通分量的端点）
#   max_cells    兜底上限（默认 5000，轨迹工具收紧）
# 输出：obs_extra.csv（monocle3_pseudotime / monocle3_cluster / monocle3_partition）
#       + 伪时间轨迹图 + 簇-伪时间表
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

# 自动选根：最大连通分量里「图距离总和最大」的端点（叶节点）
sc_monocle_pick_root <- function(principal_graph) {
  if (!requireNamespace("igraph", quietly = TRUE)) return(NULL)
  graph <- principal_graph[["UMAP"]]
  if (is.null(graph)) {
    graph <- principal_graph[[1L]]
  }
  if (is.null(graph)) return(NULL)
  # 转成 igraph（可能是 igraph 或 list 邻接）
  g <- tryCatch({
    if (inherits(graph, "igraph")) graph else igraph::graph_from_adj_list(graph)
  }, error = function(e) NULL)
  if (is.null(g)) return(NULL)
  comps <- igraph::components(g)
  keep <- which(comps$membership == which.max(comps$csize))
  sub <- igraph::induced_subgraph(g, keep)
  deg <- igraph::degree(sub)
  endpoints <- names(deg)[deg == 1]
  if (length(endpoints) == 0) endpoints <- names(deg)[order(deg)][seq_len(min(4L, length(deg)))]
  if (length(endpoints) == 1) return(endpoints)
  dist <- igraph::distances(sub, v = endpoints, to = igraph::V(sub))
  dist[!is.finite(dist)] <- 0
  # 每个端点到其它节点距离之和（无向图里叶端点在树的远端）
  totals <- rowSums(dist)
  endpoints[which.max(totals)]
}

sc_run("monocle3", "monocle3 官方轨迹", "monocle3 1.3.7（官方 new_cell_data_set/preprocess_cds/reduce_dimension/cluster_cells/learn_graph/order_cells）",
       function(ctx) {
         params <- ctx$params
         sc_require(c("monocle3", "SingleCellExperiment", "igraph"))
         resolution <- as.numeric(params$resolution %||% 0.8)
         if (!is.finite(resolution) || resolution <= 0) resolution <- 0.8
         # monocle3 的 cluster_cells resolution 量纲与 Seurat 不同（默认 1e-3）：
         # 对外参数沿用 Seurat 风格（默认 0.8），内部按 1e-3 缩放；params.monocle_resolution 可直传。
         monocle_res <- as.numeric(params$monocle_resolution %||% (resolution * 1e-3))
         if (!is.finite(monocle_res) || monocle_res <= 0) monocle_res <- 0.0008
         ctx$notes <- c(ctx$notes,
                        sprintf("cluster_cells resolution 已按 monocle3 量纲换算：用户 %.3g -> %.6g", resolution, monocle_res))
         n_top <- as.integer(params$n_top_genes %||% 2000)
         if (is.na(n_top) || n_top < 100) n_top <- 2000L
         max_cells <- as.integer(params$max_cells %||% 5000)
         if (is.na(max_cells) || max_cells < 100) max_cells <- 5000L
         root_cell <- trimws(as.character(params$root_cell %||% ""))

         b <- sc_load(ctx, need_counts = TRUE, max_cells = max_cells, max_genes = 5000)
         ctx$notes <- c(ctx$notes,
                        sprintf("轨迹工具对规模收紧：细胞 <= %d（本工具默认，可在 params.max_cells 调整）", max_cells))

         # 细胞矩阵（genes × cells）；官方接受 raw counts，preprocess 内做 log 标准化
         expr <- b$counts
         if (is.null(expr)) expr <- b$data
         expr_gc <- methods::as(Matrix::t(expr), "CsparseMatrix")
         rownames(expr_gc) <- b$genes
         colnames(expr_gc) <- b$cells

         cell_meta <- data.frame(row.names = b$cells)
         for (column in names(b$obs)) {
           cell_meta[[column]] <- as.character(b$obs[[column]])
         }
         gene_meta <- data.frame(gene_short_name = b$genes, row.names = b$genes)

         t0 <- Sys.time()
         cds <- tryCatch(
           monocle3::new_cell_data_set(expression_data = expr_gc,
                                       cell_metadata = cell_meta,
                                       gene_metadata = gene_meta),
           error = function(e) sc_prerequisite(sprintf("new_cell_data_set 失败：%s", conditionMessage(e)))
         )
         cds <- tryCatch(
           monocle3::preprocess_cds(cds, num_dim = min(30L, n_top), norm_method = "log",
                                    use_genes = NULL, verbose = FALSE),
           error = function(e) sc_prerequisite(sprintf("preprocess_cds 失败：%s", conditionMessage(e)))
         )
         cds <- tryCatch(
           monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA",
                                      umap.n_neighbors = 15L, verbose = FALSE,
                                      cores = 1L),
           error = function(e) sc_prerequisite(sprintf("reduce_dimension(UMAP) 失败：%s", conditionMessage(e)))
         )
         cds <- tryCatch(
           monocle3::cluster_cells(cds, resolution = monocle_res, random_seed = ctx$seed,
                                   verbose = FALSE),
           error = function(e) sc_prerequisite(sprintf("cluster_cells 失败：%s", conditionMessage(e)))
         )
         cds <- tryCatch(
           monocle3::learn_graph(cds, use_partition = FALSE, verbose = FALSE),
           error = function(e) sc_prerequisite(sprintf("learn_graph 失败：%s", conditionMessage(e)))
         )

         # ---- 选根并排序
         root_label <- NULL
         if (nzchar(root_cell)) {
           if (root_cell %in% b$cells) {
             cds <- tryCatch(monocle3::order_cells(cds, root_cells = root_cell, verbose = FALSE),
                             error = function(e) sc_prerequisite(sprintf("order_cells(root_cells) 失败：%s", conditionMessage(e))))
             root_label <- sprintf("用户指定根细胞 %s", root_cell)
           } else {
             cds <- tryCatch(monocle3::order_cells(cds, root_pr_nodes = root_cell, verbose = FALSE),
                             error = function(e) sc_prerequisite(sprintf("order_cells(root_pr_nodes) 失败：%s", conditionMessage(e))))
             root_label <- sprintf("用户指定根节点 %s", root_cell)
           }
         } else {
           principal <- tryCatch(monocle3::principal_graph(cds), error = function(e) NULL)
           root_node <- if (!is.null(principal)) sc_monocle_pick_root(principal) else NULL
           if (is.null(root_node)) {
             sc_prerequisite("learn_graph 未产生主图（数据可能过于离散），无法自动选根；可尝试 params.root_cell")
           }
           cds <- tryCatch(monocle3::order_cells(cds, root_pr_nodes = root_node, verbose = FALSE),
                           error = function(e) sc_prerequisite(sprintf("order_cells(root_pr_nodes) 失败：%s", conditionMessage(e))))
           root_label <- sprintf("自动选根：最大连通分量端点 %s", root_node)
         }
         elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)

         pt <- as.numeric(monocle3::pseudotime(cds))
         cl <- as.character(monocle3::clusters(cds))
         part <- as.character(monocle3::partitions(cds))
         names(pt) <- b$cells
         if (all(!is.finite(pt))) sc_prerequisite("order_cells 后没有有限伪时间（主图可能未连通）")

         # ---- obs_extra
         obs_cols <- list(
           monocle3_pseudotime = pt,
           monocle3_cluster = cl,
           monocle3_partition = part
         )
         sc_write_obs(ctx, obs_cols, b$cells)

         # ---- 簇汇总表
         ok <- is.finite(pt)
         tab <- table(cl[ok])
         cluster_df <- data.frame(
           cluster = names(tab),
           partition = vapply(names(tab), function(cc) {
             vals <- part[cl == cc]
             names(sort(table(vals), decreasing = TRUE))[1]
           }, character(1)),
           n_cells = as.integer(tab),
           pseudotime_min = as.numeric(tapply(pt[ok], cl[ok], min)[names(tab)]),
           pseudotime_median = as.numeric(tapply(pt[ok], cl[ok], stats::median)[names(tab)]),
           pseudotime_max = as.numeric(tapply(pt[ok], cl[ok], max)[names(tab)]),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         cluster_df[, 4:6] <- lapply(cluster_df[, 4:6], function(x) sc_round(x, 4))

         # ---- 图：伪时间轨迹 + cluster 轨迹
         pt_plot <- tryCatch(
           suppressWarnings(monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                                                show_trajectory_graph = TRUE,
                                                label_groups_by_cluster = FALSE,
                                                label_leaves = FALSE, label_branch_points = FALSE,
                                                cell_size = 0.7)),
           error = function(e) NULL
         )
         cl_plot <- tryCatch(
           suppressWarnings(monocle3::plot_cells(cds, color_cells_by = "cluster",
                                                show_trajectory_graph = TRUE,
                                                label_groups_by_cluster = TRUE,
                                                label_leaves = FALSE, label_branch_points = TRUE,
                                                cell_size = 0.7)),
           error = function(e) NULL
         )
         figures <- list()
         if (!is.null(pt_plot)) {
           p <- sc_png(ctx, "monocle3_pseudotime", width = 1250, height = 950, draw = function() print(pt_plot))
           figures <- c(figures, list(sc_figure("monocle3_pseudotime", "monocle3 伪时间轨迹", p)))
         }
         if (!is.null(cl_plot)) {
           p <- sc_png(ctx, "monocle3_clusters", width = 1250, height = 950, draw = function() print(cl_plot))
           figures <- c(figures, list(sc_figure("monocle3_clusters", "monocle3 簇与主图", p)))
         }
         if (length(figures) == 0) {
           sc_prerequisite("plot_cells 无法出图（monocle3 绘图失败）")
         }

         metrics <- list(
           n_cells = nrow(b$data), n_genes_used = nrow(expr_gc),
           n_clusters = length(unique(cl)), n_partitions = length(unique(part)),
           n_cells_with_pseudotime = sum(ok),
           pseudotime_min = round(min(pt, na.rm = TRUE), 4),
           pseudotime_max = round(max(pt, na.rm = TRUE), 4),
           elapsed_sec = elapsed
         )
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("monocle3 官方流程完成：%d 细胞 × %d 基因（耗时 %.1f 秒）",
                                    metrics$n_cells, metrics$n_genes_used, elapsed),
                            sprintf("cluster_cells：%d 个簇，%d 个 partition（resolution=%.3g，monocle3 内部 %.6g）",
                                    metrics$n_clusters, metrics$n_partitions, resolution, monocle_res),
                            sprintf("%s", root_label),
                            sprintf("伪时间范围：%.2f – %.2f（%d/%d 个细胞有有限伪时间）",
                                    metrics$pseudotime_min, metrics$pseudotime_max,
                                    metrics$n_cells_with_pseudotime, metrics$n_cells),
                            "obs_extra：monocle3_pseudotime / monocle3_cluster / monocle3_partition 可合并进 obs"
                          ),
                          metrics = metrics)
         res$title <- "monocle3 官方轨迹"
         res$method <- sprintf("monocle3 %s（官方 new_cell_data_set/preprocess_cds/reduce_dimension/cluster_cells/learn_graph/order_cells）",
                               sc_version("monocle3"))
         res$tables <- list(sc_table(
           "clusters", "monocle3 簇与伪时间", cluster_df,
           sc_cols(c("cluster", "Cluster", "string"), c("partition", "Partition", "string"),
                   c("n_cells", "Cells", "number"), c("pseudotime_min", "PT min", "number"),
                   c("pseudotime_median", "PT median", "number"), c("pseudotime_max", "PT max", "number")),
           max_rows = 100, note = sprintf("learn_graph(use_partition=FALSE)，cluster_cells resolution=%.6g", monocle_res)
         ))
         res$figures <- figures
         res$notes <- c(res$notes,
                        sprintf("根选择：%s", root_label),
                        "monocle3 的 UMAP 由 reduce_dimension 重新计算（不复用 bundle 的 embed.csv）",
                        "伪时间依赖 learn_graph 的主图连通性；离散/简单数据上轨迹可能退化为单线")
         res
       })
