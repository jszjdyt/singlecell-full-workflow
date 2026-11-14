#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — slingshot 官方谱系轨迹（CONTRACT §10.3）
#
# 工具：slingshot 2.14.0（官方 slingshot 谱系推断 + 平滑伪时间）
# 输入：bundle（embed.csv 优先；缺失则用 pca.csv 前两维或现场 UMAP）+ params：
#   cluster_col 聚类列（默认 clusters）
#   embed       可选：embed|pca|auto（默认 auto：embed.csv -> pca.csv -> 现场 UMAP）
# 输出：每细胞各谱系伪时间写 obs_extra.csv（sling_pseudotime_<lineage>）
#       + 嵌入图（ cluster / pseudotime + 谱系曲线）+ 谱系汇总表
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

sc_run("slingshot", "slingshot 官方谱系轨迹", "slingshot 2.14.0（官方 slingshot 谱系推断 + 平滑伪时间）",
       function(ctx) {
         params <- ctx$params
         sc_require(c("slingshot", "SingleCellExperiment", "S4Vectors"))
         embed_kind <- tolower(as.character(params$embed %||% "auto"))
         if (!embed_kind %in% c("auto", "embed", "pca", "umap")) embed_kind <- "auto"

         b <- sc_load(ctx, max_cells = 15000, max_genes = 5000)
         cluster_col <- sc_pick_column(
           b$obs,
           c(as.character(params$cluster_col %||% ""), "clusters", "cluster", "seurat_clusters", "cell_type"),
           what = "聚类列（params.cluster_col）"
         )
         grp <- sc_groups(b$obs, cluster_col)
         levels <- grp$levels
         if (length(levels) < 2) sc_prerequisite(sprintf("聚类列 %s 只有 1 个簇，无法推断轨迹", cluster_col))
         if (length(levels) > 60) sc_prerequisite("簇数过多（>60），slingshot 谱系推断不适合；请先粗聚类")

         # ---- 1. 降维坐标
         coords <- NULL
         coord_label <- NULL
         if (embed_kind %in% c("auto", "embed") && !is.null(b$embed) && ncol(b$embed) >= 2) {
           coords <- b$embed[, 1:2, drop = FALSE]
           coord_label <- "bundle 的 embed.csv（X_emb 前两维）"
         } else if (embed_kind %in% c("auto", "pca") && !is.null(b$pca) && ncol(b$pca) >= 2) {
           coords <- b$pca[, 1:2, drop = FALSE]
           coord_label <- "bundle 的 pca.csv 前两个主成分"
         }
         if (is.null(coords)) {
           if (!requireNamespace("uwot", quietly = TRUE)) {
             sc_prerequisite("bundle 无 embed.csv/pca.csv 且 uwot 不可用，无法计算降维坐标")
           }
           var <- sc_gene_variance(b$data)
           top <- order(-var, seq_along(var))[seq_len(min(2000L, length(var)))]
           dense <- as.matrix(b$data[, sort(top), drop = FALSE])
           pca <- stats::prcomp(dense, center = TRUE, scale. = FALSE, rank. = 20)$x
           set.seed(ctx$seed)
           coords <- uwot::umap(pca, n_neighbors = 15L, min_dist = 0.3, seed = ctx$seed,
                                verbose = FALSE, n_threads = 2L)
           coord_label <- "现场 UMAP（高变基因 PCA 后）"
           ctx$notes <- c(ctx$notes, "bundle 无 embed.csv/pca.csv，已现场计算 UMAP")
         }
         coords <- as.matrix(coords)
         rownames(coords) <- b$cells
         colnames(coords) <- c("dim1", "dim2")

         # ---- 2. 官方 slingshot
         t0 <- Sys.time()
         sling <- tryCatch(
           slingshot::slingshot(coords, clusterLabels = grp$labels,
                                reducedDim = "UMAP", start.clus = params$start_clus %||% NULL),
           error = function(e) sc_prerequisite(sprintf("slingshot 谱系推断失败：%s", conditionMessage(e)))
         )
         elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
         pt <- slingshot::slingPseudotime(sling)
         pt <- as.matrix(pt)
         if (is.null(pt) || ncol(pt) == 0) sc_prerequisite("slingshot 没有推断出任何谱系")
         rownames(pt) <- b$cells
         lineages <- colnames(pt)
         if (is.null(lineages) || all(is.na(lineages))) {
           lineages <- paste0("lineage", seq_len(ncol(pt)))
           colnames(pt) <- lineages
         }
         curves <- slingshot::slingCurves(sling)

         # ---- 3. obs_extra
         obs_cols <- list()
         for (j in seq_len(ncol(pt))) {
           obs_cols[[paste0("sling_pseudotime_", j)]] <- as.numeric(pt[, j])
         }
         # 平均伪时间（多谱系时方便嵌入图着色）
         obs_cols[["sling_pseudotime_mean"]] <- rowMeans(pt, na.rm = TRUE)
         sc_write_obs(ctx, obs_cols, b$cells)

         # ---- 4. 谱系汇总表
         lineage_df <- data.frame(
           lineage = lineages,
           n_cells = as.integer(colSums(is.finite(pt))),
           pseudotime_min = sc_round(apply(pt, 2, function(x) suppressWarnings(min(x, na.rm = TRUE))), 4),
           pseudotime_max = sc_round(apply(pt, 2, function(x) suppressWarnings(max(x, na.rm = TRUE))), 4),
           pseudotime_median = sc_round(apply(pt, 2, function(x) suppressWarnings(stats::median(x, na.rm = TRUE))), 4),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         # 起点/终点簇：用每条曲线首尾最近细胞判定
         cluster_vec <- grp$labels
         for (j in seq_along(curves)) {
           pts <- curves[[j]]$s
           if (is.null(pts) || nrow(pts) < 2) next
           start_cell <- which.min((coords[, 1] - pts[1, 1])^2 + (coords[, 2] - pts[1, 2])^2)
           end_cell <- which.min((coords[, 1] - pts[nrow(pts), 1])^2 + (coords[, 2] - pts[nrow(pts), 2])^2)
           lineage_df$start_cluster[j] <- cluster_vec[start_cell]
           lineage_df$end_cluster[j] <- cluster_vec[end_cell]
         }

         # ---- 5. 图：cluster + 伪时间 + 谱系曲线
         cluster_levels <- sort(unique(cluster_vec), method = "radix")
         cluster_cols <- stats::setNames(sc_palette(length(cluster_levels)), cluster_levels)
         trajectory_path <- sc_png(ctx, "slingshot_trajectory", width = 1600, height = 800, draw = function() {
           graphics::par(bg = "white", mfrow = c(1, 2), mar = c(4.0, 4.2, 3.0, 1.0), las = 1, cex = 0.8)
           # 左：cluster + 曲线
           plot(coords[, 1], coords[, 2], pch = 16, cex = 0.45,
                col = cluster_cols[cluster_vec], xlab = "UMAP 1", ylab = "UMAP 2",
                main = "slingshot lineages", axes = FALSE)
           graphics::box(col = "grey70")
           for (j in seq_along(curves)) {
             pts <- curves[[j]]$s
             if (!is.null(pts) && nrow(pts) >= 2) {
               graphics::lines(pts[, 1], pts[, 2], lwd = 3, col = sc_palette(max(3, length(curves)))[j])
             }
           }
           graphics::legend("topright", legend = cluster_levels, col = cluster_cols, pch = 16,
                            bty = "n", cex = 0.6, title = cluster_col, ncol = min(2, ceiling(length(cluster_levels) / 8)))
           # 右：平均伪时间
           pt_mean <- rowMeans(pt, na.rm = TRUE)
           pt_mean[!is.finite(pt_mean)] <- NA
           ramp <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
           rng <- range(pt_mean, na.rm = TRUE)
           idx <- if (diff(rng) > 0) round((pt_mean - rng[1]) / diff(rng) * 99) + 1 else rep(50, length(pt_mean))
           ord <- order(pt_mean)
           plot(coords[ord, 1], coords[ord, 2], pch = 16, cex = 0.45,
                col = ramp[ifelse(is.na(idx[ord]), 50, idx[ord])],
                xlab = "UMAP 1", ylab = "UMAP 2",
                main = "mean slingshot pseudotime", axes = FALSE)
           graphics::box(col = "grey70")
           for (j in seq_along(curves)) {
             pts <- curves[[j]]$s
             if (!is.null(pts) && nrow(pts) >= 2) {
               graphics::lines(pts[, 1], pts[, 2], lwd = 2, col = "grey20", lty = 2)
             }
           }
         })

         # 每条谱系单独的伪时间图（最多 4 条）
         per_lineage_path <- NULL
         n_show <- min(4L, ncol(pt))
         if (n_show >= 1) {
           per_lineage_path <- sc_png(ctx, "slingshot_lineages", width = 1600, height = 520 * ceiling(n_show / 2), draw = function() {
             graphics::par(bg = "white", mfrow = c(ceiling(n_show / 2), 2),
                           mar = c(3.6, 3.8, 2.4, 0.8), las = 1, cex = 0.72)
             ramp <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
             for (j in seq_len(n_show)) {
               vals <- as.numeric(pt[, j])
               rng <- range(vals, na.rm = TRUE)
               idx <- if (diff(rng) > 0) round((vals - rng[1]) / diff(rng) * 99) + 1 else rep(50, length(vals))
               plot(coords[, 1], coords[, 2], pch = 16, cex = 0.35,
                    col = ramp[ifelse(is.na(idx), 50, idx)],
                    xlab = "UMAP 1", ylab = "UMAP 2",
                    main = sprintf("Lineage %d pseudotime", j), axes = FALSE)
               graphics::box(col = "grey70")
               pts <- curves[[j]]$s
               if (!is.null(pts) && nrow(pts) >= 2) {
                 graphics::lines(pts[, 1], pts[, 2], lwd = 2.4, col = "grey15")
               }
             }
           })
         }

         metrics <- list(
           n_cells = nrow(pt), n_lineages = ncol(pt),
           n_cells_with_pseudotime = as.integer(sum(is.finite(rowMeans(pt, na.rm = FALSE)))),
           elapsed_sec = elapsed
         )
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("轨迹坐标：%s", coord_label),
                            sprintf("聚类列：%s（%d 个簇）", cluster_col, length(cluster_levels)),
                            sprintf("推断谱系数：%d 条（slingshot 官方实现，耗时 %.1f 秒）", ncol(pt), elapsed),
                            sprintf("伪时间范围：%.2f – %.2f（跨谱系平均）",
                                    min(pt, na.rm = TRUE), max(pt, na.rm = TRUE)),
                            sprintf("每细胞谱系伪时间已写 obs_extra（%d 列）", length(obs_cols))
                          ),
                          metrics = metrics)
         res$title <- "slingshot 官方谱系轨迹"
         res$method <- sprintf("slingshot %s（官方 slingshot：谱系推断 + 平滑伪时间）", sc_version("slingshot"))
         res$tables <- list(sc_table(
           "lineages", "谱系汇总", lineage_df,
           sc_cols(c("lineage", "Lineage", "string"), c("n_cells", "Cells", "number"),
                   c("pseudotime_min", "PT min", "number"), c("pseudotime_max", "PT max", "number"),
                   c("pseudotime_median", "PT median", "number"),
                   c("start_cluster", "Start cluster", "string"),
                   c("end_cluster", "End cluster", "string")),
           max_rows = 40, note = "伪时间按谱系分别计算；start/end 为曲线首尾最近细胞的簇"
         ))
         res$figures <- list(sc_figure("trajectory", "slingshot 谱系与伪时间", trajectory_path))
         if (!is.null(per_lineage_path)) {
           res$figures <- c(res$figures, list(sc_figure("lineages", "各谱系伪时间（前 4 条）", per_lineage_path)))
         }
         res$notes <- c(res$notes,
                        sprintf("轨迹坐标来源：%s", coord_label),
                        "slingshot 在给定降维坐标上推断谱系平滑曲线；伪时间未做根细胞生物学指定（无先验根簇）")
         res
       })
