#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — harmony 官方批次校正（CONTRACT §10.3）
#
# 工具：harmony 1.2.4（官方 RunHarmony）
# 输入：bundle（优先 pca.csv；缺失则现场 PCA）+ params：
#   batch  批次列（默认 batch）
#   dims   使用的主成分数（默认 20）
#   theta  多样性惩罚（默认 2；params.theta 可调）
#   max_cells 兜底上限（默认 15000）
# 输出：修正后的 PCA 矩阵（X_pca.csv 与 pca.csv）+ 批次混合指标表 + 校正前后 UMAP 图
#       obs_extra.csv：harmony_1..harmony_k（每细胞校正后前 k 个主成分）
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

# 批次混合分数：kNN 邻域里「同批次」占比越低越混合；1 = 完全随机混合
sc_batch_mixing <- function(coords, batches, k = 30) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n < k + 1) k <- max(2L, min(5L, n - 1L))
  nn <- FNN::get.knn(coords, k = k)$nn.index
  same <- rowMeans(matrix(batches[nn], nrow = n) == batches)
  p <- table(batches) / length(batches)
  expected_random <- sum(p^2)
  same_mean <- mean(same)
  # 1) 原始混合分数：1 - 同批次近邻占比（越高越混合）
  mixing <- 1 - same_mean
  # 2) 归一化分离指数：0 = 与随机混合一致，1 = 完全按批次分离（越低越好）
  separation <- (same_mean - expected_random) / max(1 - expected_random, 1e-9)
  list(mixing = mixing, separation = max(0, min(1, separation)),
       same_batch = same_mean, expected_random = expected_random, k = k)
}

sc_run("harmony", "harmony 官方批次校正", "harmony 1.2.4（官方 RunHarmony，kNN 批次混合指标 + UMAP 前后对比）",
       function(ctx) {
         params <- ctx$params
         sc_require("harmony")
         dims <- as.integer(params$dims %||% 20)
         if (is.na(dims) || dims < 2) dims <- 20L
         theta <- as.numeric(params$theta %||% 2)
         if (!is.finite(theta) || theta <= 0) theta <- 2
         max_cells <- as.integer(params$max_cells %||% 15000)
         if (is.na(max_cells) || max_cells < 50) max_cells <- 15000L

         b <- sc_load(ctx, max_cells = max_cells, max_genes = 5000)
         batch_col <- sc_pick_column(
           b$obs,
           c(as.character(params$batch %||% ""), "batch", "Batch", "sample", "sample_id", "orig.ident"),
           what = "批次列（params.batch）"
         )
         batches <- sc_obs_chr(b$obs, batch_col)
         batches[is.na(batches) | !nzchar(batches)] <- "NA"
         tab <- table(batches)
         if (length(tab) < 2) {
           sc_prerequisite(sprintf("批次列 %s 只有 %d 个批次，无法做批次校正", batch_col, length(tab)))
         }
         if (min(tab) < 10) ctx$warnings <- c(ctx$warnings, sprintf("最小批次只有 %d 个细胞，校正结果可能不稳定", min(tab)))

         # ---- 1. 取得 PCA
         pca <- b$pca
         pca_source <- "bundle 的 pca.csv"
         if (is.null(pca)) {
           # 现场 PCA：高变基因 + prcomp（上限 15k 细胞 × 2000 基因）
           var <- sc_gene_variance(b$data)
           top <- order(-var, seq_along(var))[seq_len(min(2000L, length(var)))]
           sub <- as.matrix(b$data[, sort(top), drop = FALSE])
           mu <- colMeans(sub)
           sdv <- apply(sub, 2, stats::sd)
           sdv[sdv == 0] <- 1
           pca <- stats::prcomp((sub - rep(mu, each = nrow(sub))) / rep(sdv, each = nrow(sub)),
                                center = FALSE, scale. = FALSE, rank. = max(dims, 20))$x
           pca_source <- sprintf("现场 PCA（%d 个高变基因）", length(top))
           ctx$notes <- c(ctx$notes, "bundle 无 pca.csv，已现场计算 PCA 后再做 harmony 校正")
         }
         pca <- as.matrix(pca)
         rownames(pca) <- b$cells
         if (ncol(pca) < dims) dims <- ncol(pca)
         if (dims < 2) sc_prerequisite("PCA 维度不足（< 2）")

         # ---- 2. harmony 官方校正
         set.seed(ctx$seed)
         t0 <- Sys.time()
         harmony_mat <- tryCatch(
           harmony::RunHarmony(
             data_mat = pca[, seq_len(dims), drop = FALSE],
             meta_data = data.frame(batch = batches, row.names = b$cells),
             vars_use = "batch",
             theta = theta, ncores = 1, verbose = FALSE
           ),
           error = function(e) sc_prerequisite(sprintf("harmony::RunHarmony 失败：%s", conditionMessage(e)))
         )
         harmony_mat <- as.matrix(harmony_mat)
         rownames(harmony_mat) <- b$cells
         elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)

         # ---- 3. 批次混合指标（校正前后）
         before <- sc_batch_mixing(pca[, seq_len(dims), drop = FALSE], batches)
         after <- sc_batch_mixing(harmony_mat, batches)
         ctx$notes <- c(ctx$notes,
                        "批次指标：mixing = 1 - 同批次近邻占比（越高越混合）；separation = 归一化批次分离指数（0=随机混合，1=完全分离）")

         # ---- 4. 校正后 PCA 落盘（X_pca.csv = 任务要求的文件名；pca.csv 供其它适配器复用）
         write_pca <- function(mat, filename) {
           df <- data.frame(cell = b$cells, check.names = FALSE)
           for (j in seq_len(ncol(mat))) df[[paste0("PC", j)]] <- as.numeric(mat[, j])
           path <- file.path(ctx$out, filename)
           utils::write.csv(df, path, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
           path
         }
         x_pca_path <- write_pca(harmony_mat, "X_pca.csv")
         write_pca(harmony_mat, "pca.csv")

         # ---- 5. obs_extra（前 5 个校正后主成分；完整矩阵在 X_pca.csv）
         k_obs <- min(5L, ncol(harmony_mat))
         obs_cols <- list()
         for (j in seq_len(k_obs)) obs_cols[[paste0("harmony_", j)]] <- as.numeric(harmony_mat[, j])
         sc_write_obs(ctx, obs_cols, b$cells)

         # ---- 6. 前后 UMAP 对比图
         umap_before <- tryCatch(
           uwot::umap(pca[, seq_len(dims), drop = FALSE], n_neighbors = 15L,
                      min_dist = 0.3, seed = ctx$seed, verbose = FALSE,
                      n_threads = 2L),
           error = function(e) NULL
         )
         umap_after <- tryCatch(
           uwot::umap(harmony_mat, n_neighbors = 15L, min_dist = 0.3,
                      seed = ctx$seed, verbose = FALSE, n_threads = 2L),
           error = function(e) NULL
         )
         before_after_path <- NULL
         if (!is.null(umap_before) && !is.null(umap_after)) {
           batch_levels <- sort(unique(batches), method = "radix")
           cols <- stats::setNames(sc_palette(length(batch_levels)), batch_levels)
           before_after_path <- sc_png(ctx, "harmony_before_after", width = 1600, height = 760, draw = function() {
             graphics::par(bg = "white", mfrow = c(1, 2), mar = c(4.0, 4.2, 3.0, 1.0), las = 1, cex = 0.8)
             for (panel in seq_along(list(umap_before, umap_after))) {
               um <- list(umap_before, umap_after)[[panel]]
               lab <- c("Before harmony", "After harmony")[panel]
               plot(um[, 1], um[, 2], pch = 16, cex = 0.45, col = cols[batches],
                    xlab = "UMAP 1", ylab = "UMAP 2",
                    main = sprintf("%s (mixing=%.2f)", lab, c(before$mixing, after$mixing)[panel]),
                    axes = FALSE)
               graphics::box(col = "grey70")
               graphics::legend("topright", legend = batch_levels, col = cols, pch = 16,
                                bty = "n", cex = 0.62, title = batch_col)
             }
           })
         }

         # ---- 7. 表
         metrics_df <- data.frame(
           metric = c("batch mixing score (1 - same-batch fraction; higher = more mixed)",
                      "normalized separation index (0 = random mixed, 1 = fully separated)",
                      "same-batch neighbor fraction",
                      "expected random same-batch fraction",
                      "n cells", "n batches", "n PCs", "theta", "elapsed sec"),
           before = c(round(before$mixing, 4), round(before$separation, 4), round(before$same_batch, 4),
                      round(before$expected_random, 4), nrow(pca), length(tab), dims, theta, NA),
           after = c(round(after$mixing, 4), round(after$separation, 4), round(after$same_batch, 4),
                     round(after$expected_random, 4), nrow(pca), length(tab), dims, theta, elapsed),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         batch_df <- data.frame(
           batch = names(tab), n_cells = as.integer(tab),
           stringsAsFactors = FALSE, check.names = FALSE
         )

         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("批次列：%s（%d 个批次，%d 个细胞）", batch_col, length(tab), nrow(pca)),
                            sprintf("PCA 来源：%s（使用前 %d 个主成分）", pca_source, dims),
                            sprintf("批次混合分数（越高越混合）：%.3f → %.3f", before$mixing, after$mixing),
                            sprintf("同批次近邻占比：%.2f → %.2f（随机期望 %.2f）",
                                    before$same_batch, after$same_batch, before$expected_random),
                            sprintf("校正后 PCA 已写出：X_pca.csv（%d × %d），耗时 %.1f 秒",
                                    nrow(harmony_mat), ncol(harmony_mat), elapsed),
                            sprintf("obs_extra：harmony_1..harmony_%d 可合并进 obs", k_obs)
                          ),
                          metrics = list(
                            n_cells = nrow(pca), n_batches = length(tab), dims = dims,
                            mixing_before = round(before$mixing, 4),
                            mixing_after = round(after$mixing, 4),
                            separation_before = round(before$separation, 4),
                            separation_after = round(after$separation, 4),
                            elapsed_sec = elapsed
                          ))
         res$title <- "harmony 官方批次校正"
         res$method <- sprintf("harmony %s（官方 RunHarmony，do_pca=FALSE，theta=%.1f）", sc_version("harmony"), theta)
         res$tables <- list(
           sc_table("mixing_metrics", "批次混合指标", metrics_df,
                    sc_cols(c("metric", "Metric", "string"), c("before", "Before", "number"),
                            c("after", "After", "number")),
                    max_rows = 20, note = sprintf("kNN (k=%d) 同批次近邻占比归一化", after$k)),
           sc_table("batch_counts", "各批次细胞数", batch_df,
                    sc_cols(c("batch", "Batch", "string"), c("n_cells", "Cells", "number")),
                    max_rows = 100)
         )
         res$notes <- c(res$notes,
                        sprintf("输入 PCA：%s", pca_source),
                        sprintf("校正后矩阵文件：%s", x_pca_path),
                        "harmony 采用随机初始化，已固定 seed；ncores=1 保证可复现")
         if (!is.null(before_after_path)) {
           res$figures <- list(sc_figure("harmony_before_after", "批次校正前后 UMAP 对比", before_after_path))
         }
         res
       })
