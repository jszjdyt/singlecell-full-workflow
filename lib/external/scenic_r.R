#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell — SCENIC 官方 R 转录调控（CONTRACT §10.3）
#
# 工具链（全部官方 R 包）：
#   GENIE3 1.28.0  共表达网络（TF → target 权重）
#   RcisTarget 1.26.0  motif 富集（mc9nr .feather rankings 数据库）
#   AUCell 1.28.0  regulon 活性打分
# 说明：这是 SCENIC 官方三段式（GENIE3 → RcisTarget → AUCell）的直接串联，
#       等价于 SCENIC 包的 runSCENIC_1/2/3，但不生成 SCENIC 中间目录（更轻量）。
# 输入：bundle（data.mtx=log-normalized）+ params：
#   species        human | mouse（默认 human）
#   n_top_genes    高变基因数（默认 1000，clamp 200..1500）
#   n_trees        GENIE3 树数（默认 100）
#   n_tfs          参与 regulon 推断的 TF 数（按 GENIE3 总权重排序，默认 50）
#   cells_subsample GENIE3/AUCell 使用的细胞上限（默认 1500）
#   nes_threshold  motif 富集 NES 阈值（默认 3）
#   db_dir         motif 数据库目录（留空=自动探测 DS_SC_SCENIC_DB / 默认用户目录）
# 输出：TF-regulon 表 + motif 表 + AUCell 活性 obs_extra + 活性热图 + regulon 条形图
# 注意：SCENIC 较重；对 demo/小数据会自动限制基因/细胞/TF 数并在 notes 说明。
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

# 找 RcisTarget 用的 R 版 mc9nr .feather 数据库（hg38/mm10）
sc_find_scenic_feather <- function(db_dir, species) {
  prefix <- if (species == "mouse") "mm10" else "hg38"
  candidates <- unique(c(
    db_dir,
    if (nzchar(db_dir)) dirname(db_dir) else "",
    if (species == "mouse") path.expand("~/scenic_databases/mouse")
    else path.expand("~/scenic_databases/human"),
    path.expand("~/cisTarget_databases")
  ))
  candidates <- candidates[nzchar(candidates)]
  for (directory in candidates) {
    if (!dir.exists(directory)) next
    files <- list.files(directory, pattern = "\\.feather$", full.names = TRUE)
    files <- files[grepl("mc9nr", basename(files)) & startsWith(basename(files), prefix)]
    if (length(files) == 0) next
    # 优先 500bp（更快、常用），其次 10kb
    ord <- order(!grepl("500bp", basename(files)), basename(files))
    return(files[ord][1L])
  }
  NULL
}

sc_load_motif_annotations <- function(species) {
  item <- if (species == "mouse") "motifAnnotations_mgi_v9" else "motifAnnotations_hgnc_v9"
  env <- new.env(parent = emptyenv())
  okay <- tryCatch({
    utils::data(list = item, package = "RcisTarget", envir = env)
    TRUE
  }, error = function(e) FALSE)
  if (!okay) {
    sc_prerequisite(sprintf("RcisTarget 缺少 motif 注释数据 %s（请检查 RcisTarget 安装）", item))
  }
  annot <- get(item, envir = env)
  # RcisTarget::cisTarget 要求 motifAnnot 是 data.table
  if (!data.table::is.data.table(annot)) annot <- data.table::as.data.table(annot)
  annot
}

sc_run("scenic_r", "SCENIC(R) 官方转录调控",
       "SCENIC 1.3.1 官方 R 链路（GENIE3 1.28.0 + RcisTarget 1.26.0 + AUCell 1.28.0，mc9nr rankings）",
       function(ctx) {
         params <- ctx$params
         sc_require(c("SCENIC", "GENIE3", "RcisTarget", "AUCell"))
         sink(stderr(), type = "output")
         on.exit(try(sink(type = "output"), silent = TRUE), add = TRUE)

         species <- tolower(as.character(params$species %||% "human"))
         if (!species %in% c("human", "mouse")) species <- "human"
         n_top <- as.integer(params$n_top_genes %||% 1000)
         if (is.na(n_top) || n_top < 100) n_top <- 1000L
         n_top <- max(200L, min(n_top, 1500L))
         n_trees <- as.integer(params$n_trees %||% 100)
         if (is.na(n_trees) || n_trees < 10) n_trees <- 100L
         n_tfs <- as.integer(params$n_tfs %||% 50)
         if (is.na(n_tfs) || n_tfs < 3) n_tfs <- 50L
         cells_sub <- as.integer(params$cells_subsample %||% 1500)
         if (is.na(cells_sub) || cells_sub < 100) cells_sub <- 1500L
         nes_thr <- as.numeric(params$nes_threshold %||% 3)
         if (!is.finite(nes_thr) || nes_thr < 1) nes_thr <- 3
         module_top <- 50L

         # ---- 1. 数据库
         db_dir <- trimws(as.character(params$db_dir %||% Sys.getenv("DS_SC_SCENIC_DB")))
         if (!nzchar(db_dir)) {
           db_dir <- path.expand("~/scenic_databases/human")
         }
         db_dir <- path.expand(db_dir)
         feather <- sc_find_scenic_feather(db_dir, species)
         if (is.null(feather)) {
           sc_prerequisite(sprintf("未找到 SCENIC RcisTarget mc9nr 数据库（%s, 目录：%s）", species, db_dir))
         }
         db_genes <- tryCatch(RcisTarget::getColumnNames(feather),
                              error = function(e) sc_prerequisite(sprintf("读取 motif 数据库失败：%s", conditionMessage(e))))
         motfi_annot <- sc_load_motif_annotations(species)
         ctx$notes <- c(ctx$notes,
                        sprintf("motif 数据库：%s（%d 个基因/motif 排名矩阵）", basename(feather), length(db_genes)),
                        sprintf("motif 注释：%d 行（RcisTarget %s）", nrow(motfi_annot), sc_version("RcisTarget")))

         # ---- 2. 表达矩阵与基因过滤
         b <- sc_load(ctx, max_cells = 3000, max_genes = 4000)
         pct <- sc_pct_expressed(b$data)
         names(pct) <- b$genes
         expressed <- b$genes[pct >= 0.10]
         expressed_db <- intersect(expressed, db_genes)
         if (length(expressed_db) < 100) {
           sc_prerequisite(sprintf("表达基因中只有 %d 个在 motif 数据库里有排名（需 >=100）；检查基因名是否为 %s 基因符号",
                                   length(expressed_db), if (species == "human") "人类" else "小鼠"))
         }
         # 高变基因 + 全部 TF（保证 TF 一定在候选里）
         tfs_all <- unique(motfi_annot$TF)
         tfs <- intersect(tfs_all, expressed_db)
         if (length(tfs) < 3) sc_prerequisite(sprintf("表达基因里只有 %d 个 TF，无法构建 regulon", length(tfs)))
         var <- sc_gene_variance(b$data[, expressed_db, drop = FALSE])
         hvg <- expressed_db[order(-var, seq_along(var))[seq_len(min(n_top, length(expressed_db)))]]
         gene_pool <- unique(c(tfs, hvg))
         if (length(gene_pool) < 100) sc_prerequisite("候选基因太少（<100）")

         # ---- 3. 细胞抽样（GENIE3/AUCell 规模控制）
         set.seed(ctx$seed)
         if (nrow(b$data) > cells_sub) {
           use_cells <- sort(sample.int(nrow(b$data), cells_sub))
           ctx$notes <- c(ctx$notes, sprintf("GENIE3/AUCell 使用 %d/%d 个细胞（确定性抽样，控制耗时）",
                                             cells_sub, nrow(b$data)))
         } else {
           use_cells <- seq_len(nrow(b$data))
         }
         # GENIE3 只接受 base matrix（genes × cells）；规模已由 n_top_genes/cells_subsample 控制
         expr <- as.matrix(Matrix::t(b$data[use_cells, gene_pool, drop = FALSE]))
         rownames(expr) <- gene_pool
         colnames(expr) <- b$cells[use_cells]

         # ---- 4. GENIE3 共表达网络
         t0 <- Sys.time()
         weights <- tryCatch(
           suppressWarnings(GENIE3::GENIE3(expr, regulators = tfs, nTrees = n_trees,
                                           nCores = 1, verbose = FALSE)),
           error = function(e) sc_prerequisite(sprintf("GENIE3 失败：%s", conditionMessage(e)))
         )
         genie_time <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
         weights <- as.matrix(weights)

         # TF 按总权重排序，取前 n_tfs
         tf_strength <- rowSums(weights)
         tf_rank <- sort(tf_strength, decreasing = TRUE)
         tfs_use <- names(tf_rank)[seq_len(min(n_tfs, length(tf_rank)))]
         modules <- list()
         for (tf in tfs_use) {
           wt <- weights[tf, ]
           ord <- order(-wt, names(wt))
           pos <- ord[wt[ord] > 0]
           if (length(pos) < 10) next
           modules[[tf]] <- names(wt)[pos][seq_len(min(module_top, length(pos)))]
         }
         if (length(modules) < 3) sc_prerequisite(sprintf("GENIE3 只得到 %d 个可用 TF 模块（<3）", length(modules)))
         ctx$notes <- c(ctx$notes,
                        sprintf("GENIE3：%d 个 TF × %d 个 target（%d 棵树，耗时 %.1f 秒），取前 %d 个 TF 模块",
                                length(tfs), ncol(weights), n_trees, genie_time, length(modules)))

         # ---- 5. RcisTarget motif 富集
         all_genes <- intersect(unique(c(unlist(modules), names(modules))), db_genes)
         modules <- lapply(modules, function(genes) intersect(genes, all_genes))
         modules <- modules[vapply(modules, length, integer(1)) >= 10]
         if (length(modules) < 3) sc_prerequisite("可用于 motif 富集的模块太少（<3）")
         rankings <- tryCatch({
           t1 <- Sys.time()
           r <- RcisTarget::importRankings(feather, columns = all_genes)
           ctx$notes <- c(ctx$notes, sprintf("importRankings 载入 %d 个基因（%.1f 秒）",
                                              length(all_genes),
                                              as.numeric(Sys.time() - t1, units = "secs")))
           r
         }, error = function(e) sc_prerequisite(sprintf("importRankings 失败：%s", conditionMessage(e))))

         t0 <- Sys.time()
         motif_res <- tryCatch(
           as.data.frame(RcisTarget::cisTarget(modules, rankings, motifAnnot = motfi_annot,
                                               nesThreshold = nes_thr, nCores = 1, verbose = FALSE)),
           error = function(e) sc_prerequisite(sprintf("RcisTarget::cisTarget 失败：%s", conditionMessage(e)))
         )
         motif_time <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
         if (nrow(motif_res) == 0) sc_prerequisite("RcisTarget 没有找到任何富集 motif（NES 阈值可下调）")

         # ---- 6. 组装 regulon（motif 注释必须直接对应本模块 TF）
         motif_res$annotTFs <- sub(" \\(.*$", "", as.character(motif_res$TF_highConf))
         keeps <- mapply(function(tf, annot) {
           if (is.na(annot) || !nzchar(annot)) return(FALSE)
           tf %in% trimws(strsplit(annot, ";", fixed = TRUE)[[1L]])
         }, as.character(motif_res$geneSet), motif_res$annotTFs)
         motif_res <- motif_res[keeps, , drop = FALSE]
         regulons <- list()
         regulon_info <- list()
         for (tf in unique(as.character(motif_res$geneSet))) {
           rows <- motif_res[motif_res$geneSet == tf, , drop = FALSE]
           rows <- rows[order(-rows$NES), , drop = FALSE]
           genes <- unique(unlist(strsplit(as.character(rows$enrichedGenes), ";", fixed = TRUE)))
           genes <- genes[nzchar(genes) & !is.na(genes)]
           if (length(genes) == 0) next
           regulons[[tf]] <- unique(c(tf, genes))
           regulon_info[[tf]] <- list(
             tf = tf,
             n_motifs = nrow(rows),
             best_motif = as.character(rows$motif[1]),
             best_NES = as.numeric(rows$NES[1])
           )
         }
         if (length(regulons) == 0) {
           sc_prerequisite("没有形成任何 regulon（motif 注释与模块 TF 不匹配）；可降低 nes_threshold 或扩大 n_tfs")
         }

         # ---- 7. AUCell 活性
         t0 <- Sys.time()
         auc <- tryCatch({
           rankings_auc <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
           AUCell::AUCell_calcAUC(regulons, rankings_auc, verbose = FALSE)
         }, error = function(e) sc_prerequisite(sprintf("AUCell 计算失败：%s", conditionMessage(e))))
         auc_time <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
         auc_mat <- t(as.matrix(AUCell::getAUC(auc)))  # cells × regulons
         colnames(auc_mat) <- names(regulons)

         # obs_extra：取活性方差最大的前 10 个 regulon
         ord <- order(-apply(auc_mat, 2, stats::var), colnames(auc_mat))
         show_sets <- colnames(auc_mat)[ord][seq_len(min(10L, ncol(auc_mat)))]
         obs_cols <- list()
         for (name in show_sets) {
           obs_cols[[paste0("regulon_", make.names(name))]] <- as.numeric(auc_mat[, name])
         }
         sc_write_obs(ctx, obs_cols, b$cells[use_cells])

         # ---- 8. 汇总表
         regulon_df <- do.call(rbind, lapply(names(regulons), function(tf) {
           info <- regulon_info[[tf]]
           targets <- setdiff(regulons[[tf]], tf)
           data.frame(
             TF = tf,
             regulon_size = length(regulons[[tf]]),
             n_motifs = info$n_motifs,
             best_motif = info$best_motif,
             best_NES = sc_round(info$best_NES, 3),
             mean_AUC = sc_round(mean(auc_mat[, tf]), 4),
             top_targets = paste(utils::head(targets, 8), collapse = ","),
             stringsAsFactors = FALSE, check.names = FALSE
           )
         }))
         regulon_df <- regulon_df[order(-regulon_df$regulon_size, -regulon_df$best_NES), , drop = FALSE]
         motif_df <- data.frame(
           TF = as.character(motif_res$geneSet),
           motif = as.character(motif_res$motif),
           NES = sc_round(as.numeric(motif_res$NES), 3),
           AUC = sc_round(as.numeric(motif_res$AUC), 5),
           nEnrGenes = as.integer(motif_res$nEnrGenes),
           enrichedGenes = vapply(as.character(motif_res$enrichedGenes),
                                  function(x) paste(utils::head(strsplit(x, ";", fixed = TRUE)[[1]], 10), collapse = ","),
                                  character(1)),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         motif_df <- motif_df[order(motif_df$TF, -motif_df$NES), , drop = FALSE]
         tf_rank_df <- data.frame(
           TF = names(tf_rank)[seq_len(min(30L, length(tf_rank)))],
           total_weight = sc_round(as.numeric(tf_rank[seq_len(min(30L, length(tf_rank)))]), 4),
           has_regulon = names(tf_rank)[seq_len(min(30L, length(tf_rank)))] %in% names(regulons),
           stringsAsFactors = FALSE, check.names = FALSE
         )

         # ---- 9. 图：AUCell 活性热图（按 cluster 均值）
         cluster_col <- sc_pick_column(
           b$obs,
           c(as.character(params$cluster_col %||% ""), "clusters", "cluster", "cell_type"),
           what = "分群列", required = FALSE
         )
         heat_path <- NULL
         heat_mat <- NULL
         if (!is.null(cluster_col)) {
           labels <- sc_obs_chr(b$obs, cluster_col)[use_cells]
           keep <- !is.na(labels)
           if (sum(keep) > 20) {
             agg <- rowsum(auc_mat[keep, , drop = FALSE], group = labels[keep], reorder = TRUE)
             n_per <- as.numeric(table(labels[keep])[rownames(agg)])
             agg <- agg / n_per
             top <- ord[seq_len(min(12L, ncol(auc_mat)))]
             heat_mat <- t(agg[, top, drop = FALSE])
             rownames(heat_mat) <- colnames(auc_mat)[top]
             colnames(heat_mat) <- rownames(agg)
             heat_path <- sc_png(ctx, "scenic_regulon_heatmap", width = 1500, height = 900, draw = function() {
               sc_heatmap(heat_mat, title = "SCENIC regulon activity (cluster mean AUC)",
                          value_label = "AUC", row_title = "TF regulon", col_title = cluster_col,
                          cex = 0.75)
             })
           }
         }
         if (is.null(heat_path)) {
           # 无分群：抽样 120 个细胞
           take <- order(rowMeans(auc_mat))
           take <- take[round(seq(1, length(take), length.out = min(120, length(take))))]
           heat_mat <- t(auc_mat[take, ord[seq_len(min(12L, ncol(auc_mat)))]])
           heat_path <- sc_png(ctx, "scenic_regulon_heatmap", width = 1500, height = 900, draw = function() {
             sc_heatmap(heat_mat, title = "SCENIC regulon activity (sampled cells)",
                        value_label = "AUC", row_title = "TF regulon", col_title = "Cells (sorted)",
                        cex = 0.72)
           })
         }

         # ---- 10. 图：regulon 大小与活性
         reg_plot <- regulon_df[order(regulon_df$regulon_size), , drop = FALSE]
         reg_plot <- utils::head(reg_plot, min(20L, nrow(reg_plot)))
         bar_path <- sc_png(ctx, "scenic_regulons", width = 1400, height = 850, draw = function() {
           graphics::par(bg = "white", mar = c(4.2, 10, 3, 1.2), las = 1, cex = 0.85)
           graphics::barplot(reg_plot$regulon_size, horiz = TRUE, names.arg = reg_plot$TF,
                             col = sc_palette(nrow(reg_plot)), border = NA,
                             xlab = "Regulon size (genes)", cex.names = 0.8,
                             main = "SCENIC regulons (GENIE3 + RcisTarget)")
         })

         metrics <- list(
           n_cells = length(use_cells), n_genes_in_db = length(expressed_db),
           n_tfs = length(tfs), n_tf_modules = length(modules),
           n_regulons = length(regulons),
           genie3_sec = genie_time, rcistarget_sec = motif_time, aucell_sec = auc_time,
           n_trees = n_trees, nes_threshold = nes_thr
         )
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("motif 数据库：%s（%d 个基因排名）", basename(feather), length(db_genes)),
                            sprintf("表达基因中 %d 个在数据库中；候选 %d 个基因 × %d 个细胞（GENIE3）",
                                    length(expressed_db), length(gene_pool), length(use_cells)),
                            sprintf("GENIE3：%d 个 TF，%d 棵树（%.1f 秒）；取前 %d 个 TF 模块",
                                    length(tfs), n_trees, genie_time, length(modules)),
                            sprintf("RcisTarget motif 富集：%d 条显著 motif（NES>=%.1f，%.1f 秒）",
                                    nrow(motif_df), nes_thr, motif_time),
                            sprintf("最终 regulon：%d 个（AUCell 打分 %.1f 秒）", length(regulons), auc_time),
                            sprintf("obs_extra：%d 个 regulon AUC 列（可嵌入图着色）", length(obs_cols))
                          ),
                          metrics = metrics)
         res$title <- "SCENIC(R) 官方转录调控"
         res$method <- sprintf("SCENIC 官方三段式：GENIE3 %s + RcisTarget %s + AUCell %s（mc9nr rankings）",
                               sc_version("GENIE3"), sc_version("RcisTarget"), sc_version("AUCell"))
         res$tables <- list(
           sc_table("regulons", "TF-regulon 表", regulon_df,
                    sc_cols(c("TF", "TF", "string"), c("regulon_size", "Regulon size", "number"),
                            c("n_motifs", "Motifs", "number"), c("best_motif", "Best motif", "string"),
                            c("best_NES", "Best NES", "number"), c("mean_AUC", "Mean AUC", "number"),
                            c("top_targets", "Top target genes", "string")),
                    max_rows = 60, note = sprintf("NES>=%.1f；motif 注释直接对应本 TF", nes_thr)),
           sc_table("motifs", "显著 motif 富集", motif_df,
                    sc_cols(c("TF", "TF", "string"), c("motif", "Motif", "string"),
                            c("NES", "NES", "number"), c("AUC", "AUC", "number"),
                            c("nEnrGenes", "Enriched genes", "number"),
                            c("enrichedGenes", "Leading edge (top 10)", "string")),
                    max_rows = 100, note = "RcisTarget::cisTarget（motif 注释 = 高置信直接注释）"),
           sc_table("tf_ranking", "GENIE3 TF 权重排序（前 30）", tf_rank_df,
                    sc_cols(c("TF", "TF", "string"), c("total_weight", "Total weight", "number"),
                            c("has_regulon", "Has regulon", "boolean")),
                    max_rows = 30)
         )
         res$figures <- list(
           sc_figure("scenic_heatmap", "AUCell regulon 活性热图", heat_path),
           sc_figure("scenic_regulons", "SCENIC regulon 大小", bar_path)
         )
         res$notes <- c(res$notes,
                        sprintf("species=%s；motif 数据库 %s", species, feather),
                        sprintf("GENIE3 nTrees=%d；模块取每个 TF 权重前 %d 个靶基因", n_trees, module_top),
                        sprintf("细胞抽样上限 cells_subsample=%d（参数可调）", cells_sub),
                        "本实现是 SCENIC 官方三段式的直接串联（未生成 SCENIC 中间目录）")
         res
       })
