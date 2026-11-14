#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — NicheNet 官方配体-靶基因（CONTRACT §10.3）
#
# 工具：nichenetr 2.2.1（官方 predict_ligand_activities + 官方先验模型）
# 先验：<nichenetPriors>/ligand_target_matrix_nsga2r_final.rds
#       <nichenetPriors>/lr_network_human_21122021.rds
#       （环境变量 DS_SC_NICHENET_PRIORS 或 params.priors_dir 可覆盖）
# 输入：bundle（data.mtx=log-normalized）+ params：
#   cell_type_col  细胞类型列（默认 cell_type）
#   sender / receiver  发送方/接收方细胞类型（留空自动：最丰富类型为 receiver，次多为 sender）
#   top_n          输出 top 配体数（默认 20）
#   n_targets      热图展示的靶基因数（默认 50）
#   max_cells      兜底上限（默认 3000，NicheNet 较重）
#   priors_dir     先验目录（默认读 DS_SC_NICHENET_PRIORS）
#   species        human | mouse（默认 human；mouse 需 mouse 版先验）
# 输出：配体活性表（AUROC/AUPR/pearson）+ L-R 表 + 配体-靶基因表 + 2 张图
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

# 简化版 Wilcoxon 秩和检验（向量化近似 + 近似正态；足够做配体排序的 receiver DE）
sc_wilcox_de <- function(mat, group1, group2, genes, top_n = 200L) {
  mat <- methods::as(mat, "CsparseMatrix")
  x1 <- mat[group1, , drop = FALSE]
  x2 <- mat[group2, , drop = FALSE]
  n1 <- nrow(x1); n2 <- nrow(x2)
  if (n1 < 3 || n2 < 3) return(NULL)
  m1 <- as.numeric(Matrix::colSums(x1)) / n1
  m2 <- as.numeric(Matrix::colSums(x2)) / n2
  logfc <- m1 - m2
  # 逐基因 rank 和（1200-5000 个基因，demo 规模足够快；真实大数据受 max_genes 限制）
  pvals <- vapply(seq_along(genes), function(j) {
    v1 <- x1[, j]; v2 <- x2[, j]
    values <- c(as.numeric(v1), as.numeric(v2))
    ranks <- rank(values, ties.method = "average")
    r1 <- sum(ranks[seq_len(n1)])
    u <- r1 - n1 * (n1 + 1) / 2
    mu <- n1 * n2 / 2
    # 结校正
    tab <- table(values)
    ties <- sum(tab^3 - tab)
    N <- n1 + n2
    sigma2 <- n1 * n2 / 12 * ((N + 1) - ties / (N * (N - 1)))
    if (!is.finite(sigma2) || sigma2 <= 0) return(1)
    z <- (u - mu) / sqrt(sigma2)
    2 * stats::pnorm(-abs(z))
  }, numeric(1))
  padj <- stats::p.adjust(pvals, method = "BH")
  df <- data.frame(gene = genes, logfc = logfc, pval = pvals, padj = padj,
                   stringsAsFactors = FALSE, check.names = FALSE)
  df <- df[order(df$pval, -df$logfc), , drop = FALSE]
  df <- df[df$logfc > 0, , drop = FALSE]
  utils::head(df, top_n)
}

sc_run("nichenet", "NicheNet 官方配体-靶基因", "nichenetr 2.2.1（官方 predict_ligand_activities + 用户本地先验模型）",
       function(ctx) {
         params <- ctx$params
         sc_require("nichenetr")
         # nichenetr 可能向 stdout 打印加载信息：重定向到 stderr
         sink(stderr(), type = "output")
         on.exit(try(sink(type = "output"), silent = TRUE), add = TRUE)

         top_n <- as.integer(params$top_n %||% 20)
         if (is.na(top_n) || top_n < 3) top_n <- 20L
         n_targets <- as.integer(params$n_targets %||% 50)
         if (is.na(n_targets) || n_targets < 5) n_targets <- 50L
         max_cells <- as.integer(params$max_cells %||% 3000)
         if (is.na(max_cells) || max_cells < 100) max_cells <- 3000L
         species <- tolower(as.character(params$species %||% "human"))
         if (!species %in% c("human", "mouse")) species <- "human"

         # ---- 1. 先验模型
         priors <- trimws(as.character(params$priors_dir %||% Sys.getenv("DS_SC_NICHENET_PRIORS")))
         if (!nzchar(priors)) {
           priors <- path.expand("~/nichenet_priors")
         }
         priors <- path.expand(priors)
         if (!dir.exists(priors)) sc_prerequisite(sprintf("NicheNet 先验模型目录不存在：%s", priors))
         suffix <- if (species == "human") "" else "_mouse"
         ltm_file <- file.path(priors, sprintf("ligand_target_matrix_nsga2r_final%s.rds", suffix))
         lr_file <- file.path(priors, sprintf("lr_network_%s_21122021.rds", species))
         if (!file.exists(ltm_file)) sc_prerequisite(sprintf("缺少 NicheNet ligand-target 先验：%s", ltm_file))
         if (!file.exists(lr_file)) sc_prerequisite(sprintf("缺少 NicheNet ligand-receptor 先验：%s", lr_file))
         t_load <- Sys.time()
         ligand_target_matrix <- readRDS(ltm_file)
         lr_network <- readRDS(lr_file)
         if (!is.matrix(ligand_target_matrix) && !is.data.frame(ligand_target_matrix)) {
           ligand_target_matrix <- as.matrix(ligand_target_matrix)
         }
         lr_network <- as.data.frame(lr_network, stringsAsFactors = FALSE)
         ctx$notes <- c(ctx$notes,
                        sprintf("先验模型：%s（%d targets × %d ligands）与 %s（%d 条 L-R）",
                                basename(ltm_file), nrow(ligand_target_matrix), ncol(ligand_target_matrix),
                                basename(lr_file), nrow(lr_network)),
                        sprintf("先验读取耗时 %.1f 秒", as.numeric(Sys.time() - t_load, units = "secs")))

         b <- sc_load(ctx, max_cells = max_cells, max_genes = 15000)
         cell_type_col <- sc_pick_column(
           b$obs,
           c(as.character(params$cell_type_col %||% ""), "cell_type", "clusters", "cluster", "seurat_clusters"),
           what = "细胞类型列（params.cell_type_col）"
         )
         labels <- sc_obs_chr(b$obs, cell_type_col)
         keep <- !is.na(labels) & nzchar(labels)
         if (sum(keep) < 50) sc_prerequisite("有效细胞类型标签不足（<50 个细胞）")
         labels <- labels[keep]
         data <- b$data[keep, , drop = FALSE]
         tab <- sort(table(labels), decreasing = TRUE)
         if (length(tab) < 2) sc_prerequisite("只有 1 个细胞类型，无法做 sender/receiver 比较")
         receiver <- trimws(as.character(params$receiver %||% ""))
         sender <- trimws(as.character(params$sender %||% ""))
         auto_pick <- !nzchar(receiver) || !nzchar(sender)
         if (!nzchar(receiver)) receiver <- names(tab)[1]
         if (!nzchar(sender)) {
           candidates <- names(tab)[names(tab) != receiver]
           sender <- candidates[1]
         }
         if (!receiver %in% names(tab)) sc_prerequisite(sprintf("receiver=%s 不在细胞类型里（可选：%s）",
                                                                receiver, paste(names(tab), collapse = ", ")))
         if (!sender %in% names(tab)) sc_prerequisite(sprintf("sender=%s 不在细胞类型里（可选：%s）",
                                                              sender, paste(names(tab), collapse = ", ")))
         if (sender == receiver) sc_prerequisite("sender 与 receiver 不能是同一细胞类型")
         if (auto_pick) ctx$notes <- c(ctx$notes, "sender/receiver 为自动选择（可显式传 params.sender/params.receiver）")

         genes <- b$genes
         pct <- sc_pct_expressed(data)
         names(pct) <- genes
         background_expressed_genes <- genes[pct >= 0.10]
         genes <- genes[genes %in% colnames(data)]

         # ---- 2. receiver 差异表达基因（官方 vignette 的与其余细胞 Wilcoxon 比较）
         is_receiver <- labels == receiver
         if (sum(is_receiver) < 10 || sum(!is_receiver) < 10) {
           sc_prerequisite("receiver 或其余细胞数太少（<10），无法做差异表达")
         }
         de <- sc_wilcox_de(data, which(is_receiver), which(!is_receiver), genes, top_n = 200L)
         if (is.null(de) || nrow(de) == 0) sc_prerequisite("receiver 没有正向差异基因，无法做配体活性分析")
         geneset_oi <- de$gene

         # ---- 3. 候选配体（sender 表达 + L-R 先验 + 有靶基因先验）
         sender_pct <- sc_pct_expressed(data[labels == sender, , drop = FALSE])
         names(sender_pct) <- genes
         expressed_ligands <- names(sender_pct)[sender_pct >= 0.10]
         potential_ligands <- intersect(intersect(unique(lr_network$from), expressed_ligands),
                                       colnames(ligand_target_matrix))
         if (length(potential_ligands) < 2) {
           sc_prerequisite(sprintf("sender=%s 中表达的候选配体太少（%d 个），无法做配体活性排序",
                                   sender, length(potential_ligands)))
         }
         geneset_oi <- intersect(geneset_oi, rownames(ligand_target_matrix))
         background_expressed_genes <- intersect(background_expressed_genes, rownames(ligand_target_matrix))
         if (length(geneset_oi) < 5) sc_prerequisite("receiver 差异基因与先验靶基因重叠太少（<5）")

         # ---- 4. 官方配体活性分析
         t0 <- Sys.time()
         activities <- tryCatch(
           suppressWarnings(nichenetr::predict_ligand_activities(
             geneset = geneset_oi,
             background_expressed_genes = background_expressed_genes,
             ligand_target_matrix = ligand_target_matrix,
             potential_ligands = potential_ligands
           )),
           error = function(e) sc_prerequisite(sprintf("predict_ligand_activities 失败：%s", conditionMessage(e)))
         )
         activities <- as.data.frame(activities)
         if (!"test_ligand" %in% names(activities)) activities$test_ligand <- potential_ligands[seq_len(nrow(activities))]
         order_col <- if ("auroc" %in% names(activities)) "auroc" else names(activities)[2]
         activities <- activities[order(-activities[[order_col]]), , drop = FALSE]
         top_ligands <- utils::head(activities$test_ligand, top_n)
         elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)

         # ---- 5. 配体-受体表
         lr_sub <- lr_network[lr_network$from %in% top_ligands, , drop = FALSE]
         lr_sub <- lr_sub[!duplicated(paste(lr_sub$from, lr_sub$to)), , drop = FALSE]
         lr_df <- data.frame(
           ligand = lr_sub$from, receptor = lr_sub$to,
           database = if ("database" %in% names(lr_sub)) lr_sub$database else NA,
           sender = sender, receiver = receiver,
           stringsAsFactors = FALSE, check.names = FALSE
         )
         lr_df <- lr_df[order(match(lr_df$ligand, top_ligands)), , drop = FALSE]

         # ---- 6. 配体-靶基因连接（官方 get_weighted_ligand_target_links）
         target_links <- list()
         for (ligand in utils::head(top_ligands, min(5L, length(top_ligands)))) {
           links <- tryCatch(
             nichenetr::get_weighted_ligand_target_links(ligand = ligand, geneset = geneset_oi,
                                                         ligand_target_matrix = ligand_target_matrix,
                                                         n = 250),
             error = function(e) NULL
           )
           if (!is.null(links) && nrow(links) > 0) target_links[[ligand]] <- links
         }
         if (length(target_links) == 0) sc_prerequisite("无法计算配体-靶基因连接（先验与差异基因重叠不足）")

         target_df <- do.call(rbind, lapply(names(target_links), function(ligand) {
           links <- as.data.frame(target_links[[ligand]], stringsAsFactors = FALSE)
           # 兼容 matrix / data.frame 两种返回
           if ("target" %in% names(links) && "weight" %in% names(links)) {
             data.frame(ligand = ligand, target = links$target, weight = links$weight,
                        stringsAsFactors = FALSE)
           } else {
             data.frame(ligand = ligand, target = rownames(links), weight = as.numeric(links[[ncol(links)]]),
                        stringsAsFactors = FALSE)
           }
         }))
         target_df <- target_df[order(target_df$ligand, -target_df$weight), , drop = FALSE]
         target_df$weight <- sc_round(target_df$weight, 4)

         # ---- 7. 图：配体活性条形图
         act_show <- utils::head(activities, min(15L, nrow(activities)))
         act_show <- act_show[order(act_show[[order_col]]), , drop = FALSE]
         activity_path <- sc_png(ctx, "nichenet_ligand_activity", width = 1400, height = 850, draw = function() {
           graphics::par(bg = "white", mar = c(4.2, 17, 3, 1.2), las = 1, cex = 0.85)
           vals <- act_show[[order_col]]
           graphics::barplot(vals, horiz = TRUE, names.arg = act_show$test_ligand,
                             col = sc_palette(nrow(act_show)), border = NA,
                             xlab = order_col, main = sprintf("NicheNet ligand activity\n(sender=%s -> receiver=%s)",
                                                              sender, receiver), cex.names = 0.78)
         })

         # ---- 8. 图：配体-靶基因热图
         heat_path <- NULL
         heat_mat <- matrix(0, nrow = 0, ncol = 0)
         target_genes_top <- utils::head(unique(target_df$target), n_targets)
         ligands_show <- names(target_links)
         if (length(ligands_show) > 0 && length(target_genes_top) > 0) {
           heat_mat <- matrix(0, nrow = length(ligands_show), ncol = length(target_genes_top),
                              dimnames = list(ligands_show, target_genes_top))
           for (ligand in ligands_show) {
             links <- target_df[target_df$ligand == ligand, , drop = FALSE]
             links <- links[match(target_genes_top, links$target), , drop = FALSE]
             heat_mat[ligand, ] <- ifelse(is.na(links$weight), 0, links$weight)
           }
           heat_path <- sc_png(ctx, "nichenet_ligand_target_heatmap", width = 1500, height = 850, draw = function() {
             nr <- nrow(heat_mat); nc <- ncol(heat_mat)
             op <- graphics::par(no.readonly = TRUE)
             on.exit(try(graphics::par(op), silent = TRUE), add = TRUE)
             graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5.6, 1))
             graphics::par(bg = "white", mar = c(6.2, 8, 3.2, 1.2), cex = 0.7)
             graphics::image(t(heat_mat[nr:1, , drop = FALSE]),
                             col = sc_heat_colors(120), axes = FALSE,
                             zlim = c(0, max(heat_mat)))
             graphics::axis(1, at = seq(0, 1, length.out = nc),
                            labels = colnames(heat_mat), las = 2, cex.axis = 0.58, tick = FALSE)
             graphics::axis(2, at = seq(0, 1, length.out = nr),
                            labels = rev(rownames(heat_mat)), las = 2, cex.axis = 0.75, tick = FALSE)
             graphics::title(main = "NicheNet ligand-target links", cex.main = 0.95)
             graphics::par(mar = c(6.2, 2.2, 3.2, 3))
             graphics::image(matrix(seq(0, max(heat_mat), length.out = 100), ncol = 1),
                             col = sc_heat_colors(120), axes = FALSE, zlim = c(0, max(heat_mat)))
             graphics::axis(4, at = seq(0, 1, length.out = 5), las = 1, cex.axis = 0.7)
             graphics::mtext("weight", side = 3, line = 0.5, cex = 0.7)
           })
         }

         # ---- 9. 表
         activity_df <- activities[, intersect(c("test_ligand", "auroc", "aupr_corrected", "pearson", "pearson_ligand", "n_targets"), names(activities)), drop = FALSE]
         names(activity_df)[1] <- "ligand"
         for (col in setdiff(names(activity_df), "ligand")) {
           activity_df[[col]] <- sc_round(as.numeric(activity_df[[col]]), 4)
         }
         de_out <- data.frame(gene = de$gene, logfc = sc_round(de$logfc, 4),
                              pval = sc_round(de$pval, 8), padj = sc_round(de$padj, 8),
                              stringsAsFactors = FALSE, check.names = FALSE)

         metrics <- list(
           n_cells = nrow(data), n_cell_types = length(tab),
           n_receiver_de_genes = nrow(de), n_potential_ligands = length(potential_ligands),
           n_top_ligands = length(top_ligands), elapsed_sec = elapsed,
           top_ligand = top_ligands[1]
         )
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("sender：%s（%d 细胞）→ receiver：%s（%d 细胞）",
                                    sender, as.integer(tab[sender]), receiver, as.integer(tab[receiver])),
                            sprintf("receiver 差异基因（Wilcoxon 正向上调）：%d 个，用于配体活性排序 %d 个",
                                    nrow(de), length(geneset_oi)),
                            sprintf("候选配体：%d 个（sender 表达 + L-R 先验 + 有靶基因先验）", length(potential_ligands)),
                            sprintf("配体活性排序（AUROC）：top ligand = %s（AUROC=%.3f，耗时 %.1f 秒）",
                                    top_ligands[1], activities[[order_col]][1], elapsed),
                            sprintf("top %d 配体与靶基因连接数：%d", length(top_ligands), nrow(target_df))
                          ),
                          metrics = metrics)
         res$title <- "NicheNet 官方配体-靶基因"
         res$method <- sprintf("nichenetr %s（官方 predict_ligand_activities + ligand_target_matrix_nsga2r_final%s）",
                               sc_version("nichenetr"), suffix)
         res$tables <- list(
           sc_table("ligand_activity", "配体活性排序", activity_df,
                    sc_cols(c("ligand", "Ligand", "string"), c("auroc", "AUROC", "number"),
                            c("aupr_corrected", "AUPR (corrected)", "number"),
                            c("pearson", "Pearson", "number"), c("n_targets", "Targets", "number")),
                    max_rows = 60, note = "NicheNet predict_ligand_activities（官方实现）"),
           sc_table("ligand_receptor", "候选配体-受体对", lr_df,
                    sc_cols(c("ligand", "Ligand", "string"), c("receptor", "Receptor", "string"),
                            c("database", "Database", "string"), c("sender", "Sender", "string"),
                            c("receiver", "Receiver", "string")),
                    max_rows = 200, note = "来自用户本地 lr_network 先验"),
           sc_table("ligand_target", "配体-靶基因连接", target_df,
                    sc_cols(c("ligand", "Ligand", "string"), c("target", "Target gene", "string"),
                            c("weight", "Weight", "number")),
                    max_rows = 200, note = "nichenetr::get_weighted_ligand_target_links（top 5 配体）"),
           sc_table("receiver_de_genes", "receiver 差异基因", de_out,
                    sc_cols(c("gene", "Gene", "string"), c("logfc", "LogFC (log data)", "number"),
                            c("pval", "P value", "number"), c("padj", "P adjust", "number")),
                    max_rows = 100, note = "receiver vs 其余细胞的 Wilcoxon 秩和检验（BH 校正）")
         )
         res$figures <- list(sc_figure("ligand_activity", "配体活性排序", activity_path))
         if (!is.null(heat_path)) {
           res$figures <- c(res$figures, list(sc_figure("ligand_target_heatmap", "配体-靶基因连接热图", heat_path)))
         }
         res$notes <- c(res$notes,
                        sprintf("先验目录：%s", priors),
                        sprintf("species=%s；sender/receiver 由 %s 列确定", species, cell_type_col),
                        "receiver 差异基因用 Wilcoxon 秩和近似（与官方 Seurat FindMarkers 的 Wilcoxon 一致）",
                        "配体活性 = NicheNet 官方 predict_ligand_activities（AUROC/AUPR/Pearson）")
         res
       })
