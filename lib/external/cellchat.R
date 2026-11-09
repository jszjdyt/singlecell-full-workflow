#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell — CellChat 官方细胞通讯（CONTRACT §10.3）
#
# 工具：CellChat 1.6.1（官方 createCellChat → setDB(CellChatDB.human) →
#       identifyOverExpressedGenes/LigandReceptor → computeCommunProb →
#       filterCommunication → computeCommunProbPathway → aggregateNet）
# 输入：bundle（data.mtx = log-normalized）+ params：
#   min_cells  每个细胞类型最小细胞数（默认 10）
#   thresh     通讯显著性 p 阈值（默认 0.05）
#   labels     标签列（默认自动 cell_type > clusters）
#   database   CellChatDB 类别，逗号分隔（默认 "Secreted Signaling,ECM-Receptor,Cell-Cell Contact"）
#   nboot      置换次数（默认 100，与官方一致；数据大可降到 20-50）
# 输出：interactions / 通路 / 细胞类型通讯强度表 + 网络图 + 泡图
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

sc_run("cellchat", "CellChat 官方细胞通讯",
       "CellChat 1.6.1（官方 createCellChat + CellChatDB.human + computeCommunProb + aggregateNet）",
       function(ctx) {
         params <- ctx$params
         sc_require("CellChat")
         # CellChat 会向 stdout 打印问候语/进度条：重定向到 stderr，保持 stdout 干净（§10.3）
         sink(stderr(), type = "output")
         on.exit(try(sink(type = "output"), silent = TRUE), add = TRUE)
         min_cells <- as.integer(params$min_cells %||% 10)
         if (is.na(min_cells) || min_cells < 3) min_cells <- 3L
         thr <- as.numeric(params$thresh %||% 0.05)
         if (!is.finite(thr) || thr <= 0 || thr > 1) thr <- 0.05
         nboot <- as.integer(params$nboot %||% 100)
         if (is.na(nboot) || nboot < 10) nboot <- 100L
         db_cats <- as.character(params$database %||% "")
         if (!nzchar(db_cats)) {
           db_cats <- "Secreted Signaling,ECM-Receptor,Cell-Cell Contact"
         }
         db_cats <- trimws(unlist(strsplit(db_cats, ",")))

         b <- sc_load(ctx, max_cells = 4000, max_genes = 15000)
         labels_col <- sc_pick_column(
           b$obs,
           c(as.character(params$labels %||% ""), "cell_type", "clusters", "cluster", "seurat_clusters"),
           what = "细胞类型/聚类列（params.labels）"
         )
         labels <- sc_obs_chr(b$obs, labels_col)
         keep <- !is.na(labels) & nzchar(labels)
         if (sum(keep) < 20) sc_prerequisite("有效标签细胞太少（<20），无法做细胞通讯")
         labels <- labels[keep]
         data <- b$data[keep, , drop = FALSE]
         # 丢弃细胞数过少的分组（CellChat 的 identifyOverExpressedGenes 需要每组有足够细胞）
         tab <- table(labels)
         keep_levels <- names(tab)[tab >= max(3L, min_cells)]
         dropped <- names(tab)[tab < max(3L, min_cells)]
         if (length(dropped) > 0) {
           ctx$warnings <- c(ctx$warnings, sprintf("细胞类型 %s 细胞数过少，已剔除",
                                                   paste(sprintf("%s(%d)", dropped, as.integer(tab[dropped])), collapse = ", ")))
         }
         keep2 <- labels %in% keep_levels
         labels <- labels[keep2]
         data <- data[keep2, , drop = FALSE]
         if (length(unique(labels)) < 2) sc_prerequisite("过滤后只剩 1 个细胞类型，无法做通讯分析")

         # CellChat 需要 genes × cells 的标准化表达矩阵
         data_gc <- methods::as(Matrix::t(data), "CsparseMatrix")
         rownames(data_gc) <- b$genes
         colnames(data_gc) <- b$cells[keep][keep2]
         meta <- data.frame(labels = factor(labels), row.names = colnames(data_gc))

         # ---- 1. 官方流程
         t0 <- Sys.time()
         cellchat <- tryCatch(
           CellChat::createCellChat(object = data_gc, meta = meta, group.by = "labels"),
           error = function(e) sc_prerequisite(sprintf("createCellChat 失败：%s", conditionMessage(e)))
         )
         db_env <- new.env(parent = emptyenv())
         utils::data("CellChatDB.human", package = "CellChat", envir = db_env)
         db <- get("CellChatDB.human", envir = db_env)
         categories <- names(table(db$interaction$annotation))
         use_cats <- intersect(db_cats, categories)
         if (length(use_cats) == 0) use_cats <- categories
         db_sub <- CellChat::subsetDB(db, search = use_cats, key = "annotation")
         cellchat@DB <- db_sub
         ctx$notes <- c(ctx$notes, sprintf("CellChatDB.human 类别：%s（共 %d 条 interaction）",
                                           paste(use_cats, collapse = ", "),
                                           nrow(db_sub$interaction)))

         cellchat <- tryCatch(CellChat::subsetData(cellchat),
                              error = function(e) sc_prerequisite(sprintf("subsetData 失败：%s", conditionMessage(e))))
         cellchat <- tryCatch(
           CellChat::identifyOverExpressedGenes(cellchat, thresh.p = 0.05),
           error = function(e) sc_prerequisite(sprintf("identifyOverExpressedGenes 失败（缺少 presto 时可能较慢）：%s",
                                                       conditionMessage(e)))
         )
         cellchat <- tryCatch(
           CellChat::identifyOverExpressedInteractions(cellchat),
           error = function(e) sc_prerequisite(sprintf("identifyOverExpressedInteractions 失败：%s",
                                                       conditionMessage(e)))
         )
         cellchat <- tryCatch(
           suppressWarnings(CellChat::computeCommunProb(cellchat, type = "triMean",
                                                        raw.use = TRUE, population.size = FALSE,
                                                        nboot = nboot, seed.use = ctx$seed)),
           error = function(e) sc_prerequisite(sprintf("computeCommunProb 失败：%s", conditionMessage(e)))
         )
         cellchat <- tryCatch(
           CellChat::filterCommunication(cellchat, min.cells = min_cells),
           error = function(e) sc_prerequisite(sprintf("filterCommunication 失败：%s", conditionMessage(e)))
         )
         cellchat <- tryCatch(
           CellChat::computeCommunProbPathway(cellchat),
           error = function(e) sc_prerequisite(sprintf("computeCommunProbPathway 失败：%s", conditionMessage(e)))
         )
         cellchat <- tryCatch(
           CellChat::aggregateNet(cellchat),
           error = function(e) sc_prerequisite(sprintf("aggregateNet 失败：%s", conditionMessage(e)))
         )
         elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)

         # ---- 2. 通讯表
         commun <- as.data.frame(CellChat::subsetCommunication(cellchat))
         if (nrow(commun) == 0) sc_prerequisite("CellChat 未预测到任何显著通讯（可尝试降低 min_cells / 换标签列）")
         # 兼容不同版本列名
         pick_col <- function(df, candidates, default = NA) {
           for (candidate in candidates) if (candidate %in% names(df)) return(df[[candidate]])
           rep(default, nrow(df))
         }
         commun_df <- data.frame(
           source = as.character(pick_col(commun, c("source"))),
           target = as.character(pick_col(commun, c("target"))),
           ligand = as.character(pick_col(commun, c("ligand"))),
           receptor = as.character(pick_col(commun, c("receptor"))),
           interaction_name = as.character(pick_col(commun, c("interaction_name"))),
           pathway = as.character(pick_col(commun, c("pathway_name", "pathway"))),
           prob = as.numeric(pick_col(commun, c("prob", "probability"), NA_real_)),
           pval = as.numeric(pick_col(commun, c("pval", "pvalue"), NA_real_)),
           annotation = as.character(pick_col(commun, c("annotation"), NA)),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         commun_df <- commun_df[order(commun_df$pval, -commun_df$prob), , drop = FALSE]
         for (col in c("prob", "pval")) commun_df[[col]] <- sc_round(commun_df[[col]], 6)
         sig <- commun_df[is.finite(commun_df$pval) & commun_df$pval <= thr, , drop = FALSE]

         # ---- 3. 通路汇总
         pathway_df <- NULL
         if (nrow(sig) > 0) {
           agg <- stats::aggregate(prob ~ pathway, data = sig, FUN = sum)
           cnt <- stats::aggregate(prob ~ pathway, data = sig, FUN = length)
           names(cnt)[2] <- "n_interactions"
           pathway_df <- merge(agg, cnt, by = "pathway")
           names(pathway_df)[2] <- "total_prob"
           pathway_df <- pathway_df[order(-pathway_df$total_prob), , drop = FALSE]
           names(pathway_df)[2] <- "total_prob"
           pathway_df$total_prob <- sc_round(pathway_df$total_prob, 4)
         }

         # ---- 4. 细胞类型通讯强度（aggregateNet 的 weight 矩阵）
         weight <- cellchat@net$weight
         ct_df <- NULL
         if (!is.null(weight) && nrow(weight) > 0) {
           ct_df <- data.frame(
             source = rep(rownames(weight), times = ncol(weight)),
             target = rep(colnames(weight), each = nrow(weight)),
             weight = as.numeric(weight),
             stringsAsFactors = FALSE, check.names = FALSE
           )
           ct_df <- ct_df[order(-ct_df$weight), , drop = FALSE]
           ct_df$weight <- sc_round(ct_df$weight, 5)
           ct_df <- ct_df[ct_df$weight > 0, , drop = FALSE]
         }

         # ---- 5. obs_extra：每个细胞的 outgoing/incoming 通讯强度（按标签汇总）
         group_size <- as.numeric(table(labels))
         names(group_size) <- names(table(labels))
         if (!is.null(weight) && nrow(weight) > 0) {
           out_strength <- rowSums(weight)
           in_strength <- colSums(weight)
           obs_cols <- list(
             cellchat_out_strength = as.numeric(out_strength[labels]),
             cellchat_in_strength = as.numeric(in_strength[labels])
           )
           sc_write_obs(ctx, obs_cols, colnames(data_gc))
         }

         # ---- 6. 图：circle 网络 + bubble
         figures <- list()
         circle_path <- tryCatch({
           sc_png(ctx, "cellchat_network", width = 1100, height = 1100, draw = function() {
             graphics::par(bg = "white", mar = c(1.2, 1.2, 3.2, 1.2))
             CellChat::netVisual_circle(weight, vertex.weight = group_size,
                                        weight.scale = TRUE, label.edge = FALSE,
                                        title.name = sprintf("CellChat interaction strength (p<=%.2f)", thr))
           })
         }, error = function(e) {
           ctx$warnings <<- c(ctx$warnings, sprintf("netVisual_circle 失败：%s", conditionMessage(e)))
           NULL
         })
         if (!is.null(circle_path)) {
           figures <- c(figures, list(sc_figure("cellchat_network", "CellChat 通讯网络（circle）", circle_path)))
         }
         bubble_path <- tryCatch({
           bubble <- suppressWarnings(CellChat::netVisual_bubble(cellchat, sources.use = NULL,
                                                                targets.use = NULL, remove.isolate = TRUE))
           sc_png(ctx, "cellchat_bubble", width = 1500, height = 1000, draw = function() print(bubble))
         }, error = function(e) {
           ctx$warnings <<- c(ctx$warnings, sprintf("netVisual_bubble 失败：%s", conditionMessage(e)))
           NULL
         })
         if (!is.null(bubble_path)) {
           figures <- c(figures, list(sc_figure("cellchat_bubble", "显著 L-R 对泡图", bubble_path)))
         }
         # 兜底图：source-target 通讯强度热图
         if (is.null(weight) || nrow(weight) < 2) {
           sc_prerequisite("aggregateNet 的通讯强度矩阵为空，无法出图")
         }
         if (length(figures) == 0) {
           heat <- sc_png(ctx, "cellchat_strength_heatmap", width = 1300, height = 900, draw = function() {
             sc_heatmap(t(weight[nrow(weight):1, , drop = FALSE]), title = "CellChat communication strength",
                        row_labels = rownames(weight), col_labels = colnames(weight),
                        value_label = "weight", row_title = "Source", col_title = "Target",
                        cex = 0.8)
           })
           figures <- c(figures, list(sc_figure("cellchat_strength_heatmap", "通讯强度热图", heat)))
         }

         metrics <- list(
           n_cell_types = nrow(weight), n_significant_interactions = nrow(sig),
           n_interactions = nrow(commun_df), n_pathways = if (is.null(pathway_df)) 0L else nrow(pathway_df),
           total_strength = round(sum(weight), 5), nboot = nboot, elapsed_sec = elapsed
         )
         top_lr <- if (nrow(sig) > 0) sprintf("%s → %s（%s-%s, p=%.2g）",
                                              sig$source[1], sig$target[1], sig$ligand[1],
                                              sig$receptor[1], sig$pval[1]) else "无"
         top_pw <- if (is.null(pathway_df)) "无" else pathway_df$pathway[1]
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("CellChatDB.human：%s（%d 条 interaction）", paste(use_cats, collapse = " / "),
                                    nrow(db_sub$interaction)),
                            sprintf("细胞类型：%d 个（标签列 %s，min_cells=%d，nboot=%d）",
                                    nrow(weight), labels_col, min_cells, nboot),
                            sprintf("显著通讯 interaction：%d / %d 条（p<=%.2f）",
                                    nrow(sig), nrow(commun_df), thr),
                            sprintf("显著通讯通路：%d 条；总通讯强度 %.3f",
                                    metrics$n_pathways, metrics$total_strength),
                            sprintf("最强 L-R：%s", top_lr),
                            sprintf("最强通路：%s（耗时 %.1f 秒）", top_pw, elapsed)
                          ),
                          metrics = metrics)
         res$title <- "CellChat 官方细胞通讯"
         res$method <- sprintf("CellChat %s（官方 computeCommunProb(triMean, nboot=%d) + aggregateNet + CellChatDB.human）",
                               sc_version("CellChat"), nboot)
         res$tables <- list()
         if (nrow(sig) > 0) {
           res$tables <- c(res$tables, list(sc_table(
             "interactions", "显著细胞通讯（L-R 对）", sig,
             sc_cols(c("source", "Source", "string"), c("target", "Target", "string"),
                     c("ligand", "Ligand", "string"), c("receptor", "Receptor", "string"),
                     c("interaction_name", "Interaction", "string"), c("pathway", "Pathway", "string"),
                     c("prob", "Prob", "number"), c("pval", "P value", "number"),
                     c("annotation", "Annotation", "string")),
             max_rows = 200, note = sprintf("subsetCommunication 的显著结果（p<=%.2f）", thr)
           )))
         }
         if (!is.null(pathway_df)) {
           res$tables <- c(res$tables, list(sc_table(
             "pathways", "通讯通路汇总", pathway_df,
             sc_cols(c("pathway", "Pathway", "string"), c("total_prob", "Total prob", "number"),
                     c("n_interactions", "Interactions", "number")),
             max_rows = 100, note = "按 pathway 汇总显著通讯概率"
           )))
         }
         if (!is.null(ct_df)) {
           res$tables <- c(res$tables, list(sc_table(
             "celltype_strength", "细胞类型通讯强度矩阵", ct_df,
             sc_cols(c("source", "Source", "string"), c("target", "Target", "string"),
                     c("weight", "Weight", "number")),
             max_rows = 300, note = "aggregateNet 的 net$weight（行=发送方，列=接收方）"
           )))
         }
         res$figures <- figures
         res$notes <- c(res$notes,
                        sprintf("标签列：%s；显著阈值 p<=%.2f", labels_col, thr),
                        "computeCommunProb 使用 triMean + 置换检验（nboot 次）；seed 固定",
                        "obs_extra：cellchat_out_strength / cellchat_in_strength（按细胞标签映射）")
         res
       })
