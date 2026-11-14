#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — WGCNA 官方共表达模块（CONTRACT §10.3）
#
# 工具：WGCNA 1.74（官方 pickSoftThreshold + blockwiseModules + moduleEigengenes）
# 输入：bundle（data.mtx=log-normalized）+ params：
#   n_top_genes     高变基因数（默认 800）
#   min_module_size 最小模块大小（默认 20）
#   power           软阈值（留空=用 pickSoftThreshold 自动选）
#   network_type    unsigned | signed（默认 signed）
# 输出：模块基因表 + hub 基因表 + 模块-簇相关热图 + 软阈值诊断图
# =============================================================================

.sc_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
  if (is.na(f)) "." else dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE), mustWork = FALSE))
})
source(file.path(.sc_dir, "_common.R"))

sc_run("wgcna", "WGCNA 官方共表达模块", "WGCNA 1.74（官方 pickSoftThreshold + blockwiseModules + moduleEigengenes）",
       function(ctx) {
         params <- ctx$params
         sc_require(c("WGCNA", "dynamicTreeCut"))
         # WGCNA 的 blockwiseModules/mergeCloseModules 内部用 match.fun("cor") 找相关函数：
         # 只 requireNamespace 时 "cor" 会解析成 stats::cor（不接受 weights.x 参数）而报错，
         # 因此这里必须 attach（library）WGCNA，使其 cor 覆盖 search path。
         suppressPackageStartupMessages(library(WGCNA))
         suppressPackageStartupMessages(library(dynamicTreeCut))
         n_top <- as.integer(params$n_top_genes %||% 800)
         if (is.na(n_top) || n_top < 100) n_top <- 800L
         min_mod <- as.integer(params$min_module_size %||% 20)
         if (is.na(min_mod) || min_mod < 5) min_mod <- 20L
         network_type <- tolower(as.character(params$network_type %||% "signed"))
         if (!network_type %in% c("signed", "unsigned")) network_type <- "signed"
         user_power <- suppressWarnings(as.numeric(params$power %||% NA))

         b <- sc_load(ctx, max_cells = 15000, max_genes = 5000)
         # 选高变基因（WGCNA 需要 genes × samples）
         var <- sc_gene_variance(b$data)
         ord <- order(-var, seq_along(var))
         genes <- b$genes[sort(ord[seq_len(min(n_top, length(ord)))])]
         # 丢弃零方差/低表达基因
         pct <- sc_pct_expressed(b$data[, genes, drop = FALSE])
         # 真实数据里很多 HVG 只在少数细胞表达；阈值从 0.05 放宽到 0.01，
         # 且保证至少留 200 个基因（否则 WGCNA 无意义/易触发内部对齐问题）
         keep <- pct >= 0.01
         if (sum(keep) < 200) {
           keep <- order(-pct)[seq_len(min(500, length(pct)))]
           keep <- sort(keep)
         }
         genes <- genes[keep]
         if (length(genes) < 50) sc_prerequisite(sprintf("可用高变基因太少（%d），无法做 WGCNA", length(genes)))
         expr <- as.matrix(b$data[, genes, drop = FALSE])  # cells × genes
         datExpr <- t(expr)  # genes × cells（WGCNA 要求）
         rownames(datExpr) <- genes
         colnames(datExpr) <- b$cells

         # debug removed
         if (FALSE) cat(sprintf("[wgcna-debug] b$data=%s dim=%s genes=%d cells=%d datExpr=%s\n",
                     class(b$data)[1], paste(dim(b$data), collapse="x"), length(b$genes), length(b$cells),
                     paste(dim(datExpr), collapse="x")), file = stderr())
         # ---- 1. 软阈值
         powers <- c(seq_len(10), seq(12, 20, by = 2))
         sft <- tryCatch(
           suppressWarnings(WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
                                                     networkType = network_type,
                                                     verbose = 0)),
           error = function(e) sc_prerequisite(sprintf("pickSoftThreshold 失败：%s", conditionMessage(e)))
         )
         if (is.finite(user_power) && user_power >= 1) {
           power <- as.integer(round(user_power))
           power_source <- sprintf("用户指定 power=%d", power)
         } else {
           fit <- sft$fitIndices
           # 选择达到 scale-free R^2 >= 0.8 的最小 power；否则取 R^2 最大者
           good <- which(fit[, "SFT.R.sq"] >= 0.8)
           idx <- if (length(good) > 0) good[1] else which.max(fit[, "SFT.R.sq"])
           power <- as.integer(fit[idx, "Power"])
           power_source <- sprintf("pickSoftThreshold 自动选择 power=%d（R^2=%.2f）",
                                   power, fit[idx, "SFT.R.sq"])
         }
         ctx$notes <- c(ctx$notes, power_source,
                        sprintf("WGCNA：%d 个高变基因 × %d 个细胞，networkType=%s",
                                nrow(datExpr), ncol(datExpr), network_type))

         # ---- 2. blockwiseModules 模块识别
         t0 <- Sys.time()
         net <- tryCatch(
           suppressWarnings(WGCNA::blockwiseModules(
             datExpr, power = power, networkType = network_type,
             TOMType = network_type, minModuleSize = min_mod,
             mergeCutHeight = 0.25, numericLabels = TRUE,
             pamRespectsDendro = FALSE, maxBlockSize = max(2000, nrow(datExpr)),
             randomSeed = ctx$seed, verbose = 0
           )),
           error = function(e) sc_prerequisite(sprintf("blockwiseModules 失败：%s", conditionMessage(e)))
         )
         elapsed <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
         module_labels <- net$colors
         if (length(module_labels) != nrow(datExpr)) {
           warning(sprintf("WGCNA 模块标签长度(%d)与基因数(%d)不一致，按基因数对齐", length(module_labels), nrow(datExpr)))
           if (length(module_labels) > nrow(datExpr)) {
             module_labels <- module_labels[seq_len(nrow(datExpr))]
           } else {
             module_labels <- c(module_labels, rep(0L, nrow(datExpr) - length(module_labels)))
           }
         }
         modules <- unique(module_labels)
         modules <- modules[modules != 0]
         if (length(modules) == 0) {
           # 默认 minModuleSize 太大（合成/小数据常见无模块）：自动降阈值重试
           for (candidate in unique(c(max(5L, as.integer(min_mod) %/% 2L), 5L, 4L, 3L))) {
             if (is.na(candidate) || candidate >= as.integer(min_mod)) next
             retry <- tryCatch(
               suppressWarnings(WGCNA::blockwiseModules(
                 datExpr, power = power, networkType = network_type, TOMType = network_type,
                 minModuleSize = candidate, mergeCutHeight = 0.25, numericLabels = TRUE,
                 pamRespectsDendro = FALSE, maxBlockSize = max(2000, nrow(datExpr)),
                 randomSeed = ctx$seed, verbose = 0)),
               error = function(e) NULL
             )
             if (is.null(retry)) next
             module_labels <- retry$colors
             if (length(module_labels) != nrow(datExpr)) {
               module_labels <- if (length(module_labels) > nrow(datExpr)) module_labels[seq_len(nrow(datExpr))] else c(module_labels, rep(0L, nrow(datExpr) - length(module_labels)))
             }
             modules <- unique(module_labels)
             modules <- modules[modules != 0]
             if (length(modules) > 0) {
               ctx$notes <- c(ctx$notes, sprintf("默认 minModuleSize=%d 无模块，自动降为 %d 重试成功", as.integer(min_mod), candidate))
               break
             }
           }
         }
         if (length(modules) == 0) sc_prerequisite("WGCNA 没有识别出任何模块（可减小 min_module_size 或增加基因数）")
         n_grey <- sum(module_labels == 0)
         module_names <- paste0("M", modules)
         names(module_names) <- as.character(modules)
         module_labels_named <- ifelse(module_labels == 0, "grey", module_names[as.character(module_labels)])

         module_df <- data.frame(
           module = c(module_names, "grey"),
           n_genes = c(as.integer(table(module_labels_named)[module_names]),
                       as.integer(n_grey)),
           stringsAsFactors = FALSE, check.names = FALSE
         )
         module_df$class <- ifelse(module_df$module == "grey", "unassigned", "module")

         # ---- 3. 模块特征基因 + hub 基因
         MEs <- NULL; hub_df <- NULL; kme_df <- NULL; corr_mat <- NULL; heat_path <- NULL
         tryCatch({
           
           
         MEs <- tryCatch(WGCNA::moduleEigengenes(datExpr, colors = module_labels,
                                                 excludeGrey = TRUE)$eigengenes,
                         error = function(e) NULL)
         hub_df <- NULL
         kme_df <- NULL
         if (!is.null(MEs)) {
           kme <- tryCatch(stats::cor(t(datExpr), MEs, use = "pairwise.complete.obs"),
                           error = function(e) NULL)
           if (!is.null(kme)) {
             module_of_gene <- module_labels_named
             kme_df <- data.frame(gene = rownames(datExpr),
                                  module = module_of_gene,
                                  kME = NA_real_,
                                  stringsAsFactors = FALSE)
             for (i in seq_len(nrow(module_df))) {
               mod <- module_df$module[i]
               if (mod == "grey") next
               me_col <- sub("^ME", "", colnames(kme))
               col_idx <- which(me_col == as.character(modules)[match(mod, module_names)])
               if (length(col_idx) == 0) next
               sel <- which(module_of_gene == mod)
               kme_df$kME[sel] <- kme[sel, col_idx]
             }
             kme_df$kME <- sc_round(kme_df$kME, 4)
             hub_df <- do.call(rbind, lapply(seq_len(nrow(module_df)), function(i) {
               mod <- module_df$module[i]
               if (mod == "grey") return(NULL)
               sel <- which(kme_df$module == mod)
               if (length(sel) == 0) return(NULL)
               hub <- kme_df[sel, , drop = FALSE]
               hub <- hub[order(-hub$kME), , drop = FALSE]
               utils::head(hub, 10)
             }))
             if (!is.null(hub_df)) {
               names(hub_df)[names(hub_df) == "module"] <- "module"
               hub_df <- hub_df[, c("module", "gene", "kME")]
               hub_df <- hub_df[order(hub_df$module, -hub_df$kME), , drop = FALSE]
             }
           }
         }

         # ---- 4. 模块-簇相关（模块特征基因 vs 簇）
         cluster_col <- sc_pick_column(
           b$obs,
           c(as.character(params$cluster_col %||% ""), "clusters", "cluster", "cell_type"),
           what = "分群列", required = FALSE
         )
         corr_mat <- NULL
         heat_path <- NULL
         if (!is.null(MEs) && !is.null(cluster_col)) {
           labels <- sc_obs_chr(b$obs, cluster_col)
           keep <- !is.na(labels)
           levels <- sort(unique(labels[keep]), method = "radix")
           # 每个簇的模块均值（module eigengene 的簇均值）
           cluster_mean <- sapply(levels, function(lv) {
             idx <- which(labels == lv & keep)
             if (length(idx) < 3) return(rep(NA_real_, ncol(MEs)))
             colMeans(MEs[idx, , drop = FALSE], na.rm = TRUE)
           })
           if (is.null(dim(cluster_mean))) cluster_mean <- matrix(cluster_mean, nrow = ncol(MEs))
           rownames(cluster_mean) <- colnames(MEs)
           colnames(cluster_mean) <- levels
           # 模块分数按行 z-score 便于热图比较
           corr_mat <- t(scale(t(cluster_mean)))
           corr_mat[!is.finite(corr_mat)] <- 0
           # 模块按最大簇均值排序
           ord <- order(-apply(corr_mat, 1, max))
           corr_mat <- corr_mat[ord, , drop = FALSE]
           heat_path <- sc_png(ctx, "wgcna_module_cluster_heatmap", width = 1450, height = 850, draw = function() {
             sc_heatmap(corr_mat, title = sprintf("WGCNA module eigengene by %s (z-score)", cluster_col),
                        value_label = "z-score", row_title = "Module", col_title = cluster_col,
                        cex = 0.8)
           })
         }

         
         }, error = function(e) {
           ctx$notes <- c(ctx$notes, sprintf("模块特征基因/热图步骤跳过：%s", conditionMessage(e)))
         })

         # ---- 5. 软阈值诊断图
         fit <- sft$fitIndices
         diag_path <- sc_png(ctx, "wgcna_soft_threshold", width = 1400, height = 700, draw = function() {
           graphics::par(bg = "white", mfrow = c(1, 2), mar = c(4.2, 4.4, 2.8, 1.0), las = 1, cex = 0.8)
           plot(fit[, "Power"], -sign(fit[, "slope"]) * fit[, "SFT.R.sq"],
                pch = 16, col = "#4E79A7", xlab = "Soft threshold power",
                ylab = "Scale-free topology R^2 (signed)",
                main = "WGCNA soft threshold")
           abline(h = 0.8, col = "#E15759", lty = 2)
           abline(v = power, col = "#59A14F", lty = 3, lwd = 2)
           plot(fit[, "Power"], fit[, "mean.k."], pch = 16, col = "#4E79A7",
                xlab = "Soft threshold power", ylab = "Mean connectivity",
                main = "Mean connectivity")
         })

         # ---- 6. 表
         all_module_df <- data.frame(
           module = module_labels_named,
           gene = rownames(datExpr),
           stringsAsFactors = FALSE
         )
         all_module_df <- all_module_df[order(all_module_df$module, all_module_df$gene), , drop = FALSE]
         if (!is.null(kme_df)) {
           all_module_df <- merge(all_module_df, kme_df[, c("gene", "kME")], by = "gene", all.x = TRUE)
           all_module_df <- all_module_df[order(all_module_df$module, -all_module_df$kME), , drop = FALSE]
           all_module_df$kME <- sc_round(all_module_df$kME, 4)
         }

         metrics <- list(
           n_genes = nrow(datExpr), n_cells = ncol(datExpr),
           n_modules = length(modules), n_grey = n_grey,
           power = power, min_module_size = min_mod,
           network_type = network_type, elapsed_sec = elapsed
         )
         top_module <- module_df$module[which.max(ifelse(module_df$module == "grey", 0, module_df$n_genes))]
         res <- sc_result(ctx,
                          bullets = c(
                            sprintf("WGCNA 输入：%d 个高变基因 × %d 个细胞（networkType=%s）",
                                    nrow(datExpr), ncol(datExpr), network_type),
                            sprintf("%s（%d 棵树的 blockwiseModules，耗时 %.1f 秒）", power_source, 1L, elapsed),
                            sprintf("识别模块：%d 个（未分配 grey：%d 个基因）", length(modules), n_grey),
                            sprintf("最大模块：%s（%d 个基因）", top_module,
                                    module_df$n_genes[match(top_module, module_df$module)]),
                            if (!is.null(heat_path)) sprintf("模块-簇热图：%s 各簇模块特征基因（z-score）", cluster_col) else "无分群列，未生成模块-簇热图"
                          ),
                          metrics = metrics)
         res$title <- "WGCNA 官方共表达模块"
         res$method <- sprintf("WGCNA %s（官方 pickSoftThreshold + blockwiseModules + moduleEigengenes）",
                               sc_version("WGCNA"))
         res$tables <- list(
           sc_table("modules", "模块汇总", module_df,
                    sc_cols(c("module", "Module", "string"), c("n_genes", "Genes", "number"),
                            c("class", "Class", "string")),
                    max_rows = 100, note = sprintf("minModuleSize=%d，mergeCutHeight=0.25", min_mod)),
           sc_table("module_genes", "模块基因", all_module_df,
                    sc_cols(c("module", "Module", "string"), c("gene", "Gene", "string"),
                            c("kME", "kME", "number")),
                    max_rows = 300, note = "kME = 基因与模块特征基因的相关性"),
           sc_table("hub_genes", "Hub 基因（每模块 top 10 kME）", hub_df %||% data.frame(module = character(), gene = character(), kME = numeric()),
                    sc_cols(c("module", "Module", "string"), c("gene", "Gene", "string"),
                            c("kME", "kME", "number")),
                    max_rows = 200, note = "按 kME 排序")
         )
         res$figures <- list(sc_figure("soft_threshold", "WGCNA 软阈值诊断", diag_path))
         if (!is.null(heat_path)) {
           res$figures <- c(res$figures, list(sc_figure("module_cluster_heatmap", "模块-簇热图", heat_path)))
         }
         res$notes <- c(res$notes,
                        "module = grey 表示未分配到任何模块",
                        "blockwiseModules 使用固定 randomSeed，结果可复现")
         res
       })
