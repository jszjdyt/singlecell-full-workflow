#!/usr/bin/env Rscript
# dsh-singlecell — CopyKAT 官方 CNV/亚克隆推断适配器（CONTRACT §10.3）
suppressPackageStartupMessages(library(copykat))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default = NULL) {
  hit <- which(args == paste0("--", name))
  if (length(hit) == 0 || hit[1] >= length(args)) return(default)
  args[hit[1] + 1]
}
bundle <- getarg("bundle"); out <- getarg("out")
params <- jsonlite::fromJSON(getarg("params", "{}"), simplifyVector = TRUE)
if (is.null(bundle) || is.null(out)) { cat("缺少 --bundle/--out\n", file = stderr()); quit(status = 1L) }
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fail <- function(msg, code = 1L) { cat(sprintf("[copykat] %s\n", msg), file = stderr()); quit(status = code) }

counts <- as(Matrix::readMM(file.path(bundle, "counts.mtx")), "CsparseMatrix")  # genes × cells
cells <- readLines(file.path(bundle, "barcodes.tsv"), warn = FALSE)
genes <- readLines(file.path(bundle, "features.tsv"), warn = FALSE)
obs <- read.csv(file.path(bundle, "obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts) != length(cells)) fail("counts.mtx 列数与 barcodes 不一致", 3L)
colnames(counts) <- cells; rownames(counts) <- genes
if (nrow(obs) != length(cells)) fail("obs.csv 行数与细胞数不一致", 3L)

n_cells <- suppressWarnings(as.integer(params$n_cells))
if (!is.finite(n_cells) || n_cells < 100) n_cells <- 2000
if (length(cells) > n_cells) {
  set.seed(0)
  keep <- sort(sample.int(length(cells), n_cells))
} else {
  keep <- seq_along(cells)
}
rawmat <- as.matrix(counts[, keep, drop = FALSE])
if (length(genes) > 3000) {
  vars <- apply(rawmat, 1, stats::var)
  rawmat <- rawmat[order(vars, decreasing = TRUE)[seq_len(3000)], , drop = FALSE]
}

cluster_col <- NULL
for (candidate in c("clusters", "cell_type", "cluster", "leiden", "seurat_clusters")) {
  if (candidate %in% colnames(obs)) { cluster_col <- candidate; break }
}
cell_group <- NULL
if (!is.null(cluster_col)) {
  cell_group <- as.character(obs[[cluster_col]])[keep]
  names(cell_group) <- cells[keep]
}

sam_name <- file.path(tempdir(), sprintf("copykat_%d", as.integer(Sys.time())))
res <- tryCatch(
  copykat(rawmat = rawmat, id.type = "S", cell.group = cell_group,
          ngene.chr = 5, win.size = 25, KS.cut = 0.1, sam.name = sam_name,
          distance = "euclidean", n.cores = 2, output.seg = FALSE),
  error = function(e) fail(sprintf("copykat 运行失败：%s", conditionMessage(e)), 1L))

pred <- res$prediction
if (is.null(pred) || nrow(pred) < 1) fail("copykat 没有返回细胞预测", 1L)
cell_names <- as.character(pred$cell.names)
labels <- as.character(pred$copykat.pred)
pred_map <- setNames(labels, cell_names)
# 肿瘤概率列（若存在）
prob_col <- grep("tumor|prob", colnames(pred), ignore.case = TRUE, value = TRUE)
prob_map <- NULL
if (length(prob_col) > 0) { prob_map <- setNames(suppressWarnings(as.numeric(pred[[prob_col[1]]])), cell_names) }

obs_out <- data.frame(cell = cells, copykat_pred = NA_character_, copykat_tumor_prob = NA_real_, check.names = FALSE)
idx <- match(cells, names(pred_map))
obs_out$copykat_pred[idx[!is.na(idx)]] <- pred_map[idx[!is.na(idx)]]
if (!is.null(prob_map)) {
  pidx <- match(cells, names(prob_map))
  obs_out$copykat_tumor_prob[pidx[!is.na(pidx)]] <- prob_map[pidx[!is.na(pidx)]]
}
write.csv(obs_out, file.path(out, "obs_extra.csv"), row.names = FALSE, quote = FALSE)

pred_tab <- table(obs_out$copykat_pred, useNA = "ifany")
table_rows <- lapply(seq_along(pred_tab), function(i) list(
  prediction = ifelse(is.na(names(pred_tab)[i]), "NA", names(pred_tab)[i]),
  n_cells = as.integer(pred_tab[[i]])))
tables <- list(list(
  id = "copykat_prediction", title = "CopyKAT cell predictions",
  columns = list(list(key = "prediction", label = "Prediction", type = "string"),
                 list(key = "n_cells", label = "Cells", type = "number")),
  rows = table_rows))

# CNA 热图（copykat 输出的 CNA 矩阵）
png(file.path(out, "copykat_heatmap.png"), width = 1400, height = 900, res = 140)
cna <- res$CNAmat
if (!is.null(cna) && nrow(cna) > 0) {
  mat <- cna[, -c(1, 2, 3), drop = FALSE]
  mat <- apply(mat, 2, as.numeric)
  lim <- max(1e-6, max(abs(mat), na.rm = TRUE))
  image(t(mat[seq_len(min(nrow(mat), 2000)), , drop = FALSE]),
        col = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(101),
        zlim = c(-lim, lim), axes = FALSE, xlab = "cells (sampled)", ylab = "genomic bins",
        main = "CopyKAT inferred CNV profile")
} else {
  plot.new(); title(main = "CopyKAT CNA profile (not available)")
}
dev.off()

result <- list(
  title = "CopyKAT 官方 CNV 推断",
  method = sprintf("CopyKAT %s 官方实现（baseline CNV + 亚克隆预测；%d 细胞子集，%d 高变基因）",
                   as.character(packageVersion("copykat")), length(keep), nrow(rawmat)),
  summary = list(
    bullets = list(
      sprintf("推断 %d 个细胞（抽样自 %d），%d 个细胞被判为 tumor/aneuploid，%d 个为 diploid/normal",
              length(keep), length(cells),
              sum(obs_out$copykat_pred == "aneuploid", na.rm = TRUE),
              sum(obs_out$copykat_pred == "diploid", na.rm = TRUE)),
      "预测结果写回 obs（copykat_pred / copykat_tumor_prob）",
      "官方方法，CNV 谱依赖基因的基因组顺序（CopyKAT 内置）"),
    metrics = list(n_cells_predicted = length(keep), n_genes_used = nrow(rawmat),
                   n_tumor = sum(obs_out$copykat_pred == "aneuploid", na.rm = TRUE),
                   n_normal = sum(obs_out$copykat_pred == "diploid", na.rm = TRUE))),
  tables = tables,
  figures = list(list(id = "cna_heatmap", title = "CopyKAT CNA heatmap", file = file.path(out, "copykat_heatmap.png"))),
  obs_file = file.path(out, "obs_extra.csv"),
  warnings = list(),
  notes = list(sprintf("抽样 %d 个细胞以控制耗时；如需全量可调 params.n_cells", length(keep))))
writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null"), file.path(out, "result.json"))
cat("[copykat] done\n", file = stderr())
