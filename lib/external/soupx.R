#!/usr/bin/env Rscript
# dsh-singlecell-workbench — SoupX 官方污染校正适配器（CONTRACT §10.3）
# 注意：只有过滤后矩阵时 SoupX 无法自动估计空液滴污染；默认用手动 rho（params.contamination）。
suppressPackageStartupMessages(library(SoupX))
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

fail <- function(msg, code = 1L) { cat(sprintf("[soupx] %s\n", msg), file = stderr()); quit(status = code) }

read_matrix <- function(path) as(Matrix::readMM(path), "CsparseMatrix")
counts <- read_matrix(file.path(bundle, "counts.mtx"))          # genes × cells
data <- tryCatch(read_matrix(file.path(bundle, "data.mtx")), error = function(e) counts)
obs <- read.csv(file.path(bundle, "obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cells <- readLines(file.path(bundle, "barcodes.tsv"), warn = FALSE)
genes <- readLines(file.path(bundle, "features.tsv"), warn = FALSE)
if (ncol(counts) != length(cells)) fail("counts.mtx 列数与 barcodes 不一致", 3L)
colnames(counts) <- cells; rownames(counts) <- genes
if (nrow(obs) != length(cells)) fail("obs.csv 行数与细胞数不一致", 3L)

cluster_col <- NULL
for (candidate in c("clusters", "cell_type", "cluster", "leiden", "seurat_clusters")) {
  if (candidate %in% colnames(obs)) { cluster_col <- candidate; break }
}
if (is.null(cluster_col)) fail("需要先聚类（obs 里没有 clusters/cell_type 列）", 3L)
clusters <- setNames(as.character(obs[[cluster_col]]), cells)
# SoupX 对全零基因/细胞敏感：先剔除，避免 adjustCounts 内部产生 NA
keep_genes <- rowSums(counts) > 0
keep_cells <- colSums(counts) > 0
if (!all(keep_cells)) {
  counts <- counts[keep_genes, keep_cells, drop = FALSE]
  cells <- cells[keep_cells]
  clusters <- clusters[keep_cells]
} else {
  counts <- counts[keep_genes, , drop = FALSE]
}

rho <- 0.1
candidate <- suppressWarnings(as.numeric(params$contamination))
if (length(candidate) == 1 && is.finite(candidate) && candidate > 0 && candidate < 1) rho <- candidate

sc <- SoupChannel(tod = counts, toc = counts)
sc <- setClusters(sc, clusters)
sc <- tryCatch(setContaminationFraction(sc, contFrac = rho, forceAccept = TRUE),
               error = function(e) fail(sprintf("setContaminationFraction 失败：%s", conditionMessage(e)), 1L))
genes <- genes[keep_genes]
out_counts <- tryCatch(adjustCounts(sc, roundToInt = TRUE), error = function(e) NULL)

if (is.null(out_counts)) {
  est <- tryCatch(calculateContaminationFraction(sc), error = function(e) NULL)
  obs_out <- data.frame(cell = cells, check.names = FALSE)
  write.csv(obs_out, file.path(out, "obs_extra.csv"), row.names = FALSE, quote = FALSE)
  result <- list(
    title = "SoupX 官方污染诊断", method = sprintf("SoupX %s 官方（adjustCounts 在过滤后矩阵上不可用，已降级为官方诊断）", as.character(packageVersion("SoupX"))),
    summary = list(bullets = list(sprintf("%d 细胞 × %d 基因，rho=%.2f 手动指定", length(cells), length(genes), rho),
                                  "adjustCounts 需要未过滤矩阵/更高覆盖度；已在真实过滤矩阵上不可用（官方内部 NA），未给校正矩阵"),
                   metrics = list(n_cells = length(cells), n_genes = length(genes), rho = rho)),
    tables = list(), figures = list(),
    warnings = list("SoupX 的 adjustCounts 在当前过滤后矩阵上不可用（官方内部 NA）；需要 raw/empty-droplet 矩阵才能做真正的去污染"),
    notes = list("建议用 raw_feature_bc_matrix 重新导入后运行"))
  writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null"), file.path(out, "result.json"))
  cat("[soupx] degraded to diagnostics\n", file = stderr()); quit(status = 0L)
}
removed <- colSums(counts) - colSums(out_counts)
total <- pmax(colSums(counts), 1)
removed_frac <- removed / total
gene_removed <- rowSums(counts) - rowSums(out_counts)
gene_frac <- gene_removed / pmax(rowSums(counts), 1)
order_gene <- order(gene_removed, decreasing = TRUE)

obs_out <- data.frame(cell = cells, soupx_removed = as.numeric(removed),
                      soupx_removed_frac = as.numeric(removed_frac), check.names = FALSE)
write.csv(obs_out, file.path(out, "obs_extra.csv"), row.names = FALSE, quote = FALSE)

table_rows <- lapply(order_gene[seq_len(min(30, length(order_gene)))], function(i) {
  list(gene = genes[i], removed_counts = as.numeric(gene_removed[i]), removed_frac = round(as.numeric(gene_frac[i]), 4))
})
table <- list(list(
  id = "soupx_genes", title = "Top genes by removed (ambient) counts",
  columns = list(list(key = "gene", label = "Gene", type = "string"),
                 list(key = "removed_counts", label = "Removed counts", type = "number"),
                 list(key = "removed_frac", label = "Removed fraction", type = "number")),
  rows = table_rows))

png(file.path(out, "soupx_removed_fraction.png"), width = 1200, height = 800, res = 140)
par(mar = c(4, 4, 3, 1))
hist(removed_frac, breaks = 40, col = "#4e79a7", border = "white",
     xlab = "Removed fraction per cell", main = sprintf("SoupX contamination correction (rho=%.2f)", rho))
abline(v = rho, col = "#e15759", lwd = 2, lty = 2)
dev.off()

result <- list(
  title = "SoupX 官方污染校正", method = sprintf("SoupX %s 官方实现（手动 rho=%.2f；过滤后矩阵无法自动估计空液滴污染）", as.character(packageVersion("SoupX")), rho),
  summary = list(
    bullets = list(
      sprintf("%d 个细胞 × %d 个基因，污染比例 rho=%.2f（手动指定）", length(cells), length(genes), rho),
      sprintf("平均每细胞移除 %.1f%% 计数；移除最多的基因：%s", mean(removed_frac) * 100, genes[order_gene[1]]),
      "校正后矩阵未写回工作台（仅统计与诊断图）；如需校正矩阵请用 export/protocol 导出"),
    metrics = list(n_cells = length(cells), n_genes = length(genes), rho = rho,
                   mean_removed_fraction = round(mean(removed_frac), 4),
                   total_removed_counts = sum(removed))),
  tables = table,
  figures = list(list(id = "removed_fraction", title = "Removed fraction per cell", file = file.path(out, "soupx_removed_fraction.png"))),
  obs_file = file.path(out, "obs_extra.csv"),
  warnings = list("SoupX 的自动污染估计需要 raw/empty-droplet 矩阵；当前用手动 rho，建议结合 AMB/背景基因检查调整"),
  notes = list("method 字段已注明为手动 rho 模式"))
writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null"), file.path(out, "result.json"))
cat("[soupx] done\n", file = stderr())
