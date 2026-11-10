#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell 官方工具适配器：Seurat 6.0 官方降维 / 聚类 / UMAP
# -----------------------------------------------------------------------------
# CONTRACT.md §10.3 / §10.5
#
# 用法（stdout 必须保持干净，日志一律走 stderr）：
#   Rscript --vanilla seurat_pipeline.R --bundle <dir> --out <dir> --params '<json>'
#
# params（JSON 对象）：
#   step        pca | neighbors | cluster | embed | full     （默认 full）
#   n_pcs       主成分数                                      （默认 30）
#   n_neighbors FindNeighbors / UMAP 的近邻数                  （默认 20）
#   resolution  FindClusters 分辨率                            （默认 0.8）
#   n_top_genes FindVariableFeatures 高变基因数                （默认 2000）
#   dims        使用前多少个 PC（数字=1:20，也接受 [1,20] / "1:20"）（默认 20）
#   umap        是否计算 UMAP（full/embed 步骤）               （默认 true）
#   seed        随机种子                                      （默认 0）
#   max_cells   细胞上限（超限按 seed 下采样）                 （默认 30000）
#   max_genes   基因上限（超限取高方差基因）                   （默认 8000）
#
# 输入 bundle（CONTRACT §10.2，由 lib/sc_external.py::export_bundle 导出）：
#   counts.mtx(genes×cells) / data.mtx(genes×cells) / features.tsv / barcodes.tsv
#   obs.csv / meta.json / pca.csv(可选) / embed.csv(可选)
#
# 输出（<out> 目录）：
#   result.json      契约 §10.3 结构
#   obs_extra.csv    第一列 cell + seurat_clusters 等
#   pca.csv          cell,PC1..PCn（当本步骤有 PCA 时）
#   embed.csv        cell,dim1,dim2（当本步骤有 UMAP 时）
#   *.png            碎石图 / UMAP 图（白底英文）
#
# 退出码：0 成功；3 缺前置/数据不满足（stderr 中文原因）；其它非 0 = 失败
# =============================================================================

options(warn = 1)  # warnings 立即输出到 stderr，便于排查

# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

.log <- function(...) {
  message(sprintf(...))
}

.prereq <- function(fmt, ...) {
  condition <- structure(
    class = c("dsh_prereq", "error", "condition"),
    list(message = sprintf(fmt, ...), call = NULL)
  )
  stop(condition)
}

.bad <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

.warn_add <- function(warnings, fmt, ...) {
  msg <- sprintf(fmt, ...)
  .log("WARN: %s", msg)
  c(warnings, msg)
}

.is_true <- function(value, default = TRUE) {
  if (is.null(value) || length(value) == 0) return(default)
  if (is.logical(value)) return(isTRUE(value[1]))
  text <- tolower(as.character(value[1]))
  if (text %in% c("true", "t", "yes", "y", "1", "on")) return(TRUE)
  if (text %in% c("false", "f", "no", "n", "0", "off")) return(FALSE)
  default
}

.as_int <- function(value, default) {
  if (is.null(value) || length(value) == 0) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(as.numeric(value[1])))
  if (is.na(parsed)) as.integer(default) else parsed
}

.as_num <- function(value, default) {
  if (is.null(value) || length(value) == 0) return(as.numeric(default))
  parsed <- suppressWarnings(as.numeric(value[1]))
  if (is.na(parsed)) as.numeric(default) else parsed
}

# ---------------------------------------------------------------------------
# 参数解析
# ---------------------------------------------------------------------------
.parse_cli <- function(argv) {
  out <- list()
  index <- 1L
  while (index <= length(argv)) {
    item <- argv[[index]]
    if (grepl("^--[^=]+=", item)) {
      key <- sub("^--", "", sub("=.*$", "", item))
      out[[key]] <- sub("^[^=]+=", "", item)
      index <- index + 1L
    } else if (grepl("^--", item)) {
      key <- sub("^--", "", item)
      if (index == length(argv)) .bad("命令行参数 --%s 缺少取值", key)
      out[[key]] <- argv[[index + 1L]]
      index <- index + 2L
    } else {
      index <- index + 1L
    }
  }
  out
}

.parse_params <- function(raw) {
  defaults <- list(
    step = "full", n_pcs = 30L, n_neighbors = 30L, resolution = NULL,
    n_top_genes = 3000L, dims = 30L, umap = TRUE, seed = 0L,
    integration = "harmony", use_harmony = TRUE, batch_key = "",
    metric = "euclidean", reduction = "",
    max_cells = 200000L, max_genes = 10000L
  )
  if (is.null(raw) || !nzchar(trimws(raw))) return(defaults)
  text <- trimws(raw)
  parsed <- NULL
  if (grepl("^[\\[{]", text)) {
    parsed <- tryCatch(jsonlite::fromJSON(text, simplifyVector = FALSE),
                       error = function(error) .bad("--params 不是合法 JSON：%s", error$message))
  } else if (file.exists(text)) {
    parsed <- tryCatch(jsonlite::fromJSON(readLines(text, warn = FALSE), simplifyVector = FALSE),
                       error = function(error) .bad("--params 文件不是合法 JSON：%s", error$message))
  } else {
    .bad("--params 既不是 JSON 文本也不是已存在文件：%s", text)
  }
  if (!is.list(parsed)) .bad("--params 必须是 JSON 对象")
  merged <- defaults
  for (name in names(parsed)) merged[[name]] <- parsed[[name]]
  merged
}

# dims：接受 20 / [1,20] / "1:20" / [1,2,3]，返回请求的 PC 序号
.requested_dims <- function(value) {
  dims <- value %||% 20L
  if (is.list(dims)) dims <- unlist(dims, use.names = FALSE)
  if (is.character(dims) && length(dims) == 1 && grepl(":", dims)) {
    parts <- strsplit(dims, ":", fixed = TRUE)[[1]]
    dims <- seq(as.integer(trimws(parts[1])), as.integer(trimws(parts[2])))
  }
  dims <- suppressWarnings(as.integer(as.numeric(dims)))
  dims <- dims[!is.na(dims)]
  if (length(dims) == 1) {
    return(seq_len(dims))
  }
  if (length(dims) == 2 && dims[2] > dims[1]) {
    return(seq(dims[1], dims[2]))
  }
  dims
}

.resolve_dims <- function(value, max_dim) {
  result <- .requested_dims(value)
  result <- result[result >= 1L & result <= max_dim]
  if (length(result) < 2L) {
    .prereq("可用主成分不足 2 个（请求 dims=%s，上限 %d），无法构图 / 计算 UMAP；请增大 n_pcs 或减少维度",
            paste(value, collapse = ","), as.integer(max_dim))
  }
  result
}

# ---------------------------------------------------------------------------
# 独立 R 库（Seurat 6.0）
# ---------------------------------------------------------------------------
.resolve_seurat_lib <- function() {
  candidates <- c(
    Sys.getenv("DS_SC_SEURAT_LIB", ""),
    Sys.getenv("DS_SC_R_LIB", ""),
    path.expand("~/.dsh/sc-plugin-rlib")
  )
  config <- Sys.getenv("DS_SC_CONFIG", "")
  if (nzchar(config)) {
    parsed <- tryCatch(jsonlite::fromJSON(config, simplifyVector = FALSE), error = function(e) NULL)
    if (is.list(parsed)) {
      candidates <- c(candidates, parsed$seuratLib, parsed$rLibs)
    }
  }
  for (item in candidates) {
    if (is.null(item) || length(item) == 0 || !nzchar(as.character(item)[1])) next
    candidate <- path.expand(as.character(item)[1])
    if (dir.exists(candidate)) return(candidate)
  }
  path.expand("~/.dsh/sc-plugin-rlib")
}

.setup_seurat <- function() {
  seurat_lib <- .resolve_seurat_lib()
  if (dir.exists(seurat_lib)) {
    .libPaths(unique(c(seurat_lib, .libPaths())))
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    .prereq("独立 R 库 %s 里没有安装 Seurat；请先在独立库安装 Seurat >= 6.0（不得覆盖系统库）", seurat_lib)
  }
  version <- tryCatch(as.character(utils::packageVersion("Seurat")), error = function(e) "0")
  if (utils::compareVersion(version, "6.0") < 0) {
    .prereq("降维要求 Seurat >= 6.0，独立库 %s 当前安装的是 %s（系统库另可能有 5.x，本适配器不会使用）",
            seurat_lib, version)
  }
  # 真实加载路径，用于 metrics 校验"用的是独立库而不是系统库"
  loaded_path <- tryCatch(dirname(system.file(package = "Seurat")), error = function(e) NA_character_)
  list(lib = seurat_lib, version = version, path = loaded_path)
}

# ---------------------------------------------------------------------------
# bundle 读取
# ---------------------------------------------------------------------------
.first_existing <- function(paths) {
  for (item in paths) if (file.exists(item)) return(item)
  NULL
}

.read_lines_nonempty <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("\r$", "", lines)
  lines[nzchar(lines)]
}

.read_features <- function(path) {
  lines <- .read_lines_nonempty(path)
  if (!length(lines)) .prereq("features.tsv 为空：%s", path)
  columns <- strsplit(lines, "\t", fixed = TRUE)
  width <- max(vapply(columns, length, integer(1)))
  if (width >= 3) {
    # 10x 风格：gene_id \t gene_symbol \t gene_type
    genes <- vapply(columns, function(parts) if (length(parts) >= 2) parts[2] else parts[1], character(1))
  } else if (width == 2) {
    first <- vapply(columns, function(parts) parts[1], character(1))
    second <- vapply(columns, function(parts) if (length(parts) >= 2) parts[2] else "", character(1))
    use_second <- all(nzchar(second)) && mean(grepl("^ENS", first)) > 0.5
    genes <- if (use_second) second else first
  } else {
    genes <- lines
  }
  genes <- trimws(genes)
  if (any(!nzchar(genes))) .prereq("features.tsv 存在空基因名：%s", path)
  genes
}

.read_obs_csv <- function(path, barcodes) {
  if (!file.exists(path)) return(NULL)
  frame <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    colClasses = "character", quote = "", comment.char = ""),
    error = function(error) NULL
  )
  if (is.null(frame) || ncol(frame) < 1) return(NULL)
  cells <- as.character(frame[[1]])
  frame <- frame[, -1, drop = FALSE]
  rownames(frame) <- cells
  if (nrow(frame) != length(barcodes) || !setequal(cells, barcodes)) {
    .log("WARN: obs.csv 细胞名与 barcodes.tsv 不一致，按 barcodes.tsv 顺序对齐")
    frame <- frame[match(barcodes, cells), , drop = FALSE]
    rownames(frame) <- barcodes
  } else {
    frame <- frame[match(barcodes, cells), , drop = FALSE]
    rownames(frame) <- barcodes
  }
  frame
}

.read_csv_matrix <- function(path, barcodes, prefix) {
  if (!file.exists(path)) return(NULL)
  frame <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(error) NULL
  )
  if (is.null(frame) || ncol(frame) < 2) return(NULL)
  cells <- as.character(frame[[1]])
  matrix <- suppressWarnings(as.matrix(frame[, -1, drop = FALSE]))
  storage.mode(matrix) <- "double"
  if (nrow(matrix) != length(barcodes) || !setequal(cells, barcodes)) {
    .log("WARN: %s 细胞名与 barcodes.tsv 不一致，忽略该文件", basename(path))
    return(NULL)
  }
  matrix <- matrix[match(barcodes, cells), , drop = FALSE]
  rownames(matrix) <- barcodes
  colnames(matrix) <- paste0(prefix, seq_len(ncol(matrix)))
  if (any(!is.finite(matrix))) {
    .log("WARN: %s 含非有限值，忽略该文件", basename(path))
    return(NULL)
  }
  matrix
}

.read_bundle <- function(bundle_dir, warnings) {
  if (is.null(bundle_dir) || !nzchar(bundle_dir)) .bad("缺少 --bundle 参数")
  if (!dir.exists(bundle_dir)) .prereq("bundle 目录不存在：%s", bundle_dir)

  meta <- list()
  meta_path <- file.path(bundle_dir, "meta.json")
  if (file.exists(meta_path)) {
    meta <- tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = TRUE), error = function(error) {
      .log("WARN: meta.json 解析失败（%s），忽略", error$message)
      list()
    })
    if (!is.list(meta)) meta <- list()
  }

  mtx_path <- .first_existing(c(
    file.path(bundle_dir, "counts.mtx"), file.path(bundle_dir, "counts.mtx.gz"),
    file.path(bundle_dir, "matrix.mtx"), file.path(bundle_dir, "matrix.mtx.gz"),
    file.path(bundle_dir, "data.mtx"), file.path(bundle_dir, "data.mtx.gz")
  ))
  if (is.null(mtx_path)) {
    .prereq("bundle 缺少表达矩阵（counts.mtx / data.mtx）：%s", bundle_dir)
  }
  is_counts_file <- grepl("^counts\\.|^matrix\\.", basename(mtx_path))

  genes_path <- .first_existing(c(file.path(bundle_dir, "features.tsv"),
                                  file.path(bundle_dir, "features.tsv.gz"),
                                  file.path(bundle_dir, "genes.tsv")))
  barcodes_path <- .first_existing(c(file.path(bundle_dir, "barcodes.tsv"),
                                     file.path(bundle_dir, "barcodes.tsv.gz")))
  if (is.null(genes_path) || is.null(barcodes_path)) {
    .prereq("bundle 缺少 features.tsv / barcodes.tsv：%s", bundle_dir)
  }
  genes <- .read_features(genes_path)
  barcodes <- .read_lines_nonempty(barcodes_path)
  if (!length(barcodes)) .prereq("barcodes.tsv 为空：%s", barcodes_path)
  if (anyDuplicated(barcodes)) {
    barcodes <- make.unique(barcodes, sep = "-")
    warnings <- .warn_add(warnings, "barcodes.tsv 存在重复细胞名，已自动去重")
  }
  n_cells <- length(barcodes)
  n_genes <- length(genes)
  if (n_cells < 10L) .prereq("细胞数太少（%d），无法完成 PCA / 聚类（至少需要 10 个细胞）", n_cells)
  if (n_genes < 10L) .prereq("基因数太少（%d），无法完成 PCA / 聚类（至少需要 10 个基因）", n_genes)

  .log("读取矩阵 %s …", basename(mtx_path))
  matrix <- tryCatch(Matrix::readMM(mtx_path), error = function(error) {
    .prereq("MatrixMarket 矩阵读取失败（%s）：%s", basename(mtx_path), error$message)
  })
  matrix <- methods::as(matrix, "CsparseMatrix")
  # 契约：genes × cells；容忍 cells × genes 的转置矩阵
  if (nrow(matrix) == n_genes && ncol(matrix) == n_cells) {
    matrix <- Matrix::t(matrix)
  } else if (nrow(matrix) == n_cells && ncol(matrix) == n_genes) {
    warnings <- .warn_add(warnings, "矩阵是 cells × genes 方向，已自动转置（契约要求 genes × cells）")
  } else {
    .prereq("矩阵维度 %d×%d 与 features(%d) / barcodes(%d) 不匹配",
            nrow(matrix), ncol(matrix), n_genes, n_cells)
  }
  matrix <- methods::as(matrix, "CsparseMatrix")

  if (anyDuplicated(genes)) {
    genes <- make.unique(genes, sep = "-")
    warnings <- .warn_add(warnings, "features.tsv 存在重复基因名，已自动去重")
  }
  rownames(matrix) <- barcodes
  colnames(matrix) <- genes

  obs <- .read_obs_csv(file.path(bundle_dir, "obs.csv"), barcodes)
  pca <- .read_csv_matrix(file.path(bundle_dir, "pca.csv"), barcodes, "PC")
  if (!is.null(pca)) colnames(pca) <- paste0("PC", seq_len(ncol(pca)))
  embed <- .read_csv_matrix(file.path(bundle_dir, "embed.csv"), barcodes, "dim")

  # X_kind：meta.json 里的 "lognorm" / "counts"；缺失时按数值特征推断
  x_kind <- ""
  for (key in c("X_kind", "x_kind", "matrix_kind", "kind")) {
    if (!is.null(meta[[key]]) && length(meta[[key]]) >= 1 && nzchar(as.character(meta[[key]])[1])) {
      x_kind <- tolower(as.character(meta[[key]])[1])
      break
    }
  }
  if (!nzchar(x_kind)) {
    if (is_counts_file) {
      x_kind <- "counts"
    } else {
      values <- if (methods::is(matrix, "sparseMatrix")) matrix@x else as.numeric(matrix)
      looks_integer <- length(values) > 0 && all(abs(values - round(values)) < 1e-4) && max(values) > 50
      x_kind <- if (looks_integer) "counts" else "lognorm"
      .log("meta.json 未提供 X_kind，按数值特征推断为 %s", x_kind)
    }
  }
  if (!(x_kind %in% c("counts", "lognorm"))) x_kind <- if (is_counts_file) "counts" else "lognorm"

  list(
    dir = bundle_dir, meta = meta, warnings = warnings,
    matrix = matrix, genes = genes, barcodes = barcodes,
    obs = obs, pca = pca, embed = embed, x_kind = x_kind,
    from_counts_file = is_counts_file
  )
}

# ---------------------------------------------------------------------------
# 下采样 / 基因上限
# ---------------------------------------------------------------------------
.limit_bundle <- function(bundle, params, warnings) {
  matrix <- bundle$matrix
  barcodes <- bundle$barcodes
  genes <- bundle$genes
  seed <- as.integer(params$seed %||% 0L)
  max_cells <- max(10L, .as_int(params$max_cells, 200000L))
  max_genes <- max(10L, .as_int(params$max_genes, 10000L))
  pca <- bundle$pca
  obs <- bundle$obs

  if (length(barcodes) > max_cells) {
    set.seed(seed)
    keep <- sort(sample.int(length(barcodes), max_cells))
    warnings <- .warn_add(warnings,
      "细胞数 %d 超过上限 %d，按 seed=%d 随机下采样到 %d 个细胞",
      length(barcodes), max_cells, seed, max_cells)
    matrix <- matrix[keep, , drop = FALSE]
    barcodes <- barcodes[keep]
    if (!is.null(pca)) pca <- pca[keep, , drop = FALSE]
    if (!is.null(obs)) obs <- obs[keep, , drop = FALSE]
    if (!is.null(bundle$embed)) bundle$embed <- bundle$embed[keep, , drop = FALSE]
  }

  if (length(genes) > max_genes) {
    warnings <- .warn_add(warnings,
      "基因数 %d 超过上限 %d，按跨细胞方差取前 %d 个高方差基因（随后仍在其内部做 HVG）",
      length(genes), max_genes, max_genes)
    mu <- Matrix::colMeans(matrix)
    mu2 <- Matrix::colMeans(matrix * matrix)
    variance <- mu2 - mu * mu
    order_index <- order(variance, decreasing = TRUE)[seq_len(max_genes)]
    order_index <- sort(order_index)
    matrix <- matrix[, order_index, drop = FALSE]
    genes <- genes[order_index]
    rownames(matrix) <- barcodes
    colnames(matrix) <- genes
    if (!is.null(pca)) pca <- NULL  # PCA 与基因子集不对应，必须重算
  }

  bundle$matrix <- matrix
  bundle$barcodes <- barcodes
  bundle$genes <- genes
  bundle$pca <- pca
  bundle$obs <- obs
  bundle$warnings <- warnings
  bundle
}

# ---------------------------------------------------------------------------
# 构建 Seurat 对象
# ---------------------------------------------------------------------------
.obs_batch_key <- function(obs, key = "") {
  if (is.null(obs) || !is.data.frame(obs) || ncol(obs) == 0) return("")
  key <- trimws(as.character(key %||% ""))
  if (nzchar(key) && key %in% colnames(obs)) return(key)
  candidates <- c("sample", "orig.ident", "batch", "group", "donor", "patient", "library", "lane")
  hit <- candidates[candidates %in% colnames(obs)]
  if (length(hit)) return(hit[1])
  ""
}

.build_seurat <- function(bundle, warnings, batch_key = "") {
  # 内部统一用 cells × genes（bundle 契约方向转换后）；Seurat 的 counts 要求 genes × cells
  n_cells <- nrow(bundle$matrix)
  n_genes <- ncol(bundle$matrix)
  counts_matrix <- Matrix::t(bundle$matrix)  # genes × cells
  .log("CreateSeuratObject：%d 细胞 × %d 基因", n_cells, n_genes)
  object <- tryCatch(
    Seurat::CreateSeuratObject(
      counts = counts_matrix,
      project = "dsh_singlecell",
      min.cells = 0, min.features = 0
    ),
    error = function(error) .bad("CreateSeuratObject 失败：%s", error$message)
  )
  object <- tryCatch(Seurat::RenameCells(object, new.names = bundle$barcodes), error = function(error) object)
  .finish_seurat_build(object, bundle, n_cells, warnings)
}

.finish_seurat_build <- function(object, bundle, n_cells, warnings) {
  # 关键：把 bundle 里的细胞元数据（sample/group/batch/...）塞进 Seurat meta.data，
  # 否则 Harmony 找不到批次列（用户软件按样本做 Harmony 整合）。
  obs <- bundle$obs
  if (!is.null(obs) && nrow(obs) == n_cells && ncol(obs) > 0) {
    meta <- as.data.frame(obs, stringsAsFactors = FALSE)
    rownames(meta) <- bundle$barcodes
    existing <- colnames(object[[]])
    for (column in setdiff(colnames(meta), existing)) {
      values <- meta[[column]]
      if (is.list(values)) {
        values <- vapply(values, function(item) paste(as.character(item), collapse = ";"), character(1))
      }
      object[[column]] <- values
    }
    .log("附加细胞元数据：%s", paste(setdiff(colnames(meta), existing), collapse = ", "))
  }
  already_normalized <- FALSE
  if (!bundle$from_counts_file && identical(bundle$x_kind, "lognorm")) {
    already_normalized <- TRUE
    object <- .set_data_layer(object, Matrix::t(as(bundle$matrix, "CsparseMatrix")))
    warnings <- .warn_add(warnings,
      "bundle 只有 data.mtx（meta.X_kind=lognorm），直接作为 data 层使用，跳过 NormalizeData 以免二次标准化")
  }
  list(object = object, warnings = warnings, already_normalized = already_normalized)
}

.set_data_layer <- function(object, matrix) {
  assay <- SeuratObject::DefaultAssay(object)
  assay_object <- object[[assay]]
  if (inherits(assay_object, "Assay5")) {
    SeuratObject::LayerData(assay_object, layer = "data") <- matrix
    object[[assay]] <- assay_object
  } else {
    assay_object@data <- matrix
    object[[assay]] <- assay_object
  }
  object
}

# ---------------------------------------------------------------------------
# 上游标准化 → HVG → scale → PCA
# ---------------------------------------------------------------------------
.compute_pca <- function(object, params, warnings, already_normalized = FALSE) {
  n_top_genes <- max(50L, .as_int(params$n_top_genes, 2000L))
  n_pcs <- max(2L, .as_int(params$n_pcs, 30L))
  seed <- as.integer(params$seed %||% 0L)
  n_cells <- ncol(object)
  n_genes <- nrow(object)

  if (!already_normalized) {
    .log("NormalizeData(LogNormalize, scale.factor=1e4) …")
    object <- tryCatch(
      Seurat::NormalizeData(object, normalization.method = "LogNormalize",
                            scale.factor = 1e4, verbose = FALSE),
      error = function(error) .bad("NormalizeData 失败：%s", error$message)
    )
  }

  nfeatures <- max(10L, min(n_top_genes, n_genes))
  if (n_top_genes > n_genes) {
    warnings <- .warn_add(warnings,
      "请求 n_top_genes=%d 超过基因总数 %d，FindVariableFeatures 使用全部基因", n_top_genes, n_genes)
  }
  .log("FindVariableFeatures(vst, nfeatures=%d) …", nfeatures)
  object <- tryCatch(
    Seurat::FindVariableFeatures(object, selection.method = "vst",
                                 nfeatures = nfeatures, verbose = FALSE),
    error = function(error) .bad("FindVariableFeatures 失败：%s", error$message)
  )
  hvg <- Seurat::VariableFeatures(object)
  if (length(hvg) < 10L) .prereq("高变基因不足 10 个（%d），无法完成 PCA", length(hvg))

  .log("ScaleData(features=%d) …", length(hvg))
  object <- tryCatch(
    Seurat::ScaleData(object, features = hvg, verbose = FALSE),
    error = function(error) .bad("ScaleData 失败：%s", error$message)
  )

  npcs <- max(2L, min(n_pcs, n_cells - 1L, length(hvg) - 1L))
  if (npcs < 2L) .prereq("可用主成分不足 2 个（cells=%d, HVG=%d），无法完成 PCA", n_cells, length(hvg))
  if (npcs < n_pcs) {
    warnings <- .warn_add(warnings, "请求 n_pcs=%d 超过数据上限，已收敛为 %d 个主成分", n_pcs, npcs)
  }
  .log("RunPCA(npcs=%d, seed.use=%d) …", npcs, seed)
  object <- tryCatch(
    Seurat::RunPCA(object, npcs = npcs, features = hvg, seed.use = seed, verbose = FALSE),
    error = function(error) .bad("RunPCA 失败：%s", error$message)
  )
  list(object = object, warnings = warnings, npcs = npcs, hvg = hvg)
}

.pca_variance <- function(object) {
  stdev <- tryCatch(SeuratObject::Stdev(object, reduction = "pca"), error = function(e) NULL)
  if (is.null(stdev) || !length(stdev)) {
    return(list(ratio = numeric(0), stdev = numeric(0)))
  }
  total <- sum(stdev^2)
  if (!is.finite(total) || total <= 0) total <- 1
  list(ratio = as.numeric(stdev^2 / total), stdev = as.numeric(stdev))
}

.attach_bundle_pca <- function(object, pca_matrix, seed) {
  assay <- SeuratObject::DefaultAssay(object)
  reduction <- SeuratObject::CreateDimReducObject(
    embeddings = pca_matrix,
    key = "PC_",
    assay = assay
  )
  object[["pca"]] <- reduction
  object
}

.detect_batch_key <- function(object, batch_key) {
  meta <- tryCatch(object[[]], error = function(e) NULL)
  if (is.null(meta) || ncol(meta) == 0) return("")
  key <- trimws(as.character(batch_key %||% ""))
  if (nzchar(key) && key %in% colnames(meta)) return(key)
  candidates <- c("sample", "orig.ident", "batch", "group", "donor", "patient", "library", "lane")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit)) return(hit[1])
  ""
}

# 按用户 pro.R 的官方写法：Seurat v5/v6 IntegrateLayers(method=HarmonyIntegration / RPCAIntegration)
.run_integration <- function(object, integration, batch_key, warnings, seed) {
  method <- tolower(as.character(integration %||% "none"))
  info <- list(used = FALSE, method = method, batch_key = NULL, reduction = "pca")
  if (!(method %in% c("harmony", "rpca"))) {
    return(list(object = object, warnings = warnings, reduction = "pca", info = info))
  }
  key <- .detect_batch_key(object, batch_key)
  if (!nzchar(key)) {
    warnings <- .warn_add(warnings,
      "未找到批次列（sample/orig.ident/group），跳过 %s 整合（可用 params.batch_key 指定）", method)
    return(list(object = object, warnings = warnings, reduction = "pca", info = info))
  }
  meta <- tryCatch(object[[]], error = function(e) NULL)
  if (is.null(meta) || length(unique(meta[[key]][!is.na(meta[[key]])])) < 2L) {
    warnings <- .warn_add(warnings, "批次列 %s 只有 1 个取值，跳过 %s 整合", key, method)
    return(list(object = object, warnings = warnings, reduction = "pca", info = info))
  }
  reduction_name <- if (method == "harmony") "harmony" else "integrated.rpca"
  .log("%s 整合（IntegrateLayers，批次列 %s）…", toupper(method), key)
  # Seurat v5 需要按样本拆层后 IntegrateLayers 才有 layer 可整合；已拆过的（merge 出来的多层对象）直接用
  already_split <- tryCatch({
    assay_name <- SeuratObject::DefaultAssay(object)
    layer_names <- SeuratObject::Layers(object[[assay_name]])
    length(layer_names) >= 2L && any(grepl("\\.[^.]+$", layer_names))
  }, error = function(error) FALSE)
  split_ok <- TRUE
  if (!already_split) {
    split_ok <- tryCatch({
      object[["RNA"]] <- split(object[["RNA"]], f = factor(meta[[key]]))
      TRUE
    }, error = function(error) {
      .log("WARN: 拆分 RNA 层失败：%s", conditionMessage(error))
      FALSE
    })
  } else {
    .log("RNA 层已按样本拆分（%s），直接 IntegrateLayers", paste(SeuratObject::Layers(object[[SeuratObject::DefaultAssay(object)]]), collapse = ", "))
  }
  if (!split_ok) {
    warnings <- .warn_add(warnings, "RNA 层拆分失败（%s），跳过整合", method)
    return(list(object = object, warnings = warnings, reduction = "pca", info = info))
  }
  # IntegrateLayers 需要按层分组的 data 层（merge 出来的对象只有 counts.* 层）：
  # 缺 data.* 时按 counts 层做 LogNormalize + HVG + Scale，补齐后再整合。
  assay_name <- tryCatch(SeuratObject::DefaultAssay(object), error = function(e) "RNA")
  data_missing <- tryCatch({
    layer_names <- SeuratObject::Layers(object[[assay_name]])
    sum(grepl("^data", layer_names)) < 2L
  }, error = function(e) TRUE)
  if (data_missing) {
    object <- tryCatch({
      assay_object <- object[[assay_name]]
      count_layers <- grep("^counts", SeuratObject::Layers(assay_object), value = TRUE)
      for (layer_name in count_layers) {
        counts <- SeuratObject::LayerData(assay_object, layer = layer_name)
        totals <- Matrix::colSums(counts)
        totals[totals <= 0] <- 1
        normalized <- log1p(Matrix::t(Matrix::t(counts) / totals * 10000))
        SeuratObject::LayerData(assay_object, layer = sub("^counts", "data", layer_name)) <- normalized
      }
      object[[assay_name]] <- assay_object
      object <- Seurat::FindVariableFeatures(object, nfeatures = 3000, verbose = FALSE)
      object <- Seurat::ScaleData(object, verbose = FALSE)
      .log("补齐 data/scale 层（%d 个批次层）", length(count_layers))
      object
    }, error = function(error) {
      .log("WARN: 补齐 data/scale 层失败：%s", conditionMessage(error))
      object
    })
  }
  integrate_fn <- if (method == "harmony") Seurat::HarmonyIntegration else Seurat::RPCAIntegration
  integrated <- if (method == "rpca") {
    # RPCA 的 k.weight 不能超过最小样本的细胞数，否则小样本数据会直接报错
    layer_sizes <- tryCatch({
      a <- object[[assay_name]]
      vapply(SeuratObject::Layers(a), function(ly) ncol(SeuratObject::LayerData(a, layer = ly)), integer(1))
    }, error = function(e) NULL)
    min_cells <- if (is.null(layer_sizes) || !length(layer_sizes)) 0L else as.integer(min(layer_sizes))
    candidates <- unique(as.integer(c(
      min(100L, max(20L, min_cells - 2L)),
      min(100L, max(20L, floor(min_cells / 2))),
      30L, 20L
    )))
    candidates <- candidates[candidates >= 2L & candidates < min_cells]
    result_integrated <- NULL
    for (k_weight in candidates) {
      .log("RPCA 尝试 k.weight=%d（最小样本 %d 细胞）", k_weight, min_cells)
      attempt <- tryCatch(
        Seurat::IntegrateLayers(object = object, method = integrate_fn,
                                orig.reduction = "pca", new.reduction = reduction_name,
                                k.weight = k_weight, verbose = FALSE),
        error = function(error) error
      )
      if (!inherits(attempt, "error")) {
        result_integrated <- attempt
        break
      }
      .log("RPCA k.weight=%d 失败：%s", k_weight, conditionMessage(attempt))
    }
    if (is.null(result_integrated)) {
      warnings <- .warn_add(warnings, "RPCA 整合失败（小样本/锚点不足），回退 Harmony 整合")
      .log("RPCA 失败，回退 Harmony…")
      result_integrated <- tryCatch(
        Seurat::IntegrateLayers(object = object, method = Seurat::HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony",
                                verbose = FALSE),
        error = function(error) error
      )
      if (!inherits(result_integrated, "error")) {
        reduction_name <- "harmony"
        integrate_fn <- Seurat::HarmonyIntegration
        info$fallback <- TRUE
      } else {
        warnings <- .warn_add(warnings, "Harmony 兜底也失败（%s），回退 PCA", conditionMessage(result_integrated))
        return(list(object = object, warnings = warnings, reduction = "pca", info = info))
      }
    }
    result_integrated
  } else {
    tryCatch(
      Seurat::IntegrateLayers(object = object, method = integrate_fn,
                              orig.reduction = "pca", new.reduction = reduction_name,
                              verbose = FALSE),
      error = function(error) error
    )
  }
  if (inherits(integrated, "error")) {
    warnings <- .warn_add(warnings, "%s 整合失败（%s），回退 PCA", toupper(method), conditionMessage(integrated))
    return(list(object = object, warnings = warnings, reduction = "pca", info = info))
  }
  info$used <- TRUE
  info$batch_key <- key
  info$reduction <- reduction_name
  list(object = integrated, warnings = warnings, reduction = reduction_name, info = info)
}

.graph_names <- function(object) {
  tryCatch(SeuratObject::Graphs(object), error = function(e) character(0))
}

.find_neighbors <- function(object, dims, n_neighbors, warnings, reduction = "pca") {
  k <- max(2L, .as_int(n_neighbors, 20L))
  k <- min(k, ncol(object) - 1L)
  if (k < 2L) .prereq("细胞数不足，无法构建 kNN 图")
  .log("FindNeighbors(dims=%d:%d, k=%d) …", min(dims), max(dims), k)
  # 注意：这里必须显式调用而不是 do.call(…, list(object = object, …))——
  # 把 S4 Seurat 对象作为 do.call 的实参会触发 Seurat 内部对巨大对象的 deparse，
  # 2208 细胞时 FindNeighbors 会从 0.1 秒退化到 5 秒、UMAP 退化到 60 秒。
  object <- tryCatch(
    Seurat::FindNeighbors(object = object, reduction = reduction, dims = dims,
                          k.param = as.integer(k), annoy.metric = "euclidean",
                          verbose = FALSE),
    error = function(error) .bad("FindNeighbors 失败：%s", error$message)
  )
  graphs <- .graph_names(object)
  snn <- graphs[grepl("(^|_)snn$", graphs)]
  list(object = object, warnings = warnings, k = as.integer(k), snn = if (length(snn)) snn[1] else NULL)
}

.find_clusters <- function(object, resolution, seed, snn_name, warnings) {
  resolution <- max(0.01, .as_num(resolution, 0.8))
  .log("FindClusters(resolution=%.2f, seed=%d) …", resolution, seed)
  object <- tryCatch(
    if (!is.null(snn_name)) {
      Seurat::FindClusters(object = object, graph.name = snn_name, resolution = resolution,
                           random.seed = as.integer(seed), verbose = FALSE)
    } else {
      Seurat::FindClusters(object = object, resolution = resolution,
                           random.seed = as.integer(seed), verbose = FALSE)
    },
    error = function(error) .bad("FindClusters 失败：%s", error$message)
  )
  labels <- NULL
  idents <- tryCatch(SeuratObject::Idents(object), error = function(error) NULL)
  if (!is.null(idents)) {
    labels <- as.character(idents)
    names(labels) <- names(idents)
    labels <- unname(labels[colnames(object)])
  }
  if (is.null(labels) || length(labels) != ncol(object) || all(is.na(labels))) {
    metadata <- object[[]]
    candidates <- c("seurat_clusters", grep("_res\\.", colnames(metadata), value = TRUE), "clusters")
    column <- NULL
    for (candidate in candidates) {
      if (candidate %in% colnames(metadata)) { column <- candidate; break }
    }
    if (is.null(column)) .bad("FindClusters 完成但找不到簇标签列")
    labels <- as.character(metadata[[column]])
  }
  labels[is.na(labels)] <- "-1"
  list(object = object, warnings = warnings, labels = labels, resolution = resolution)
}

.run_umap <- function(object, dims, n_neighbors, seed, warnings, reduction = "pca", metric = "euclidean") {
  k <- max(2L, .as_int(n_neighbors, 20L))
  k <- min(k, max(2L, ncol(object) - 2L))
  .log("RunUMAP(dims=%d:%d, n.neighbors=%d, seed.use=%d) …", min(dims), max(dims), k, seed)
  object <- tryCatch(
    Seurat::RunUMAP(object = object, reduction = reduction, dims = dims,
                    n.neighbors = as.integer(k), seed.use = as.integer(seed),
                    umap.method = "uwot", min.dist = 0.3, metric = metric,
                    verbose = FALSE),
    error = function(error) .bad("RunUMAP 失败：%s", error$message)
  )
  embedding <- tryCatch(Seurat::Embeddings(object, reduction = "umap"), error = function(e) NULL)
  if (is.null(embedding) || nrow(embedding) != ncol(object)) .bad("RunUMAP 完成但取不到 UMAP 坐标")
  embedding <- embedding[, seq_len(min(2L, ncol(embedding))), drop = FALSE]
  colnames(embedding) <- c("UMAP_1", "UMAP_2")[seq_len(ncol(embedding))]
  rownames(embedding) <- colnames(object)
  list(object = object, warnings = warnings, embedding = embedding, k = as.integer(k))
}

# ---------------------------------------------------------------------------
# 出图（白底英文）
# ---------------------------------------------------------------------------
.new_png <- function(path, width = 1400L, height = 1000L) {
  grDevices::png(filename = path, width = width, height = height, res = 140, bg = "white")
}

.plot_pca_variance <- function(ratio, path) {
  if (!length(ratio)) return(NULL)
  ratio <- ratio[seq_len(min(30L, length(ratio)))]
  tryCatch({
    .new_png(path)
    graphics::par(mar = c(5, 5.2, 4, 2))
    labels <- if (length(ratio) > 20) seq_along(ratio) else seq_along(ratio)
    graphics::barplot(
      ratio * 100,
      names.arg = as.character(labels),
      col = "#4C78A8", border = NA,
      xlab = "Principal component", ylab = "Variance explained (%)",
      main = "Seurat 6.0 PCA variance", cex.names = 0.8, las = 2
    )
    graphics::abline(h = pretty(c(0, max(ratio) * 100)), col = "grey85", lty = 3)
    grDevices::dev.off()
    file.exists(path)
  }, error = function(error) {
    try(grDevices::dev.off(), silent = TRUE)
    .log("WARN: PCA 碎石图生成失败：%s", error$message)
    FALSE
  })
}

.cluster_colors <- function(labels) {
  levels <- sort(unique(as.character(labels)), na.last = TRUE)
  ordered_levels <- levels[order(suppressWarnings(as.numeric(levels)), levels, na.last = TRUE)]
  n <- length(ordered_levels)
  colors <- if (n <= 1) "#4C78A8" else grDevices::hcl.colors(n, palette = "Dark 3")
  list(levels = ordered_levels, colors = colors)
}

.plot_umap <- function(embedding, labels, path, title = "Seurat 6.0 UMAP") {
  tryCatch({
    .new_png(path)
    graphics::par(mar = c(5, 5.2, 4, 2))
    x <- embedding[, 1]
    y <- embedding[, 2]
    if (is.null(labels)) {
      graphics::plot(x, y, pch = 19, cex = 0.55, col = "#4C78A8",
                     xlab = "UMAP 1", ylab = "UMAP 2", main = title)
    } else {
      palette <- .cluster_colors(labels)
      codes <- match(as.character(labels), palette$levels)
      graphics::plot(x, y, pch = 19, cex = 0.55,
                     col = palette$colors[codes],
                     xlab = "UMAP 1", ylab = "UMAP 2", main = title)
      if (length(palette$levels) <= 25) {
        graphics::legend("topright", legend = palette$levels, col = palette$colors,
                         pch = 19, pt.cex = 1.1, cex = 0.78, bty = "n",
                         title = "Cluster", title.cex = 0.85)
      }
    }
    grDevices::dev.off()
    file.exists(path)
  }, error = function(error) {
    try(grDevices::dev.off(), silent = TRUE)
    .log("WARN: UMAP 图生成失败：%s", error$message)
    FALSE
  })
}

.plot_cluster_sizes <- function(labels, path) {
  tryCatch({
    counts <- table(labels)
    if (!length(counts)) return(NULL)
    .new_png(path, 1200L, 900L)
    graphics::par(mar = c(5, 6.5, 4, 2))
    graphics::barplot(
      as.numeric(counts),
      names.arg = names(counts),
      horiz = TRUE, las = 1,
      col = "#4C78A8", border = NA,
      xlab = "Number of cells", ylab = "Cluster",
      main = "Seurat 6.0 cluster sizes", cex.names = 0.85
    )
    grDevices::dev.off()
    file.exists(path)
  }, error = function(error) {
    try(grDevices::dev.off(), silent = TRUE)
    .log("WARN: 簇大小图生成失败：%s", error$message)
    FALSE
  })
}

# ---------------------------------------------------------------------------
# 写输出
# ---------------------------------------------------------------------------
.write_obs_extra <- function(path, barcodes, columns) {
  frame <- data.frame(cell = barcodes, stringsAsFactors = FALSE, check.names = FALSE)
  for (name in names(columns)) {
    frame[[name]] <- columns[[name]]
  }
  utils::write.csv(frame, path, row.names = FALSE, quote = FALSE, na = "")
  path
}

.write_pca <- function(path, barcodes, pca_matrix) {
  frame <- data.frame(cell = barcodes, stringsAsFactors = FALSE, check.names = FALSE)
  frame <- cbind(frame, as.data.frame(pca_matrix, stringsAsFactors = FALSE, check.names = FALSE))
  utils::write.csv(frame, path, row.names = FALSE, quote = FALSE)
  path
}

.write_embed <- function(path, barcodes, embedding) {
  frame <- data.frame(cell = barcodes, dim1 = embedding[, 1], dim2 = embedding[, 2],
                      stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.csv(frame, path, row.names = FALSE, quote = FALSE)
  path
}

.write_result <- function(path, payload) {
  jsonlite::write_json(payload, path, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null", na = "null")
  path
}

.png_size <- function(path) {
  tryCatch({
    size <- file.info(path)$size
    # 从 PNG IHDR 读宽高
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    header <- readBin(connection, "raw", n = 24)
    width <- sum(as.integer(header[17:20]) * 256^(3:0))
    height <- sum(as.integer(header[21:24]) * 256^(3:0))
    list(width = width, height = height, bytes = as.numeric(size))
  }, error = function(e) list(width = 0, height = 0, bytes = 0))
}

.figure_entry <- function(id, title, path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  size <- .png_size(path)
  list(id = id, title = title, file = normalizePath(path, mustWork = TRUE),
       width = size$width, height = size$height)
}

# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
.main <- function() {
  started <- Sys.time()
  cli <- .parse_cli(commandArgs(trailingOnly = TRUE))
  params <- .parse_params(cli$params)
  step <- tolower(as.character(params$step %||% "full")[1])
  allowed_steps <- c("pca", "neighbors", "cluster", "embed", "full")
  if (!(step %in% allowed_steps)) {
    .bad("不支持的 step：%s（可选：%s）", step, paste(allowed_steps, collapse = "|"))
  }
  out_dir <- cli$out %||% ""
  if (!nzchar(out_dir)) .bad("缺少 --out 参数")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(out_dir)) .bad("无法创建输出目录：%s", out_dir)
  bundle_dir <- cli$bundle %||% ""

  seed <- as.integer(params$seed %||% 0L)
  n_pcs <- max(2L, .as_int(params$n_pcs, 30L))
  n_neighbors <- max(2L, .as_int(params$n_neighbors, 20L))
  integration <- tolower(trimws(as.character(params$integration %||%
    (if (.is_true(params$use_harmony, TRUE)) "harmony" else "none"))))
  if (!(integration %in% c("none", "harmony", "rpca", "pca"))) integration <- "harmony"
  if (integration == "pca") integration <- "none"
  resolution_default <- if (integration == "none") 0.5 else 0.8
  resolution <- if (is.null(params$resolution)) resolution_default else max(0.01, .as_num(params$resolution, resolution_default))
  metric <- as.character(params$metric %||% "euclidean")
  n_top_genes <- max(50L, .as_int(params$n_top_genes, 3000L))
  do_umap <- .is_true(params$umap, TRUE)
  use_harmony <- .is_true(params$use_harmony, TRUE)
  batch_key <- as.character(params$batch_key %||% "")

  .log("dsh-singlecell Seurat 适配器：step=%s seed=%d n_pcs=%d resolution=%.2f", step, seed, n_pcs, resolution)

  seurat_info <- .setup_seurat()
  .log("使用 Seurat %s（独立库 %s）", seurat_info$version, seurat_info$lib)

  suppressPackageStartupMessages({
    library(Matrix)     # noqa
    library(Seurat)     # noqa
  })
  options(future.globals.maxSize = 4 * 1024^3)
  if (requireNamespace("future", quietly = TRUE)) {
    try(future::plan("sequential"), silent = TRUE)
  }
  set.seed(seed)

  warnings <- character(0)
  bundle <- .read_bundle(bundle_dir, warnings)
  warnings <- bundle$warnings
  bundle <- .limit_bundle(bundle, params, warnings)
  warnings <- bundle$warnings
  barcodes <- bundle$barcodes

  built <- .build_seurat(bundle, warnings)
  object <- built$object
  warnings <- built$warnings
  already_normalized <- built$already_normalized

  final_labels <- NULL
  pca_matrix <- NULL
  variance_ratio <- numeric(0)
  embedding <- NULL
  hvg_count <- NULL
  npcs <- 0L
  k_used <- n_neighbors
  used_bundle_pca <- FALSE
  neighbor_out <- NULL

  need_pca <- TRUE
  if (step %in% c("neighbors", "cluster", "embed") && !is.null(bundle$pca)) {
    # 上游步骤已经把 PCA 写回 bundle 时直接复用（保证与 pca 步骤一致）
    n_avail <- ncol(bundle$pca)
    n_use <- min(n_pcs, n_avail)
    pca_matrix <- bundle$pca[, seq_len(n_use), drop = FALSE]
    colnames(pca_matrix) <- paste0("PC", seq_len(ncol(pca_matrix)))
    if (n_use < n_pcs) {
      warnings <- .warn_add(warnings,
        "bundle 里只有 %d 个 PC（请求 %d），本步骤按 %d 个 PC 计算", n_avail, n_pcs, n_use)
    }
    dims <- .resolve_dims(params$dims, n_use)
    object <- .attach_bundle_pca(object, pca_matrix, seed)
    npcs <- n_use
    used_bundle_pca <- TRUE
    .log("复用 bundle pca.csv：%d 个 PC", n_use)
  }

  if (step == "pca" || step == "full" || is.null(pca_matrix)) {
    pca_out <- .compute_pca(object, params, warnings, already_normalized)
    object <- pca_out$object
    warnings <- pca_out$warnings
    npcs <- pca_out$npcs
    hvg_count <- length(pca_out$hvg)
    variance <- .pca_variance(object)
    variance_ratio <- variance$ratio
    pca_matrix <- Seurat::Embeddings(object, reduction = "pca")[, seq_len(min(n_pcs, npcs)), drop = FALSE]
    colnames(pca_matrix) <- paste0("PC", seq_len(ncol(pca_matrix)))
  }
  # ---- 批次整合（pro.R：IntegrateLayers Harmony/RPCA；none = 未整合 PCA）
  harmony_reduction <- "pca"
  harmony_info <- list(used = FALSE, method = integration, batch_key = NULL, reduction = "pca")
  requested_reduction <- trimws(as.character(params$reduction %||% ""))
  if (nzchar(requested_reduction)) {
    harmony_reduction <- requested_reduction
  }
  if (step %in% c("neighbors", "cluster", "embed", "full")) {
    # 集群/嵌入步骤可能只拿到 bundle pca：若显式要求 harmony/rpca 而没有对应 reduction，就现场整合
    if (integration %in% c("harmony", "rpca") && (step != "pca")) {
      integration_out <- .run_integration(object, integration, batch_key, warnings, seed)
      object <- integration_out$object
      warnings <- integration_out$warnings
      harmony_reduction <- integration_out$reduction
      harmony_info <- integration_out$info
    }
  }

  requested_dims <- .requested_dims(params$dims)
  dims <- .resolve_dims(params$dims, npcs)
  if (length(requested_dims) && max(requested_dims) > npcs && step %in% c("neighbors", "cluster", "embed", "full")) {
    warnings <- .warn_add(warnings,
      "请求使用前 %d 个主成分，但本次只有 %d 个 PC，构图/UMAP 已收敛为 1:%d",
      max(requested_dims), npcs, npcs)
  }

  tables <- list()
  figures <- list()

  pca_path <- NULL
  if (step %in% c("pca", "full")) {
    pca_path <- .write_pca(file.path(out_dir, "pca.csv"), barcodes, pca_matrix)
    .log("写出 %s（%d × %d）", pca_path, nrow(pca_matrix), ncol(pca_matrix))
    if (length(variance_ratio)) {
      ratio_top <- variance_ratio[seq_len(min(20L, length(variance_ratio)))]
      rows <- lapply(seq_along(ratio_top), function(index) {
        list(PC = paste0("PC", index), variance_ratio = round(ratio_top[index], 6))
      })
      tables <- c(tables, list(list(
        id = "pca_variance", title = "Seurat PCA variance explained",
        columns = list(
          list(key = "PC", label = "PC", type = "string"),
          list(key = "variance_ratio", label = "Variance ratio", type = "number")
        ),
        rows = rows
      )))
      elbow_path <- file.path(out_dir, "seurat_pca_variance.png")
      elbow_ok <- .plot_pca_variance(variance_ratio, elbow_path)
      figures <- c(figures, list(.figure_entry(
        "pca_elbow", "Seurat PCA variance", if (isTRUE(elbow_ok)) elbow_path else NULL
      )))
    }
  } else {
    # 非 pca 步骤也写回 pca.csv，方便引擎把 PCA 存回 obsm["X_pca"]
    pca_path <- .write_pca(file.path(out_dir, "pca.csv"), barcodes, pca_matrix)
  }

  if (step %in% c("neighbors", "cluster", "full")) {
    neighbor_out <- .find_neighbors(object, dims, n_neighbors, warnings, reduction = harmony_reduction)
    object <- neighbor_out$object
    warnings <- neighbor_out$warnings
    k_used <- neighbor_out$k
  }

  if (step %in% c("cluster", "full")) {
    cluster_out <- .find_clusters(object, resolution, seed, if (!is.null(neighbor_out)) neighbor_out$snn else NULL, warnings)
    object <- cluster_out$object
    warnings <- cluster_out$warnings
    final_labels <- cluster_out$labels
    cluster_counts <- table(final_labels)
    rows <- lapply(names(cluster_counts), function(name) {
      list(cluster = as.character(name), n_cells = as.integer(cluster_counts[[name]]))
    })
    tables <- c(tables, list(list(
      id = "cluster_sizes", title = "Seurat cluster sizes",
      columns = list(
        list(key = "cluster", label = "Cluster", type = "string"),
        list(key = "n_cells", label = "Cells", type = "number")
      ),
      rows = rows
    )))
    n_cluster_now <- length(unique(final_labels))
    if (n_cluster_now > 30L) {
      warnings <- .warn_add(warnings,
        "聚类数 %d 偏多（参考 20–30）：建议把 resolution 调低（当前 %.2f）", n_cluster_now, resolution)
    } else if (n_cluster_now < 10L) {
      warnings <- .warn_add(warnings,
        "聚类数 %d 偏少（参考 20–30）：建议把 resolution 调高（当前 %.2f）", n_cluster_now, resolution)
    }
    if (step == "cluster") {
      sizes_path <- file.path(out_dir, "seurat_cluster_sizes.png")
      sizes_ok <- .plot_cluster_sizes(final_labels, sizes_path)
      figures <- c(figures, list(.figure_entry(
        "cluster_sizes", "Seurat cluster sizes", if (isTRUE(sizes_ok)) sizes_path else NULL
      )))
    }
  }

  embed_path <- NULL
  if (step == "embed" || (step == "full" && do_umap)) {
    umap_out <- .run_umap(object, dims, n_neighbors, seed, warnings, reduction = harmony_reduction, metric = metric)
    object <- umap_out$object
    warnings <- umap_out$warnings
    embedding <- umap_out$embedding
    embed_path <- .write_embed(file.path(out_dir, "embed.csv"), barcodes, embedding)
    .log("写出 %s（%d × 2）", embed_path, nrow(embedding))
    labels_for_plot <- final_labels
    if (is.null(labels_for_plot) && !is.null(bundle$obs) && "clusters" %in% colnames(bundle$obs)) {
      labels_for_plot <- as.character(bundle$obs[["clusters"]])
    }
    umap_path <- file.path(out_dir, "seurat_umap.png")
    umap_ok <- .plot_umap(embedding, labels_for_plot, umap_path)
    figures <- c(figures, list(.figure_entry(
      "umap", "Seurat 6.0 UMAP", if (isTRUE(umap_ok)) umap_path else NULL
    )))
  } else if (step == "full" && !do_umap) {
    warnings <- .warn_add(warnings, "params.umap=false，full 流程未计算 UMAP")
  }

  obs_columns <- list()
  if (!is.null(final_labels)) obs_columns[["seurat_clusters"]] <- final_labels
  if (!is.null(embedding) && step == "full") {
    obs_columns[["seurat_umap1"]] <- as.numeric(embedding[, 1])
    obs_columns[["seurat_umap2"]] <- as.numeric(embedding[, 2])
  }
  obs_path <- .write_obs_extra(file.path(out_dir, "obs_extra.csv"), barcodes, obs_columns)

  n_clusters <- if (!is.null(final_labels)) length(unique(final_labels)) else NULL
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  titles <- list(
    pca = "Seurat 6.0 官方 PCA",
    neighbors = "Seurat 6.0 官方 kNN 图",
    cluster = "Seurat 6.0 官方聚类（Louvain）",
    embed = "Seurat 6.0 官方 UMAP",
    full = "Seurat 6.0 官方降维全流程"
  )
  methods <- list(
    pca = "Seurat 6.0 官方 R 包：NormalizeData(LogNormalize) → FindVariableFeatures(vst) → ScaleData → RunPCA",
    neighbors = "Seurat 6.0 官方 R 包：RunPCA（或复用上游 PCA）→ FindNeighbors",
    cluster = "Seurat 6.0 官方 R 包：FindNeighbors → FindClusters（Louvain）",
    embed = "Seurat 6.0 官方 R 包：RunPCA（或复用上游 PCA）→ RunUMAP（uwot）",
    full = "Seurat 6.0 官方 R 包：NormalizeData(LogNormalize) → FindVariableFeatures(vst) → ScaleData → RunPCA → Harmony 整合（harmony 包）→ FindNeighbors → FindClusters（Louvain）→ RunUMAP"
  )

  metrics <- list(
    step = step,
    seurat_version = seurat_info$version,
    seurat_object_version = tryCatch(as.character(utils::packageVersion("SeuratObject")), error = function(e) NA_character_),
    seurat_lib = seurat_info$lib,
    seurat_path = seurat_info$path,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    n_cells = length(barcodes),
    n_genes = length(bundle$genes),
    n_hvg = if (is.null(hvg_count)) NULL else as.integer(hvg_count),
    n_pcs = as.integer(npcs),
    dims = paste(min(dims), max(dims), sep = ":"),
    n_neighbors = as.integer(k_used),
    resolution = resolution,
    seed = seed,
    n_clusters = if (is.null(n_clusters)) NULL else as.integer(n_clusters),
    used_bundle_pca = used_bundle_pca,
    integration = harmony_info$method,
    integration_used = isTRUE(harmony_info$used),
    harmony_batch_key = if (isTRUE(harmony_info$used)) harmony_info$batch_key else NULL,
    reduction = harmony_reduction,
    elapsed_seconds = round(elapsed, 3)
  )
  metrics <- metrics[!vapply(metrics, is.null, logical(1))]

  bullets <- c(
    sprintf("Seurat %s（独立 R 库 %s，R %s）：%d 个细胞 × %d 个基因",
            seurat_info$version, seurat_info$lib, paste(R.version$major, R.version$minor, sep = "."),
            length(barcodes), length(bundle$genes))
  )
  if (step %in% c("pca", "full")) {
    bullets <- c(bullets, sprintf("HVG（vst）= %d，RunPCA npcs = %d（seed=%d）",
                                  as.integer(hvg_count), as.integer(npcs), seed))
  } else if (used_bundle_pca) {
    bullets <- c(bullets, sprintf("复用上游 PCA（%d 个 PC），未重新标准化", npcs))
  }
  if (step %in% c("neighbors", "cluster", "full")) {
    bullets <- c(bullets, sprintf("FindNeighbors k = %d，使用 %d 个主成分", as.integer(k_used), length(dims)))
  }
  if (!is.null(final_labels)) {
    counts <- table(final_labels)
    bullets <- c(bullets, sprintf("FindClusters（Louvain，resolution=%.2f）得到 %d 个簇，最大簇 %d 个细胞",
                                  resolution, length(counts), max(as.integer(counts))))
  }
  if (!is.null(embedding)) {
    bullets <- c(bullets, sprintf("RunUMAP（uwot，n.neighbors=%d，前 %d 个 PC）输出 2 维坐标",
                                  as.integer(k_used), length(dims)))
  }
  if (isTRUE(harmony_info$used)) {
    label <- if (identical(harmony_info$method, "rpca")) "RPCA" else "Harmony"
    if (isTRUE(harmony_info$fallback)) label <- sprintf("%s（失败，已回退 Harmony）", label)
    bullets <- c(bullets, sprintf("%s 批次整合（Seurat IntegrateLayers）：按 %s 分组，使用 %s 嵌入构图/聚类/UMAP",
                                  label, harmony_info$batch_key, harmony_info$reduction))
  } else if (step %in% c("neighbors", "cluster", "embed", "full")) {
    bullets <- c(bullets, "未做批次整合（无批次列或 integration=none），使用未整合 PCA 构图")
  }
  bullets <- c(bullets, sprintf("Seurat 步骤耗时 %.2f 秒", elapsed))

  payload <- list(
    id = "ext_seurat6",
    title = titles[[step]],
    method = sprintf("%s（Seurat %s）", methods[[step]], seurat_info$version),
    summary = list(bullets = as.list(bullets), metrics = metrics),
    tables = tables,
    figures = figures[vapply(figures, function(item) !is.null(item), logical(1))],
    obs_file = normalizePath(obs_path, mustWork = TRUE),
    files = list(
      obs_extra = normalizePath(obs_path, mustWork = TRUE),
      pca = if (!is.null(pca_path)) normalizePath(pca_path, mustWork = TRUE) else NULL,
      embed = if (!is.null(embed_path)) normalizePath(embed_path, mustWork = TRUE) else NULL
    ),
    warnings = as.list(warnings),
    notes = list("Seurat >= 6.0 安装在独立 R 库，系统库里的 Seurat 5.4 未被改动")
  )
  payload$files <- payload$files[!vapply(payload$files, is.null, logical(1))]

  result_path <- .write_result(file.path(out_dir, "result.json"), payload)
  .log("完成：%s", result_path)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 入口：统一错误处理 / 退出码
# ---------------------------------------------------------------------------
exit_code <- 0L
tryCatch(
  {
    .main()
    exit_code <- 0L
  },
  dsh_prereq = function(condition) {
    message(sprintf("前置条件不满足：%s", conditionMessage(condition)))
    exit_code <<- 3L
  },
  error = function(condition) {
    message(sprintf("Seurat 6.0 适配器失败：%s", conditionMessage(condition)))
    exit_code <<- 1L
  }
)
quit(save = "no", status = exit_code)
