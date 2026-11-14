#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell-workbench — 外部 R 适配器公共库（CONTRACT §10）
#
# 所有 lib/external/*.R 适配器共用本文件：
#   1. CLI 解析：--bundle DIR --out DIR --params '<json>'
#   2. bundle 读取（§10.2：counts.mtx/data.mtx 是 genes × cells，R 侧转成 cells × genes）
#   3. 统一的错误语义：退出码 3 = 缺前置（中文 stderr）；其它非 0 = 运行失败
#   4. result.json 组装（§10.3 严格键）+ obs_extra.csv 写出
#   5. PNG 白底 dpi=144、图上文字英文的绘图助手
#   6. 规模上限与确定性（固定 seed）
#
# 适配器写法：
#   .sc_dir <- local({
#     a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")][1L])
#     if (is.na(f)) "." else dirname(normalizePath(f))
#   })
#   source(file.path(.sc_dir, "_common.R"))
#   sc_run("cellchat", "CellChat 官方细胞通讯", "CellChat 1.6.1（官方 R 包）",
#          function(ctx) { ... return(list(summary = ..., tables = ..., figures = ...)) })
#
# 约束：只用 base R + Matrix + jsonlite（jsonlite 2.x 系统已装）；其它包按需
#       requireNamespace 并在缺失时抛 sc_prerequisite（退出码 3），不硬崩。
# =============================================================================

options(warn = 1)
options(encoding = "UTF-8")
options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 0. 基础工具
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) == 1 && is.na(a)) return(b)
  a
}

sc_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- args[startsWith(args, "--file=")]
  if (length(hit) == 0) return(getwd())
  # Rscript 前端会把 --file 路径里的空格编码成 ~+~
  path <- gsub("~+~", " ", sub("^--file=", "", hit[[1L]]), fixed = TRUE)
  dirname(normalizePath(path, mustWork = FALSE))
}

sc_log <- function(...) {
  cat(paste0(..., collapse = ""), "\n", sep = "", file = stderr())
}

# 退出码：2 参数/内部错误；3 缺前置（由 sc_prerequisite 触发）
sc_die <- function(msg, code = 2L) {
  sc_log("[sc-external] ERROR: ", msg)
  quit(save = "no", status = code, runLast = FALSE)
}

sc_prerequisite <- function(msg) {
  stop(structure(class = c("sc_prerequisite", "error", "condition"),
                 list(message = msg, call = NULL)))
}

sc_is_prereq <- function(e) inherits(e, "sc_prerequisite")

sc_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "unknown")
}

sc_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    sc_prerequisite(sprintf("缺少 R 包：%s（Rscript 环境未安装；请在 R 里 install.packages/BiocManager::install 后重试）",
                            paste(missing, collapse = "、")))
  }
  invisible(TRUE)
}

sc_have <- function(pkgs) vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)

# ---------------------------------------------------------------------------
# 1. CLI / 上下文
# ---------------------------------------------------------------------------

sc_flag <- function(argv, name) {
  eq <- paste0("--", name, "=")
  hit <- which(startsWith(argv, eq))
  if (length(hit) > 0) return(substring(argv[hit[1L]], nchar(eq) + 1L))
  hit <- which(argv == paste0("--", name))
  if (length(hit) > 0) {
    idx <- hit[1L]
    if (idx < length(argv) && !startsWith(argv[idx + 1L], "--")) return(argv[idx + 1L])
  }
  NULL
}

sc_init <- function(tool) {
  argv <- commandArgs(trailingOnly = TRUE)
  bundle <- sc_flag(argv, "bundle")
  out <- sc_flag(argv, "out")
  params_raw <- sc_flag(argv, "params") %||% "{}"
  if (is.null(bundle) || is.null(out)) {
    sc_die(sprintf("[%s] 用法：Rscript --vanilla %s.R --bundle <dir> --out <dir> --params '<json>'", tool, tool), 2L)
  }
  if (!dir.exists(bundle)) sc_die(sprintf("[%s] bundle 目录不存在：%s", tool, bundle), 2L)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  params <- tryCatch(jsonlite::fromJSON(params_raw, simplifyVector = FALSE),
                     error = function(e) {
                       sc_die(sprintf("[%s] --params 不是合法 JSON：%s", tool, conditionMessage(e)), 2L)
                     })
  if (!is.list(params)) params <- list()

  seed <- suppressWarnings(as.integer(params$seed %||% 0L))
  if (is.na(seed)) seed <- 0L
  set.seed(seed)

  ctx <- list(tool = tool, bundle = normalizePath(bundle), out = normalizePath(out),
              params = params, seed = seed, warnings = character(), notes = character())
  sc_log(sprintf("[%s] R %s.%s | bundle=%s", tool, R.version$major, R.version$minor, ctx$bundle))
  ctx
}

sc_warn <- function(ctx, msg) {
  ctx$warnings <- c(ctx$warnings, as.character(msg))
  sc_log("[warn] ", msg)
  invisible(ctx)
}

# 注意：R 值语义 —— 适配器里请用 ctx <- sc_warn(ctx, "…") 或直接收集 warnings 向量。
sc_note <- function(ctx, msg) {
  ctx$notes <- c(ctx$notes, as.character(msg))
  invisible(ctx)
}

# ---------------------------------------------------------------------------
# 2. bundle 读取（§10.2）
# ---------------------------------------------------------------------------

sc_read_mtx <- function(path, label) {
  if (!file.exists(path)) sc_prerequisite(sprintf("bundle 缺少 %s（%s）", label, basename(path)))
  m <- tryCatch(Matrix::readMM(path), error = function(e) {
    sc_prerequisite(sprintf("读取 %s 失败：%s", basename(path), conditionMessage(e)))
  })
  # genes × cells -> cells × genes
  m <- methods::as(Matrix::t(m), "CsparseMatrix")
  m
}

sc_read_tsv <- function(path, label) {
  if (!file.exists(path)) sc_prerequisite(sprintf("bundle 缺少 %s：%s", label, basename(path)))
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- x[nzchar(trimws(x))]
  if (length(x) == 0) sc_prerequisite(sprintf("%s 为空：%s", label, basename(path)))
  x
}

sc_lognormalize <- function(m, target = 1e4) {
  if (inherits(m, "sparseMatrix")) {
    libs <- Matrix::rowSums(m)
    scale <- target / pmax(libs, 1e-9)
    scale[!is.finite(scale)] <- 0
    out <- methods::as(Matrix::Diagonal(x = scale) %*% m, "CsparseMatrix")
    out@x <- log1p(out@x)
    dimnames(out) <- dimnames(m)
    return(out)
  }
  m <- as.matrix(m)
  libs <- rowSums(m)
  m <- m / pmax(libs, 1e-9) * target
  log1p(m)
}

# 基因检出比例（稀疏友好；用于「表达基因」过滤）
sc_pct_expressed <- function(m) {
  if (inherits(m, "sparseMatrix")) {
    as.numeric(Matrix::colSums(m > 0)) / max(nrow(m), 1)
  } else {
    colMeans(as.matrix(m) > 0)
  }
}

sc_gene_variance <- function(m) {
  if (inherits(m, "sparseMatrix")) {
    n <- nrow(m)
    mu <- as.numeric(Matrix::colSums(m)) / max(n, 1)
    mu2 <- as.numeric(Matrix::colSums(m * m)) / max(n, 1)
    var <- mu2 - mu^2
  } else {
    var <- apply(as.matrix(m), 2, stats::var)
  }
  var[!is.finite(var)] <- 0
  pmax(var, 0)
}

sc_cell_variance_genes <- function(m) {
  # 行（细胞）方向不用于选基因；保留给 WGCNA 需要的高变基因选择（列方向）
  sc_gene_variance(m)
}

# 读 bundle；返回 list(data, counts, genes, cells, obs, pca, embed, meta, warnings, X_kind)
#   need_data=TRUE 时保证 data 是 log-normalized（若 meta X_kind=counts 则从 counts 现场标准化）
#   max_cells/max_genes 超限确定性下采样，并写入 warnings
sc_load <- function(ctx, need_counts = FALSE, max_cells = NULL, max_genes = NULL) {
  b <- ctx$bundle
  meta_path <- file.path(b, "meta.json")
  meta <- if (file.exists(meta_path)) {
    tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = TRUE), error = function(e) list())
  } else list()

  cells <- sc_read_tsv(file.path(b, "barcodes.tsv"), "barcodes.tsv")
  genes <- sc_read_tsv(file.path(b, "features.tsv"), "features.tsv")

  has_data <- file.exists(file.path(b, "data.mtx"))
  has_counts <- file.exists(file.path(b, "counts.mtx"))
  if (!has_data && !has_counts) sc_prerequisite("bundle 缺少表达矩阵（data.mtx / counts.mtx）")

  data <- NULL
  counts <- NULL
  if (has_data) data <- sc_read_mtx(file.path(b, "data.mtx"), "data.mtx")
  if (has_counts) counts <- sc_read_mtx(file.path(b, "counts.mtx"), "counts.mtx")

  for (nm in c("data", "counts")) {
    m <- get(nm)
    if (is.null(m)) next
    if (nrow(m) != length(cells)) {
      sc_prerequisite(sprintf("%s 维度与 barcodes.tsv 不一致（矩阵行=%d，细胞=%d）",
                              nm, nrow(m), length(cells)))
    }
    if (ncol(m) != length(genes)) {
      sc_prerequisite(sprintf("%s 维度与 features.tsv 不一致（矩阵列=%d，基因=%d）",
                              nm, ncol(m), length(genes)))
    }
  }
  if (is.null(data)) {
    # 只有 counts：为需要 lognorm 的工具现场标准化
    data <- sc_lognormalize(counts)
  }
  x_kind <- as.character(meta$X_kind %||% if (has_data) "lognorm" else "counts")
  if (is.null(counts) && need_counts) {
    # 无原始 counts：用 data 兜底并明确标注（不伪装成 raw counts）
    counts <- data
    ctx$warnings <- c(ctx$warnings, "bundle 无 counts.mtx，含原始计数假设的分析用 data.mtx 兜底")
  }

  dimnames(data) <- list(cells, genes)
  if (!is.null(counts)) dimnames(counts) <- list(cells, genes)

  # ---- 细胞下采样（确定性）
  if (!is.null(max_cells) && length(cells) > max_cells) {
    n_before <- length(cells)
    idx <- sort(sample.int(n_before, size = as.integer(max_cells)))
    data <- data[idx, , drop = FALSE]
    if (!is.null(counts)) counts <- counts[idx, , drop = FALSE]
    cells <- cells[idx]
    ctx$warnings <- c(ctx$warnings, sprintf("%s：为控制耗时，细胞已从 %d 确定性下采样到 %d",
                                            ctx$tool, n_before, length(cells)))
  }

  # ---- 基因筛选（按方差，确定性）
  if (!is.null(max_genes) && length(genes) > max_genes) {
    var <- sc_gene_variance(data)
    idx <- sort(order(-var, seq_along(var))[seq_len(as.integer(max_genes))])
    data <- data[, idx, drop = FALSE]
    if (!is.null(counts)) counts <- counts[, idx, drop = FALSE]
    genes <- genes[idx]
    ctx$warnings <- c(ctx$warnings, sprintf("%s：基因数超过上限，已取 %d 个高变基因（按方差）",
                                            ctx$tool, length(genes)))
  }

  # ---- obs
  obs_path <- file.path(b, "obs.csv")
  if (file.exists(obs_path)) {
    obs <- tryCatch(utils::read.csv(obs_path, header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE, colClasses = "character",
                                    fileEncoding = "UTF-8", na.strings = character()),
                    error = function(e) sc_prerequisite(sprintf("读取 obs.csv 失败：%s", conditionMessage(e))))
    if (ncol(obs) >= 1) {
      m <- match(cells, as.character(obs[[1L]]))
      if (anyNA(m)) {
        sc_prerequisite(sprintf("obs.csv 与 barcodes.tsv 细胞名不一致（缺失 %d 个）", sum(is.na(m))))
      }
      obs <- obs[m, , drop = FALSE]
      rownames(obs) <- cells
      obs <- obs[, -1L, drop = FALSE]
    }
  } else {
    obs <- data.frame(row.names = cells)
  }

  # ---- pca / embed（可选）
  read_coord <- function(name) {
    path <- file.path(b, name)
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(utils::read.csv(path, header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE, fileEncoding = "UTF-8"),
                   error = function(e) NULL)
    if (is.null(df) || ncol(df) < 2) return(NULL)
    m <- match(cells, as.character(df[[1L]]))
    if (anyNA(m)) return(NULL)
    mat <- as.matrix(df[m, -1L, drop = FALSE])
    storage.mode(mat) <- "double"
    mat
  }
  pca <- read_coord("pca.csv")
  embed <- read_coord("embed.csv")

  list(data = data, counts = counts, genes = genes, cells = cells, obs = obs,
       pca = pca, embed = embed, meta = meta, X_kind = x_kind,
       warnings = ctx$warnings)
}

# ---------------------------------------------------------------------------
# 3. obs / 分组
# ---------------------------------------------------------------------------

sc_obs_chr <- function(obs, col) {
  if (is.null(col) || !(col %in% names(obs))) return(NULL)
  as.character(obs[[col]])
}

sc_obs_num <- function(obs, col) {
  if (is.null(col) || !(col %in% names(obs))) return(NULL)
  suppressWarnings(as.numeric(obs[[col]]))
}

sc_pick_column <- function(obs, candidates, what = "分组列", required = TRUE) {
  candidates <- as.character(candidates %||% character())
  candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
  for (candidate in candidates) {
    if (candidate %in% names(obs)) return(candidate)
  }
  if (!required) return(NULL)
  sc_prerequisite(sprintf("obs 里找不到%s（候选：%s）；现有列：%s",
                          what, paste(candidates, collapse = " / "),
                          paste(names(obs), collapse = ", ")))
}

# 返回 list(labels, levels, keep_mask)
sc_groups <- function(obs, col, drop_na = TRUE) {
  labels <- sc_obs_chr(obs, col)
  if (is.null(labels)) sc_prerequisite(sprintf("obs 缺少分组列 %s", col))
  if (drop_na) labels[is.na(labels) | !nzchar(labels)] <- NA_character_
  levels <- sort(unique(labels[!is.na(labels)]), method = "radix")
  list(labels = labels, levels = levels)
}

sc_drop_small_groups <- function(ctx, labels, min_cells) {
  min_cells <- as.integer(min_cells %||% 3)
  if (is.null(min_cells) || is.na(min_cells) || min_cells < 1) min_cells <- 1L
  tab <- table(labels, useNA = "no")
  keep_levels <- names(tab)[tab >= min_cells]
  dropped <- names(tab)[tab < min_cells]
  if (length(dropped) > 0) {
    ctx$warnings <- c(ctx$warnings, sprintf("%s：分组 %s 细胞数 < %d，已剔除",
                                            ctx$tool, paste(dropped, collapse = ", "), min_cells))
    sc_log("[info] 剔除小分组：", paste(sprintf("%s(%d)", dropped, as.integer(tab[dropped])), collapse = ", "))
  }
  keep <- labels %in% keep_levels
  list(labels = labels, keep = keep, levels = sort(keep_levels, method = "radix"), tab = tab,
       dropped = dropped, ctx = ctx)
}

# ---------------------------------------------------------------------------
# 4. 表 / 图 / 结果
# ---------------------------------------------------------------------------

sc_round <- function(x, digits = 6) {
  if (is.numeric(x)) {
    out <- round(as.numeric(x), digits)
    out[!is.finite(out)] <- NA_real_
    return(out)
  }
  x
}

.sc_json_scalar <- function(x) {
  if (length(x) == 0) return(NULL)
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x)) {
    x <- as.numeric(x)
    if (!is.finite(x)) return(NULL)
    return(x)
  }
  if (is.logical(x)) {
    if (is.na(x)) return(NULL)
    return(x)
  }
  x <- as.character(x)
  if (is.na(x)) return(NULL)
  x
}

# columns: list(list(key="a", label="A", type="string|number"), ...)；df 的列名与 key 对应
sc_table <- function(id, title, df, columns, max_rows = 500, note = NULL) {
  total <- nrow(df)
  if (total > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
    note <- paste0(note %||% "", if (nzchar(note %||% "")) "；" else "",
                   sprintf("共 %d 行，仅展示前 %d 行", total, max_rows))
  }
  rows <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    row <- list()
    for (spec in columns) {
      key <- spec$key
      value <- if (key %in% names(df)) df[[key]][i] else NA
      row[[key]] <- .sc_json_scalar(value)
    }
    rows[[i]] <- row
  }
  list(id = id, title = title, columns = columns, rows = rows, note = note %||% "")
}

sc_cols <- function(...) {
  specs <- list(...)
  lapply(specs, function(spec) {
    list(key = spec[[1L]], label = spec[[2L]], type = spec[[3L]] %||% "string")
  })
}

sc_png <- function(ctx, name, width = 1200, height = 900, res = 144, draw) {
  path <- file.path(ctx$out, paste0(name, ".png"))
  # 设备回退：cairo（Linux）→ quartz（macOS，无需 X11）→ 默认 png
  types <- character()
  if (grepl("darwin", R.version$platform)) {
    types <- c(types, "quartz")
  } else if (isTRUE(capabilities("cairo"))) {
    types <- c(types, "cairo")
  }
  types <- c(types, NA)
  opened <- FALSE
  for (tp in types) {
    okay <- tryCatch({
      if (is.na(tp)) {
        grDevices::png(path, width = width, height = height, res = res, bg = "white")
      } else {
        grDevices::png(path, width = width, height = height, res = res, bg = "white", type = tp)
      }
      TRUE
    }, error = function(e) FALSE)
    if (okay) { opened <- TRUE; break }
  }
  if (!opened) stop("无法打开 PNG 绘图设备（cairo/quartz 均不可用）")
  ok <- FALSE
  tryCatch({
    graphics::par(bg = "white")
    draw()
    ok <- TRUE
  }, finally = {
    try(grDevices::dev.off(), silent = TRUE)
  })
  if (!ok || !file.exists(path) || file.info(path)$size < 800) {
    stop(sprintf("绘图失败：%s", paste0(name, ".png")))
  }
  path
}

sc_figure <- function(id, title, path) {
  list(id = id, title = title, file = normalizePath(path, mustWork = FALSE))
}

sc_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (n <= 2) return(c("#4E79A7", "#E15759")[seq_len(n)])
  grDevices::hcl.colors(n, "Dark 3")
}

sc_heat_colors <- function(n = 100) {
  grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(n)
}

# 通用 base-R 热图（白底、英文标签；不依赖 ComplexHeatmap）
#   mat: rows = 行标签（如通路），cols = 列标签（如簇）
sc_heatmap <- function(mat, row_labels = rownames(mat), col_labels = colnames(mat),
                       title = "", value_label = "score", cex = 0.75,
                       row_title = "", col_title = "") {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0 || ncol(mat) == 0) stop("热图矩阵为空")
  nr <- nrow(mat); nc <- ncol(mat)
  zlim <- range(mat, na.rm = TRUE)
  if (!all(is.finite(zlim)) || diff(zlim) == 0) zlim <- c(zlim[1] - 1, zlim[1] + 1)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(op), silent = TRUE), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4.6, 1))
  graphics::par(bg = "white", mar = c(6, 19, 3.2, 1.2), cex = cex)
  graphics::image(t(mat[nr:1, , drop = FALSE]), col = sc_heat_colors(120),
                  axes = FALSE, zlim = zlim)
  graphics::axis(1, at = seq(0, 1, length.out = nc), labels = col_labels,
                 las = 2, cex.axis = cex, tick = FALSE, line = -0.6)
  graphics::axis(2, at = seq(0, 1, length.out = nr), labels = rev(row_labels),
                 las = 2, cex.axis = cex, tick = FALSE, line = -0.6)
  graphics::title(main = title, cex.main = 1.05, line = 1.4)
  graphics::mtext(col_title, side = 1, line = 4.6, cex = 0.8)
  graphics::mtext(row_title, side = 2, line = 17.5, cex = 0.8)
  graphics::par(mar = c(6, 2.2, 3.2, 3.2))
  legend_vals <- seq(zlim[1], zlim[2], length.out = 121)
  graphics::image(matrix(legend_vals, ncol = 1), col = sc_heat_colors(120),
                  axes = FALSE, zlim = zlim)
  tick_vals <- seq(zlim[1], zlim[2], length.out = 5)
  graphics::axis(4, at = seq(0, 1, length.out = 5),
                 labels = formatC(tick_vals, format = "g", digits = 3),
                 las = 1, cex.axis = 0.7)
  graphics::mtext(value_label, side = 3, line = 0.6, cex = 0.72)
  invisible(NULL)
}

# 基因集读取
sc_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  sets <- list()
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 3) next
    genes <- parts[-c(1L, 2L)]
    genes <- genes[nzchar(genes)]
    if (length(genes) > 0) sets[[parts[1L]]] <- unique(genes)
  }
  sets[!duplicated(names(sets))]
}

# 选基因集：params.gmt_file 优先，否则内置精选集
sc_gene_sets <- function(ctx, genes, max_sets = NULL, min_size = 5L, max_size = 500L) {
  gmt <- as.character(ctx$params$gmt_file %||% "")
  source_label <- "dsh-singlecell-workbench 内置精选基因集（非 MSigDB 全量）"
  sets <- NULL
  if (nzchar(gmt)) {
    gmt_path <- path.expand(gmt)
    if (!file.exists(gmt_path)) sc_prerequisite(sprintf("params.gmt_file 不存在：%s", gmt))
    sets <- sc_read_gmt(gmt_path)
    source_label <- sprintf("用户 GMT：%s", basename(gmt_path))
    if (length(sets) == 0) sc_prerequisite(sprintf("GMT 文件里没有可用基因集：%s", gmt_path))
  } else {
    sets <- get("SC_BUILTIN_GENE_SETS", envir = globalenv())
  }
  genes <- unique(as.character(genes))
  sets <- lapply(sets, function(gs) intersect(unique(as.character(gs)), genes))
  sizes <- vapply(sets, length, integer(1))
  keep <- sizes >= min_size & sizes <= max_size
  sets <- sets[keep]
  if (length(sets) == 0) {
    sc_prerequisite("没有基因集与数据基因名有足够重叠（>=5 个基因）；请提供 params.gmt_file 或使用真实全转录组数据")
  }
  if (!is.null(max_sets) && length(sets) > max_sets) {
    # 按重叠基因数降序（确定性），再按名字稳定排序
    sizes2 <- vapply(sets, length, integer(1))
    ord <- order(-sizes2, names(sets))
    sets <- sets[ord][seq_len(max_sets)]
  }
  list(sets = sets, source = source_label,
       sizes = vapply(sets, length, integer(1)))
}

# ---------------------------------------------------------------------------
# 5. 写 result.json / obs_extra.csv
# ---------------------------------------------------------------------------

sc_write_result <- function(ctx, res) {
  res$title <- as.character(res$title %||% ctx$tool)
  res$method <- as.character(res$method %||% "官方 R 实现")
  summary <- res$summary %||% list()
  bullets <- as.character(summary$bullets %||% character())
  metrics <- summary$metrics
  if (is.null(metrics) || length(metrics) == 0) {
    metrics <- structure(list(), names = character(0))
  }
  res$summary <- list(bullets = I(bullets), metrics = metrics)
  res$tables <- res$tables %||% list()
  res$figures <- res$figures %||% list()
  warnings <- as.character(res$warnings %||% character())
  warnings <- c(warnings, ctx$warnings)
  warnings <- unique(warnings[nzchar(warnings)])
  notes <- as.character(res$notes %||% character())
  notes <- c(notes, ctx$notes)
  notes <- unique(notes[nzchar(notes)])
  res$warnings <- I(warnings)
  res$notes <- I(notes)
  res <- res[c("title", "method", "summary", "tables", "figures", "warnings", "notes")]

  json <- jsonlite::toJSON(res, auto_unbox = TRUE, null = "null", na = "null",
                           digits = NA, pretty = TRUE, force = TRUE)
  final <- file.path(ctx$out, "result.json")
  tmp <- paste0(final, ".tmp")
  con <- file(tmp, open = "wt", encoding = "UTF-8")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  writeLines(as.character(json), con, useBytes = TRUE)
  close(con)
  on.exit(NULL, add = FALSE)
  if (!file.rename(tmp, final)) {
    file.copy(tmp, final, overwrite = TRUE)
    unlink(tmp)
  }
  sc_log(sprintf("[%s] 写出 %s（%d tables, %d figures）", ctx$tool, final,
                 length(res$tables), length(res$figures)))
  invisible(final)
}

sc_write_obs <- function(ctx, columns, cells = NULL) {
  if (is.null(columns) || length(columns) == 0) return(NULL)
  cells <- cells %||% ctx$cells
  if (is.null(cells) || length(cells) == 0) {
    stop("sc_write_obs 需要细胞名向量（cells）")
  }
  n <- length(cells)
  df <- data.frame(cell = cells, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(columns)) {
    v <- columns[[nm]]
    if (length(v) != n) {
      stop(sprintf("obs_extra 列 %s 长度 %d 与细胞数 %d 不一致", nm, length(v), n))
    }
    if (is.factor(v)) v <- as.character(v)
    if (is.character(v)) v <- gsub("[,\r\n]", " ", v)
    if (is.numeric(v)) {
      v[!is.finite(v)] <- NA
      # 保留多位小数，不用科学计数法
      v <- signif(v, 8)
    }
    df[[nm]] <- v
  }
  path <- file.path(ctx$out, "obs_extra.csv")
  utils::write.csv(df, path, row.names = FALSE, quote = FALSE, na = "",
                   fileEncoding = "UTF-8")
  sc_log(sprintf("[%s] 写出 %s（%d 列 × %d 细胞）", ctx$tool, path,
                 length(columns), nrow(df)))
  path
}

# ---------------------------------------------------------------------------
# 6. 主流程包装
# ---------------------------------------------------------------------------

sc_result <- function(ctx, bullets = character(), metrics = NULL) {
  list(title = NULL, method = NULL,
       summary = list(bullets = as.character(bullets),
                      metrics = metrics %||% structure(list(), names = character(0))),
       tables = list(), figures = list(),
       warnings = character(), notes = character())
}

# fn(ctx) 返回 list(title?, method?, summary?, tables?, figures?, warnings?, notes?) 或空
sc_run <- function(tool, title, method, fn) {
  ctx <- tryCatch(sc_init(tool), error = function(e) {
    sc_die(sprintf("[%s] 初始化失败：%s", tool, conditionMessage(e)), 2L)
  })
  fn_ctx <- ctx
  debug <- identical(Sys.getenv("DS_SC_DEBUG"), "1")
  run_fn <- function() {
    if (!debug) return(fn(fn_ctx))
    # 调试模式：error 发生时立刻捕获完整调用栈
    withCallingHandlers(
      fn(fn_ctx),
      error = function(e) {
        calls <- sys.calls()
        txt <- vapply(calls, function(cl) paste(deparse(cl, nlines = 1), collapse = ""), character(1))
        sc_log("[debug] 出错调用栈：\n", paste(utils::tail(txt, 30), collapse = "\n"))
      }
    )
  }
  result <- tryCatch(
    run_fn(),
    sc_prerequisite = function(e) {
      sc_log(sprintf("[%s] 前置条件不满足：%s", tool, conditionMessage(e)))
      quit(save = "no", status = 3L, runLast = FALSE)
    },
    error = function(e) {
      sc_log(sprintf("[%s] 运行失败：%s", tool, conditionMessage(e)))
      if (identical(Sys.getenv("DS_SC_DEBUG"), "1")) {
        calls <- sys.calls()
        txt <- vapply(calls, function(cl) paste(deparse(cl, nlines = 1), collapse = ""), character(1))
        sc_log("[debug] 调用栈：\n", paste(utils::tail(txt, 12), collapse = "\n"))
      }
      quit(save = "no", status = 1L, runLast = FALSE)
    }
  )
  if (!is.list(result)) result <- list()
  result$title <- result$title %||% title
  result$method <- result$method %||% method
  result$warnings <- c(result$warnings %||% character(), fn_ctx$warnings)
  result$notes <- c(result$notes %||% character(), fn_ctx$notes)
  result <- sc_write_result(fn_ctx, result)
  invisible(result)
}

# ---------------------------------------------------------------------------
# 7. 加载内置基因集（生成的 R 文件）
# ---------------------------------------------------------------------------

.sc_gene_sets_file <- file.path(sc_script_dir(), "_gene_sets.R")
if (file.exists(.sc_gene_sets_file)) {
  source(.sc_gene_sets_file, local = FALSE)
}

sc_log("[sc-common] 公共库已加载（R ", R.version$major, ".", R.version$minor, "）")
