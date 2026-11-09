#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# dsh-singlecell — R bridge (CONTRACT §8.3)
#
# 把 R 的单细胞对象（.rds / .RData / .rda）导出成 Python 侧能直接读的中立文件：
#   <out>/matrix.mtx     MatrixMarket，cells × genes（orientation 见 meta.json）
#   <out>/features.tsv   每行一个基因名
#   <out>/barcodes.tsv   每行一个细胞名
#   <out>/obs.csv        row.names = 细胞名，列 = 细胞元数据
#   <out>/meta.json      {status,class,assay,assays,layers,orientation,dim,n_obs,n_vars,
#                         has_counts,warnings,notes[,candidates]}
#
# 支持对象：Seurat v3/v4/v5、SingleCellExperiment / SummarizedExperiment、
#           matrix / Matrix / dgCMatrix、data.frame（数值列）。
# 多对象歧义：meta.status="ambiguous" + candidates（object/class/dim/label），
#            Python 侧转成 EngineError("ambiguous")；--object 可指定对象名。
#
# 调用（由 lib/sc_io.py 负责）：
#   Rscript --vanilla r_bridge.R --in <file> --out <dir> [--object <name>] [--inspect]
#
# 约束：只用 base R + 已安装的 Seurat/SeuratObject/SingleCellExperiment/
#       SummarizedExperiment/Matrix；不依赖 jsonlite（手写 JSON，控制转义）。
#       所有文本输出 UTF-8、quote=FALSE（obs 里会破坏 CSV 的逗号/换行先替换成空格）。
# =============================================================================

options(warn = 1)
options(encoding = "UTF-8")
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || (is.character(a) && length(a) == 1 && is.na(a))) b else a

# ---------------------------------------------------------------------------
# 0. CLI
# ---------------------------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)

flag <- function(name, default = NULL) {
  eq <- paste0("--", name, "=")
  hit <- which(startsWith(argv, eq))
  if (length(hit) > 0) return(substring(argv[hit[1]], nchar(eq) + 1L))
  hit <- which(argv == paste0("--", name))
  if (length(hit) > 0) {
    idx <- hit[1]
    if (idx < length(argv) && !startsWith(argv[idx + 1], "--")) return(argv[idx + 1])
    return("1")
  }
  default
}

IN <- flag("in")
OUT <- flag("out")
OBJ <- flag("object")
INSPECT <- !is.null(flag("inspect"))

if (is.null(IN) || is.null(OUT)) {
  cat("[r-bridge] 缺少参数：--in <file> --out <dir>\n", file = stderr())
  quit(status = 2L)
}

if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. JSON 输出（手写；不依赖 jsonlite）
# ---------------------------------------------------------------------------

json_escape <- function(s) {
  s <- enc2utf8(as.character(s))
  if (is.na(s)) return("null")
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  paste0("\"", s, "\"")
}

json_num <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1 || is.na(x) || !is.finite(x)) return("null")
  if (x == round(x) && abs(x) < 1e15) return(format(x, scientific = FALSE, trim = TRUE))
  format(x, scientific = FALSE, digits = 15, trim = TRUE)
}

json_value <- function(x) {
  if (is.null(x)) return("null")
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms) && length(nms) == length(x) && all(nzchar(nms))) {
      keep <- !vapply(x, is.null, logical(1))
      x <- x[keep]
      nms <- nms[keep]
      parts <- vapply(seq_along(x),
                      function(i) paste0(json_escape(nms[i]), ":", json_value(x[[i]])),
                      character(1))
      return(paste0("{", paste(parts, collapse = ","), "}"))
    }
    parts <- vapply(x, json_value, character(1))
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    if (length(x) == 0) return("[]")
    parts <- vapply(x, json_escape, character(1))
    if (length(parts) == 1) return(parts)
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.logical(x)) {
    if (length(x) == 0) return("[]")
    parts <- ifelse(is.na(x), "null", ifelse(x, "true", "false"))
    if (length(parts) == 1) return(parts)
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.numeric(x)) {
    if (length(x) == 0) return("[]")
    parts <- vapply(x, json_num, character(1))
    if (length(parts) == 1) return(parts)
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  json_escape(paste(as.character(x), collapse = " "))
}

write_meta <- function(meta, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(json_value(meta), con, useBytes = TRUE)
}

bridge_error <- function(msg) {
  cat(sprintf("[r-bridge] ERROR: %s\n", msg), file = stderr())
  try(write_meta(list(status = "error", error = msg), file.path(OUT, "meta.json")), silent = TRUE)
  quit(save = "no", status = 1L, runLast = FALSE)
}

# ---------------------------------------------------------------------------
# 2. 对象识别工具
# ---------------------------------------------------------------------------

class_label <- function(x) paste(class(x), collapse = ",")

is_supported_object <- function(x) {
  inherits(x, "Seurat") ||
    inherits(x, "SingleCellExperiment") ||
    inherits(x, "SummarizedExperiment") ||
    is.matrix(x) ||
    inherits(x, "Matrix") ||
    inherits(x, "DelayedMatrix") ||
    is.data.frame(x)
}

is_high_level <- function(x) {
  inherits(x, "Seurat") || inherits(x, "SingleCellExperiment") || inherits(x, "SummarizedExperiment")
}

obj_dim <- function(x) {
  d <- tryCatch(dim(x), error = function(e) NULL)
  if (is.null(d) || length(d) != 2L) return(NULL)
  as.integer(d)
}

# 统一表达成 c(cells, genes)；data.frame 按「行=细胞」估计（只是候选提示用）。
cell_gene_dim <- function(x) {
  d <- obj_dim(x)
  if (is.null(d)) return(NULL)
  if (is.data.frame(x)) return(c(cells = d[1], genes = d[2]))
  c(cells = d[2], genes = d[1])
}

candidate_entry <- function(name, x) {
  d <- cell_gene_dim(x)
  edim <- if (is.null(d)) list(NA_integer_, NA_integer_) else list(unname(d["cells"]), unname(d["genes"]))
  lab <- if (is.null(d)) {
    sprintf("%s · %s · 维度未知", name, class_label(x))
  } else {
    sprintf("%s · %s · %d 细胞 × %d 基因", name, class_label(x), d["cells"], d["genes"])
  }
  list(object = name, class = class_label(x), dim = edim, label = lab)
}

is_barcode_vec <- function(x) {
  x <- as.character(x)
  grepl("^[ACGTN]{10,}", x) | grepl("-[0-9]+$", x)
}

looks_barcode <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(FALSE)
  mean(is_barcode_vec(x), na.rm = TRUE) >= 0.6
}

looks_gene_symbol <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(FALSE)
  mask <- grepl("^[A-Za-z][A-Za-z0-9._-]{0,40}$", x) &
    !is_barcode_vec(x) &
    !grepl("^(GSM|GSE|SRR|ERR|DRR|SAM|SAMPLE)[-_]?[0-9]+", x, ignore.case = TRUE)
  mean(mask, na.rm = TRUE) >= 0.6
}

looks_sample_id <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) return(FALSE)
  mask <- grepl("^(GSM|GSE|SRR|ERR|DRR|SAM|SAMPLE)[-_]?[0-9]+", x, ignore.case = TRUE) |
    grepl("^[0-9]+$", x)
  mean(mask, na.rm = TRUE) >= 0.6
}

# 返回当前矩阵方向："genes_cells" 或 "cells_genes"。
guess_orientation <- function(rows, cols, default = "genes_cells") {
  if (is.null(rows) || is.null(cols) || length(rows) == 0 || length(cols) == 0) return(default)
  rb <- looks_barcode(rows)
  cb <- looks_barcode(cols)
  if (rb && !cb) return("cells_genes")
  if (cb && !rb) return("genes_cells")
  rg <- looks_gene_symbol(rows)
  cg <- looks_gene_symbol(cols)
  if (rg && !cg) return("genes_cells")
  if (cg && !rg) return("cells_genes")
  rs <- looks_sample_id(rows)
  cs <- looks_sample_id(cols)
  if (rs && !cs) return("cells_genes")
  if (cs && !rs) return("genes_cells")
  default
}

is_default_names <- function(x, n) {
  if (is.null(x) || length(x) != n) return(TRUE)
  x <- as.character(x)
  identical(x, as.character(seq_len(n)))
}

axis_names <- function(x, n, prefix, warnings) {
  bad <- is.null(x) || length(x) != n || any(!nzchar(as.character(x))) || is_default_names(x, n)
  if (bad) {
    warnings <- c(warnings, sprintf("%s 名缺失或无效，已自动生成 %s1..%sN", prefix, prefix, prefix))
    return(list(names = paste0(prefix, seq_len(n)), warnings = warnings))
  }
  list(names = enc2utf8(as.character(x)), warnings = warnings)
}

# ---------------------------------------------------------------------------
# 3. 读入 .rds / .RData
# ---------------------------------------------------------------------------

read_any_r_object <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    return(list(kind = "rds", objects = list(object = obj)))
  }
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (length(loaded) == 0) return(list(kind = "rdata", objects = list()))
  objs <- lapply(loaded, function(nm) get(nm, envir = env))
  names(objs) <- loaded
  list(kind = "rdata", objects = objs)
}

# rds 里常见 list(seurat_a=..., seurat_b=...) 的打包方式，摊平成候选对象。
# 只展开「纯 list」（class(x) 恰好是 list），避免把 lm 之类的有类 list 拆坏。
expand_objects <- function(objs) {
  out <- list()
  for (nm in names(objs)) {
    x <- objs[[nm]]
    if (!is_supported_object(x) && is.list(x) && identical(class(x), "list")) {
      inner <- names(x)
      if (is.null(inner) || any(!nzchar(inner))) inner <- paste0(nm, "_", seq_along(x))
      supported_inner <- vapply(x, is_supported_object, logical(1))
      if (any(supported_inner)) {
        for (i in which(supported_inner)) out[[inner[i]]] <- x[[i]]
        next
      }
    }
    out[[nm]] <- x
  }
  out
}

# 返回 list(status="ok", object=x) 或 list(status="ambiguous", candidates=[...])
select_object <- function(objs) {
  names_all <- names(objs)
  if (length(names_all) == 0) bridge_error("R 文件里没有找到任何对象")
  supported <- names_all[vapply(objs, is_supported_object, logical(1))]

  if (!is.null(OBJ) && nzchar(OBJ) && OBJ != "1") {
    if (!(OBJ %in% names_all)) {
      bridge_error(sprintf("对象 %s 不存在；可用对象：%s", OBJ, paste(names_all, collapse = ", ")))
    }
    if (!is_supported_object(objs[[OBJ]])) {
      bridge_error(sprintf("对象 %s 的类型 %s 不支持（支持 Seurat/SCE/matrix/Matrix/data.frame）",
                           OBJ, class_label(objs[[OBJ]])))
    }
    return(list(status = "ok", object = objs[[OBJ]]))
  }

  if (length(supported) == 0) {
    kinds <- paste(sprintf("%s(%s)", names_all, vapply(objs, class_label, character(1))), collapse = ", ")
    bridge_error(sprintf("R 文件里没有可识别的单细胞对象（找到：%s）", kinds))
  }

  high <- supported[vapply(objs[supported], is_high_level, logical(1))]
  pool <- if (length(high) > 0) high else supported
  if (length(pool) == 1) return(list(status = "ok", object = objs[[pool]]))

  list(status = "ambiguous",
       candidates = lapply(pool, function(nm) candidate_entry(nm, objs[[nm]])))
}

# ---------------------------------------------------------------------------
# 4. 提取：统一返回 genes × cells 的稀疏矩阵 + obs
# ---------------------------------------------------------------------------

seurat_assay_layers <- function(obj, assay) {
  layers <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
  if (length(layers) == 0) {
    # v3/v4 的 Assay 没有 Layers()，退化成固定槽名
    layers <- c("counts", "data", "scale.data")
  }
  layers
}

seurat_get_layer <- function(obj, assay, layer) {
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      "LayerData" %in% getNamespaceExports("SeuratObject")) {
    m <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = layer),
                  error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  }
  m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = layer),
                error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  m <- tryCatch({
    a <- obj@assays[[assay]]
    switch(layer,
           counts = a@counts,
           data = a@data,
           scale.data = a@scale.data,
           NULL)
  }, error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) return(m)
  NULL
}

extract_seurat <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  assays <- tryCatch(names(obj@assays), error = function(e) character(0))
  assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
  if (is.null(assay) || !nzchar(assay) || !(assay %in% assays)) {
    assay <- if (length(assays) > 0) assays[1] else NULL
  }
  if (is.null(assay)) bridge_error("Seurat 对象里没有 assay")
  if (length(assays) > 1) {
    notes <- c(notes, sprintf("Seurat 对象含 %d 个 assay（%s），已导出默认 assay %s",
                              length(assays), paste(assays, collapse = ","), assay))
  }
  layers <- seurat_assay_layers(obj, assay)

  m <- seurat_get_layer(obj, assay, "counts")
  has_counts <- !is.null(m)
  layer_used <- if (has_counts) "counts" else "data"
  if (is.null(m)) {
    m <- seurat_get_layer(obj, assay, "data")
    warnings <- c(warnings, "Seurat 对象没有 counts 层，已使用 data 层（X 可能已标准化）")
  }

  # Seurat v5 多 layer：counts.1 / counts.2 ... 按列拼接
  if (is.null(m)) {
    count_layers <- grep("^counts(\\.|$)", layers, value = TRUE)
    if (length(count_layers) == 0) count_layers <- grep("^((data)|(scale\\.data))(\\.|$)", layers, value = TRUE)
    if (length(count_layers) > 0) {
      mats <- lapply(count_layers, function(l) seurat_get_layer(obj, assay, l))
      mats <- mats[!vapply(mats, is.null, logical(1))]
      if (length(mats) > 0) {
        m <- do.call(cbind, mats)
        layer_used <- paste(count_layers, collapse = "+")
        has_counts <- grepl("^counts", layer_used)
        warnings <- c(warnings, sprintf("Seurat v5 多 layer（%s），已按列拼接", paste(count_layers, collapse = ",")))
      }
    }
  }

  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    bridge_error("Seurat 对象里没有可导出的表达矩阵（counts/data 都为空）")
  }

  genes <- rownames(m)
  cells <- colnames(m)
  obs <- tryCatch(obj@meta.data, error = function(e) NULL)

  list(matrix = m, genes = genes, cells = cells, obs = obs,
       class = class_label(obj), assay = assay, assays = assays, layers = layers,
       has_counts = has_counts, layer_used = layer_used,
       warnings = warnings, notes = notes)
}

extract_sce <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  assays <- tryCatch(SummarizedExperiment::assayNames(obj), error = function(e) character(0))
  m <- NULL
  layer_used <- NULL
  has_counts <- FALSE
  if ("counts" %in% assays) {
    m <- tryCatch(SummarizedExperiment::assay(obj, "counts"), error = function(e) NULL,
                  warning = function(w) NULL)
    if (!is.null(m)) {
      has_counts <- TRUE
      layer_used <- "counts"
    }
  }
  if (is.null(m)) {
    first <- if (length(assays) > 0) assays[1] else 1
    m <- tryCatch(SummarizedExperiment::assay(obj, first), error = function(e) NULL,
                  warning = function(w) NULL)
    layer_used <- if (length(assays) > 0) first else "assay1"
    warnings <- c(warnings, "SingleCellExperiment 没有 counts assay，已使用第一个 assay（X 可能已标准化）")
  }
  if (is.null(m) || nrow(m) == 0 || ncol(m) == 0) {
    bridge_error("SingleCellExperiment 对象里没有可导出的 assay")
  }
  if (length(assays) > 1) {
    notes <- c(notes, sprintf("SingleCellExperiment 含 %d 个 assay（%s），已导出 %s",
                              length(assays), paste(assays, collapse = ","), layer_used))
  }
  obs <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)), error = function(e) NULL)
  list(matrix = m, genes = rownames(m), cells = colnames(m), obs = obs,
       class = class_label(obj), assay = layer_used, assays = assays, layers = assays,
       has_counts = has_counts, layer_used = layer_used,
       warnings = warnings, notes = notes)
}

extract_dataframe <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  numeric_cols <- vapply(obj, is.numeric, logical(1))
  dropped <- names(obj)[!numeric_cols]
  if (!any(numeric_cols)) bridge_error("data.frame 里没有数值列，无法作为表达矩阵")
  if (length(dropped) > 0) {
    warnings <- c(warnings, sprintf("data.frame 有 %d 个非数值列已忽略（%s）",
                                    length(dropped), paste(utils::head(dropped, 8), collapse = ",")))
  }
  raw <- as.matrix(obj[, numeric_cols, drop = FALSE])
  rownames(raw) <- rownames(obj)
  if (is.null(colnames(raw))) colnames(raw) <- names(obj)[numeric_cols]
  row_names <- rownames(raw)
  col_names <- colnames(raw)
  orient <- guess_orientation(row_names, col_names, default = "cells_genes")
  warnings <- c(warnings, sprintf("data.frame 方向启发式判定：%s",
                                  if (orient == "cells_genes") "行=细胞（默认）" else "行=基因，已转置"))
  if (orient == "cells_genes") {
    cells <- row_names
    genes <- col_names
    m <- Matrix::t(Matrix::Matrix(raw, sparse = TRUE))   # 统一成 genes × cells
  } else {
    genes <- row_names
    cells <- col_names
    m <- Matrix::Matrix(raw, sparse = TRUE)
  }
  list(matrix = m, genes = genes, cells = cells,
       obs = data.frame(row.names = cells),
       class = class_label(obj), assay = "", assays = character(0), layers = character(0),
       has_counts = TRUE, layer_used = "matrix",
       warnings = c(warnings, "矩阵没有独立的 counts 层，已把唯一矩阵同时作为 counts 与 X"),
       notes = notes)
}

extract_matrix_obj <- function(obj) {
  warnings <- character(0)
  notes <- character(0)
  m <- obj
  if (!is.matrix(m) && is.data.frame(m)) {
    return(extract_dataframe(m))
  }
  if (!is.matrix(m)) {
    m <- tryCatch(methods::as(m, "CsparseMatrix"),
                  error = function(e) tryCatch(as.matrix(m), error = function(e2) NULL))
  }
  if (is.null(m)) bridge_error(sprintf("无法把 %s 转成矩阵", class_label(obj)))
  row_names <- rownames(m)
  col_names <- colnames(m)
  orient <- guess_orientation(row_names, col_names, default = "genes_cells")
  warnings <- c(warnings, sprintf("矩阵方向启发式判定：%s",
                                  if (orient == "genes_cells") "行=基因（默认）" else "行=细胞，已转置"))
  if (is.null(row_names) || is.null(col_names)) {
    warnings <- c(warnings, "矩阵 dimnames 缺失，已自动生成基因名/细胞名")
  }
  if (orient == "cells_genes") {
    cells <- row_names
    genes <- col_names
    m <- Matrix::t(if (is.matrix(m)) Matrix::Matrix(m, sparse = TRUE) else m)  # genes × cells
  } else {
    genes <- row_names
    cells <- col_names
  }
  list(matrix = m, genes = genes, cells = cells,
       obs = data.frame(row.names = cells),
       class = class_label(obj), assay = "", assays = character(0), layers = character(0),
       has_counts = TRUE, layer_used = "matrix",
       warnings = c(warnings, "矩阵没有独立的 counts 层，已把唯一矩阵同时作为 counts 与 X"),
       notes = notes)
}

extract_any <- function(obj) {
  if (inherits(obj, "Seurat")) return(extract_seurat(obj))
  if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) return(extract_sce(obj))
  if (is.data.frame(obj)) return(extract_dataframe(obj))
  if (is.matrix(obj) || inherits(obj, "Matrix") || inherits(obj, "DelayedMatrix")) return(extract_matrix_obj(obj))
  bridge_error(sprintf("不支持的 R 对象类型：%s", class_label(obj)))
}

# 统一补基因名/细胞名、obs 对齐
finalize_extract <- function(result) {
  warnings <- result$warnings %||% character(0)
  n_genes <- nrow(result$matrix)
  n_cells <- ncol(result$matrix)

  genes <- axis_names(result$genes, n_genes, "gene_", warnings)
  warnings <- genes$warnings
  genes <- genes$names

  cells <- axis_names(result$cells, n_cells, "cell_", warnings)
  warnings <- cells$warnings
  cells <- cells$names

  obs <- result$obs
  if (is.null(obs) || nrow(obs) != n_cells) {
    if (!is.null(obs)) warnings <- c(warnings, "obs 行数与矩阵细胞数不一致，已只保留细胞名")
    obs <- data.frame(row.names = cells)
  } else {
    obs <- as.data.frame(obs, stringsAsFactors = FALSE)
    rn <- rownames(obs)
    if (is.null(rn) || any(!nzchar(rn))) {
      rownames(obs) <- cells
    } else if (setequal(rn, cells)) {
      obs <- obs[cells, , drop = FALSE]
    } else {
      warnings <- c(warnings, "obs 行名与细胞名不完全一致，已按顺序对齐")
      rownames(obs) <- cells
    }
  }

  result$matrix <- result$matrix
  result$genes <- genes
  result$cells <- cells
  result$obs <- obs
  result$warnings <- warnings
  result
}

# ---------------------------------------------------------------------------
# 5. 写输出文件
# ---------------------------------------------------------------------------

write_lines_utf8 <- function(x, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(as.character(x)), con, useBytes = TRUE)
}

to_csparse <- function(m) {
  if (is.matrix(m)) m <- Matrix::Matrix(m, sparse = TRUE)
  if (!inherits(m, "sparseMatrix")) m <- methods::as(m, "CsparseMatrix")
  if (inherits(m, "nMatrix") || inherits(m, "lMatrix") || !is.numeric(m@x)) m <- m * 1
  if (!inherits(m, "dgCMatrix")) {
    m <- methods::as(methods::as(m, "CsparseMatrix"), "dgCMatrix")
  }
  m
}

sanitize_matrix <- function(m) {
  warnings <- character(0)
  if (inherits(m, "sparseMatrix")) {
    if (!is.numeric(m@x)) m <- m * 1
    bad <- !is.finite(m@x)
    if (any(bad)) {
      m@x[bad] <- 0
      warnings <- c(warnings, sprintf("矩阵中 %d 个非有限值已替换为 0", sum(bad)))
    }
  } else {
    bad <- !is.finite(m)
    if (any(bad)) {
      m[bad] <- 0
      warnings <- c(warnings, sprintf("矩阵中 %d 个非有限值已替换为 0", sum(bad)))
    }
  }
  list(matrix = m, warnings = warnings)
}

write_obs_csv <- function(obs, cells, path) {
  warnings <- character(0)
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  if (nrow(obs) != length(cells)) obs <- data.frame(row.names = cells)
  rownames(obs) <- cells
  if (ncol(obs) == 0) {
    file.create(path)
    return(warnings)
  }
  # 只保留可写成一列标量的列（list 列、matrix 列丢掉）
  keep <- vapply(obs, function(col) is.atomic(col) || is.factor(col), logical(1))
  dropped <- names(obs)[!keep]
  if (length(dropped) > 0) {
    warnings <- c(warnings, sprintf("obs 有 %d 个非标量列已忽略（%s）",
                                    length(dropped), paste(utils::head(dropped, 8), collapse = ",")))
    obs <- obs[, keep, drop = FALSE]
  }
  for (col in names(obs)) if (is.factor(obs[[col]])) obs[[col]] <- as.character(obs[[col]])
  bad_values <- 0
  for (col in names(obs)) {
    if (is.character(obs[[col]])) {
      hit <- grepl("[,\r\n]", obs[[col]])
      bad_values <- bad_values + sum(hit, na.rm = TRUE)
      obs[[col]] <- gsub("[,\r\n]", " ", obs[[col]])
    }
  }
  if (bad_values > 0) {
    warnings <- c(warnings, sprintf("obs 中 %d 个值含逗号/换行，已替换为空格以保证 CSV 可解析", bad_values))
  }
  names(obs) <- gsub("[,\r\n]", "_", names(obs))
  names(obs) <- make.unique(ifelse(nzchar(names(obs)), names(obs), "obs"))
  utils::write.table(obs, path, sep = ",", row.names = TRUE, col.names = NA,
                     quote = FALSE, na = "", fileEncoding = "UTF-8")
  warnings
}

write_outputs <- function(result, out_dir) {
  warnings <- character(0)
  san <- sanitize_matrix(result$matrix)
  warnings <- c(warnings, san$warnings)
  mm <- to_csparse(san$matrix)
  mm <- methods::as(Matrix::t(mm), "CsparseMatrix")   # cells × genes
  mm <- to_csparse(mm)
  Matrix::writeMM(mm, file.path(out_dir, "matrix.mtx"))

  write_lines_utf8(result$genes, file.path(out_dir, "features.tsv"))
  write_lines_utf8(result$cells, file.path(out_dir, "barcodes.tsv"))
  warnings <- c(warnings, write_obs_csv(result$obs, result$cells, file.path(out_dir, "obs.csv")))
  warnings
}

# ---------------------------------------------------------------------------
# 6. 轻量 inspect（只报 class/dim，不导出矩阵）
# ---------------------------------------------------------------------------

inspect_object <- function(obj) {
  d <- cell_gene_dim(obj)
  assays <- character(0)
  assay <- ""
  layers <- character(0)
  if (inherits(obj, "Seurat")) {
    assays <- tryCatch(names(obj@assays), error = function(e) character(0))
    assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) "") %||% ""
    if (nzchar(assay) && assay %in% assays) layers <- seurat_assay_layers(obj, assay)
  } else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SummarizedExperiment")) {
    assays <- tryCatch(SummarizedExperiment::assayNames(obj), error = function(e) character(0))
    assay <- if (length(assays) > 0) assays[1] else ""
    layers <- assays
  }
  list(class = class_label(obj),
       assay = assay,
       assays = as.list(assays),
       layers = as.list(layers),
       dim = if (is.null(d)) list(NA_integer_, NA_integer_) else list(unname(d["cells"]), unname(d["genes"])),
       n_obs = if (is.null(d)) NA_integer_ else unname(d["cells"]),
       n_vars = if (is.null(d)) NA_integer_ else unname(d["genes"]))
}

# ---------------------------------------------------------------------------
# 7. main
# ---------------------------------------------------------------------------

main <- function() {
  if (!file.exists(IN)) bridge_error(sprintf("文件不存在：%s", IN))
  loaded <- tryCatch(read_any_r_object(IN),
                     error = function(e) bridge_error(sprintf("读取 R 文件失败：%s", conditionMessage(e))))
  objs <- expand_objects(loaded$objects)
  if (length(objs) == 0) bridge_error(sprintf("文件里没有对象：%s", IN))

  sel <- select_object(objs)
  if (identical(sel$status, "ambiguous")) {
    meta <- list(status = "ambiguous",
                 kind = loaded$kind,
                 candidates = sel$candidates,
                 warnings = list("R 文件里存在多个可用的单细胞对象，请用 object 参数指定"),
                 notes = list())
    write_meta(meta, file.path(OUT, "meta.json"))
    quit(save = "no", status = 0L, runLast = FALSE)
  }

  obj <- sel$object
  if (INSPECT) {
    info <- inspect_object(obj)
    meta <- c(list(status = "ok", kind = loaded$kind, inspect = TRUE,
                   orientation = "cells_genes", warnings = list(), notes = list()),
              info)
    write_meta(meta, file.path(OUT, "meta.json"))
    quit(save = "no", status = 0L, runLast = FALSE)
  }

  result <- finalize_extract(extract_any(obj))
  extra_warnings <- write_outputs(result, OUT)
  warnings <- unique(c(result$warnings, extra_warnings))
  meta <- list(
    status = "ok",
    kind = loaded$kind,
    class = result$class,
    assay = result$assay %||% "",
    assays = as.list(result$assays %||% character(0)),
    layers = as.list(result$layers %||% character(0)),
    layer_used = result$layer_used %||% "",
    orientation = "cells_genes",
    dim = list(ncol(result$matrix), nrow(result$matrix)),
    n_obs = ncol(result$matrix),
    n_vars = nrow(result$matrix),
    has_counts = isTRUE(result$has_counts),
    warnings = as.list(warnings),
    notes = as.list(result$notes %||% character(0))
  )
  write_meta(meta, file.path(OUT, "meta.json"))
  quit(save = "no", status = 0L, runLast = FALSE)
}

main()
