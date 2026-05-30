options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(copykat)
})

set.seed(20260527)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pre_root <- file.path(root, "preanalysis")
analysis_root <- file.path(root, "analysis")
geo_dir <- file.path(pre_root, "data", "geo")
result_dir <- file.path(analysis_root, "results", "copykat")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

max_epi <- as.integer(Sys.getenv("COPYKAT_MAX_EPI", "1200"))
max_ref <- as.integer(Sys.getenv("COPYKAT_MAX_REF", "600"))
min_epi <- as.integer(Sys.getenv("COPYKAT_MIN_EPI", "50"))
min_ref <- as.integer(Sys.getenv("COPYKAT_MIN_REF", "100"))
n_cores <- as.integer(Sys.getenv("COPYKAT_N_CORES", "1"))
parse_only <- tolower(Sys.getenv("COPYKAT_PARSE_ONLY", "false")) %in% c("true", "1", "yes")
sample_filter <- Sys.getenv("COPYKAT_SAMPLES", "")
if (nzchar(sample_filter)) {
  sample_filter <- trimws(unlist(strsplit(sample_filter, ",")))
} else {
  sample_filter <- character()
}

message("copykat settings: max_epi=", max_epi, ", max_ref=", max_ref,
        ", min_epi=", min_epi, ", min_ref=", min_ref,
        ", n_cores=", n_cores, ", parse_only=", parse_only)
if (length(sample_filter)) message("Sample filter: ", paste(sample_filter, collapse = ", "))

cell_meta <- as.data.table(readRDS(file.path(pre_root, "results", "single_cell_per_cell_scores.rds")))
normal_types_primary <- c("T_NK", "B_Plasma", "Myeloid", "Endothelial")
normal_types_fallback <- c(normal_types_primary, "CAF")

select_cells <- function(meta) {
  epi <- meta[cell_type == "Epithelial_Thyroid"]
  ref <- meta[cell_type %in% normal_types_primary]
  if (nrow(ref) < min_ref) ref <- meta[cell_type %in% normal_types_fallback]
  if (nrow(epi) < min_epi || nrow(ref) < min_ref) {
    return(list(skip = TRUE, reason = paste0("n_epi=", nrow(epi), "; n_ref=", nrow(ref))))
  }
  if (nrow(epi) > max_epi) {
    q <- quantile(epi$RAIR_like, 0.75, na.rm = TRUE)
    high <- epi[RAIR_like >= q]
    n_high_keep <- min(nrow(high), ceiling(max_epi * 0.45))
    keep_high <- if (nrow(high) > n_high_keep) high[sample(.N, n_high_keep)] else high
    rest <- epi[!cell %in% keep_high$cell]
    n_rest <- max_epi - nrow(keep_high)
    keep_rest <- if (nrow(rest) > n_rest) rest[sample(.N, n_rest)] else rest
    epi <- rbindlist(list(keep_high, keep_rest), fill = TRUE)
  }
  if (nrow(ref) > max_ref) ref <- ref[sample(.N, max_ref)]
  list(skip = FALSE, epi = epi, ref = ref,
       selected_cells = unique(c(epi$cell, ref$cell)),
       norm_cells = ref$cell)
}

dedupe_gene_matrix <- function(mat, genes) {
  keep <- !is.na(genes) & genes != "" & !duplicated(genes)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- genes[keep]
  mat
}

read_gse184_counts <- function(sample_id, selected_cells) {
  files <- list.files(file.path(geo_dir, "GSE184362_RAW"),
                      pattern = paste0("_", sample_id, "_matrix\\.mtx\\.gz$"),
                      full.names = TRUE)
  if (length(files) == 0) stop("No GSE184362 matrix for ", sample_id)
  matrix_file <- files[1]
  prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
  feature_file <- file.path(dirname(matrix_file), paste0(prefix, "_features.tsv.gz"))
  barcode_file <- file.path(dirname(matrix_file), paste0(prefix, "_barcodes.tsv.gz"))
  features <- fread(feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
  genes <- features[[2]]
  barcodes <- fread(barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE)[[1]]
  all_cells <- paste(sample_id, barcodes, sep = "_")
  keep_cols <- match(selected_cells, all_cells)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  mat <- Matrix::readMM(gzfile(matrix_file))
  mat <- mat[, keep_cols, drop = FALSE]
  colnames(mat) <- all_cells[keep_cols]
  dedupe_gene_matrix(as(mat, "dgCMatrix"), genes)
}

read_gse232_counts <- function(sample_id, selected_cells) {
  files <- list.files(file.path(geo_dir, "GSE232237_RAW"),
                      pattern = paste0("_", sample_id, "\\.count\\.tsv\\.gz$"),
                      full.names = TRUE)
  if (length(files) == 0) stop("No GSE232237 count table for ", sample_id)
  header <- names(fread(files[1], nrows = 0, showProgress = FALSE))
  selected_cells <- intersect(selected_cells, header)
  dt <- fread(files[1], select = c("gene_name", selected_cells),
              data.table = FALSE, showProgress = FALSE)
  genes <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  dedupe_gene_matrix(mat, genes)
}

summarize_prediction <- function(pred, meta, norm_cells, dataset_id, sample_id, out_dir, pred_file) {
  if (!"cell" %in% names(pred)) {
    if ("cell.names" %in% names(pred)) setnames(pred, "cell.names", "cell")
  }
  pred <- merge(pred, meta[, .(cell, disease_group, tissue_group, cell_type,
                               RAIR_like, TDS_raw, OneCarbon_raw, MAPK_raw,
                               Dediff_EMT_raw)],
                by = "cell", all.x = TRUE)
  pred[, dataset := dataset_id]
  pred[, sample := sample_id]
  pred_class_col <- intersect(c("copykat.pred", "prediction", "copykat_prediction"), names(pred))[1]
  if (is.na(pred_class_col)) pred_class_col <- setdiff(names(pred), "cell")[1]
  pred[, copykat_class := as.character(get(pred_class_col))]
  fwrite(pred, pred_file)
  sample_summary <- pred[, .(
    n_cells = .N,
    n_epithelial = sum(cell_type == "Epithelial_Thyroid", na.rm = TRUE),
    n_normal_ref = sum(cell %in% norm_cells),
    n_aneuploid = sum(copykat_class == "aneuploid", na.rm = TRUE),
    n_diploid = sum(copykat_class == "diploid", na.rm = TRUE),
    n_not_defined = sum(!copykat_class %in% c("aneuploid", "diploid"), na.rm = TRUE),
    epithelial_aneuploid_fraction = mean(copykat_class[cell_type == "Epithelial_Thyroid"] == "aneuploid", na.rm = TRUE),
    epithelial_median_RAIR_aneuploid = median(RAIR_like[cell_type == "Epithelial_Thyroid" & copykat_class == "aneuploid"], na.rm = TRUE),
    epithelial_median_RAIR_diploid = median(RAIR_like[cell_type == "Epithelial_Thyroid" & copykat_class == "diploid"], na.rm = TRUE)
  ), by = .(dataset, sample, disease_group, tissue_group)]
  fwrite(sample_summary, file.path(out_dir, paste0(sample_id, "_copykat_summary.csv")))
  invisible(pred)
}

run_one_sample <- function(dataset_id, sample_id) {
  key <- paste(dataset_id, sample_id, sep = ":")
  if (length(sample_filter) && !key %in% sample_filter && !sample_id %in% sample_filter) {
    return(NULL)
  }
  meta <- cell_meta[cell_meta$dataset == dataset_id & cell_meta$sample == sample_id]
  sel <- select_cells(meta)
  if (isTRUE(sel$skip)) {
    message("Skipping ", key, ": ", sel$reason)
    return(data.table(dataset = dataset_id, sample = sample_id, status = "skipped",
                      reason = sel$reason))
  }

  out_dir <- file.path(result_dir, dataset_id, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pred_file <- file.path(out_dir, paste0(sample_id, "_copykat_prediction.csv"))
  native_pred_file <- file.path(out_dir, paste0(sample_id, "_copykat_prediction.txt"))
  if (file.exists(pred_file) && file.info(pred_file)$size > 0) {
    pred <- fread(pred_file)
    return(data.table(dataset = dataset_id, sample = sample_id, status = "cached",
                      n_cells_run = nrow(pred),
                      n_epi_input = nrow(sel$epi), n_ref_input = nrow(sel$ref)))
  }
  if (file.exists(native_pred_file) && file.info(native_pred_file)$size > 0) {
    pred <- fread(native_pred_file)
    summarize_prediction(pred, meta, sel$norm_cells, dataset_id, sample_id, out_dir, pred_file)
    return(data.table(dataset = dataset_id, sample = sample_id, status = "cached_native_prediction",
                      n_cells_run = nrow(pred),
                      n_epi_input = nrow(sel$epi), n_ref_input = nrow(sel$ref)))
  }
  if (parse_only) {
    return(data.table(dataset = dataset_id, sample = sample_id, status = "no_native_prediction",
                      n_epi_input = nrow(sel$epi), n_ref_input = nrow(sel$ref)))
  }

  message("Running copykat: ", key, " epi=", nrow(sel$epi), " ref=", nrow(sel$ref))
  rawmat <- if (dataset_id == "GSE184362") {
    read_gse184_counts(sample_id, sel$selected_cells)
  } else {
    read_gse232_counts(sample_id, sel$selected_cells)
  }
  rawmat <- rawmat[, intersect(colnames(rawmat), sel$selected_cells), drop = FALSE]
  norm_cells <- intersect(sel$norm_cells, colnames(rawmat))
  if (length(norm_cells) < min_ref) {
    return(data.table(dataset = dataset_id, sample = sample_id, status = "skipped",
                      reason = paste0("norm_cells_after_read=", length(norm_cells))))
  }
  if (inherits(rawmat, "Matrix")) {
    rawmat <- as.matrix(rawmat)
    storage.mode(rawmat) <- "double"
  }

  old_wd <- getwd()
  setwd(out_dir)
  ck <- try(copykat(
    rawmat = rawmat,
    id.type = "S",
    cell.line = "no",
    ngene.chr = 5,
    min.gene.per.cell = 200,
    LOW.DR = 0.05,
    UP.DR = 0.10,
    win.size = 25,
    norm.cell.names = norm_cells,
    KS.cut = 0.10,
    sam.name = sample_id,
    distance = "euclidean",
    output.seg = "FALSE",
    plot.genes = "FALSE",
    genome = "hg20",
    n.cores = n_cores
  ), silent = TRUE)
  setwd(old_wd)
  rm(rawmat)
  gc()

  if (inherits(ck, "try-error") || is.null(ck$prediction)) {
    err <- if (inherits(ck, "try-error")) as.character(ck) else "copykat returned no prediction"
    writeLines(err, file.path(out_dir, paste0(sample_id, "_copykat_error.txt")))
    message("copykat failed: ", key, ": ", substr(err, 1, 180))
    return(data.table(dataset = dataset_id, sample = sample_id, status = "failed",
                      reason = substr(err, 1, 500),
                      n_epi_input = nrow(sel$epi), n_ref_input = nrow(sel$ref)))
  }

  pred <- as.data.table(ck$prediction)
  if (!"cell.names" %in% names(pred)) pred[, cell.names := rownames(ck$prediction)]
  summarize_prediction(pred, meta, norm_cells, dataset_id, sample_id, out_dir, pred_file)
  saveRDS(ck, file.path(out_dir, paste0(sample_id, "_copykat_result.rds")))
  data.table(dataset = dataset_id, sample = sample_id, status = "completed",
             n_cells_run = nrow(pred),
             n_epi_input = nrow(sel$epi), n_ref_input = nrow(sel$ref))
}

samples <- unique(cell_meta[, .(dataset, sample)])
setorder(samples, dataset, sample)
run_status <- rbindlist(lapply(seq_len(nrow(samples)), function(i) {
  run_one_sample(samples$dataset[i], samples$sample[i])
}), fill = TRUE)
fwrite(run_status, file.path(result_dir, "copykat_run_status.csv"))

pred_files <- list.files(result_dir, pattern = "_copykat_prediction\\.csv$", recursive = TRUE, full.names = TRUE)
all_pred <- rbindlist(lapply(pred_files, fread), fill = TRUE)
if (nrow(all_pred) > 0) {
  class_col <- intersect(c("copykat_class", "copykat.pred", "prediction"), names(all_pred))[1]
  if (is.na(class_col)) class_col <- "copykat_class"
  all_pred[, copykat_class := as.character(get(class_col))]
  fwrite(all_pred, file.path(result_dir, "copykat_all_predictions.csv"))
  sample_summary <- all_pred[, .(
    n_cells = .N,
    n_epithelial = sum(cell_type == "Epithelial_Thyroid", na.rm = TRUE),
    epithelial_aneuploid_fraction = mean(copykat_class[cell_type == "Epithelial_Thyroid"] == "aneuploid", na.rm = TRUE),
    epithelial_median_RAIR_aneuploid = median(RAIR_like[cell_type == "Epithelial_Thyroid" & copykat_class == "aneuploid"], na.rm = TRUE),
    epithelial_median_RAIR_diploid = median(RAIR_like[cell_type == "Epithelial_Thyroid" & copykat_class == "diploid"], na.rm = TRUE),
    n_aneuploid = sum(copykat_class == "aneuploid", na.rm = TRUE),
    n_diploid = sum(copykat_class == "diploid", na.rm = TRUE)
  ), by = .(dataset, sample, disease_group, tissue_group)]
  fwrite(sample_summary, file.path(result_dir, "copykat_sample_summary.csv"))

  epithelial <- all_pred[cell_type == "Epithelial_Thyroid" & copykat_class %in% c("aneuploid", "diploid")]
  tests <- epithelial[, {
    if (length(unique(copykat_class)) < 2) {
      data.table(n_aneuploid = sum(copykat_class == "aneuploid"),
                 n_diploid = sum(copykat_class == "diploid"),
                 median_RAIR_aneuploid = NA_real_,
                 median_RAIR_diploid = NA_real_,
                 delta_RAIR_aneuploid_minus_diploid = NA_real_,
                 p_wilcox = NA_real_)
    } else {
      wt <- suppressWarnings(wilcox.test(RAIR_like ~ copykat_class, exact = FALSE))
      data.table(n_aneuploid = sum(copykat_class == "aneuploid"),
                 n_diploid = sum(copykat_class == "diploid"),
                 median_RAIR_aneuploid = median(RAIR_like[copykat_class == "aneuploid"], na.rm = TRUE),
                 median_RAIR_diploid = median(RAIR_like[copykat_class == "diploid"], na.rm = TRUE),
                 delta_RAIR_aneuploid_minus_diploid =
                   median(RAIR_like[copykat_class == "aneuploid"], na.rm = TRUE) -
                   median(RAIR_like[copykat_class == "diploid"], na.rm = TRUE),
                 p_wilcox = wt$p.value)
    }
  }, by = dataset]
  fwrite(tests, file.path(result_dir, "copykat_epithelial_RAIR_tests.csv"))
}

sink(file.path(result_dir, "copykat_summary.txt"))
cat("copyKAT malignant/diploid classification batch\n\n")
cat("Settings:\n")
cat("max_epi=", max_epi, "; max_ref=", max_ref, "; min_epi=", min_epi,
    "; min_ref=", min_ref, "; n_cores=", n_cores, "\n\n", sep = "")
cat("Run status:\n")
print(run_status, row.names = FALSE)
if (exists("sample_summary")) {
  cat("\nSample summary:\n")
  print(sample_summary, row.names = FALSE)
}
if (exists("tests")) {
  cat("\nEpithelial RAIR-like by copyKAT class:\n")
  print(tests, row.names = FALSE)
}
sink()

message("copykat batch complete: ", result_dir)
