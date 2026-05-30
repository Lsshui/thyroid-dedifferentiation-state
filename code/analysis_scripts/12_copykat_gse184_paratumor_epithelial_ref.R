options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(copykat)
})

set.seed(20260527)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pre_root <- file.path(root, "preanalysis")
geo_dir <- file.path(pre_root, "data", "geo", "GSE184362_RAW")
result_dir <- file.path(root, "analysis", "results", "copykat_paratumor_ref", "GSE184362")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

max_target <- as.integer(Sys.getenv("COPYKAT_MAX_TARGET_EPI", "180"))
max_ref <- as.integer(Sys.getenv("COPYKAT_MAX_PARATUMOR_REF", "180"))
min_target <- as.integer(Sys.getenv("COPYKAT_MIN_TARGET_EPI", "50"))
min_ref <- as.integer(Sys.getenv("COPYKAT_MIN_PARATUMOR_REF", "100"))
n_cores <- as.integer(Sys.getenv("COPYKAT_N_CORES", "1"))
parse_only <- tolower(Sys.getenv("COPYKAT_PARSE_ONLY", "false")) %in% c("true", "1", "yes")
sample_filter <- Sys.getenv("COPYKAT_SAMPLES", "PTC10_T,PTC10_RightLN,PTC11_SC")
sample_filter <- trimws(unlist(strsplit(sample_filter, ",")))

message("GSE184362 paratumor epithelial-reference copyKAT settings: max_target=", max_target,
        ", max_ref=", max_ref, ", parse_only=", parse_only)

cell_meta <- as.data.table(readRDS(file.path(pre_root, "results", "single_cell_per_cell_scores.rds")))
cell_meta <- cell_meta[dataset == "GSE184362"]

matrix_files <- list.files(geo_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
sample_to_file <- setNames(matrix_files, sub("^GSM[0-9]+_", "", sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files))))

select_target_epi <- function(sample_id) {
  epi <- cell_meta[sample == sample_id & cell_type == "Epithelial_Thyroid"]
  if (nrow(epi) < min_target) return(NULL)
  if (nrow(epi) > max_target) {
    q <- quantile(epi$RAIR_like, 0.75, na.rm = TRUE)
    high <- epi[RAIR_like >= q]
    n_high_keep <- min(nrow(high), ceiling(max_target * 0.50))
    keep_high <- if (nrow(high) > n_high_keep) high[sample(.N, n_high_keep)] else high
    rest <- epi[!cell %in% keep_high$cell]
    keep_rest <- if (nrow(rest) > max_target - nrow(keep_high)) rest[sample(.N, max_target - nrow(keep_high))] else rest
    epi <- rbindlist(list(keep_high, keep_rest), fill = TRUE)
  }
  epi
}

select_paratumor_ref <- function() {
  ref <- cell_meta[tissue_group == "Paratumor" & cell_type == "Epithelial_Thyroid" &
                     nchar(cell) > 0]
  ref <- ref[sample %in% names(sample_to_file)]
  if (nrow(ref) < min_ref) stop("Too few paratumor epithelial reference cells.")
  ## Avoid a single paratumor sample dominating the normal reference.
  per_sample <- max(20, ceiling(max_ref / max(1, length(unique(ref$sample)))))
  ref <- ref[, if (.N > per_sample) .SD[sample(.N, per_sample)] else .SD, by = sample]
  if (nrow(ref) > max_ref) ref <- ref[sample(.N, max_ref)]
  ref
}

read_gse184_counts <- function(sample_id, selected_cells) {
  matrix_file <- sample_to_file[[sample_id]]
  if (is.null(matrix_file) || is.na(matrix_file)) stop("No matrix file for ", sample_id)
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
  keep <- !is.na(genes) & genes != "" & !duplicated(genes)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- genes[keep]
  as.matrix(mat)
}

build_matrix <- function(target_epi, ref_epi) {
  selected <- rbindlist(list(target_epi[, .(sample, cell)],
                            ref_epi[, .(sample, cell)]), fill = TRUE)
  mats <- lapply(split(selected, selected$sample), function(dt) {
    read_gse184_counts(dt$sample[1], dt$cell)
  })
  common_genes <- Reduce(intersect, lapply(mats, rownames))
  mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
  mat <- do.call(cbind, mats)
  storage.mode(mat) <- "double"
  mat
}

parse_prediction <- function(pred, target_epi, ref_epi, sample_id, out_dir, csv_file) {
  if ("cell.names" %in% names(pred)) setnames(pred, "cell.names", "cell")
  meta <- rbindlist(list(target_epi, ref_epi), fill = TRUE)
  pred <- merge(pred, meta[, .(cell, disease_group, tissue_group, cell_type,
                               RAIR_like, TDS_raw, OneCarbon_raw, MAPK_raw,
                               Dediff_EMT_raw)],
                by = "cell", all.x = TRUE)
  pred[, dataset := "GSE184362"]
  pred[, sample_target := sample_id]
  pred_class_col <- intersect(c("copykat.pred", "prediction", "copykat_prediction"), names(pred))[1]
  if (is.na(pred_class_col)) pred_class_col <- setdiff(names(pred), "cell")[1]
  pred[, copykat_class := as.character(get(pred_class_col))]
  fwrite(pred, csv_file)
  summary <- pred[cell %in% target_epi$cell, .(
    n_target_epi = .N,
    n_aneuploid = sum(copykat_class == "aneuploid", na.rm = TRUE),
    n_diploid = sum(copykat_class == "diploid", na.rm = TRUE),
    aneuploid_fraction = mean(copykat_class == "aneuploid", na.rm = TRUE),
    median_RAIR_aneuploid = median(RAIR_like[copykat_class == "aneuploid"], na.rm = TRUE),
    median_RAIR_diploid = median(RAIR_like[copykat_class == "diploid"], na.rm = TRUE)
  ), by = .(dataset, sample_target)]
  fwrite(summary, file.path(out_dir, paste0(sample_id, "_paratumor_ref_copykat_summary.csv")))
  invisible(pred)
}

run_one <- function(sample_id) {
  target_epi <- select_target_epi(sample_id)
  if (is.null(target_epi)) {
    message("Skipping ", sample_id, ": too few target epithelial cells")
    return(data.table(sample_target = sample_id, status = "skipped_too_few_target"))
  }
  ref_epi <- select_paratumor_ref()
  out_dir <- file.path(result_dir, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  csv_file <- file.path(out_dir, paste0(sample_id, "_paratumor_ref_copykat_prediction.csv"))
  native_file <- file.path(out_dir, paste0(sample_id, "_paratumor_ref_copykat_prediction.txt"))
  if (file.exists(csv_file) && file.info(csv_file)$size > 0) {
    pred <- fread(csv_file)
    return(data.table(sample_target = sample_id, status = "cached", n_cells = nrow(pred)))
  }
  if (file.exists(native_file) && file.info(native_file)$size > 0) {
    pred <- fread(native_file)
    parse_prediction(pred, target_epi, ref_epi, sample_id, out_dir, csv_file)
    return(data.table(sample_target = sample_id, status = "cached_native_prediction",
                      n_cells = nrow(pred), n_target_epi = nrow(target_epi), n_ref_epi = nrow(ref_epi)))
  }
  if (parse_only) {
    return(data.table(sample_target = sample_id, status = "no_native_prediction",
                      n_target_epi = nrow(target_epi), n_ref_epi = nrow(ref_epi)))
  }

  message("Running paratumor-ref copyKAT for ", sample_id,
          ": target_epi=", nrow(target_epi), ", ref_epi=", nrow(ref_epi))
  rawmat <- build_matrix(target_epi, ref_epi)
  norm_cells <- intersect(ref_epi$cell, colnames(rawmat))
  old_wd <- getwd()
  setwd(out_dir)
  ck <- try(copykat(rawmat = rawmat, id.type = "S", cell.line = "no",
                    ngene.chr = 5, min.gene.per.cell = 200,
                    LOW.DR = 0.05, UP.DR = 0.10, win.size = 25,
                    norm.cell.names = norm_cells, KS.cut = 0.10,
                    sam.name = paste0(sample_id, "_paratumor_ref"),
                    distance = "euclidean", output.seg = "FALSE",
                    plot.genes = "FALSE", genome = "hg20", n.cores = n_cores),
            silent = TRUE)
  setwd(old_wd)
  if (inherits(ck, "try-error") || is.null(ck$prediction)) {
    err <- if (inherits(ck, "try-error")) as.character(ck) else "copykat returned no prediction"
    writeLines(err, file.path(out_dir, paste0(sample_id, "_paratumor_ref_copykat_error.txt")))
    return(data.table(sample_target = sample_id, status = "failed",
                      reason = substr(err, 1, 500),
                      n_target_epi = nrow(target_epi), n_ref_epi = nrow(ref_epi)))
  }
  pred <- as.data.table(ck$prediction)
  if (!"cell.names" %in% names(pred)) pred[, cell.names := rownames(ck$prediction)]
  parse_prediction(pred, target_epi, ref_epi, sample_id, out_dir, csv_file)
  saveRDS(ck, file.path(out_dir, paste0(sample_id, "_paratumor_ref_copykat_result.rds")))
  data.table(sample_target = sample_id, status = "completed",
             n_cells = ncol(rawmat), n_target_epi = nrow(target_epi), n_ref_epi = nrow(ref_epi))
}

run_status <- rbindlist(lapply(sample_filter, run_one), fill = TRUE)
fwrite(run_status, file.path(result_dir, "paratumor_ref_copykat_run_status.csv"))

pred_files <- list.files(result_dir, pattern = "_paratumor_ref_copykat_prediction\\.csv$",
                         recursive = TRUE, full.names = TRUE)
if (length(pred_files) > 0) {
  all_pred <- rbindlist(lapply(pred_files, fread), fill = TRUE)
  fwrite(all_pred, file.path(result_dir, "paratumor_ref_copykat_all_predictions.csv"))
  target_pred <- all_pred[tissue_group != "Paratumor" & cell_type == "Epithelial_Thyroid" &
                            copykat_class %in% c("aneuploid", "diploid")]
  sample_summary <- target_pred[, .(
    n_target_epi = .N,
    aneuploid_fraction = mean(copykat_class == "aneuploid", na.rm = TRUE),
    median_RAIR_aneuploid = median(RAIR_like[copykat_class == "aneuploid"], na.rm = TRUE),
    median_RAIR_diploid = median(RAIR_like[copykat_class == "diploid"], na.rm = TRUE),
    n_aneuploid = sum(copykat_class == "aneuploid"),
    n_diploid = sum(copykat_class == "diploid")
  ), by = .(dataset, sample_target)]
  fwrite(sample_summary, file.path(result_dir, "paratumor_ref_copykat_sample_summary.csv"))
  tests <- target_pred[, {
    if (length(unique(copykat_class)) < 2) {
      data.table(n_aneuploid = sum(copykat_class == "aneuploid"),
                 n_diploid = sum(copykat_class == "diploid"),
                 median_RAIR_aneuploid = NA_real_, median_RAIR_diploid = NA_real_,
                 delta_RAIR = NA_real_, p_wilcox = NA_real_)
    } else {
      wt <- suppressWarnings(wilcox.test(RAIR_like ~ copykat_class, exact = FALSE))
      data.table(n_aneuploid = sum(copykat_class == "aneuploid"),
                 n_diploid = sum(copykat_class == "diploid"),
                 median_RAIR_aneuploid = median(RAIR_like[copykat_class == "aneuploid"], na.rm = TRUE),
                 median_RAIR_diploid = median(RAIR_like[copykat_class == "diploid"], na.rm = TRUE),
                 delta_RAIR = median(RAIR_like[copykat_class == "aneuploid"], na.rm = TRUE) -
                   median(RAIR_like[copykat_class == "diploid"], na.rm = TRUE),
                 p_wilcox = wt$p.value)
    }
  }]
  fwrite(tests, file.path(result_dir, "paratumor_ref_copykat_RAIR_test.csv"))
}

sink(file.path(result_dir, "paratumor_ref_copykat_summary.txt"))
cat("GSE184362 copyKAT with paratumor epithelial reference\n\n")
cat("Settings: max_target=", max_target, "; max_ref=", max_ref, "; parse_only=", parse_only, "\n", sep = "")
cat("Run status:\n")
print(run_status, row.names = FALSE)
if (exists("sample_summary")) {
  cat("\nTarget epithelial sample summary:\n")
  print(sample_summary, row.names = FALSE)
}
if (exists("tests")) {
  cat("\nTarget epithelial RAIR-like by copyKAT class:\n")
  print(tests, row.names = FALSE)
}
sink()

message("GSE184362 paratumor-ref copyKAT complete: ", result_dir)
