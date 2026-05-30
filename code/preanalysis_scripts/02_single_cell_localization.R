options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

root <- "D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/preanalysis"
data_dir <- file.path(root, "data", "geo")
result_dir <- file.path(root, "results")
figure_dir <- file.path(root, "figures")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

gene_sets <- list(
  TDS = c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1",
          "PAX8", "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB",
          "TPO", "TSHR"),
  OneCarbon = c("SHMT2", "MTHFD2", "MTHFD1L", "PHGDH", "PSAT1", "PSPH",
                "SLC1A5", "SLC7A5", "TYMS", "DHFR", "GART"),
  MAPK = c("DUSP4", "DUSP5", "DUSP6", "ETV4", "ETV5", "SPRY1", "SPRY2",
           "SPRY4", "CCND1", "FOSL1"),
  Dediff_EMT = c("HMGA2", "AXL", "VIM", "FN1", "CD44", "ITGA6", "KRT19",
                 "EPCAM", "LAMC2", "ITGB1")
)

marker_sets <- list(
  Epithelial_Thyroid = c("EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "TG",
                         "TPO", "SLC5A5", "TSHR", "SLC26A4"),
  T_NK = c("PTPRC", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "GZMB"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1"),
  Myeloid = c("LYZ", "LST1", "CD68", "CD14", "FCGR3A", "S100A8", "S100A9"),
  CAF = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "PDGFRA", "TAGLN"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "RAMP2", "CLDN5")
)

all_needed_genes <- unique(c(unlist(gene_sets), unlist(marker_sets)))

module_mean <- function(log_norm, genes, gene_index) {
  idx <- gene_index[intersect(genes, names(gene_index))]
  if (length(idx) == 0) return(rep(NA_real_, ncol(log_norm)))
  Matrix::colMeans(log_norm[idx, , drop = FALSE])
}

classify_cells <- function(score_df, min_score = 0.12) {
  marker_cols <- names(marker_sets)
  marker_mat <- as.matrix(score_df[, marker_cols, drop = FALSE])
  top <- max.col(marker_mat, ties.method = "first")
  cell_type <- marker_cols[top]
  cell_type[marker_mat[cbind(seq_len(nrow(marker_mat)), top)] < min_score] <- "Other"
  cell_type
}

add_rair_like <- function(df) {
  df$TDS_z <- as.numeric(scale(df$TDS_raw))
  df$OneCarbon_z <- as.numeric(scale(df$OneCarbon_raw))
  df$MAPK_z <- as.numeric(scale(df$MAPK_raw))
  df$Dediff_EMT_z <- as.numeric(scale(df$Dediff_EMT_raw))
  df$RAIR_like <- df$OneCarbon_z + df$MAPK_z + df$Dediff_EMT_z - df$TDS_z
  df
}

score_log_norm <- function(log_norm, genes, dataset, sample_id, disease_group, tissue_group) {
  gene_index <- seq_along(genes)
  names(gene_index) <- genes
  out <- data.frame(
    dataset = dataset,
    sample = sample_id,
    disease_group = disease_group,
    tissue_group = tissue_group,
    cell = colnames(log_norm),
    TDS_raw = module_mean(log_norm, gene_sets$TDS, gene_index),
    OneCarbon_raw = module_mean(log_norm, gene_sets$OneCarbon, gene_index),
    MAPK_raw = module_mean(log_norm, gene_sets$MAPK, gene_index),
    Dediff_EMT_raw = module_mean(log_norm, gene_sets$Dediff_EMT, gene_index),
    stringsAsFactors = FALSE
  )
  for (nm in names(marker_sets)) {
    out[[nm]] <- module_mean(log_norm, marker_sets[[nm]], gene_index)
  }
  out$cell_type <- classify_cells(out)
  out
}

process_gse232237_file <- function(file) {
  sample_id <- sub("^GSM[0-9]+_", "", basename(file))
  sample_id <- sub("\\.count\\.tsv\\.gz$", "", sample_id)
  disease_group <- ifelse(grepl("^AT", sample_id), "ATC", "PTC")
  message("GSE232237: ", sample_id)
  con <- gzfile(file, "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1)
  kept <- character()
  repeat {
    chunk <- readLines(con, n = 1000)
    if (length(chunk) == 0) break
    g <- sub("\t.*$", "", chunk)
    kept <- c(kept, chunk[g %in% all_needed_genes])
  }
  if (length(kept) == 0) stop("No target genes found in ", file)
  dt <- data.table::fread(text = paste(c(header, kept), collapse = "\n"),
                          data.table = FALSE, showProgress = FALSE)
  genes <- make.unique(dt[[1]])
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  ## These TSVs are dense count tables. To keep the preanalysis tractable on
  ## Windows, only target genes are streamed from gzip; use log1p counts for
  ## localization rather than exact library-size normalization.
  log_norm <- log1p(mat)
  rm(dt, mat)
  gc()
  score_log_norm(Matrix(log_norm, sparse = TRUE), rownames(log_norm),
                 "GSE232237", sample_id, disease_group, disease_group)
}

parse_gse184_sample <- function(prefix) {
  x <- sub("^GSM[0-9]+_", "", prefix)
  patient <- sub("^(PTC[0-9]+).*", "\\1", x)
  tissue <- sub("^PTC[0-9]+_", "", x)
  tissue_group <- ifelse(tissue == "P", "Paratumor",
                         ifelse(tissue == "T", "Primary_tumor",
                                ifelse(grepl("LN", tissue), "Lymph_node_metastasis",
                                       ifelse(tissue == "SC", "Subcutaneous_metastasis", tissue))))
  data.frame(sample_id = x, patient = patient, tissue_group = tissue_group,
             disease_group = ifelse(tissue_group == "Paratumor", "Paratumor", "PTC_related"),
             stringsAsFactors = FALSE)
}

process_gse184362_sample <- function(matrix_file) {
  prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
  meta <- parse_gse184_sample(prefix)
  base <- file.path(dirname(matrix_file), prefix)
  feature_file <- paste0(base, "_features.tsv.gz")
  barcode_file <- paste0(base, "_barcodes.tsv.gz")
  message("GSE184362: ", meta$sample_id)
  features <- data.table::fread(feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
  genes <- make.unique(features[[2]])
  barcodes <- data.table::fread(barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE)[[1]]
  mat <- Matrix::readMM(gzfile(matrix_file))
  rownames(mat) <- genes
  colnames(mat) <- paste(meta$sample_id, barcodes, sep = "_")
  lib <- Matrix::colSums(mat)
  keep <- intersect(all_needed_genes, rownames(mat))
  log_norm <- log1p(t(t(mat[keep, , drop = FALSE]) / lib) * 10000)
  out <- score_log_norm(log_norm, rownames(log_norm),
                        "GSE184362", meta$sample_id, meta$disease_group, meta$tissue_group)
  out$patient <- meta$patient
  out
}

summarize_localization <- function(df, dataset) {
  df <- add_rair_like(df)
  cutoff <- quantile(df$RAIR_like, 0.9, na.rm = TRUE)
  df$RAIR_top10 <- df$RAIR_like >= cutoff
  med <- aggregate(cbind(TDS_raw, OneCarbon_raw, MAPK_raw, Dediff_EMT_raw, RAIR_like) ~
                     dataset + disease_group + tissue_group + cell_type,
                   data = df, FUN = median, na.rm = TRUE)
  n_by_type <- aggregate(cell ~ dataset + disease_group + tissue_group + cell_type,
                         data = df, FUN = length)
  names(n_by_type)[names(n_by_type) == "cell"] <- "n_cells"
  med <- merge(med, n_by_type, by = c("dataset", "disease_group", "tissue_group", "cell_type"))

  top_dist <- as.data.frame(table(df$dataset, df$cell_type, df$RAIR_top10), stringsAsFactors = FALSE)
  names(top_dist) <- c("dataset", "cell_type", "RAIR_top10", "n_cells")
  top_dist <- top_dist[top_dist$RAIR_top10 == "TRUE", c("dataset", "cell_type", "n_cells")]
  top_dist$fraction_of_top10 <- top_dist$n_cells / sum(top_dist$n_cells)

  sample_dist <- aggregate(RAIR_like ~ dataset + sample + disease_group + tissue_group + cell_type,
                           data = df, FUN = median, na.rm = TRUE)
  n_sample <- aggregate(cell ~ dataset + sample + disease_group + tissue_group + cell_type,
                        data = df, FUN = length)
  names(n_sample)[names(n_sample) == "cell"] <- "n_cells"
  sample_dist <- merge(sample_dist, n_sample,
                       by = c("dataset", "sample", "disease_group", "tissue_group", "cell_type"))

  list(scores = df, medians = med, top_distribution = top_dist, sample_medians = sample_dist)
}

## GSE232237
g232_dir <- file.path(data_dir, "GSE232237_RAW")
g232_files <- list.files(g232_dir, pattern = "\\.count\\.tsv\\.gz$", full.names = TRUE)
g232_list <- lapply(g232_files, process_gse232237_file)
g232_df <- do.call(rbind, g232_list)
g232_res <- summarize_localization(g232_df, "GSE232237")

## GSE184362
g184_dir <- file.path(data_dir, "GSE184362_RAW")
g184_files <- list.files(g184_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
g184_list <- lapply(g184_files, process_gse184362_sample)
g184_df <- do.call(rbind, g184_list)
g184_res <- summarize_localization(g184_df, "GSE184362")

if (!"patient" %in% colnames(g232_res$scores)) g232_res$scores$patient <- NA_character_
if (!"patient" %in% colnames(g184_res$scores)) g184_res$scores$patient <- NA_character_
cell_scores <- rbind(g232_res$scores[, sort(colnames(g232_res$scores))],
                     g184_res$scores[, sort(colnames(g184_res$scores))])
cell_medians <- rbind(g232_res$medians, g184_res$medians)
top_dist <- rbind(g232_res$top_distribution, g184_res$top_distribution)
sample_medians <- rbind(g232_res$sample_medians, g184_res$sample_medians)

write.csv(cell_medians, file.path(result_dir, "single_cell_module_medians_by_celltype.csv"), row.names = FALSE)
write.csv(top_dist, file.path(result_dir, "single_cell_top_RAIR_like_celltype_distribution.csv"), row.names = FALSE)
write.csv(sample_medians, file.path(result_dir, "single_cell_sample_celltype_medians.csv"), row.names = FALSE)

## Store per-cell scores in compressed RDS to avoid huge CSV.
saveRDS(cell_scores, file.path(result_dir, "single_cell_per_cell_scores.rds"))

plot_medians <- function(med, dataset) {
  sub <- med[med$dataset == dataset & med$n_cells >= 50, ]
  p <- ggplot(sub, aes(x = cell_type, y = RAIR_like, fill = cell_type)) +
    geom_col(width = 0.7) +
    facet_grid(disease_group ~ tissue_group, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
    labs(x = NULL, y = "Median RAIR-like score", title = paste(dataset, "RAIR-like localization"))
  ggsave(file.path(figure_dir, paste0(dataset, "_single_cell_RAIR_like_by_celltype.png")),
         p, width = 9, height = 5.5, dpi = 300)
}

plot_top <- function(td, dataset) {
  sub <- td[td$dataset == dataset, ]
  p <- ggplot(sub, aes(x = reorder(cell_type, fraction_of_top10), y = fraction_of_top10, fill = cell_type)) +
    geom_col(width = 0.7) +
    coord_flip() +
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Fraction among top 10% RAIR-like cells",
         title = paste(dataset, "top RAIR-like cells"))
  ggsave(file.path(figure_dir, paste0(dataset, "_top_RAIR_like_celltype_fraction.png")),
         p, width = 5.5, height = 4.2, dpi = 300)
}

plot_medians(cell_medians, "GSE232237")
plot_medians(cell_medians, "GSE184362")
plot_top(top_dist, "GSE232237")
plot_top(top_dist, "GSE184362")

sink(file.path(result_dir, "single_cell_localization_summary.txt"))
cat("Single-cell localization preanalysis completed\n\n")
cat("Cell counts by dataset and assigned broad cell type:\n")
print(table(cell_scores$dataset, cell_scores$cell_type))
cat("\nTop 10% RAIR-like cell-type distribution:\n")
print(top_dist)
cat("\nMedian module scores by cell type/tissue group:\n")
print(cell_medians[order(cell_medians$dataset, cell_medians$tissue_group,
                         -cell_medians$RAIR_like), ])
cat("\nSample-level medians for epithelial/thyroid cells:\n")
print(sample_medians[sample_medians$cell_type == "Epithelial_Thyroid",
                     c("dataset", "sample", "disease_group", "tissue_group",
                       "cell_type", "n_cells", "RAIR_like")])
sink()

message("Single-cell localization complete. Results written to: ", result_dir)
