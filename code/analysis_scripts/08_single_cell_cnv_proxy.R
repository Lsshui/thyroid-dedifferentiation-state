options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(ggplot2)
})

set.seed(20260527)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pre_root <- file.path(root, "preanalysis")
analysis_root <- file.path(root, "analysis")
data_dir <- file.path(pre_root, "data", "geo")
result_dir <- file.path(analysis_root, "results")
figure_dir <- file.path(analysis_root, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

window_genes <- 100
max_reference_cells <- 2000
min_reference_cells <- 100
reference_types_primary <- c("T_NK", "B_Plasma", "Myeloid", "Endothelial")
reference_types_fallback <- c(reference_types_primary, "CAF")

message("Loading per-cell score metadata...")
cell_scores <- readRDS(file.path(pre_root, "results", "single_cell_per_cell_scores.rds"))
cell_scores <- as.data.table(cell_scores)

message("Loading gene genomic positions from cached TCGA object...")
se <- readRDS(file.path(analysis_root, "data", "tcga", "TCGA_THCA_STAR_counts_primary_tumor.rds"))
rd <- as.data.frame(rowData(se))
rr <- rowRanges(se)
gene_pos <- data.table(
  gene = rd$gene_name,
  chr = as.character(GenomeInfoDb::seqnames(rr)),
  start = GenomicRanges::start(rr)
)
gene_pos <- gene_pos[!is.na(gene) & gene != ""]
gene_pos[, chr := sub("^chr", "", chr)]
gene_pos <- gene_pos[chr %in% as.character(1:22)]
gene_pos[, chr_num := as.integer(chr)]
setorder(gene_pos, gene, chr_num, start)
gene_pos <- gene_pos[!duplicated(gene)]
setorder(gene_pos, chr_num, start)
gene_pos[, gene_rank_chr := seq_len(.N), by = chr_num]
gene_pos[, bin := paste0("chr", chr_num, "_bin", ceiling(gene_rank_chr / window_genes))]
gene_pos[, bin_order := .GRP, by = bin]
gene_pos <- gene_pos[!is.na(bin_order)]

make_unique_gene_matrix <- function(mat, genes) {
  keep <- !is.na(genes) & genes != "" & genes %in% gene_pos$gene & !duplicated(genes)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- genes[keep]
  mat
}

select_reference_cells <- function(meta) {
  ref <- meta[cell_type %in% reference_types_primary, cell]
  if (length(ref) < min_reference_cells) {
    ref <- meta[cell_type %in% reference_types_fallback, cell]
  }
  if (length(ref) > max_reference_cells) {
    ref <- sample(ref, max_reference_cells)
  }
  ref
}

compute_cnv_proxy <- function(counts, sample_meta, sample_id, dataset) {
  if (!inherits(counts, "Matrix")) counts <- Matrix(counts, sparse = FALSE)
  common_cells <- intersect(colnames(counts), sample_meta$cell)
  sample_meta <- sample_meta[match(common_cells, cell)]
  counts <- counts[, common_cells, drop = FALSE]

  epi_cells <- sample_meta[cell_type == "Epithelial_Thyroid", cell]
  ref_cells <- select_reference_cells(sample_meta)
  if (length(epi_cells) == 0 || length(ref_cells) < min_reference_cells) {
    message("Skipping ", dataset, " ", sample_id, ": epithelial cells = ", length(epi_cells),
            ", reference cells = ", length(ref_cells))
    return(NULL)
  }

  selected_cells <- unique(c(epi_cells, ref_cells))
  counts <- counts[, selected_cells, drop = FALSE]
  lib <- Matrix::colSums(counts)
  keep_cell <- lib > 0
  counts <- counts[, keep_cell, drop = FALSE]
  selected_cells <- colnames(counts)
  epi_cells <- intersect(epi_cells, selected_cells)
  ref_cells <- intersect(ref_cells, selected_cells)
  if (length(epi_cells) == 0 || length(ref_cells) < min_reference_cells) return(NULL)

  norm <- t(t(counts) / Matrix::colSums(counts)) * 10000
  norm@x <- log1p(norm@x)
  rm(counts)
  gc()

  ref_idx <- match(ref_cells, colnames(norm))
  ref_mean <- Matrix::rowMeans(norm[, ref_idx, drop = FALSE])

  genes <- rownames(norm)
  pos <- gene_pos[match(genes, gene)]
  valid <- !is.na(pos$bin_order)
  norm <- norm[valid, , drop = FALSE]
  ref_mean <- ref_mean[valid]
  pos <- pos[valid]
  bins <- split(seq_len(nrow(norm)), pos$bin_order)

  bin_means <- vapply(bins, function(idx) {
    Matrix::colMeans(norm[idx, , drop = FALSE]) - mean(ref_mean[idx], na.rm = TRUE)
  }, FUN.VALUE = numeric(ncol(norm)))
  bin_means <- t(bin_means)
  rownames(bin_means) <- names(bins)
  colnames(bin_means) <- colnames(norm)

  amp <- Matrix::colMeans(abs(bin_means))
  ref_amp <- amp[ref_cells]
  ref_med <- median(ref_amp, na.rm = TRUE)
  ref_mad <- mad(ref_amp, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(ref_mad) || ref_mad <= 0) ref_mad <- sd(ref_amp, na.rm = TRUE)
  if (!is.finite(ref_mad) || ref_mad <= 0) ref_mad <- 1e-6

  cell_amp <- data.table(
    dataset = dataset,
    sample = sample_id,
    cell = epi_cells,
    cnv_proxy_amp = as.numeric(amp[epi_cells]),
    cnv_proxy_z_vs_ref = as.numeric((amp[epi_cells] - ref_med) / ref_mad),
    ref_cells_used = length(ref_cells),
    ref_amp_median = ref_med,
    ref_amp_mad = ref_mad
  )
  cell_amp[, cnv_proxy_high := cnv_proxy_z_vs_ref > 3]
  cell_amp <- merge(cell_amp, sample_meta[, .(cell, disease_group, tissue_group, cell_type,
                                              RAIR_like, TDS_raw, OneCarbon_raw, MAPK_raw,
                                              Dediff_EMT_raw)],
                    by = "cell", all.x = TRUE)
  cell_amp
}

process_gse184_sample <- function(matrix_file) {
  prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
  sample_id <- sub("^GSM[0-9]+_", "", prefix)
  meta <- cell_scores[dataset == "GSE184362" & sample == sample_id]
  if (nrow(meta) == 0) return(NULL)
  message("CNV proxy GSE184362: ", sample_id)

  feature_file <- file.path(dirname(matrix_file), paste0(prefix, "_features.tsv.gz"))
  barcode_file <- file.path(dirname(matrix_file), paste0(prefix, "_barcodes.tsv.gz"))
  features <- fread(feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
  genes <- features[[2]]
  barcodes <- fread(barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE)[[1]]
  mat <- Matrix::readMM(gzfile(matrix_file))
  colnames(mat) <- paste(sample_id, barcodes, sep = "_")
  mat <- make_unique_gene_matrix(mat, genes)
  compute_cnv_proxy(mat, meta, sample_id, "GSE184362")
}

process_gse232_file <- function(file) {
  sample_id <- sub("^GSM[0-9]+_", "", basename(file))
  sample_id <- sub("\\.count\\.tsv\\.gz$", "", sample_id)
  meta <- cell_scores[dataset == "GSE232237" & sample == sample_id]
  if (nrow(meta) == 0) return(NULL)
  message("CNV proxy GSE232237: ", sample_id)

  dt <- fread(file, data.table = FALSE, showProgress = FALSE)
  genes <- dt[[1]]
  keep <- !is.na(genes) & genes != "" & genes %in% gene_pos$gene & !duplicated(genes)
  mat <- as.matrix(dt[keep, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes[keep]
  rm(dt)
  gc()
  compute_cnv_proxy(mat, meta, sample_id, "GSE232237")
}

cnv_file <- file.path(result_dir, "single_cell_cnv_proxy_epithelial_cells.csv")
if (file.exists(cnv_file) && file.info(cnv_file)$size > 0) {
  message("Loading cached per-cell CNV proxy results...")
  cnv <- fread(cnv_file)
} else {
  message("Processing GSE184362 sparse matrices...")
  g184_files <- list.files(file.path(data_dir, "GSE184362_RAW"),
                           pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
  g184_res <- rbindlist(lapply(g184_files, process_gse184_sample), fill = TRUE)

  message("Processing GSE232237 dense count tables...")
  g232_files <- list.files(file.path(data_dir, "GSE232237_RAW"),
                           pattern = "\\.count\\.tsv\\.gz$", full.names = TRUE)
  g232_res <- rbindlist(lapply(g232_files, process_gse232_file), fill = TRUE)

  cnv <- rbindlist(list(g184_res, g232_res), fill = TRUE)
  if (nrow(cnv) == 0) stop("No CNV proxy results produced.")
  fwrite(cnv, cnv_file)
}

sample_summary <- cnv[, .(
  n_epithelial = as.integer(.N),
  median_RAIR_like = as.numeric(median(RAIR_like, na.rm = TRUE)),
  median_cnv_proxy_amp = as.numeric(median(cnv_proxy_amp, na.rm = TRUE)),
  median_cnv_proxy_z = as.numeric(median(cnv_proxy_z_vs_ref, na.rm = TRUE)),
  fraction_cnv_proxy_high = as.numeric(mean(cnv_proxy_high, na.rm = TRUE)),
  fraction_RAIR_top_quartile_cnv_high = {
    q <- quantile(RAIR_like, 0.75, na.rm = TRUE)
    as.numeric(mean(cnv_proxy_high[RAIR_like >= q], na.rm = TRUE))
  },
  ref_cells_used = as.numeric(median(ref_cells_used, na.rm = TRUE))
), by = .(dataset, sample, disease_group, tissue_group)]
fwrite(sample_summary, file.path(result_dir, "single_cell_cnv_proxy_sample_summary.csv"))

cor_tests <- cnv[, {
  ok <- is.finite(RAIR_like) & is.finite(cnv_proxy_amp)
  ct <- if (sum(ok) >= 10) suppressWarnings(cor.test(RAIR_like[ok], cnv_proxy_amp[ok], method = "spearman")) else NULL
  data.table(n_cells = sum(ok),
             spearman_rho = if (is.null(ct)) NA_real_ else unname(ct$estimate),
             p_value = if (is.null(ct)) NA_real_ else ct$p.value)
}, by = .(dataset)]

status_tests <- cnv[, {
  if (length(unique(cnv_proxy_high)) < 2) {
    data.table(n_high = sum(cnv_proxy_high), n_low = sum(!cnv_proxy_high),
               median_RAIR_high = NA_real_, median_RAIR_low = NA_real_,
               delta_RAIR_high_minus_low = NA_real_, p_wilcox = NA_real_)
  } else {
    wt <- suppressWarnings(wilcox.test(RAIR_like ~ cnv_proxy_high, exact = FALSE))
    data.table(n_high = sum(cnv_proxy_high), n_low = sum(!cnv_proxy_high),
               median_RAIR_high = median(RAIR_like[cnv_proxy_high], na.rm = TRUE),
               median_RAIR_low = median(RAIR_like[!cnv_proxy_high], na.rm = TRUE),
               delta_RAIR_high_minus_low = median(RAIR_like[cnv_proxy_high], na.rm = TRUE) -
                 median(RAIR_like[!cnv_proxy_high], na.rm = TRUE),
               p_wilcox = wt$p.value)
  }
}, by = .(dataset)]

g184_epi <- sample_summary[dataset == "GSE184362" & n_epithelial >= 50]
if (nrow(g184_epi) > 0) {
  g184_epi[, tumor_met_group := ifelse(tissue_group == "Paratumor", "Paratumor", "Tumor_met")]
}
g184_tests <- if (nrow(g184_epi) >= 4 && length(unique(g184_epi$tumor_met_group)) == 2) {
  wt1 <- suppressWarnings(wilcox.test(median_cnv_proxy_z ~ tumor_met_group, data = g184_epi, exact = FALSE))
  wt2 <- suppressWarnings(wilcox.test(fraction_cnv_proxy_high ~ tumor_met_group, data = g184_epi, exact = FALSE))
  data.table(
    comparison = c("GSE184362 tumor/met vs paratumor sample median CNV-z",
                   "GSE184362 tumor/met vs paratumor CNV-high fraction"),
    n_tumor_met = sum(g184_epi$tumor_met_group == "Tumor_met"),
    n_paratumor = sum(g184_epi$tumor_met_group == "Paratumor"),
    median_tumor_met = c(median(g184_epi$median_cnv_proxy_z[g184_epi$tumor_met_group == "Tumor_met"], na.rm = TRUE),
                         median(g184_epi$fraction_cnv_proxy_high[g184_epi$tumor_met_group == "Tumor_met"], na.rm = TRUE)),
    median_paratumor = c(median(g184_epi$median_cnv_proxy_z[g184_epi$tumor_met_group == "Paratumor"], na.rm = TRUE),
                         median(g184_epi$fraction_cnv_proxy_high[g184_epi$tumor_met_group == "Paratumor"], na.rm = TRUE)),
    p_wilcox = c(wt1$p.value, wt2$p.value)
  )
} else data.table()

g232_epi <- sample_summary[dataset == "GSE232237" & n_epithelial >= 50]
g232_tests <- if (nrow(g232_epi) >= 4 && length(unique(g232_epi$disease_group)) == 2) {
  wt1 <- suppressWarnings(wilcox.test(median_cnv_proxy_z ~ disease_group, data = g232_epi, exact = FALSE))
  wt2 <- suppressWarnings(wilcox.test(fraction_cnv_proxy_high ~ disease_group, data = g232_epi, exact = FALSE))
  data.table(
    comparison = c("GSE232237 ATC vs PTC sample median CNV-z",
                   "GSE232237 ATC vs PTC CNV-high fraction"),
    n_atc = sum(g232_epi$disease_group == "ATC"),
    n_ptc = sum(g232_epi$disease_group == "PTC"),
    median_atc = c(median(g232_epi$median_cnv_proxy_z[g232_epi$disease_group == "ATC"], na.rm = TRUE),
                   median(g232_epi$fraction_cnv_proxy_high[g232_epi$disease_group == "ATC"], na.rm = TRUE)),
    median_ptc = c(median(g232_epi$median_cnv_proxy_z[g232_epi$disease_group == "PTC"], na.rm = TRUE),
                   median(g232_epi$fraction_cnv_proxy_high[g232_epi$disease_group == "PTC"], na.rm = TRUE)),
    p_wilcox = c(wt1$p.value, wt2$p.value)
  )
} else data.table()

all_tests <- list(cor_tests = cor_tests, status_tests = status_tests,
                  g184_sample_tests = g184_tests, g232_sample_tests = g232_tests)
saveRDS(all_tests, file.path(result_dir, "single_cell_cnv_proxy_tests.rds"))
fwrite(cor_tests, file.path(result_dir, "single_cell_cnv_proxy_cor_tests.csv"))
fwrite(status_tests, file.path(result_dir, "single_cell_cnv_proxy_status_tests.csv"))
fwrite(rbindlist(list(g184_tests, g232_tests), fill = TRUE),
       file.path(result_dir, "single_cell_cnv_proxy_group_tests.csv"))

p1 <- ggplot(cnv, aes(x = cnv_proxy_amp, y = RAIR_like, color = cnv_proxy_high)) +
  geom_point(size = 0.35, alpha = 0.25) +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_manual(values = c("FALSE" = "#6f7a85", "TRUE" = "#b73e3e")) +
  theme_bw(base_size = 10) +
  labs(x = "Lightweight CNV proxy amplitude", y = "RAIR-like score",
       color = "CNV-high proxy",
       title = "RAIR-like state across epithelial/thyroid cells by CNV proxy")
ggsave(file.path(figure_dir, "single_cell_cnv_proxy_vs_RAIR_epithelial.png"),
       p1, width = 7.2, height = 4.6, dpi = 300)

p2 <- ggplot(sample_summary[n_epithelial >= 50],
             aes(x = tissue_group, y = fraction_cnv_proxy_high, color = disease_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.15, height = 0, size = 2.0, alpha = 0.8) +
  facet_wrap(~ dataset, scales = "free_x") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = "Fraction of epithelial/thyroid cells with CNV-high proxy",
       color = "Disease group",
       title = "Sample-level CNV proxy burden in epithelial/thyroid cells")
ggsave(file.path(figure_dir, "single_cell_cnv_proxy_high_fraction_by_sample.png"),
       p2, width = 7.4, height = 4.5, dpi = 300)

p3 <- ggplot(cnv, aes(x = cnv_proxy_high, y = RAIR_like, fill = cnv_proxy_high)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.18, size = 0.25, alpha = 0.18) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "#6f7a85", "TRUE" = "#b73e3e")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = "CNV-high proxy", y = "RAIR-like score",
       title = "RAIR-like score in epithelial/thyroid cells stratified by CNV proxy")
ggsave(file.path(figure_dir, "single_cell_RAIR_by_cnv_proxy_status.png"),
       p3, width = 6.2, height = 4.4, dpi = 300)

sink(file.path(result_dir, "single_cell_cnv_proxy_summary.txt"))
cat("Single-cell lightweight CNV proxy analysis\n\n")
cat("Method guardrail:\n")
cat("This is not inferCNV/copyKAT. It orders genes by TCGA/GENCODE genomic coordinates, centers log-normalized expression on non-epithelial reference cells, smooths across 100-gene genomic windows, and reports per-cell absolute deviation amplitude.\n")
cat("Use it as a reproducible proxy to test whether RAIR-like epithelial/thyroid cells are enriched for malignant-like genome-wide expression imbalance.\n\n")
cat("Epithelial/thyroid cells scored:\n")
print(cnv[, .N, by = dataset], row.names = FALSE)
cat("\nCorrelation between RAIR-like and CNV proxy amplitude:\n")
print(cor_tests, row.names = FALSE)
cat("\nRAIR-like score by CNV proxy high/low status:\n")
print(status_tests, row.names = FALSE)
cat("\nSample-level group tests:\n")
print(rbindlist(list(g184_tests, g232_tests), fill = TRUE), row.names = FALSE)
cat("\nSample-level summary, first rows:\n")
print(head(sample_summary[order(dataset, sample)], 20), row.names = FALSE)
sink()

message("Saved CNV proxy analysis outputs.")
