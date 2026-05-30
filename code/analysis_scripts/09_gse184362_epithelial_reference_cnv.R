options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(ggplot2)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pre_root <- file.path(root, "preanalysis")
analysis_root <- file.path(root, "analysis")
data_dir <- file.path(pre_root, "data", "geo", "GSE184362_RAW")
result_dir <- file.path(analysis_root, "results")
figure_dir <- file.path(analysis_root, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

window_genes <- 100

message("Loading per-cell metadata...")
cell_scores <- as.data.table(readRDS(file.path(pre_root, "results", "single_cell_per_cell_scores.rds")))
cell_scores <- cell_scores[dataset == "GSE184362"]

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
ordered_genes <- gene_pos$gene

read_gse184_matrix <- function(matrix_file) {
  prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
  sample_id <- sub("^GSM[0-9]+_", "", prefix)
  feature_file <- file.path(dirname(matrix_file), paste0(prefix, "_features.tsv.gz"))
  barcode_file <- file.path(dirname(matrix_file), paste0(prefix, "_barcodes.tsv.gz"))
  features <- fread(feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
  genes <- features[[2]]
  barcodes <- fread(barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE)[[1]]
  mat <- Matrix::readMM(gzfile(matrix_file))
  colnames(mat) <- paste(sample_id, barcodes, sep = "_")
  keep <- !is.na(genes) & genes != "" & genes %in% ordered_genes & !duplicated(genes)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- genes[keep]
  list(sample_id = sample_id, mat = as(mat, "dgCMatrix"))
}

log_normalize <- function(mat) {
  lib <- Matrix::colSums(mat)
  keep <- lib > 0
  mat <- mat[, keep, drop = FALSE]
  norm <- t(t(mat) / Matrix::colSums(mat)) * 10000
  norm@x <- log1p(norm@x)
  norm
}

message("Building paratumor epithelial/thyroid reference...")
matrix_files <- list.files(data_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
ref_sum <- numeric(length(ordered_genes))
ref_n <- numeric(length(ordered_genes))
names(ref_sum) <- ordered_genes
names(ref_n) <- ordered_genes

for (mf in matrix_files) {
  m <- read_gse184_matrix(mf)
  meta <- cell_scores[sample == m$sample_id & tissue_group == "Paratumor" &
                        cell_type == "Epithelial_Thyroid" & cell %in% colnames(m$mat)]
  if (nrow(meta) == 0) next
  message("Reference cells: ", m$sample_id, " n=", nrow(meta))
  norm <- log_normalize(m$mat[, meta$cell, drop = FALSE])
  idx <- match(rownames(norm), ordered_genes)
  ref_sum[idx] <- ref_sum[idx] + Matrix::rowSums(norm)
  ref_n[idx] <- ref_n[idx] + ncol(norm)
  rm(m, norm)
  gc()
}

ref_mean <- ref_sum / ref_n
ref_valid <- is.finite(ref_mean) & ref_n >= 100
message("Reference genes retained: ", sum(ref_valid), "; reference epithelial cells represented per retained gene median n = ",
        median(ref_n[ref_valid]))
if (sum(ref_valid) < 5000) stop("Too few reference genes for epithelial-reference CNV proxy.")

score_sample <- function(matrix_file) {
  m <- read_gse184_matrix(matrix_file)
  meta <- cell_scores[sample == m$sample_id & cell_type == "Epithelial_Thyroid" &
                        cell %in% colnames(m$mat)]
  if (nrow(meta) == 0) return(NULL)
  message("Scoring epithelial-reference CNV: ", m$sample_id, " n=", nrow(meta))
  norm <- log_normalize(m$mat[, meta$cell, drop = FALSE])
  common <- intersect(rownames(norm), names(ref_mean)[ref_valid])
  norm <- norm[common, , drop = FALSE]
  pos <- gene_pos[match(common, gene)]
  ref_vec <- ref_mean[common]
  bins <- split(seq_len(nrow(norm)), pos$bin_order)
  bin_means <- vapply(bins, function(idx) {
    Matrix::colMeans(norm[idx, , drop = FALSE]) - mean(ref_vec[idx], na.rm = TRUE)
  }, FUN.VALUE = numeric(ncol(norm)))
  bin_means <- t(bin_means)
  amp <- colMeans(abs(bin_means), na.rm = TRUE)
  out <- data.table(
    dataset = "GSE184362",
    sample = m$sample_id,
    cell = colnames(norm),
    cnv_epithelial_ref_amp = as.numeric(amp)
  )
  out <- merge(out, meta[, .(cell, disease_group, tissue_group, cell_type, RAIR_like,
                             TDS_raw, OneCarbon_raw, MAPK_raw, Dediff_EMT_raw)],
               by = "cell", all.x = TRUE)
  rm(m, norm)
  gc()
  out
}

cnv <- rbindlist(lapply(matrix_files, score_sample), fill = TRUE)
if (nrow(cnv) == 0) stop("No epithelial-reference CNV proxy results produced.")

ref_amp <- cnv[tissue_group == "Paratumor", cnv_epithelial_ref_amp]
ref_med <- median(ref_amp, na.rm = TRUE)
ref_mad <- mad(ref_amp, constant = 1.4826, na.rm = TRUE)
if (!is.finite(ref_mad) || ref_mad <= 0) ref_mad <- sd(ref_amp, na.rm = TRUE)
if (!is.finite(ref_mad) || ref_mad <= 0) ref_mad <- 1e-6
cnv[, cnv_epithelial_ref_z := (cnv_epithelial_ref_amp - ref_med) / ref_mad]
cnv[, cnv_epithelial_ref_high := cnv_epithelial_ref_z > 3]

fwrite(cnv, file.path(result_dir, "GSE184362_epithelial_reference_cnv_proxy_cells.csv"))

sample_summary <- cnv[, .(
  n_epithelial = .N,
  median_RAIR_like = median(RAIR_like, na.rm = TRUE),
  median_cnv_epithelial_ref_z = median(cnv_epithelial_ref_z, na.rm = TRUE),
  fraction_cnv_epithelial_ref_high = mean(cnv_epithelial_ref_high, na.rm = TRUE),
  fraction_RAIR_top_quartile_cnv_high = {
    q <- quantile(RAIR_like, 0.75, na.rm = TRUE)
    mean(cnv_epithelial_ref_high[RAIR_like >= q], na.rm = TRUE)
  }
), by = .(dataset, sample, disease_group, tissue_group)]
sample_summary[, tumor_met_group := ifelse(tissue_group == "Paratumor", "Paratumor", "Tumor_met")]
fwrite(sample_summary, file.path(result_dir, "GSE184362_epithelial_reference_cnv_proxy_sample_summary.csv"))

cell_cor <- cnv[, {
  ok <- is.finite(RAIR_like) & is.finite(cnv_epithelial_ref_amp)
  ct <- suppressWarnings(cor.test(RAIR_like[ok], cnv_epithelial_ref_amp[ok], method = "spearman"))
  data.table(n_cells = sum(ok), spearman_rho = unname(ct$estimate), p_value = ct$p.value)
}]

status_test <- if (length(unique(cnv$cnv_epithelial_ref_high)) == 2) {
  wt <- suppressWarnings(wilcox.test(RAIR_like ~ cnv_epithelial_ref_high, data = cnv, exact = FALSE))
  data.table(
    n_high = sum(cnv$cnv_epithelial_ref_high),
    n_low = sum(!cnv$cnv_epithelial_ref_high),
    median_RAIR_high = median(cnv$RAIR_like[cnv$cnv_epithelial_ref_high], na.rm = TRUE),
    median_RAIR_low = median(cnv$RAIR_like[!cnv$cnv_epithelial_ref_high], na.rm = TRUE),
    delta_RAIR_high_minus_low = median(cnv$RAIR_like[cnv$cnv_epithelial_ref_high], na.rm = TRUE) -
      median(cnv$RAIR_like[!cnv$cnv_epithelial_ref_high], na.rm = TRUE),
    p_wilcox = wt$p.value
  )
} else data.table()

sample_for_tests <- sample_summary[n_epithelial >= 50]
group_tests <- if (length(unique(sample_for_tests$tumor_met_group)) == 2) {
  wt1 <- suppressWarnings(wilcox.test(median_cnv_epithelial_ref_z ~ tumor_met_group,
                                      data = sample_for_tests, exact = FALSE))
  wt2 <- suppressWarnings(wilcox.test(fraction_cnv_epithelial_ref_high ~ tumor_met_group,
                                      data = sample_for_tests, exact = FALSE))
  wt3 <- suppressWarnings(wilcox.test(median_RAIR_like ~ tumor_met_group,
                                      data = sample_for_tests, exact = FALSE))
  data.table(
    comparison = c("Tumor/met vs paratumor median CNV-z",
                   "Tumor/met vs paratumor CNV-high fraction",
                   "Tumor/met vs paratumor median RAIR-like"),
    n_tumor_met = sum(sample_for_tests$tumor_met_group == "Tumor_met"),
    n_paratumor = sum(sample_for_tests$tumor_met_group == "Paratumor"),
    median_tumor_met = c(median(sample_for_tests$median_cnv_epithelial_ref_z[sample_for_tests$tumor_met_group == "Tumor_met"], na.rm = TRUE),
                         median(sample_for_tests$fraction_cnv_epithelial_ref_high[sample_for_tests$tumor_met_group == "Tumor_met"], na.rm = TRUE),
                         median(sample_for_tests$median_RAIR_like[sample_for_tests$tumor_met_group == "Tumor_met"], na.rm = TRUE)),
    median_paratumor = c(median(sample_for_tests$median_cnv_epithelial_ref_z[sample_for_tests$tumor_met_group == "Paratumor"], na.rm = TRUE),
                         median(sample_for_tests$fraction_cnv_epithelial_ref_high[sample_for_tests$tumor_met_group == "Paratumor"], na.rm = TRUE),
                         median(sample_for_tests$median_RAIR_like[sample_for_tests$tumor_met_group == "Paratumor"], na.rm = TRUE)),
    p_wilcox = c(wt1$p.value, wt2$p.value, wt3$p.value)
  )
} else data.table()

fwrite(cell_cor, file.path(result_dir, "GSE184362_epithelial_reference_cnv_proxy_cor.csv"))
fwrite(status_test, file.path(result_dir, "GSE184362_epithelial_reference_cnv_proxy_status_test.csv"))
fwrite(group_tests, file.path(result_dir, "GSE184362_epithelial_reference_cnv_proxy_group_tests.csv"))

p1 <- ggplot(cnv, aes(x = cnv_epithelial_ref_amp, y = RAIR_like,
                      color = tissue_group == "Paratumor")) +
  geom_point(size = 0.35, alpha = 0.25) +
  scale_color_manual(values = c("FALSE" = "#b73e3e", "TRUE" = "#3274a1"),
                     labels = c("FALSE" = "Tumor/met", "TRUE" = "Paratumor")) +
  theme_bw(base_size = 10) +
  labs(x = "CNV proxy amplitude vs paratumor epithelial reference",
       y = "RAIR-like score", color = "Group",
       title = "GSE184362 epithelial-reference CNV proxy")
ggsave(file.path(figure_dir, "GSE184362_epithelial_reference_CNV_vs_RAIR.png"),
       p1, width = 6.6, height = 4.6, dpi = 300)

p2 <- ggplot(sample_for_tests,
             aes(x = tissue_group, y = fraction_cnv_epithelial_ref_high,
                 color = tumor_met_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.14, size = 2.0, alpha = 0.8) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Fraction CNV-high epithelial/thyroid cells",
       title = "GSE184362 CNV-high burden using paratumor epithelial reference")
ggsave(file.path(figure_dir, "GSE184362_epithelial_reference_CNV_high_fraction.png"),
       p2, width = 6.4, height = 4.4, dpi = 300)

p3 <- ggplot(sample_for_tests, aes(x = median_cnv_epithelial_ref_z, y = median_RAIR_like,
                                   color = tissue_group)) +
  geom_point(size = 2.4, alpha = 0.85) +
  theme_bw(base_size = 10) +
  labs(x = "Sample median CNV proxy z",
       y = "Sample median RAIR-like",
       color = "Tissue group",
       title = "GSE184362 sample-level RAIR-like and epithelial-reference CNV proxy")
ggsave(file.path(figure_dir, "GSE184362_sample_CNV_z_vs_RAIR.png"),
       p3, width = 5.8, height = 4.4, dpi = 300)

sink(file.path(result_dir, "GSE184362_epithelial_reference_cnv_proxy_summary.txt"))
cat("GSE184362 epithelial-reference lightweight CNV proxy\n\n")
cat("Reference: paratumor Epithelial_Thyroid cells only.\n")
cat("Reference genes retained:", sum(ref_valid), "\n")
cat("Reference epithelial cell representation per retained gene median n:", median(ref_n[ref_valid]), "\n")
cat("Cells scored:\n")
print(cnv[, .N, by = tissue_group], row.names = FALSE)
cat("\nCell-level correlation between RAIR-like and epithelial-reference CNV amplitude:\n")
print(cell_cor, row.names = FALSE)
cat("\nRAIR-like by CNV-high proxy status:\n")
print(status_test, row.names = FALSE)
cat("\nSample-level tests, requiring >=50 epithelial/thyroid cells:\n")
print(group_tests, row.names = FALSE)
cat("\nSample-level summary:\n")
print(sample_summary[order(tissue_group, sample)], row.names = FALSE)
sink()

message("Saved GSE184362 epithelial-reference CNV proxy outputs.")
