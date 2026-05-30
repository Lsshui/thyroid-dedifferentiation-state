options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(clariomshumantranscriptcluster.db)
  library(ggplot2)
})

analysis_root <- "D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/analysis"
pre_root <- "D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/preanalysis"
data_dir <- file.path(pre_root, "data", "geo")
result_dir <- file.path(analysis_root, "results")
figure_dir <- file.path(analysis_root, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

target_priority <- read.csv(file.path(result_dir, "TCGA_THCA_target_priority_table.csv"))
candidate_genes <- unique(c(head(target_priority$gene, 20),
                            "SHMT2", "MTHFD2", "MTHFD1L", "SLC1A5", "SLC7A5",
                            "MET", "AXL", "DUSP5", "DUSP6", "FN1", "HMGA2", "KRT19"))

load_eset <- function(gse_id) {
  x <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
  if (is.list(x)) x <- x[[1]]
  x
}

maybe_log2 <- function(x) {
  q <- quantile(as.numeric(x), c(0.01, 0.99), na.rm = TRUE)
  if (q[[2]] > 100) log2(x + 1) else x
}

collapse_to_gene <- function(expr, platform) {
  db <- switch(platform,
               GPL23159 = clariomshumantranscriptcluster.db,
               GPL570 = hgu133plus2.db,
               stop("Unsupported platform"))
  ann <- AnnotationDbi::select(db, keys = rownames(expr), keytype = "PROBEID", columns = "SYMBOL")
  ann <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", c("PROBEID", "SYMBOL")]
  ann <- ann[ann$PROBEID %in% rownames(expr), ]
  ann <- ann[!duplicated(ann), ]
  v <- apply(expr[unique(ann$PROBEID), , drop = FALSE], 1, var, na.rm = TRUE)
  ann$probe_var <- v[ann$PROBEID]
  ann <- ann[order(ann$SYMBOL, -ann$probe_var), ]
  ann <- ann[!duplicated(ann$SYMBOL), ]
  out <- expr[ann$PROBEID, , drop = FALSE]
  rownames(out) <- ann$SYMBOL
  out
}

test_gene <- function(expr_gene, meta, group_col, positive, negative, dataset, comparison) {
  genes <- intersect(candidate_genes, rownames(expr_gene))
  do.call(rbind, lapply(genes, function(g) {
    sub <- meta[meta[[group_col]] %in% c(positive, negative), ]
    if (length(unique(sub[[group_col]])) < 2) return(NULL)
    e <- as.numeric(expr_gene[g, sub$sample])
    grp <- factor(sub[[group_col]], levels = c(negative, positive))
    data.frame(
      dataset = dataset,
      comparison = comparison,
      gene = g,
      n_positive = sum(grp == positive),
      n_negative = sum(grp == negative),
      median_positive = median(e[grp == positive], na.rm = TRUE),
      median_negative = median(e[grp == negative], na.rm = TRUE),
      delta_positive_negative = median(e[grp == positive], na.rm = TRUE) -
        median(e[grp == negative], na.rm = TRUE),
      p_wilcox = suppressWarnings(wilcox.test(e ~ grp)$p.value),
      stringsAsFactors = FALSE
    )
  }))
}

## GSE33630
g336 <- load_eset("GSE33630")
expr336 <- collapse_to_gene(maybe_log2(exprs(g336)), "GPL570")
pd336 <- pData(g336)
grp336 <- ifelse(grepl("anaplastic", pd336$characteristics_ch1, ignore.case = TRUE), "ATC",
                 ifelse(grepl("papillary", pd336$characteristics_ch1, ignore.case = TRUE), "PTC",
                        ifelse(grepl("normal", paste(pd336$characteristics_ch1, pd336$source_name_ch1), ignore.case = TRUE),
                               "Normal", NA)))
meta336 <- data.frame(sample = rownames(pd336), group = grp336, stringsAsFactors = FALSE)
res336_atc_ptc <- test_gene(expr336, meta336, "group", "ATC", "PTC", "GSE33630", "ATC_vs_PTC")
res336_ptc_normal <- test_gene(expr336, meta336, "group", "PTC", "Normal", "GSE33630", "PTC_vs_Normal")

## GSE76039
g760 <- load_eset("GSE76039")
expr760 <- collapse_to_gene(maybe_log2(exprs(g760)), "GPL570")
pd760 <- pData(g760)
grp760 <- ifelse(grepl("Anaplastic", pd760$source_name_ch1, ignore.case = TRUE), "ATC",
                 ifelse(grepl("Poorly", pd760$source_name_ch1, ignore.case = TRUE), "PDTC", NA))
meta760 <- data.frame(sample = rownames(pd760), group = grp760, stringsAsFactors = FALSE)
res760 <- test_gene(expr760, meta760, "group", "ATC", "PDTC", "GSE76039", "ATC_vs_PDTC")

validation <- rbind(res336_atc_ptc, res336_ptc_normal, res760)
validation$p_adj <- ave(validation$p_wilcox, validation$dataset, validation$comparison,
                        FUN = function(p) p.adjust(p, method = "BH"))

tcga <- target_priority[, c("gene", "delta_high_low", "rho_RAIR_like", "rho_TDS", "priority_score")]
names(tcga) <- c("gene", "TCGA_delta_RAIR_high_low", "TCGA_rho_RAIR_like", "TCGA_rho_TDS", "TCGA_priority_score")
validation_merged <- merge(validation, tcga, by = "gene", all.x = TRUE)

wide_delta <- reshape(validation_merged[, c("gene", "dataset", "comparison", "delta_positive_negative")],
                      idvar = "gene", timevar = c("dataset", "comparison"),
                      direction = "wide")
names(wide_delta) <- sub("^delta_positive_negative\\.", "", names(wide_delta))
wide_delta <- merge(tcga, wide_delta, by = "gene", all.x = TRUE)

score_consistency <- validation_merged
score_consistency$up_validated <- score_consistency$delta_positive_negative > 0 & score_consistency$p_adj < 0.05
consistent <- aggregate(up_validated ~ gene, data = score_consistency, FUN = sum)
names(consistent)[2] <- "n_external_up_validations"
consistent <- merge(tcga, consistent, by = "gene", all.x = TRUE)
consistent <- consistent[order(-consistent$n_external_up_validations, consistent$TCGA_priority_score), ]

write.csv(validation_merged, file.path(result_dir, "external_target_validation_long.csv"), row.names = FALSE)
write.csv(wide_delta, file.path(result_dir, "external_target_validation_delta_wide.csv"), row.names = FALSE)
write.csv(consistent, file.path(result_dir, "cross_cohort_target_consistency.csv"), row.names = FALSE)

plot_genes <- head(consistent$gene, 18)
plot_df <- validation_merged[validation_merged$gene %in% plot_genes &
                               validation_merged$comparison %in% c("ATC_vs_PTC", "ATC_vs_PDTC"), ]
p <- ggplot(plot_df, aes(x = reorder(gene, delta_positive_negative), y = delta_positive_negative,
                         fill = comparison)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Median expression difference", title = "External validation of TCGA-prioritized target/state genes")
ggsave(file.path(figure_dir, "external_target_validation_delta.png"), p, width = 7.0, height = 5.4, dpi = 300)

sink(file.path(result_dir, "external_target_validation_summary.txt"))
cat("External target/state validation\n\n")
cat("Most consistent genes:\n")
print(head(consistent, 20))
cat("\nValidation long table, selected one-carbon and therapeutic genes:\n")
print(validation_merged[validation_merged$gene %in% c("SHMT2", "MTHFD2", "MTHFD1L", "SLC1A5", "SLC7A5", "MET", "AXL", "DUSP6", "FN1"), ])
sink()

message("External target validation complete.")
