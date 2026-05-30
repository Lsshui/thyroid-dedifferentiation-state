options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(clariomshumantranscriptcluster.db)
  library(hgu133plus2.db)
  library(ggplot2)
  library(limma)
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

load_eset <- function(gse_id) {
  message("Loading ", gse_id)
  eset <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
  if (is.list(eset)) eset <- eset[[1]]
  eset
}

maybe_log2 <- function(x) {
  q <- quantile(as.numeric(x), c(0.01, 0.5, 0.99), na.rm = TRUE)
  if (q[[3]] > 100) log2(x + 1) else x
}

collapse_to_gene <- function(expr, platform) {
  db <- switch(platform,
               GPL23159 = clariomshumantranscriptcluster.db,
               GPL570 = hgu133plus2.db,
               stop("Unsupported platform: ", platform))
  ann <- AnnotationDbi::select(
    db,
    keys = rownames(expr),
    keytype = "PROBEID",
    columns = c("SYMBOL")
  )
  ann <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", c("PROBEID", "SYMBOL")]
  ann <- ann[ann$PROBEID %in% rownames(expr), ]
  ann <- ann[!duplicated(ann), ]
  probe_var <- apply(expr[unique(ann$PROBEID), , drop = FALSE], 1, var, na.rm = TRUE)
  ann$probe_var <- probe_var[ann$PROBEID]
  ann <- ann[order(ann$SYMBOL, -ann$probe_var), ]
  ann <- ann[!duplicated(ann$SYMBOL), ]
  out <- expr[ann$PROBEID, , drop = FALSE]
  rownames(out) <- ann$SYMBOL
  out
}

module_score <- function(expr_gene, genes) {
  genes <- intersect(genes, rownames(expr_gene))
  if (length(genes) < 2) {
    return(rep(NA_real_, ncol(expr_gene)))
  }
  z <- t(scale(t(expr_gene[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

score_dataset <- function(eset, platform, dataset) {
  expr <- maybe_log2(exprs(eset))
  expr_gene <- collapse_to_gene(expr, platform)
  scores <- data.frame(
    sample = colnames(expr_gene),
    dataset = dataset,
    TDS = module_score(expr_gene, gene_sets$TDS),
    OneCarbon = module_score(expr_gene, gene_sets$OneCarbon),
    MAPK = module_score(expr_gene, gene_sets$MAPK),
    Dediff_EMT = module_score(expr_gene, gene_sets$Dediff_EMT),
    stringsAsFactors = FALSE
  )
  scores$RAIR_like <- as.numeric(scale(scores$OneCarbon)) +
    as.numeric(scale(scores$MAPK)) +
    as.numeric(scale(scores$Dediff_EMT)) -
    as.numeric(scale(scores$TDS))

  coverage <- do.call(rbind, lapply(names(gene_sets), function(nm) {
    data.frame(
      dataset = dataset,
      module = nm,
      n_available = length(intersect(gene_sets[[nm]], rownames(expr_gene))),
      available_genes = paste(intersect(gene_sets[[nm]], rownames(expr_gene)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  list(expr_gene = expr_gene, scores = scores, coverage = coverage)
}

auc_binary <- function(score, group, positive) {
  ok <- !is.na(score) & !is.na(group)
  score <- score[ok]
  group <- group[ok]
  pos <- group == positive
  n_pos <- sum(pos)
  n_neg <- sum(!pos)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

compare_two <- function(df, group_col, score_col, positive, negative, expected = "higher") {
  sub <- df[df[[group_col]] %in% c(positive, negative), ]
  g <- factor(sub[[group_col]], levels = c(negative, positive))
  s <- sub[[score_col]]
  if (length(unique(g[!is.na(s)])) < 2) {
    return(data.frame(
      score = score_col, positive = positive, negative = negative,
      n_positive = sum(g == positive, na.rm = TRUE),
      n_negative = sum(g == negative, na.rm = TRUE),
      median_positive = NA_real_, median_negative = NA_real_,
      delta_median = NA_real_, p_wilcox = NA_real_, auc = NA_real_,
      expected = expected
    ))
  }
  direction_score <- if (expected == "lower") -s else s
  data.frame(
    score = score_col,
    positive = positive,
    negative = negative,
    n_positive = sum(g == positive, na.rm = TRUE),
    n_negative = sum(g == negative, na.rm = TRUE),
    median_positive = median(s[g == positive], na.rm = TRUE),
    median_negative = median(s[g == negative], na.rm = TRUE),
    delta_median = median(s[g == positive], na.rm = TRUE) - median(s[g == negative], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(s ~ g)$p.value),
    auc = auc_binary(direction_score, as.character(g), positive),
    expected = expected,
    stringsAsFactors = FALSE
  )
}

plot_score <- function(df, x, y, filename, title) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.65) +
    geom_jitter(width = 0.12, size = 1.7, alpha = 0.75) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(x = NULL, y = y, title = title)
  ggsave(file.path(figure_dir, filename), p, width = 5.5, height = 4.1, dpi = 300)
}

## GSE151179: RAI-refractory vs RAI-avid/remission
g151 <- load_eset("GSE151179")
s151 <- score_dataset(g151, "GPL23159", "GSE151179")
pd151 <- pData(g151)
meta151 <- data.frame(
  sample = rownames(pd151),
  title = pd151$title,
  source = pd151$source_name_ch1,
  tissue_type = pd151$`tissue type:ch1`,
  collection_rai = pd151$`collection before/after rai:ch1`,
  disease = pd151$`disease:ch1`,
  rai_response = pd151$`patient rai responce:ch1`,
  rai_uptake = pd151$`rai uptake at the metastatic site:ch1`,
  purity = pd151$`tumor purity class by cibersort:ch1`,
  stringsAsFactors = FALSE
)
df151 <- merge(meta151, s151$scores, by = "sample")
df151$rai_group <- ifelse(df151$source == "papillary thyroid carcinoma" &
                            df151$rai_response == "Refractory", "RAI_refractory",
                          ifelse(df151$source == "papillary thyroid carcinoma" &
                                   df151$rai_response == "Avid", "RAI_avid", NA))
df151$rai_group_remission <- ifelse(df151$source == "papillary thyroid carcinoma" &
                                      df151$rai_response == "Refractory", "RAI_refractory",
                                    ifelse(df151$source == "papillary thyroid carcinoma" &
                                             df151$rai_response == "Avid" &
                                             df151$disease == "Remission",
                                           "RAI_avid_remission", NA))

tests151 <- rbind(
  compare_two(df151[!is.na(df151$rai_group), ], "rai_group", "TDS",
              "RAI_refractory", "RAI_avid", expected = "lower"),
  compare_two(df151[!is.na(df151$rai_group), ], "rai_group", "RAIR_like",
              "RAI_refractory", "RAI_avid", expected = "higher"),
  compare_two(df151[!is.na(df151$rai_group), ], "rai_group", "OneCarbon",
              "RAI_refractory", "RAI_avid", expected = "higher"),
  compare_two(df151[!is.na(df151$rai_group), ], "rai_group", "MAPK",
              "RAI_refractory", "RAI_avid", expected = "higher"),
  compare_two(df151[!is.na(df151$rai_group), ], "rai_group", "Dediff_EMT",
              "RAI_refractory", "RAI_avid", expected = "higher"),
  compare_two(df151[!is.na(df151$rai_group_remission), ], "rai_group_remission", "TDS",
              "RAI_refractory", "RAI_avid_remission", expected = "lower"),
  compare_two(df151[!is.na(df151$rai_group_remission), ], "rai_group_remission", "RAIR_like",
              "RAI_refractory", "RAI_avid_remission", expected = "higher")
)

## GSE33630: normal/PTC/ATC
g336 <- load_eset("GSE33630")
s336 <- score_dataset(g336, "GPL570", "GSE33630")
pd336 <- pData(g336)
diag336 <- pd336$characteristics_ch1
grp336 <- ifelse(grepl("anaplastic", diag336, ignore.case = TRUE), "ATC",
                 ifelse(grepl("papillary", diag336, ignore.case = TRUE), "PTC",
                        ifelse(grepl("normal|non-tumor|non tumour", paste(diag336, pd336$source_name_ch1), ignore.case = TRUE),
                               "Normal", NA)))
df336 <- merge(data.frame(sample = rownames(pd336), title = pd336$title,
                          diagnostic = diag336, histology_group = grp336,
                          stringsAsFactors = FALSE),
               s336$scores, by = "sample")
df336$histology_group <- factor(df336$histology_group, levels = c("Normal", "PTC", "ATC"))

tests336 <- rbind(
  compare_two(df336, "histology_group", "TDS", "ATC", "PTC", expected = "lower"),
  compare_two(df336, "histology_group", "RAIR_like", "ATC", "PTC", expected = "higher"),
  compare_two(df336, "histology_group", "OneCarbon", "ATC", "PTC", expected = "higher"),
  compare_two(df336, "histology_group", "MAPK", "ATC", "PTC", expected = "higher"),
  compare_two(df336, "histology_group", "Dediff_EMT", "ATC", "PTC", expected = "higher")
)
ord336 <- as.numeric(df336$histology_group)
trend336 <- data.frame(
  dataset = "GSE33630",
  score = c("TDS", "RAIR_like", "OneCarbon", "MAPK", "Dediff_EMT"),
  spearman_rho = sapply(c("TDS", "RAIR_like", "OneCarbon", "MAPK", "Dediff_EMT"),
                        function(v) suppressWarnings(cor(ord336, df336[[v]], method = "spearman", use = "complete.obs"))),
  p_spearman = sapply(c("TDS", "RAIR_like", "OneCarbon", "MAPK", "Dediff_EMT"),
                      function(v) suppressWarnings(cor.test(ord336, df336[[v]], method = "spearman")$p.value)),
  stringsAsFactors = FALSE
)

## GSE76039: PDTC/ATC
g760 <- load_eset("GSE76039")
s760 <- score_dataset(g760, "GPL570", "GSE76039")
pd760 <- pData(g760)
src760 <- pd760$source_name_ch1
grp760 <- ifelse(grepl("Anaplastic", src760, ignore.case = TRUE), "ATC",
                 ifelse(grepl("Poorly", src760, ignore.case = TRUE), "PDTC",
                        ifelse(grepl("Papillary", src760, ignore.case = TRUE), "PTC", NA)))
df760 <- merge(data.frame(sample = rownames(pd760), title = pd760$title,
                          source = src760, histology_group = grp760,
                          stringsAsFactors = FALSE),
               s760$scores, by = "sample")
df760$histology_group <- factor(df760$histology_group, levels = c("PTC", "PDTC", "ATC"))

tests760 <- rbind(
  compare_two(df760, "histology_group", "TDS", "ATC", "PDTC", expected = "lower"),
  compare_two(df760, "histology_group", "RAIR_like", "ATC", "PDTC", expected = "higher"),
  compare_two(df760, "histology_group", "OneCarbon", "ATC", "PDTC", expected = "higher"),
  compare_two(df760, "histology_group", "MAPK", "ATC", "PDTC", expected = "higher"),
  compare_two(df760, "histology_group", "Dediff_EMT", "ATC", "PDTC", expected = "higher")
)

bulk_tests <- rbind(
  cbind(dataset_test = "GSE151179_RAI", tests151),
  cbind(dataset_test = "GSE33630_ATC_vs_PTC", tests336),
  cbind(dataset_test = "GSE76039_ATC_vs_PDTC", tests760)
)
coverage <- rbind(s151$coverage, s336$coverage, s760$coverage)

write.csv(df151, file.path(result_dir, "GSE151179_scores_metadata.csv"), row.names = FALSE)
write.csv(df336, file.path(result_dir, "GSE33630_scores_metadata.csv"), row.names = FALSE)
write.csv(df760, file.path(result_dir, "GSE76039_scores_metadata.csv"), row.names = FALSE)
write.csv(bulk_tests, file.path(result_dir, "bulk_score_group_tests.csv"), row.names = FALSE)
write.csv(trend336, file.path(result_dir, "GSE33630_histology_trend_tests.csv"), row.names = FALSE)
write.csv(coverage, file.path(result_dir, "bulk_gene_set_coverage.csv"), row.names = FALSE)

plot_score(df151[!is.na(df151$rai_group), ], "rai_group", "TDS",
           "GSE151179_TDS_by_RAI.png", "GSE151179 TDS by RAI response")
plot_score(df151[!is.na(df151$rai_group), ], "rai_group", "RAIR_like",
           "GSE151179_RAIR_like_by_RAI.png", "GSE151179 RAIR-like score by RAI response")
plot_score(df336[!is.na(df336$histology_group), ], "histology_group", "TDS",
           "GSE33630_TDS_by_histology.png", "GSE33630 TDS by histology")
plot_score(df336[!is.na(df336$histology_group), ], "histology_group", "RAIR_like",
           "GSE33630_RAIR_like_by_histology.png", "GSE33630 RAIR-like score by histology")
plot_score(df760[!is.na(df760$histology_group), ], "histology_group", "TDS",
           "GSE76039_TDS_by_histology.png", "GSE76039 TDS by histology")
plot_score(df760[!is.na(df760$histology_group), ], "histology_group", "RAIR_like",
           "GSE76039_RAIR_like_by_histology.png", "GSE76039 RAIR-like score by histology")

sink(file.path(result_dir, "bulk_preanalysis_summary.txt"))
cat("Bulk preanalysis completed\n\n")
cat("GSE151179 RAI groups:\n")
print(table(df151$rai_group, useNA = "ifany"))
cat("\nGSE151179 tissue types by RAI group:\n")
print(table(df151$rai_group, df151$tissue_type, useNA = "ifany"))
cat("\nGSE33630 histology groups:\n")
print(table(df336$histology_group, useNA = "ifany"))
cat("\nGSE76039 histology groups:\n")
print(table(df760$histology_group, useNA = "ifany"))
cat("\nGene-set coverage:\n")
print(coverage)
cat("\nGroup tests:\n")
print(bulk_tests)
cat("\nGSE33630 ordinal trend tests:\n")
print(trend336)
sink()

message("Bulk preanalysis complete. Results written to: ", result_dir)
