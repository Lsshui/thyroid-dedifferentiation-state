options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pre_result_dir <- file.path(root, "preanalysis", "results")
result_dir <- file.path(root, "analysis", "results")
figure_dir <- file.path(root, "analysis", "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sample_medians <- fread(file.path(pre_result_dir, "single_cell_sample_celltype_medians.csv"))
top_dist <- fread(file.path(pre_result_dir, "single_cell_top_RAIR_like_celltype_distribution.csv"))

top_summary <- top_dist[, .(
  top10_fraction_epithelial_thyroid = sum(fraction_of_top10[cell_type == "Epithelial_Thyroid"], na.rm = TRUE),
  top10_fraction_myeloid = sum(fraction_of_top10[cell_type == "Myeloid"], na.rm = TRUE),
  top10_fraction_caf = sum(fraction_of_top10[cell_type == "CAF"], na.rm = TRUE),
  top10_fraction_t_nk = sum(fraction_of_top10[cell_type == "T_NK"], na.rm = TRUE),
  top10_total_cells = sum(n_cells, na.rm = TRUE)
), by = dataset]

test_one <- function(df, dataset_name, group_col, group_a, group_b, label) {
  sub <- df[dataset == dataset_name & cell_type == "Epithelial_Thyroid" &
              n_cells >= 50 & get(group_col) %in% c(group_a, group_b)]
  if (nrow(sub) < 4 || length(unique(sub[[group_col]])) < 2) return(NULL)
  wt <- suppressWarnings(wilcox.test(sub$RAIR_like ~ factor(sub[[group_col]], levels = c(group_b, group_a)),
                                     exact = FALSE))
  data.table(
    comparison = label,
    dataset = dataset_name,
    cell_type = "Epithelial_Thyroid",
    n_group_a = sum(sub[[group_col]] == group_a),
    n_group_b = sum(sub[[group_col]] == group_b),
    group_a = group_a,
    group_b = group_b,
    median_group_a = median(sub$RAIR_like[sub[[group_col]] == group_a], na.rm = TRUE),
    median_group_b = median(sub$RAIR_like[sub[[group_col]] == group_b], na.rm = TRUE),
    delta_a_minus_b = median(sub$RAIR_like[sub[[group_col]] == group_a], na.rm = TRUE) -
      median(sub$RAIR_like[sub[[group_col]] == group_b], na.rm = TRUE),
    p_wilcox = wt$p.value
  )
}

sample_medians[, tumor_met_group := ifelse(tissue_group == "Paratumor", "Paratumor", "Tumor_met")]
tests <- rbindlist(list(
  test_one(sample_medians, "GSE184362", "tumor_met_group", "Tumor_met", "Paratumor",
           "Epithelial/thyroid RAIR-like in tumor/metastasis vs paratumor"),
  test_one(sample_medians, "GSE232237", "disease_group", "ATC", "PTC",
           "Epithelial/thyroid RAIR-like in ATC vs PTC")
), fill = TRUE)

fwrite(top_summary, file.path(result_dir, "single_cell_top_decile_ecology_summary.csv"))
fwrite(tests, file.path(result_dir, "single_cell_sample_balanced_tests.csv"))

plot_df <- sample_medians[n_cells >= 50]
plot_df[, display_group := ifelse(dataset == "GSE184362", tissue_group, disease_group)]
plot_df$cell_type <- factor(plot_df$cell_type,
                            levels = c("Epithelial_Thyroid", "CAF", "Myeloid", "Endothelial",
                                       "T_NK", "B_Plasma", "Other"))
p <- ggplot(plot_df, aes(x = cell_type, y = RAIR_like, color = display_group)) +
  geom_hline(yintercept = 0, color = "grey82", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2, linewidth = 0.35) +
  geom_jitter(width = 0.18, height = 0, size = 1.5, alpha = 0.75) +
  facet_wrap(~ dataset, scales = "free_x") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = "Sample-level median RAIR-like score",
       color = "Sample group",
       title = "Sample-balanced single-cell localization of the RAIR-like state")
ggsave(file.path(figure_dir, "single_cell_sample_balanced_RAIR_by_celltype.png"),
       p, width = 8.0, height = 4.8, dpi = 300)

sink(file.path(result_dir, "single_cell_ecology_summary.txt"))
cat("Single-cell ecology localization summary\n\n")
cat("Top 10% RAIR-like cells by broad marker-based cell class:\n")
print(top_summary, row.names = FALSE)
cat("\nSample-balanced epithelial/thyroid lineage tests, requiring >=50 cells per sample-celltype:\n")
print(tests, row.names = FALSE)
cat("\nInterpretation guardrail:\n")
cat("Cell-type labels are marker-based because raw GEO matrices lack author malignant labels/CNV calls in this local analysis.\n")
cat("The current result supports epithelial/thyroid-lineage enrichment plus microenvironmental participation, not definitive malignant-cell exclusivity.\n")
sink()

message("Saved single-cell ecology summary.")
