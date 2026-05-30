options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
copykat_dir <- file.path(root, "analysis", "results", "copykat")
para_dir <- file.path(root, "analysis", "results", "copykat_paratumor_ref", "GSE184362")
figure_dir <- file.path(root, "analysis", "figures")
result_dir <- file.path(root, "analysis", "results")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

general_sample <- fread(file.path(copykat_dir, "copykat_sample_summary.csv"))
general_tests <- fread(file.path(copykat_dir, "copykat_epithelial_RAIR_tests.csv"))
general_pred <- fread(file.path(copykat_dir, "copykat_all_predictions.csv"))

para_sample <- fread(file.path(para_dir, "paratumor_ref_copykat_sample_summary.csv"))
para_tests <- fread(file.path(para_dir, "paratumor_ref_copykat_RAIR_test.csv"))
para_pred <- fread(file.path(para_dir, "paratumor_ref_copykat_all_predictions.csv"))

general_sample[, analysis := "copyKAT immune/stromal reference"]
para_sample[, analysis := "copyKAT paratumor epithelial reference"]
setnames(para_sample, "sample_target", "sample", skip_absent = TRUE)
sample_plot <- rbindlist(list(
  general_sample[, .(analysis, dataset, sample, disease_group, tissue_group,
                     n_epithelial, epithelial_aneuploid_fraction)],
  para_sample[, .(analysis, dataset, sample, disease_group = "PTC_related",
                  tissue_group = "Tumor/met target",
                  n_epithelial = n_target_epi,
                  epithelial_aneuploid_fraction = aneuploid_fraction)]
), fill = TRUE)

p1 <- ggplot(sample_plot,
             aes(x = reorder(sample, epithelial_aneuploid_fraction),
                 y = epithelial_aneuploid_fraction, fill = tissue_group)) +
  geom_col(width = 0.72) +
  coord_flip() +
  facet_wrap(~ analysis, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "Aneuploid fraction among epithelial/thyroid cells",
       fill = "Tissue group",
       title = "copyKAT malignant-like classification in epithelial/thyroid cells")
ggsave(file.path(figure_dir, "copykat_epithelial_aneuploid_fraction.png"),
       p1, width = 8.0, height = 5.6, dpi = 300)

plot_pred <- general_pred[cell_type == "Epithelial_Thyroid" &
                            copykat_class %in% c("aneuploid", "diploid")]
p2 <- ggplot(plot_pred,
             aes(x = copykat_class, y = RAIR_like, fill = copykat_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.18, size = 0.45, alpha = 0.22) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("aneuploid" = "#b73e3e", "diploid" = "#3274a1")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = "copyKAT class", y = "RAIR-like score",
       title = "RAIR-like score in copyKAT-classified epithelial/thyroid cells")
ggsave(file.path(figure_dir, "copykat_epithelial_RAIR_by_class.png"),
       p2, width = 6.8, height = 4.6, dpi = 300)

para_target <- para_pred[tissue_group != "Paratumor" &
                           cell_type == "Epithelial_Thyroid" &
                           copykat_class %in% c("aneuploid", "diploid")]
p3 <- ggplot(para_target, aes(x = sample_target, y = RAIR_like, fill = copykat_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.18, size = 0.45, alpha = 0.22) +
  scale_fill_manual(values = c("aneuploid" = "#b73e3e", "diploid" = "#3274a1")) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = NULL, y = "RAIR-like score", fill = "copyKAT class",
       title = "GSE184362 tumor/met epithelial cells using paratumor epithelial reference")
ggsave(file.path(figure_dir, "copykat_GSE184362_paratumor_ref_RAIR_by_class.png"),
       p3, width = 6.4, height = 4.4, dpi = 300)

combined_summary <- list(
  immune_stromal_ref_sample_summary = general_sample,
  immune_stromal_ref_tests = general_tests,
  paratumor_epithelial_ref_sample_summary = para_sample,
  paratumor_epithelial_ref_tests = para_tests
)
saveRDS(combined_summary, file.path(result_dir, "copykat_combined_summary.rds"))

sink(file.path(result_dir, "copykat_interpretation_summary.txt"))
cat("copyKAT malignant-cell classification summary\n\n")
cat("Important method note:\n")
cat("copyKAT was run on subsampled epithelial/thyroid cells plus reference cells because full datasets are too large for per-sample Windows execution. Native copyKAT prediction files were parsed after classification; heatmap generation was terminated when necessary.\n\n")
cat("Immune/stromal-reference copyKAT sample summary:\n")
print(general_sample, row.names = FALSE)
cat("\nImmune/stromal-reference epithelial RAIR-like by copyKAT class:\n")
print(general_tests, row.names = FALSE)
cat("\nGSE184362 paratumor-epithelial-reference copyKAT sample summary:\n")
print(para_sample, row.names = FALSE)
cat("\nGSE184362 paratumor-epithelial-reference RAIR-like by copyKAT class:\n")
print(para_tests, row.names = FALSE)
cat("\nInterpretation:\n")
cat("- GSE232237 supports higher RAIR-like score in copyKAT-aneuploid epithelial/thyroid cells than diploid epithelial/thyroid cells.\n")
cat("- GSE184362 immune/stromal-reference runs overcall some paratumor epithelial cells as aneuploid, so that reference design is not sufficient for malignant specificity claims.\n")
cat("- GSE184362 paratumor-epithelial-reference runs classify nearly all sampled tumor/metastatic epithelial cells as aneuploid, supporting malignant-like compartment localization, but diploid target epithelial cells are too few for a strong within-target RAIR comparison.\n")
sink()

message("Saved copyKAT summary figures and interpretation.")
