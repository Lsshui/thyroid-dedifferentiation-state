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

z <- function(x) as.numeric(scale(x))

coef_table_lm <- function(fit, term, label) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NULL)
  ci <- suppressWarnings(confint(fit, term))
  data.table(
    model = label,
    term = term,
    estimate = sm[term, "Estimate"],
    ci_low = ci[1],
    ci_high = ci[2],
    p_value = sm[term, "Pr(>|t|)"],
    n = nobs(fit),
    r_squared = summary(fit)$r.squared
  )
}

coef_table_glm <- function(fit, term, label) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NULL)
  ci <- suppressWarnings(confint.default(fit, term))
  data.table(
    model = label,
    term = term,
    log_or = sm[term, "Estimate"],
    OR = exp(sm[term, "Estimate"]),
    ci_low = exp(ci[1]),
    ci_high = exp(ci[2]),
    p_value = sm[term, "Pr(>|z|)"],
    n = nobs(fit)
  )
}

message("TCGA multivariable robustness...")
tcga <- fread(file.path(result_dir, "TCGA_THCA_scores_clinical_mutation.csv"))
for (v in c("RAIR_like", "OneCarbon", "MAPK", "Dediff_EMT", "TDS",
            "IodineHandling", "CAF", "ImmuneCheckpoint")) {
  tcga[[paste0(v, "_z")]] <- z(tcga[[v]])
}
tcga[, BRAF_mut := as.logical(BRAF_mut)]
tcga[, RAS_mut := as.logical(RAS_mut)]
tcga[, N_positive_bin := fifelse(N_positive == "N1", 1, fifelse(N_positive == "N0", 0, NA_real_))]
tcga[, T_high_bin := fifelse(T_high == "T3_T4", 1, fifelse(T_high == "T1_T2", 0, NA_real_))]
tcga[, Stage_high_bin := fifelse(stage_high == "Stage_III_IV", 1,
                                 fifelse(stage_high == "Stage_I_II", 0, NA_real_))]

tcga_lm_results <- rbindlist(list(
  coef_table_lm(lm(IodineHandling_z ~ RAIR_like_z + CAF_z + ImmuneCheckpoint_z +
                     BRAF_mut + RAS_mut, data = tcga),
                "RAIR_like_z", "IodineHandling ~ RAIR-like + CAF/immune + BRAF/RAS"),
  coef_table_lm(lm(OneCarbon_z ~ MAPK_z + Dediff_EMT_z + CAF_z + ImmuneCheckpoint_z +
                     BRAF_mut + RAS_mut + TDS_z, data = tcga),
                "MAPK_z", "OneCarbon ~ MAPK + dediff + CAF/immune + drivers + TDS"),
  coef_table_lm(lm(OneCarbon_z ~ MAPK_z + Dediff_EMT_z + CAF_z + ImmuneCheckpoint_z +
                     BRAF_mut + RAS_mut + TDS_z, data = tcga),
                "Dediff_EMT_z", "OneCarbon ~ MAPK + dediff + CAF/immune + drivers + TDS"),
  coef_table_lm(lm(RAIR_like_z ~ CAF_z + ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = tcga),
                "CAF_z", "RAIR-like ~ CAF/immune + BRAF/RAS"),
  coef_table_lm(lm(RAIR_like_z ~ CAF_z + ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = tcga),
                "ImmuneCheckpoint_z", "RAIR-like ~ CAF/immune + BRAF/RAS"),
  coef_table_lm(lm(RAIR_like_z ~ CAF_z + ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = tcga),
                "BRAF_mutTRUE", "RAIR-like ~ CAF/immune + BRAF/RAS"),
  coef_table_lm(lm(RAIR_like_z ~ CAF_z + ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = tcga),
                "RAS_mutTRUE", "RAIR-like ~ CAF/immune + BRAF/RAS")
), fill = TRUE)

tcga_glm_results <- rbindlist(list(
  coef_table_glm(glm(N_positive_bin ~ RAIR_like_z + CAF_z + ImmuneCheckpoint_z +
                       BRAF_mut + RAS_mut, data = tcga, family = binomial()),
                 "RAIR_like_z", "N1 vs N0 ~ RAIR-like + CAF/immune + BRAF/RAS"),
  coef_table_glm(glm(T_high_bin ~ RAIR_like_z + CAF_z + ImmuneCheckpoint_z +
                       BRAF_mut + RAS_mut, data = tcga, family = binomial()),
                 "RAIR_like_z", "T3/T4 vs T1/T2 ~ RAIR-like + CAF/immune + BRAF/RAS"),
  coef_table_glm(glm(Stage_high_bin ~ RAIR_like_z + CAF_z + ImmuneCheckpoint_z +
                       BRAF_mut + RAS_mut, data = tcga, family = binomial()),
                 "RAIR_like_z", "Stage III/IV vs I/II ~ RAIR-like + CAF/immune + BRAF/RAS")
), fill = TRUE)

candidate_genes <- c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1")
target_tbl <- fread(file.path(result_dir, "TCGA_THCA_scores_clinical_mutation.csv"))
## Candidate expression is stored in the target priority table for high-vs-low testing,
## but not per sample. Use correlations and tier outputs for targets; multivariable
## robustness here focuses on state-level scores.

fwrite(tcga_lm_results, file.path(result_dir, "TCGA_multivariable_lm_robustness.csv"))
fwrite(tcga_glm_results, file.path(result_dir, "TCGA_multivariable_logistic_robustness.csv"))

message("Single-cell sample-level composition robustness...")
sc <- as.data.table(readRDS(file.path(pre_result_dir, "single_cell_per_cell_scores.rds")))
sample_counts <- sc[, .N, by = .(dataset, sample, disease_group, tissue_group, cell_type)]
totals <- sample_counts[, .(total_cells = sum(N)), by = .(dataset, sample)]
sample_counts <- merge(sample_counts, totals, by = c("dataset", "sample"))
sample_counts[, fraction := N / total_cells]
frac_wide <- dcast(sample_counts, dataset + sample + disease_group + tissue_group + total_cells ~
                     cell_type, value.var = "fraction", fill = 0)
for (nm in c("Epithelial_Thyroid", "Myeloid", "CAF", "T_NK", "B_Plasma", "Endothelial")) {
  if (!nm %in% colnames(frac_wide)) frac_wide[[nm]] <- 0
}

sample_medians <- fread(file.path(pre_result_dir, "single_cell_sample_celltype_medians.csv"))
epi <- sample_medians[cell_type == "Epithelial_Thyroid" & n_cells >= 50,
                      .(dataset, sample, epithelial_RAIR = RAIR_like,
                        epithelial_n = n_cells)]
sc_sample <- merge(frac_wide, epi, by = c("dataset", "sample"), all.x = FALSE)
sc_sample[, log_total_cells := log10(total_cells + 1)]
sc_sample[, tumor_met_group := fifelse(tissue_group == "Paratumor", "Paratumor", "Tumor_met")]
fwrite(sc_sample, file.path(result_dir, "single_cell_sample_composition_epithelial_RAIR.csv"))

fit_sc <- function(data, formula, term, label) {
  if (nrow(data) < 6 || length(unique(data[[all.vars(formula)[2]]])) < 2) return(NULL)
  fit <- lm(formula, data = data)
  coef_table_lm(fit, term, label)
}

g184 <- sc_sample[dataset == "GSE184362"]
g232 <- sc_sample[dataset == "GSE232237"]
sc_lm_results <- rbindlist(list(
  fit_sc(g184, epithelial_RAIR ~ tumor_met_group + Myeloid + CAF + T_NK +
           Epithelial_Thyroid + log_total_cells,
         "tumor_met_groupTumor_met",
         "GSE184362 epithelial RAIR ~ tumor/met + cell fractions + total cells"),
  fit_sc(g232, epithelial_RAIR ~ disease_group + Myeloid + CAF + T_NK +
           Epithelial_Thyroid + log_total_cells,
         "disease_groupPTC",
         "GSE232237 epithelial RAIR ~ PTC(vs ATC) + cell fractions + total cells")
), fill = TRUE)
fwrite(sc_lm_results, file.path(result_dir, "single_cell_sample_composition_lm_robustness.csv"))

top_sample <- fread(file.path(pre_result_dir, "single_cell_top_RAIR_like_by_sample.csv"))
top_by_sample <- top_sample[, .(n_top_cells = sum(n_top_cells),
                                fraction_dataset_top10 = sum(fraction_dataset_top10)),
                            by = .(dataset, sample, tissue_group)]
setorder(top_by_sample, dataset, -fraction_dataset_top10)
dominance <- top_by_sample[, .(
  n_samples_with_top_cells = .N,
  top1_sample = sample[1],
  top1_fraction = fraction_dataset_top10[1],
  top3_fraction = sum(head(fraction_dataset_top10, 3)),
  n_samples_over_1pct = sum(fraction_dataset_top10 >= 0.01)
), by = dataset]
fwrite(top_by_sample, file.path(result_dir, "single_cell_top_RAIR_sample_contribution.csv"))
fwrite(dominance, file.path(result_dir, "single_cell_top_RAIR_sample_dominance.csv"))

p1 <- ggplot(sc_sample, aes(x = Myeloid + CAF + T_NK + B_Plasma + Endothelial,
                            y = epithelial_RAIR, color = tissue_group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  facet_wrap(~ dataset, scales = "free") +
  theme_bw(base_size = 10) +
  labs(x = "Non-epithelial cell fraction", y = "Sample epithelial/thyroid median RAIR-like",
       color = "Tissue group",
       title = "Epithelial RAIR-like score versus microenvironment fraction")
ggsave(file.path(figure_dir, "single_cell_epithelial_RAIR_vs_microenvironment_fraction.png"),
       p1, width = 7.2, height = 4.6, dpi = 300)

p2 <- ggplot(top_by_sample,
             aes(x = reorder(sample, fraction_dataset_top10),
                 y = fraction_dataset_top10, fill = tissue_group)) +
  geom_col(width = 0.72) +
  coord_flip() +
  facet_wrap(~ dataset, scales = "free_y") +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Fraction of dataset top 10% RAIR-like cells",
       fill = "Tissue group",
       title = "Sample contribution to top RAIR-like cells")
ggsave(file.path(figure_dir, "single_cell_top_RAIR_sample_contribution.png"),
       p2, width = 7.6, height = 7.2, dpi = 300)

sink(file.path(result_dir, "confounder_robustness_summary.txt"))
cat("Confounder robustness analysis\n\n")
cat("TCGA multivariable linear models:\n")
print(tcga_lm_results, row.names = FALSE)
cat("\nTCGA multivariable logistic models:\n")
print(tcga_glm_results, row.names = FALSE)
cat("\nSingle-cell sample-level composition-adjusted models:\n")
print(sc_lm_results, row.names = FALSE)
cat("\nSingle-cell top RAIR-like sample dominance:\n")
print(dominance, row.names = FALSE)
cat("\nInterpretation notes:\n")
cat("- TCGA models test whether RAIR-like/state associations remain after CAF/immune checkpoint scores and BRAF/RAS mutation status.\n")
cat("- Single-cell models are small-n and exploratory; they test whether epithelial RAIR-like signal is merely a sample-level microenvironment fraction artifact.\n")
cat("- Candidate genes assessed in target tiering remain the primary target-level evidence.\n")
sink()

message("Saved confounder robustness outputs.")
