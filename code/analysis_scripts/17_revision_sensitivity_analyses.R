options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
result_dir <- file.path(root, "analysis", "results")
pre_dir <- file.path(root, "preanalysis", "results")
out_dir <- file.path(result_dir, "revision_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## GSE151179 power context for the exploratory RAI-response contrast.
g151 <- fread(file.path(pre_dir, "GSE151179_scores_metadata.csv"))
rai <- g151[!is.na(rai_group) & rai_group %in% c("RAI_refractory", "RAI_avid")]
n_refractory <- rai[rai_group == "RAI_refractory", .N]
n_avid <- rai[rai_group == "RAI_avid", .N]

power_auc <- data.table(
  contrast = "GSE151179 RAI-refractory vs RAI-avid",
  n_RAI_refractory = n_refractory,
  n_RAI_avid = n_avid,
  target_AUC = c(0.70, 0.75, 0.80, 0.85)
)
power_auc[, power_alpha_0_05 := vapply(
  target_AUC,
  function(a) pROC::power.roc.test(
    ncases = n_RAI_refractory[1],
    ncontrols = n_RAI_avid[1],
    auc = a,
    sig.level = 0.05
  )$power,
  numeric(1)
)]
fwrite(power_auc, file.path(out_dir, "GSE151179_RAI_AUC_power_context.csv"))

## GSE232237 unadjusted sample-level epithelial RAIR comparison.
sc <- fread(file.path(result_dir, "single_cell_sample_composition_epithelial_RAIR.csv"))
g232 <- sc[dataset == "GSE232237" & disease_group %in% c("ATC", "PTC") &
             !is.na(epithelial_RAIR)]
g232[, disease_group := factor(disease_group, levels = c("PTC", "ATC"))]
wt <- wilcox.test(epithelial_RAIR ~ disease_group, data = g232, exact = FALSE)

g232_summary <- g232[, .(
  n_samples = .N,
  median_epithelial_RAIR = median(epithelial_RAIR, na.rm = TRUE),
  min_epithelial_RAIR = min(epithelial_RAIR, na.rm = TRUE),
  max_epithelial_RAIR = max(epithelial_RAIR, na.rm = TRUE),
  median_epithelial_fraction = median(Epithelial_Thyroid, na.rm = TRUE),
  min_epithelial_fraction = min(Epithelial_Thyroid, na.rm = TRUE),
  max_epithelial_fraction = max(Epithelial_Thyroid, na.rm = TRUE)
), by = disease_group]

g232_test <- data.table(
  comparison = "GSE232237 unadjusted sample-level epithelial RAIR: ATC vs PTC",
  n_ATC = g232[disease_group == "ATC", .N],
  n_PTC = g232[disease_group == "PTC", .N],
  median_ATC = g232[disease_group == "ATC", median(epithelial_RAIR, na.rm = TRUE)],
  median_PTC = g232[disease_group == "PTC", median(epithelial_RAIR, na.rm = TRUE)],
  delta_ATC_minus_PTC = g232[disease_group == "ATC", median(epithelial_RAIR, na.rm = TRUE)] -
    g232[disease_group == "PTC", median(epithelial_RAIR, na.rm = TRUE)],
  p_wilcox = wt$p.value
)

fwrite(g232[, .(dataset, sample, disease_group, total_cells, Epithelial_Thyroid,
                Myeloid, CAF, T_NK, epithelial_n, epithelial_RAIR)],
       file.path(out_dir, "GSE232237_sample_level_epithelial_RAIR.csv"))
fwrite(g232_summary, file.path(out_dir, "GSE232237_sample_level_epithelial_RAIR_summary.csv"))
fwrite(g232_test, file.path(out_dir, "GSE232237_unadjusted_epithelial_RAIR_test.csv"))

sink(file.path(out_dir, "revision_sensitivity_summary.txt"))
cat("Revision sensitivity analyses\n\n")
cat("GSE151179 RAI endpoint power context:\n")
print(power_auc)
cat("\nGSE232237 unadjusted epithelial/thyroid sample-level RAIR comparison:\n")
print(g232_test)
cat("\nGSE232237 disease-group summary:\n")
print(g232_summary)
sink()

message("Revision sensitivity analyses complete: ", out_dir)
