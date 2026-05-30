options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_root <- file.path(root, "analysis")
data_dir <- file.path(analysis_root, "data", "hpa")
result_dir <- file.path(analysis_root, "results")
figure_dir <- file.path(analysis_root, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

targets <- c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1")
tier1 <- c("MTHFD2", "SLC1A5", "SLC7A5", "AXL")

download_if_missing <- function(url, file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    message("Downloading ", url)
    download.file(url, file, mode = "wb", quiet = TRUE)
  }
}

download_if_missing("https://www.proteinatlas.org/download/proteinatlas.tsv.zip",
                    file.path(data_dir, "proteinatlas.tsv.zip"))
download_if_missing("https://v23.proteinatlas.org/download/pathology.tsv.zip",
                    file.path(data_dir, "pathology_v23.tsv.zip"))
download_if_missing("https://v23.proteinatlas.org/download/normal_tissue.tsv.zip",
                    file.path(data_dir, "normal_tissue_v23.tsv.zip"))

message("Parsing Human Protein Atlas tables...")
pa <- fread(cmd = paste("powershell -NoProfile -Command",
                        shQuote(sprintf("Expand-Archive -Path '%s' -DestinationPath '%s/tmp_pa' -Force; Get-Content '%s/tmp_pa/proteinatlas.tsv'",
                                        file.path(data_dir, "proteinatlas.tsv.zip"),
                                        data_dir, data_dir))),
            sep = "\t", data.table = TRUE, showProgress = FALSE)
pathology <- fread(cmd = paste("powershell -NoProfile -Command",
                               shQuote(sprintf("Expand-Archive -Path '%s' -DestinationPath '%s/tmp_pathology' -Force; Get-Content '%s/tmp_pathology/pathology.tsv'",
                                               file.path(data_dir, "pathology_v23.tsv.zip"),
                                               data_dir, data_dir))),
                   sep = "\t", data.table = TRUE, showProgress = FALSE)
normal <- fread(cmd = paste("powershell -NoProfile -Command",
                            shQuote(sprintf("Expand-Archive -Path '%s' -DestinationPath '%s/tmp_normal' -Force; Get-Content '%s/tmp_normal/normal_tissue.tsv'",
                                            file.path(data_dir, "normal_tissue_v23.tsv.zip"),
                                            data_dir, data_dir))),
                sep = "\t", data.table = TRUE, showProgress = FALSE)

hpa_general <- pa[Gene %in% targets, .(
  gene = Gene,
  protein_class = `Protein class`,
  evidence = Evidence,
  hpa_evidence = `HPA evidence`,
  reliability_IH = `Reliability (IH)`,
  antibody = Antibody,
  subcellular_location = `Subcellular location`,
  rna_cancer_specificity = `RNA cancer specificity`,
  rna_cancer_distribution = `RNA cancer distribution`,
  rna_cancer_specific_pTPM = `RNA cancer specific pTPM`,
  thyroid_carcinoma_prognostics_TCGA = `Cancer prognostics - Thyroid Carcinoma (TCGA)`
)]

thyroid_normal <- normal[`Gene name` %in% targets & Tissue == "thyroid gland",
                         .(gene = `Gene name`, tissue = Tissue, cell_type = `Cell type`,
                           normal_thyroid_IHC_level = Level, normal_reliability = Reliability)]

thyroid_pathology <- pathology[`Gene name` %in% targets & Cancer == "thyroid cancer",
                               .(gene = `Gene name`, cancer = Cancer,
                                 pathology_high = as.numeric(High),
                                 pathology_medium = as.numeric(Medium),
                                 pathology_low = as.numeric(Low),
                                 pathology_not_detected = as.numeric(`Not detected`),
                                 favorable_p = `prognostic - favorable`,
                                 unfavorable_p = `prognostic - unfavorable`)]
thyroid_pathology[, pathology_total := pathology_high + pathology_medium + pathology_low + pathology_not_detected]
thyroid_pathology[, pathology_detected_fraction := (pathology_high + pathology_medium + pathology_low) / pathology_total]
thyroid_pathology[, pathology_high_medium_fraction := (pathology_high + pathology_medium) / pathology_total]

fwrite(hpa_general, file.path(result_dir, "HPA_target_general_annotation.csv"))
fwrite(thyroid_normal, file.path(result_dir, "HPA_target_normal_thyroid_IHC.csv"))
fwrite(thyroid_pathology, file.path(result_dir, "HPA_target_thyroid_cancer_pathology.csv"))

message("Querying ChEMBL target/actionability records...")
chembl_targets <- data.table(
  gene = c("MTHFD2", "SLC1A5", "SLC7A5", "AXL"),
  target_chembl_id = c("CHEMBL4105963", "CHEMBL3562162", "CHEMBL4459", "CHEMBL4895")
)

safe_json <- function(url) {
  tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
}

chembl <- rbindlist(lapply(seq_len(nrow(chembl_targets)), function(i) {
  gene <- chembl_targets$gene[i]
  tid <- chembl_targets$target_chembl_id[i]
  activity <- safe_json(sprintf("https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=%s&limit=1", tid))
  mechanism <- safe_json(sprintf("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id=%s&limit=100", tid))
  target <- safe_json(sprintf("https://www.ebi.ac.uk/chembl/api/data/target/%s.json", tid))
  mech_count <- if (is.null(mechanism)) NA_integer_ else mechanism$page_meta$total_count
  mech_tbl <- if (!is.null(mechanism) && length(mechanism$mechanisms) > 0) as.data.table(mechanism$mechanisms) else data.table()
  molecule_ids <- if (nrow(mech_tbl) > 0) unique(head(mech_tbl$molecule_chembl_id, 8)) else character()
  max_phase <- if (nrow(mech_tbl) > 0 && "max_phase" %in% names(mech_tbl)) max(mech_tbl$max_phase, na.rm = TRUE) else NA_real_
  if (!is.finite(max_phase)) max_phase <- NA_real_
  data.table(
    gene = gene,
    target_chembl_id = tid,
    target_pref_name = if (is.null(target)) NA_character_ else target$pref_name,
    target_type = if (is.null(target)) NA_character_ else target$target_type,
    activity_total_count = if (is.null(activity)) NA_integer_ else activity$page_meta$total_count,
    mechanism_total_count = mech_count,
    max_mechanism_phase = max_phase,
    mechanism_action_types = if (nrow(mech_tbl) > 0 && "action_type" %in% names(mech_tbl)) paste(unique(mech_tbl$action_type), collapse = ";") else "",
    mechanism_molecule_ids = paste(molecule_ids, collapse = ";")
  )
}), fill = TRUE)

chembl[, actionability_bin := fifelse(!is.na(max_mechanism_phase) & max_mechanism_phase >= 4,
                                      "Clinical/approved mechanism record",
                                      fifelse(!is.na(max_mechanism_phase) & max_mechanism_phase >= 2,
                                              "Clinical-stage mechanism record",
                                              fifelse(!is.na(activity_total_count) & activity_total_count >= 100,
                                                      "Preclinical chemical tractability",
                                                      fifelse(!is.na(activity_total_count) & activity_total_count > 0,
                                                              "Sparse chemical evidence", "No ChEMBL activity"))))]
fwrite(chembl, file.path(result_dir, "ChEMBL_target_actionability.csv"))

tiers <- fread(file.path(result_dir, "candidate_target_tiers.csv"))
trans <- merge(tiers[gene %in% targets], hpa_general, by = "gene", all.x = TRUE)
path_slim <- thyroid_pathology[, .(gene, pathology_high, pathology_medium, pathology_low,
                                   pathology_not_detected, pathology_total,
                                   pathology_detected_fraction,
                                   pathology_high_medium_fraction)]
normal_slim <- thyroid_normal[, .(normal_thyroid_IHC_summary =
                                    paste(cell_type, normal_thyroid_IHC_level,
                                          paste0("(", normal_reliability, ")"),
                                          collapse = "; ")), by = gene]
trans <- merge(trans, path_slim, by = "gene", all.x = TRUE)
trans <- merge(trans, normal_slim, by = "gene", all.x = TRUE)
trans <- merge(trans, chembl[, .(gene, target_chembl_id, activity_total_count,
                                mechanism_total_count, max_mechanism_phase,
                                actionability_bin)], by = "gene", all.x = TRUE)
fwrite(trans, file.path(result_dir, "translational_target_annotation_HPA_ChEMBL.csv"))

path_long <- melt(thyroid_pathology[gene %in% targets],
                  id.vars = c("gene", "cancer"),
                  measure.vars = c("pathology_high", "pathology_medium", "pathology_low", "pathology_not_detected"),
                  variable.name = "IHC_category", value.name = "n_cases")
path_long[, IHC_category := factor(IHC_category,
                                   levels = c("pathology_high", "pathology_medium", "pathology_low", "pathology_not_detected"),
                                   labels = c("High", "Medium", "Low", "Not detected"))]
p1 <- ggplot(path_long, aes(x = gene, y = n_cases, fill = IHC_category)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("High" = "#8c1d1d", "Medium" = "#c65f38",
                               "Low" = "#d5a253", "Not detected" = "#a9b0b8")) +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "HPA pathology cases in thyroid cancer",
       fill = "IHC level",
       title = "Human Protein Atlas thyroid cancer IHC evidence")
ggsave(file.path(figure_dir, "HPA_thyroid_cancer_IHC_targets.png"),
       p1, width = 6.4, height = 4.2, dpi = 300)

p2 <- ggplot(chembl, aes(x = reorder(gene, activity_total_count), y = activity_total_count,
                         fill = actionability_bin)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_log10() +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "ChEMBL activity records, log10 scale",
       fill = "Actionability bin",
       title = "ChEMBL chemical tractability of Tier 1A targets")
ggsave(file.path(figure_dir, "ChEMBL_Tier1A_actionability.png"),
       p2, width = 6.4, height = 4.0, dpi = 300)

sink(file.path(result_dir, "HPA_ChEMBL_translational_annotation_summary.txt"))
cat("HPA and ChEMBL translational annotation for candidate targets\n\n")
cat("Sources:\n")
cat("- Human Protein Atlas current proteinatlas.tsv.zip: https://www.proteinatlas.org/download/proteinatlas.tsv.zip\n")
cat("- Human Protein Atlas v23 pathology.tsv.zip: https://v23.proteinatlas.org/download/pathology.tsv.zip\n")
cat("- Human Protein Atlas v23 normal_tissue.tsv.zip: https://v23.proteinatlas.org/download/normal_tissue.tsv.zip\n")
cat("- ChEMBL API: https://www.ebi.ac.uk/chembl/api/data\n\n")
cat("HPA thyroid cancer pathology:\n")
print(thyroid_pathology[gene %in% targets], row.names = FALSE)
cat("\nHPA normal thyroid IHC:\n")
print(thyroid_normal[gene %in% targets], row.names = FALSE)
cat("\nChEMBL target actionability:\n")
print(chembl, row.names = FALSE)
cat("\nIntegrated target annotation:\n")
print(trans[, .(gene, evidence_tier, n_external_up_validations,
                pathology_high_medium_fraction, normal_thyroid_IHC_summary,
                reliability_IH, activity_total_count, max_mechanism_phase,
                actionability_bin, actionability_note)], row.names = FALSE)
cat("\nInterpretation:\n")
cat("AXL has the strongest clinical actionability in ChEMBL. SLC1A5/SLC7A5 have substantial preclinical chemical activity records. MTHFD2 is biologically strong in omics validation but has sparse ChEMBL mechanism/activity evidence, so it should be framed as an experimental metabolic vulnerability rather than a near-clinic target.\n")
sink()

message("Saved HPA/ChEMBL translational annotation outputs.")
