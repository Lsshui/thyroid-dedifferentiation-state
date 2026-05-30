options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(openxlsx)
  library(grid)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
pre_dir <- file.path(root, "preanalysis", "results")
res_dir <- file.path(root, "analysis", "results")
fig_src_dir <- file.path(root, "analysis", "figures")
asset_dir <- file.path(root, "analysis", "manuscript_assets")
main_fig_dir <- file.path(asset_dir, "figures", "main")
supp_fig_dir <- file.path(asset_dir, "figures", "supplementary")
table_dir <- file.path(asset_dir, "tables")
source_dir <- file.path(asset_dir, "source_data")
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
  navy = "#26364A",
  blue = "#3978A8",
  sky = "#8EBAD9",
  teal = "#2E9F8F",
  green = "#5DA35A",
  red = "#B6403A",
  orange = "#D9923B",
  gold = "#D8B365",
  purple = "#7A5EA8",
  grey1 = "#2F3437",
  grey2 = "#6F7A82",
  grey3 = "#B9C0C7",
  grey4 = "#E8EAED"
)
`[.palette_drop_names` <- function(x, i, ...) unname(NextMethod("["))
class(pal) <- c("palette_drop_names", class(pal))

theme_jtm <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.30, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.7, colour = "black"),
      legend.title = element_text(size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.9),
      legend.key.size = unit(3.5, "mm"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
      strip.text = element_text(size = base_size - 0.3, face = "bold"),
      plot.title = element_text(size = base_size + 0.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = pal["grey2"], hjust = 0),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}
theme_set(theme_jtm())

pathway_label <- function(x, width = 42) {
  x <- gsub("^(REACTOME|HALLMARK)_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("\\bRho\\b", "RHO", x)
  x <- gsub("\\bGtpase\\b", "GTPase", x)
  x <- gsub("\\bMhc\\b", "MHC", x)
  x <- gsub("\\bTp53\\b", "TP53", x)
  x <- gsub("\\bSars Cov\\b", "SARS-CoV", x)
  x <- gsub("\\bKras\\b", "KRAS", x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

save_pub <- function(plot, file_base, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(file_base, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(file_base, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(file_base, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(file_base, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot)
  dev.off()
}

pval_fmt <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 3)))
}

score_long <- function(dt, score_cols = c("TDS", "OneCarbon", "RAIR_like")) {
  melt(dt, id.vars = intersect(c("dataset", "sample", "histology_group", "rai_group",
                                 "rai_group_remission", "disease_group", "tissue_group"),
                               names(dt)),
       measure.vars = score_cols, variable.name = "score", value.name = "value")
}

message("Loading analysis tables...")
tcga <- fread(file.path(res_dir, "TCGA_THCA_scores_clinical_mutation.csv"))
tcga_resid <- fread(file.path(res_dir, "TCGA_MAPK_residual_models.csv"))
tcga_braf <- fread(file.path(res_dir, "TCGA_BRAF_mutant_internal_models.csv"))
tcga_target_cor <- fread(file.path(res_dir, "TCGA_BRAF_mutant_target_gene_correlations.csv"))
hallmark <- fread(file.path(res_dir, "TCGA_RAIR_hallmark_fgsea.csv"))
reactome <- fread(file.path(res_dir, "TCGA_RAIR_reactome_fgsea.csv"))
gse33630 <- fread(file.path(pre_dir, "GSE33630_scores_metadata.csv"))
gse76039 <- fread(file.path(pre_dir, "GSE76039_scores_metadata.csv"))
gse151179 <- fread(file.path(pre_dir, "GSE151179_scores_metadata.csv"))
bulk_tests <- fread(file.path(pre_dir, "bulk_score_group_tests.csv"))
external_targets <- fread(file.path(res_dir, "external_target_validation_long.csv"))
target_tiers <- fread(file.path(res_dir, "candidate_target_tiers.csv"))
single_top <- fread(file.path(pre_dir, "single_cell_top_RAIR_like_celltype_distribution.csv"))
single_bal <- fread(file.path(res_dir, "single_cell_sample_composition_epithelial_RAIR.csv"))
copykat_summary <- fread(file.path(res_dir, "copykat", "copykat_sample_summary.csv"))
copykat_tests <- fread(file.path(res_dir, "copykat", "copykat_epithelial_RAIR_tests.csv"))
copykat_pred <- fread(file.path(res_dir, "copykat", "copykat_all_predictions.csv"))
copykat_para <- fread(file.path(res_dir, "copykat_paratumor_ref", "GSE184362", "paratumor_ref_copykat_sample_summary.csv"))
copykat_para_pred <- fread(file.path(res_dir, "copykat_paratumor_ref", "GSE184362", "paratumor_ref_copykat_all_predictions.csv"))
hpa_path <- fread(file.path(res_dir, "HPA_target_thyroid_cancer_pathology.csv"))
chembl <- fread(file.path(res_dir, "ChEMBL_target_actionability.csv"))
trans_target <- fread(file.path(res_dir, "translational_target_annotation_HPA_ChEMBL.csv"))
cnv_summary <- fread(file.path(res_dir, "GSE184362_epithelial_reference_cnv_proxy_sample_summary.csv"))
sc_ecology <- fread(file.path(res_dir, "single_cell_sample_balanced_tests.csv"))

## -----------------------------
## Main Figure 1: study design.
## -----------------------------
message("Drawing Figure 1...")
workflow <- data.table(
  x = 1:5,
  step = c("Discovery", "Dedifferentiation\nvalidation", "Single-cell\nlocalization",
           "copyKAT malignant\nclassification", "Target\ntranslation"),
  detail = c("TCGA-THCA\n505 tumors", "GSE33630/GSE76039\nPTC-PDTC-ATC",
             "GSE184362/GSE232237\n282,758 cells", "Aneuploid epithelial\ncompartment",
             "HPA + ChEMBL\nTiered nodes"),
  fill = c(pal["navy"], pal["blue"], pal["teal"], pal["orange"], pal["red"])
)
p1a <- ggplot(workflow, aes(x = x, y = 1)) +
  geom_segment(data = workflow[x < 5], aes(x = x + 0.32, xend = x + 0.68, y = 1, yend = 1),
               linewidth = 0.45, arrow = arrow(length = unit(2.2, "mm")), colour = pal["grey2"]) +
  geom_label(aes(label = step, fill = fill), colour = "white", size = 2.4, fontface = "bold",
             linewidth = 0, label.r = unit(1.2, "mm"), label.padding = unit(2.5, "mm")) +
  geom_text(aes(y = 0.67, label = detail), size = 2.1, lineheight = 0.9, colour = pal["grey1"]) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.55, 5.45), ylim = c(0.45, 1.25), clip = "off") +
  theme_void(base_family = "Arial") +
  labs(title = "Study workflow")

module_dt <- data.table(
  module = c("TDS", "One-carbon", "MAPK", "Dediff/EMT"),
  n_genes = c(16, 11, 10, 10),
  direction = c("subtracted", "added", "added", "added")
)
module_dt[, module := factor(module, levels = module)]
p1b <- ggplot(module_dt, aes(x = module, y = n_genes, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_genes), vjust = -0.35, size = 2.2) +
  scale_fill_manual(values = c("added" = pal["teal"], "subtracted" = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "RAIR-like score definition",
       subtitle = "z(One-carbon) + z(MAPK) + z(Dediff/EMT) - z(TDS)",
       x = NULL, y = "Genes in module", fill = "Score term")

cohort_dt <- data.table(
  cohort = c("TCGA-THCA", "GSE151179", "GSE33630", "GSE76039", "GSE184362", "GSE232237"),
  role = c("Discovery", "RAI exploratory", "PTC/ATC validation", "PDTC/ATC validation",
           "PTC scRNA ecology", "PTC/ATC scRNA ecology"),
  n = c(505, 49, 105, 37, 197955, 84803),
  modality = c("Bulk RNA-seq", "Bulk array", "Bulk array", "Bulk array", "scRNA-seq", "scRNA-seq")
)
cohort_dt[, cohort_label := paste0(cohort, "\n", role)]
cohort_dt[, cohort_label := factor(cohort_label, levels = rev(cohort_label))]
p1c <- ggplot(cohort_dt, aes(x = n, y = cohort_label, fill = modality)) +
  geom_col(width = 0.72) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = c("Bulk RNA-seq" = pal["navy"], "Bulk array" = pal["blue"], "scRNA-seq" = pal["teal"])) +
  theme_jtm() +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.8, lineheight = 0.86)) +
  labs(title = "Cohorts used for orthogonal evidence", x = "Samples/cells, log scale", y = NULL, fill = NULL)

claim_dt <- data.table(
  evidence = c("RAIR-like vs iodine-handling", "PTC-PDTC-ATC validation",
               "copyKAT malignant localization", "Target actionability"),
  status = c("Strong", "Strong", "Moderate", "Moderate"),
  note = c("P=4.46e-37", "AUC 0.87-0.89", "copyKAT P=4.18e-11", "SLC7A5/AXL prioritized")
)
claim_dt[, evidence := factor(evidence, levels = rev(evidence))]
claim_dt[, strength := fifelse(status == "Strong", 2, 1)]
p1d <- ggplot(claim_dt, aes(x = strength, y = evidence, fill = status)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.35) +
  geom_text(aes(x = strength + 0.06, label = note), size = 2.0, colour = pal["grey1"], hjust = 0) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Moderate", "Strong"),
                     limits = c(0, 3.55), expand = expansion(mult = c(0.01, 0.02))) +
  scale_fill_manual(values = c("Strong" = pal["red"], "Moderate" = pal["orange"])) +
  theme_jtm() +
  theme(axis.title = element_blank(), legend.position = "none", axis.text.y = element_text(size = 6.0)) +
  labs(title = "Evidence strength for the manuscript claim")

bottom_row <- (p1b | p1c | p1d) + plot_layout(widths = c(0.82, 1.33, 1.45))
fig1 <- (p1a / bottom_row) +
  plot_layout(heights = c(0.85, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig1, file.path(main_fig_dir, "Figure_1_study_design"), 190, 126)
fwrite(workflow, file.path(source_dir, "Figure_1A_workflow.csv"))
fwrite(module_dt, file.path(source_dir, "Figure_1B_score_modules.csv"))
fwrite(cohort_dt, file.path(source_dir, "Figure_1C_cohorts.csv"))
fwrite(claim_dt, file.path(source_dir, "Figure_1D_evidence_strength.csv"))

## -----------------------------
## Main Figure 2: TCGA discovery.
## -----------------------------
message("Drawing Figure 2...")
tcga[, driver_group := fifelse(BRAF_mut == TRUE, "BRAF-mut", fifelse(RAS_mut == TRUE, "RAS-mut", "Other/WT"))]
tcga[, driver_group := factor(driver_group, levels = c("Other/WT", "RAS-mut", "BRAF-mut"))]
tcga[, RAIR_z := as.numeric(scale(RAIR_like))]
tcga[, Iodine_z := as.numeric(scale(IodineHandling))]
tcga[, MAPK_z := as.numeric(scale(MAPK))]
tcga[, CAF_z := as.numeric(scale(CAF))]
tcga[, Immune_z := as.numeric(scale(ImmuneCheckpoint))]
tcga[, RAIR_resid_MAPK := residuals(lm(RAIR_z ~ MAPK_z + CAF_z + Immune_z + BRAF_mut + RAS_mut, data = tcga))]

p2a <- ggplot(tcga, aes(x = RAIR_like, y = IodineHandling, colour = driver_group)) +
  geom_point(size = 1.2, alpha = 0.70) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.45) +
  scale_colour_manual(values = c("Other/WT" = pal["grey2"], "RAS-mut" = pal["blue"], "BRAF-mut" = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(title = "RAIR-like vs iodine-handling", x = "RAIR-like score", y = "Iodine-handling score", colour = "Driver")

p2b <- ggplot(tcga, aes(x = driver_group, y = RAIR_like, fill = driver_group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.12, size = 0.55, alpha = 0.35) +
  scale_fill_manual(values = c("Other/WT" = pal["grey3"], "RAS-mut" = pal["blue"], "BRAF-mut" = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "none") +
  labs(title = "Driver-associated RAIR-like shift", x = NULL, y = "RAIR-like score")

resid_group <- copy(tcga)
resid_group[, residual_group := ifelse(RAIR_resid_MAPK >= median(RAIR_resid_MAPK, na.rm = TRUE), "High residual RAIR", "Low residual RAIR")]
p2c <- ggplot(resid_group, aes(x = residual_group, y = IodineHandling, fill = residual_group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.12, size = 0.55, alpha = 0.35) +
  scale_fill_manual(values = c("High residual RAIR" = pal["red"], "Low residual RAIR" = pal["grey3"])) +
  theme_jtm() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 18, hjust = 1)) +
  labs(title = "MAPK-adjusted residual retains iodine signal", x = NULL, y = "Iodine-handling score")

forest_lm <- tcga_resid[model %like% "IodineHandling" & term %in% c("RAIR_resid_MAPK", "OneCarbon_resid_MAPK"),
                        .(analysis = ifelse(term == "RAIR_resid_MAPK", "RAIR residual", "One-carbon residual"),
                          effect = estimate, ci_low, ci_high, p_value, type = "beta")]
forest_glm <- tcga_resid[model %like% "N1" & term %in% c("RAIR_resid_MAPK", "OneCarbon_resid_MAPK"),
                         .(analysis = ifelse(term == "RAIR_resid_MAPK", "RAIR residual", "One-carbon residual"),
                           effect = log(OR), ci_low = log(ci_low), ci_high = log(ci_high), p_value, type = "log OR")]
forest <- rbindlist(list(forest_lm, forest_glm), fill = TRUE)
forest[, label := paste0(analysis, "\n", type)]
forest[, label := factor(label, levels = rev(unique(label)))]
p2d <- ggplot(forest, aes(x = effect, y = label, colour = analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = pal["grey3"]) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.45) +
  geom_point(size = 2.0) +
  geom_text(aes(label = paste0("P=", pval_fmt(p_value))), hjust = -0.05, size = 2.0, colour = pal["grey1"]) +
  scale_colour_manual(values = c("RAIR residual" = pal["red"], "One-carbon residual" = pal["teal"])) +
  coord_cartesian(xlim = c(min(forest$ci_low, na.rm = TRUE) - 0.2, max(forest$ci_high, na.rm = TRUE) + 1.0)) +
  theme_jtm() +
  theme(legend.position = "none") +
  labs(title = "Residual models beyond MAPK/CAF/immune", x = "Standardized effect", y = NULL)

hall_plot <- hallmark[!is.na(padj) & NES > 0][order(padj)][1:10]
hall_plot[, pathway_clean := gsub("HALLMARK_|_", " ", pathway)]
hall_plot[, pathway_clean := factor(pathway_clean, levels = pathway_clean[order(NES)])]
p2e <- ggplot(hall_plot, aes(x = pathway_clean, y = NES, fill = -log10(padj))) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_gradient(low = pal["sky"], high = pal["red"]) +
  theme_jtm() +
  labs(title = "RAIR-high tumors show EMT/inflammatory programs", x = NULL, y = "Hallmark NES", fill = "-log10 FDR")

fig2 <- ((p2a | p2b) / (p2c | p2d) / p2e) +
  plot_layout(heights = c(1, 1, 0.9)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig2, file.path(main_fig_dir, "Figure_2_TCGA_discovery"), 190, 180)
fwrite(tcga[, .(patient, driver_group, RAIR_like, IodineHandling, RAIR_resid_MAPK)], file.path(source_dir, "Figure_2ABC_TCGA_scores.csv"))
fwrite(forest, file.path(source_dir, "Figure_2D_residual_models.csv"))
fwrite(hall_plot, file.path(source_dir, "Figure_2E_hallmark.csv"))

## -----------------------------
## Main Figure 3: GEO validation.
## -----------------------------
message("Drawing Figure 3...")
geo <- rbindlist(list(
  gse33630[, .(dataset, sample, histology_group, TDS, OneCarbon, RAIR_like)],
  gse76039[, .(dataset, sample, histology_group, TDS, OneCarbon, RAIR_like)]
), fill = TRUE)
geo[, histology_group := factor(histology_group, levels = c("Normal", "PTC", "PDTC", "ATC"))]
geo_long <- score_long(geo)
geo_long[, score := factor(score, levels = c("TDS", "OneCarbon", "RAIR_like"))]
p3a <- ggplot(geo_long[score %in% c("TDS", "OneCarbon", "RAIR_like")],
              aes(x = histology_group, y = value, fill = histology_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.76) +
  geom_jitter(width = 0.12, size = 0.65, alpha = 0.45) +
  facet_grid(score ~ dataset, scales = "free", space = "free_x") +
  scale_fill_manual(values = c("Normal" = pal["grey3"], "PTC" = pal["sky"], "PDTC" = pal["orange"], "ATC" = pal["red"]),
                    na.value = pal["grey3"]) +
  theme_jtm() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "The same state score follows the PTC-PDTC-ATC continuum", x = NULL, y = "Module score", fill = "Histology")

auc_dt <- bulk_tests[dataset_test %in% c("GSE33630_ATC_vs_PTC", "GSE76039_ATC_vs_PDTC") &
                       score %in% c("TDS", "OneCarbon", "RAIR_like")]
auc_dt[, dataset_test := factor(dataset_test, levels = c("GSE33630_ATC_vs_PTC", "GSE76039_ATC_vs_PDTC"),
                                labels = c("GSE33630\nATC vs PTC", "GSE76039\nATC vs PDTC"))]
auc_dt[, score := factor(score, levels = c("TDS", "OneCarbon", "RAIR_like"))]
p3b <- ggplot(auc_dt, aes(x = score, y = auc, fill = score)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = pal["grey2"], linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", auc)), vjust = -0.35, size = 2.0) +
  facet_wrap(~ dataset_test) +
  scale_fill_manual(values = c("TDS" = pal["blue"], "OneCarbon" = pal["teal"], "RAIR_like" = pal["red"])) +
  coord_cartesian(ylim = c(0, 1.08)) +
  theme_jtm() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Validation performance", x = NULL, y = "AUC")

target_genes <- c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1", "MTHFD1L", "SHMT2")
target_heat <- external_targets[gene %in% target_genes,
                                .(gene, dataset, comparison, delta = delta_positive_negative, p_adj)]
target_heat[, sig := p_adj < 0.05]
target_heat[, gene := factor(gene, levels = rev(target_genes))]
target_heat[, cohort := paste(dataset, comparison, sep = "\n")]
p3c <- ggplot(target_heat, aes(x = cohort, y = gene, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = ifelse(sig, "*", "")), size = 3.0) +
  scale_fill_gradient2(low = pal["blue"], mid = "white", high = pal["red"], midpoint = 0) +
  theme_jtm() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "External validation of candidate nodes", fill = "Median delta")

fig3 <- (p3a / (p3b | p3c)) +
  plot_layout(heights = c(1.35, 1), widths = c(0.8, 1.2), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"), legend.position = "right")
save_pub(fig3, file.path(main_fig_dir, "Figure_3_GEO_validation"), 183, 150)
fwrite(geo_long, file.path(source_dir, "Figure_3A_GEO_score_long.csv"))
fwrite(auc_dt, file.path(source_dir, "Figure_3B_AUC.csv"))
fwrite(target_heat, file.path(source_dir, "Figure_3C_external_target_validation.csv"))

## -----------------------------
## Main Figure 4: single-cell/copyKAT.
## -----------------------------
message("Drawing Figure 4...")
single_top[, cell_type := factor(cell_type, levels = c("Epithelial_Thyroid", "Myeloid", "CAF", "Endothelial", "T_NK", "B_Plasma", "Other"))]
p4a <- ggplot(single_top, aes(x = dataset, y = fraction_of_top10, fill = cell_type)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Epithelial_Thyroid" = pal["red"], "Myeloid" = pal["orange"],
                               "CAF" = pal["purple"], "Endothelial" = pal["teal"],
                               "T_NK" = pal["blue"], "B_Plasma" = pal["sky"], "Other" = pal["grey3"]),
                    labels = c("Epithelial_Thyroid" = "Epithelial/thyroid", "Myeloid" = "Myeloid",
                               "CAF" = "CAF", "Endothelial" = "Endothelial", "T_NK" = "T/NK",
                               "B_Plasma" = "B/plasma", "Other" = "Other")) +
  theme_jtm() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Top RAIR-like cells localize to epithelial/thyroid lineage",
       x = NULL, y = "Fraction among top 10% RAIR-like cells", fill = "Cell type")

sc_plot <- copy(single_bal)
sc_plot[, non_epithelial_fraction := Myeloid + CAF + T_NK + B_Plasma + Endothelial]
p4b <- ggplot(sc_plot, aes(x = non_epithelial_fraction, y = epithelial_RAIR, colour = tissue_group)) +
  geom_point(size = 2.0, alpha = 0.85) +
  facet_wrap(~ dataset, scales = "free") +
  scale_colour_manual(values = c("Paratumor" = pal["blue"], "Primary_tumor" = pal["teal"],
                                 "Lymph_node_metastasis" = pal["orange"], "Subcutaneous_metastasis" = pal["red"],
                                 "PTC" = pal["sky"], "ATC" = pal["red"]),
                      labels = c("ATC" = "ATC", "Lymph_node_metastasis" = "LN metastasis",
                                 "Paratumor" = "Paratumor", "Primary_tumor" = "Primary tumor",
                                 "PTC" = "PTC", "Subcutaneous_metastasis" = "Subcutaneous metastasis")) +
  theme_jtm() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Sample epithelial RAIR-like vs non-epithelial fraction",
       x = "Non-epithelial fraction", y = "Sample epithelial RAIR-like", colour = "Tissue")

copykat_plot <- rbindlist(list(
  copykat_summary[, .(analysis = "Immune/stromal reference", dataset, sample, tissue_group,
                      n_epithelial, aneuploid_fraction = epithelial_aneuploid_fraction)],
  copykat_para[, .(analysis = "Paratumor epithelial reference", dataset, sample = sample_target,
                   tissue_group = "Tumor/met target", n_epithelial = n_target_epi,
                   aneuploid_fraction)]
), fill = TRUE)
copykat_plot[, analysis := fifelse(analysis == "Immune/stromal reference", "Immune/stromal ref.",
                                   "Paratumor epithelial ref.")]
p4c <- ggplot(copykat_plot, aes(x = reorder(sample, aneuploid_fraction), y = aneuploid_fraction, fill = analysis)) +
  geom_col(width = 0.70) +
  coord_flip() +
  facet_wrap(~ analysis, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("Immune/stromal ref." = pal["grey2"], "Paratumor epithelial ref." = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "none", strip.text = element_text(size = 5.8, face = "bold")) +
  labs(title = "copyKAT aneuploid fraction in epithelial cells",
       x = NULL, y = "Aneuploid fraction", fill = "Tissue")

pred_plot <- copykat_pred[cell_type == "Epithelial_Thyroid" & copykat_class %in% c("aneuploid", "diploid")]
pred_plot[, copykat_class := factor(copykat_class, levels = c("diploid", "aneuploid"))]
p4d <- ggplot(pred_plot, aes(x = copykat_class, y = RAIR_like, fill = copykat_class)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.78) +
  geom_jitter(width = 0.16, size = 0.38, alpha = 0.22) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("diploid" = pal["blue"], "aneuploid" = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "none") +
  labs(title = "RAIR-like by copyKAT class",
       x = "copyKAT class", y = "RAIR-like score")

fig4 <- (p4a | p4b) / (p4c | p4d) +
  plot_layout(widths = c(1, 1.08), heights = c(1, 1.02)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig4, file.path(main_fig_dir, "Figure_4_single_cell_copyKAT"), 190, 165)
fwrite(single_top, file.path(source_dir, "Figure_4A_top_RAIR_celltype_fraction.csv"))
fwrite(sc_plot, file.path(source_dir, "Figure_4B_sample_composition.csv"))
fwrite(copykat_plot, file.path(source_dir, "Figure_4C_copykat_aneuploid_fraction.csv"))
fwrite(pred_plot[, .(dataset, sample, cell, copykat_class, RAIR_like, tissue_group)], file.path(source_dir, "Figure_4D_copykat_RAIR_by_class.csv"))

## -----------------------------
## Main Figure 5: target translation.
## -----------------------------
message("Drawing Figure 5...")
tier_plot <- target_tiers[gene %in% c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1", "MTHFD1L", "SHMT2")]
tier_plot[, validation := factor(n_external_up_validations, levels = c(0, 1, 2))]
p5a <- ggplot(tier_plot, aes(x = rho_RAIR_like, y = delta_high_low,
                             size = n_external_up_validations, colour = evidence_tier)) +
  geom_hline(yintercept = 0, colour = pal["grey3"], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = pal["grey3"], linewidth = 0.3) +
  geom_point(alpha = 0.88) +
  ggrepel::geom_text_repel(aes(label = gene), size = 2.2, max.overlaps = 20, min.segment.length = 0.1) +
  scale_colour_manual(values = c("Tier 1A_core_actionable_axis" = pal["red"],
                                 "Tier 1B_ecology_state_axis" = pal["orange"],
                                 "Tier 2_supportive_marker_axis" = pal["blue"]),
                      labels = c("Tier 1A_core_actionable_axis" = "Core actionable",
                                 "Tier 1B_ecology_state_axis" = "Ecology/state",
                                 "Tier 2_supportive_marker_axis" = "Supportive")) +
  scale_size_continuous(range = c(1.8, 4.2), breaks = c(0, 1, 2)) +
  theme_jtm() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), size = guide_legend(nrow = 1)) +
  labs(title = "Prioritized nodes within the RAIR-like state",
       x = "TCGA rho with RAIR-like", y = "RAIR-high vs low expression delta",
       colour = "Tier", size = "GEO validations")

hpa_long <- melt(hpa_path[gene %in% c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1")],
                 id.vars = "gene",
                 measure.vars = c("pathology_high", "pathology_medium", "pathology_low", "pathology_not_detected"),
                 variable.name = "category", value.name = "n_cases")
hpa_long[, category := factor(category,
                              levels = c("pathology_high", "pathology_medium", "pathology_low", "pathology_not_detected"),
                              labels = c("High", "Medium", "Low", "Not detected"))]
p5b <- ggplot(hpa_long, aes(x = gene, y = n_cases, fill = category)) +
  geom_col(width = 0.70) +
  scale_fill_manual(values = c("High" = "#8C1D1D", "Medium" = "#C65F38", "Low" = "#D5A253", "Not detected" = pal["grey3"])) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(title = "HPA thyroid cancer IHC",
       x = NULL, y = "Cases", fill = "IHC")

chem_plot <- chembl
chem_plot[, gene := factor(gene, levels = gene[order(activity_total_count)])]
p5c <- ggplot(chem_plot, aes(x = gene, y = activity_total_count, fill = actionability_bin)) +
  geom_col(width = 0.70) +
  coord_flip() +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c("Clinical/approved mechanism record" = pal["red"],
                               "Preclinical chemical tractability" = pal["teal"],
                               "Sparse chemical evidence" = pal["grey2"]),
                    labels = c("Clinical/approved mechanism record" = "Clinical/approved",
                               "Preclinical chemical tractability" = "Preclinical",
                               "Sparse chemical evidence" = "Sparse")) +
  theme_jtm() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "ChEMBL actionability",
       x = NULL, y = "Activity records, log scale", fill = "Actionability")

matrix_dt <- trans_target[gene %in% c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1"),
                          .(gene,
                            `TCGA RAIR rho` = rho_RAIR_like,
                            `GEO validations` = n_external_up_validations / 2,
                            `HPA high/medium` = pathology_high_medium_fraction,
                            `ChEMBL tractability` = fifelse(actionability_bin == "Clinical/approved mechanism record", 1,
                                                          fifelse(actionability_bin == "Preclinical chemical tractability", 0.65,
                                                                  fifelse(actionability_bin == "Sparse chemical evidence", 0.25, 0))))]
matrix_long <- melt(matrix_dt, id.vars = "gene", variable.name = "evidence", value.name = "scaled_value")
matrix_long[, gene := factor(gene, levels = rev(c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1")))]
p5d <- ggplot(matrix_long, aes(x = evidence, y = gene, fill = scaled_value)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", scaled_value)), size = 1.9) +
  scale_fill_gradient(low = "white", high = pal["red"], limits = c(0, 1), oob = squish) +
  theme_jtm() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  labs(title = "Integrated translational evidence matrix", fill = "Scaled")

fig5 <- (p5a | p5d) / (p5b | p5c) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig5, file.path(main_fig_dir, "Figure_5_target_translation"), 190, 165)
fwrite(tier_plot, file.path(source_dir, "Figure_5A_target_tiers.csv"))
fwrite(hpa_long, file.path(source_dir, "Figure_5B_HPA_IHC.csv"))
fwrite(chem_plot, file.path(source_dir, "Figure_5C_ChEMBL.csv"))
fwrite(matrix_long, file.path(source_dir, "Figure_5D_integrated_matrix.csv"))

## -----------------------------
## Supplementary figures.
## -----------------------------
message("Drawing supplementary figures...")
g151_long <- score_long(gse151179, c("TDS", "OneCarbon", "RAIR_like"))
g151_long <- g151_long[!is.na(rai_group) & rai_group %in% c("RAI_refractory", "RAI_avid")]
p_s1a <- ggplot(g151_long[score %in% c("TDS", "RAIR_like")],
                aes(x = rai_group, y = value, fill = rai_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.76) +
  geom_jitter(width = 0.12, size = 0.75, alpha = 0.6) +
  facet_wrap(~ score, scales = "free_y") +
  scale_fill_manual(values = c("RAI_avid" = pal["blue"], "RAI_refractory" = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "GSE151179 direct RAI endpoint is underpowered/negative", x = NULL, y = "Score")
s1_auc <- bulk_tests[dataset_test == "GSE151179_RAI" & score %in% c("TDS", "OneCarbon", "RAIR_like")]
p_s1b <- ggplot(s1_auc, aes(x = score, y = auc, fill = score)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = pal["grey2"]) +
  geom_text(aes(label = sprintf("AUC %.2f\nP=%s", auc, pval_fmt(p_wilcox))), vjust = -0.25, size = 2.0) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_fill_manual(values = c("TDS" = pal["blue"], "OneCarbon" = pal["teal"], "RAIR_like" = pal["red"])) +
  theme_jtm() +
  theme(legend.position = "none") +
  labs(title = "Direct RAI discrimination did not pass go/no-go", x = NULL, y = "AUC")
fig_s1 <- (p_s1a | p_s1b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig_s1, file.path(supp_fig_dir, "Figure_S1_GSE151179_RAI_endpoint"), 183, 85)
fwrite(g151_long, file.path(source_dir, "Figure_S1A_GSE151179_scores.csv"))
fwrite(s1_auc, file.path(source_dir, "Figure_S1B_GSE151179_auc.csv"))

p_s2a <- ggplot(tcga[driver_group == "BRAF-mut"], aes(x = OneCarbon, y = IodineHandling)) +
  geom_point(size = 1.1, alpha = 0.65, colour = pal["grey1"]) +
  geom_smooth(method = "lm", se = TRUE, colour = pal["red"]) +
  theme_jtm() +
  labs(title = "BRAF-mutant: one-carbon alone is not a low-iodine surrogate", x = "One-carbon score", y = "Iodine-handling")
target_cor_plot <- tcga_target_cor[gene %in% c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1")]
target_cor_plot[, gene := factor(gene, levels = gene[order(rho_IodineHandling)])]
p_s2b <- ggplot(target_cor_plot, aes(x = gene, y = rho_IodineHandling, fill = p_adj_IodineHandling < 0.05)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = pal["red"], "FALSE" = pal["grey3"])) +
  theme_jtm() +
  theme(legend.position = "none") +
  labs(title = "Target genes versus iodine-handling within BRAF-mutant tumors", x = NULL, y = "Spearman rho")
fig_s2 <- (p_s2a | p_s2b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig_s2, file.path(supp_fig_dir, "Figure_S2_BRAF_MAPK_residual_detail"), 183, 85)
fwrite(tcga_target_cor, file.path(source_dir, "Figure_S2B_BRAF_target_correlations.csv"))

react_plot <- reactome[!is.na(padj) & NES > 0][order(padj)][1:12]
react_plot[, pathway_clean := pathway_label(pathway, width = 46)]
react_plot[, pathway_clean := factor(pathway_clean, levels = pathway_clean[order(NES)])]
react_plot[, neg_log10_fdr := -log10(padj)]
p_s3a <- ggplot(react_plot, aes(x = NES, y = pathway_clean)) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_clean),
               linewidth = 0.42, colour = pal["grey3"], lineend = "round") +
  geom_point(aes(fill = neg_log10_fdr), shape = 21, size = 2.4,
             colour = "white", stroke = 0.25) +
  scale_x_continuous(limits = c(0, 2.05), breaks = seq(0, 2, 0.5),
                     expand = expansion(mult = c(0, 0.04))) +
  scale_fill_gradient(low = pal["sky"], high = pal["red"],
                      guide = guide_colorbar(title.position = "top",
                                             barwidth = unit(28, "mm"),
                                             barheight = unit(2.4, "mm"))) +
  theme_jtm(base_size = 6.6) +
  theme(axis.text.y = element_text(size = 5.6, lineheight = 0.90),
        legend.position = "bottom",
        legend.justification = "left",
        plot.margin = margin(4, 7, 4, 4)) +
  labs(title = "Reactome pathways enriched in high-score tumors",
       x = "Normalized enrichment score", y = NULL, fill = "-log10 FDR")
hall_neg <- hallmark[!is.na(padj)][order(NES)][1:10]
hall_neg[, pathway_clean := pathway_label(pathway, width = 46)]
hall_neg[, pathway_clean := factor(pathway_clean, levels = pathway_clean[order(NES)])]
hall_neg[, neg_log10_fdr := -log10(padj)]
p_s3b <- ggplot(hall_neg, aes(x = NES, y = pathway_clean)) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_clean),
               linewidth = 0.42, colour = pal["grey3"], lineend = "round") +
  geom_point(aes(fill = neg_log10_fdr), shape = 21, size = 2.4,
             colour = "white", stroke = 0.25) +
  scale_x_continuous(limits = c(0, 1.50), breaks = seq(0, 1.5, 0.5),
                     expand = expansion(mult = c(0, 0.04))) +
  scale_fill_gradient(low = pal["sky"], high = pal["red"],
                      guide = guide_colorbar(title.position = "top",
                                             barwidth = unit(28, "mm"),
                                             barheight = unit(2.4, "mm"))) +
  theme_jtm(base_size = 6.6) +
  theme(axis.text.y = element_text(size = 5.6, lineheight = 0.90),
        legend.position = "bottom",
        legend.justification = "left",
        plot.margin = margin(4, 7, 4, 4)) +
  labs(title = "Least-enriched Hallmark programs",
       subtitle = "All Hallmark gene sets had positive NES values in this contrast",
       x = "Normalized enrichment score", y = NULL, fill = "-log10 FDR")
fig_s3 <- (p_s3a / p_s3b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8.5, face = "bold"))
save_pub(fig_s3, file.path(supp_fig_dir, "Figure_S3_GSEA_detail"), 183, 165)
fwrite(react_plot, file.path(source_dir, "Figure_S3A_reactome.csv"))
fwrite(hall_neg, file.path(source_dir, "Figure_S3B_hallmark_lowest_NES.csv"))

cnv_plot <- cnv_summary[n_epithelial >= 50]
cnv_plot[, tissue_label := fcase(
  tissue_group == "Lymph_node_metastasis", "LN metastasis",
  tissue_group == "Primary_tumor", "Primary tumor",
  tissue_group == "Subcutaneous_metastasis", "SC metastasis",
  default = tissue_group
)]
p_s4a <- ggplot(cnv_plot, aes(x = tissue_label, y = fraction_cnv_epithelial_ref_high, colour = tumor_met_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) +
  geom_jitter(width = 0.14, size = 1.8, alpha = 0.8) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
  labs(title = "Lightweight CNV proxy does not explain RAIR-like signal", x = NULL, y = "CNV-high proxy fraction")
p_s4b <- ggplot(copykat_para_pred[tissue_group != "Paratumor" & cell_type == "Epithelial_Thyroid"],
                aes(x = sample_target, y = RAIR_like, fill = copykat_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.15, size = 0.40, alpha = 0.25) +
  scale_fill_manual(values = c("aneuploid" = pal["red"], "diploid" = pal["blue"])) +
  scale_x_discrete(labels = c("PTC10_RightLN" = "PTC10 LN", "PTC10_T" = "PTC10 T", "PTC11_SC" = "PTC11 SC")) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom") +
  labs(title = "GSE184362 copyKAT with paratumor epithelial reference", x = NULL, y = "RAIR-like", fill = "copyKAT")
fig_s4 <- (p_s4a | p_s4b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig_s4, file.path(supp_fig_dir, "Figure_S4_CNV_copyKAT_sensitivity"), 183, 85)
fwrite(cnv_plot, file.path(source_dir, "Figure_S4A_CNV_proxy.csv"))
fwrite(copykat_para_pred, file.path(source_dir, "Figure_S4B_copykat_paratumor_reference.csv"))

normal_hpa <- fread(file.path(res_dir, "HPA_target_normal_thyroid_IHC.csv"))
normal_hpa[, gene := factor(gene, levels = c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1"))]
p_s5a <- ggplot(normal_hpa, aes(x = gene, y = normal_thyroid_IHC_level, colour = normal_reliability)) +
  geom_point(size = 2.4) +
  theme_jtm() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Normal thyroid HPA IHC context", x = NULL, y = "Normal glandular cell level", colour = "Reliability")
p_s5b <- ggplot(matrix_long, aes(x = evidence, y = gene, fill = scaled_value)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", scaled_value)), size = 1.9) +
  scale_fill_gradient(low = "white", high = pal["red"], limits = c(0, 1), oob = squish) +
  theme_jtm() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Detailed target evidence scoring", fill = "Scaled")
fig_s5 <- (p_s5a | p_s5b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_pub(fig_s5, file.path(supp_fig_dir, "Figure_S5_target_annotation_detail"), 183, 85)
fwrite(normal_hpa, file.path(source_dir, "Figure_S5A_HPA_normal_thyroid.csv"))
fwrite(matrix_long, file.path(source_dir, "Figure_S5B_target_evidence_matrix.csv"))

## -----------------------------
## Tables.
## -----------------------------
message("Writing publication tables...")
table1 <- cohort_dt[, .(Cohort = as.character(cohort), Modality = modality, Role = role, `N samples/cells` = n)]
table2 <- rbindlist(list(
  data.table(Evidence = "TCGA RAIR residual -> iodine-handling",
             Cohort = "TCGA-THCA", Statistic = "beta", Effect = -1.016, P = 4.459886e-37,
             Interpretation = "RAIR-like retains iodine-handling loss beyond MAPK/CAF/immune/drivers"),
  data.table(Evidence = "TCGA RAIR residual -> N1",
             Cohort = "TCGA-THCA", Statistic = "OR", Effect = 2.158, P = 1.423127e-03,
             Interpretation = "Residual RAIR-like associated with nodal disease"),
  data.table(Evidence = "GSE33630 ATC vs PTC",
             Cohort = "GSE33630", Statistic = "RAIR-like AUC", Effect = auc_dt[dataset_test == "GSE33630\nATC vs PTC" & score == "RAIR_like", auc], P = auc_dt[dataset_test == "GSE33630\nATC vs PTC" & score == "RAIR_like", p_wilcox],
             Interpretation = "External dedifferentiation validation"),
  data.table(Evidence = "GSE76039 ATC vs PDTC",
             Cohort = "GSE76039", Statistic = "RAIR-like AUC", Effect = auc_dt[dataset_test == "GSE76039\nATC vs PDTC" & score == "RAIR_like", auc], P = auc_dt[dataset_test == "GSE76039\nATC vs PDTC" & score == "RAIR_like", p_wilcox],
             Interpretation = "External high-grade validation"),
  data.table(Evidence = "GSE232237 copyKAT aneuploid vs diploid epithelial",
             Cohort = "GSE232237", Statistic = "Median delta", Effect = copykat_tests[dataset == "GSE232237", delta_RAIR_aneuploid_minus_diploid], P = copykat_tests[dataset == "GSE232237", p_wilcox],
             Interpretation = "RAIR-like enriched in copyKAT-aneuploid epithelial cells")
), fill = TRUE)

table3 <- trans_target[, .(
  Gene = gene,
  Tier = evidence_tier,
  Axis = biological_axis,
  `External validations` = n_external_up_validations,
  `TCGA rho RAIR-like` = rho_RAIR_like,
  `HPA thyroid cancer high/medium fraction` = pathology_high_medium_fraction,
  `ChEMBL activity records` = activity_total_count,
  `ChEMBL max phase` = max_mechanism_phase,
  `Actionability bin` = actionability_bin,
  `Recommended interpretation` = actionability_note
)]

supp_tables <- list(
  Table_1_Cohorts = table1,
  Table_2_Key_Evidence = table2,
  Table_3_Target_Prioritization = table3,
  Sup_Bulk_Tests = bulk_tests,
  Sup_TCGA_Residual = tcga_resid,
  Sup_BRAF_Internal = tcga_braf,
  Sup_CopyKAT = copykat_summary,
  Sup_HPA_ChEMBL = trans_target
)

fwrite(table1, file.path(table_dir, "Table_1_cohorts.csv"))
fwrite(table2, file.path(table_dir, "Table_2_key_evidence.csv"))
fwrite(table3, file.path(table_dir, "Table_3_target_prioritization.csv"))
wb <- createWorkbook()
for (nm in names(supp_tables)) {
  addWorksheet(wb, nm)
  writeDataTable(wb, nm, supp_tables[[nm]])
  setColWidths(wb, nm, cols = 1:ncol(supp_tables[[nm]]), widths = "auto")
}
saveWorkbook(wb, file.path(table_dir, "JTM_tables_and_supplementary_tables.xlsx"), overwrite = TRUE)

## Figure manifest.
manifest <- data.table(
  item = c(paste0("Figure ", 1:5), paste0("Figure S", 1:5)),
  file_base = c("Figure_1_study_design", "Figure_2_TCGA_discovery",
                "Figure_3_GEO_validation", "Figure_4_single_cell_copyKAT",
                "Figure_5_target_translation",
                "Figure_S1_GSE151179_RAI_endpoint", "Figure_S2_BRAF_MAPK_residual_detail",
                "Figure_S3_GSEA_detail", "Figure_S4_CNV_copyKAT_sensitivity",
                "Figure_S5_target_annotation_detail"),
  conclusion = c(
    "Multi-cohort design tests a composite RAIR-like dedifferentiation state and targetable nodes.",
    "TCGA identifies a MAPK/BRAF-related but nonredundant RAIR-like state linked to iodine-handling loss.",
    "Independent GEO cohorts validate the RAIR-like dedifferentiation continuum and candidate nodes.",
    "Single-cell and copyKAT analyses localize RAIR-like signal to malignant-like epithelial/thyroid compartments.",
    "SLC7A5/SLC1A5/MTHFD2/AXL form a prioritized translational target map.",
    "Direct RAI endpoint in GSE151179 is underpowered/negative.",
    "BRAF-mutant internal analysis refines the role of one-carbon metabolism.",
    "GSEA details support EMT, inflammatory, RTK and cell-cycle programs.",
    "CNV proxy and copyKAT reference sensitivity define localization limits.",
    "HPA and ChEMBL details support revised target prioritization."
  )
)
fwrite(manifest, file.path(asset_dir, "figure_manifest.csv"))

summary_file <- file(file.path(asset_dir, "publication_asset_summary.txt"), open = "wt", encoding = "UTF-8")
sink(summary_file)
cat("JTM manuscript figure/table asset package\n\n")
cat("Main figures exported to:\nanalysis/manuscript_assets/figures/main\n\n", sep = "")
cat("Supplementary figures exported to:\nanalysis/manuscript_assets/figures/supplementary\n\n", sep = "")
cat("Tables exported to:\nanalysis/manuscript_assets/tables\n\n", sep = "")
cat("Source data exported to:\nanalysis/manuscript_assets/source_data\n\n", sep = "")
cat("Figure manifest:\n")
print(manifest, row.names = FALSE)
cat("\nTables:\n")
cat("- Table_1_cohorts.csv\n- Table_2_key_evidence.csv\n- Table_3_target_prioritization.csv\n- JTM_tables_and_supplementary_tables.xlsx\n")
sink()
close(summary_file)

message("Publication figure/table package complete: ", asset_dir)
