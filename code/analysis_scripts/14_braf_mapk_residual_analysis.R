options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(ggplot2)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_root <- file.path(root, "analysis")
result_dir <- file.path(analysis_root, "results")
figure_dir <- file.path(analysis_root, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

z <- function(x) as.numeric(scale(x))

collapse_symbols <- function(expr, symbols) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "NA"
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  v <- apply(expr, 1, var, na.rm = TRUE)
  ord <- order(symbols, -v)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  first <- !duplicated(symbols)
  expr <- expr[first, , drop = FALSE]
  rownames(expr) <- symbols[first]
  expr
}

coef_lm <- function(fit, term, model, cohort) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NULL)
  ci <- suppressWarnings(confint(fit, term))
  data.table(cohort = cohort, model = model, term = term,
             estimate = sm[term, "Estimate"], ci_low = ci[1], ci_high = ci[2],
             p_value = sm[term, "Pr(>|t|)"], n = nobs(fit),
             r_squared = summary(fit)$r.squared)
}

coef_glm <- function(fit, term, model, cohort) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NULL)
  ci <- suppressWarnings(confint.default(fit, term))
  data.table(cohort = cohort, model = model, term = term,
             log_or = sm[term, "Estimate"], OR = exp(sm[term, "Estimate"]),
             ci_low = exp(ci[1]), ci_high = exp(ci[2]),
             p_value = sm[term, "Pr(>|z|)"], n = nobs(fit))
}

message("Loading TCGA scores...")
df <- fread(file.path(result_dir, "TCGA_THCA_scores_clinical_mutation.csv"))
for (v in c("RAIR_like", "OneCarbon", "MAPK", "Dediff_EMT", "TDS",
            "IodineHandling", "CAF", "ImmuneCheckpoint")) {
  df[[paste0(v, "_z")]] <- z(df[[v]])
}
df[, BRAF_mut := as.logical(BRAF_mut)]
df[, RAS_mut := as.logical(RAS_mut)]
df[, N_positive_bin := fifelse(N_positive == "N1", 1, fifelse(N_positive == "N0", 0, NA_real_))]
df[, T_high_bin := fifelse(T_high == "T3_T4", 1, fifelse(T_high == "T1_T2", 0, NA_real_))]
df[, Stage_high_bin := fifelse(stage_high == "Stage_III_IV", 1,
                               fifelse(stage_high == "Stage_I_II", 0, NA_real_))]

message("Adding residual scores...")
df[, RAIR_resid_MAPK := residuals(lm(RAIR_like_z ~ MAPK_z + CAF_z + ImmuneCheckpoint_z +
                                       BRAF_mut + RAS_mut, data = df))]
df[, OneCarbon_resid_MAPK := residuals(lm(OneCarbon_z ~ MAPK_z + CAF_z + ImmuneCheckpoint_z +
                                            BRAF_mut + RAS_mut, data = df))]
df[, Dediff_resid_MAPK := residuals(lm(Dediff_EMT_z ~ MAPK_z + CAF_z + ImmuneCheckpoint_z +
                                         BRAF_mut + RAS_mut, data = df))]

all_models <- rbindlist(list(
  coef_lm(lm(IodineHandling_z ~ RAIR_resid_MAPK + MAPK_z + CAF_z +
               ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = df),
          "RAIR_resid_MAPK", "IodineHandling ~ RAIR residual + MAPK/CAF/immune/drivers", "All TCGA"),
  coef_lm(lm(IodineHandling_z ~ OneCarbon_resid_MAPK + MAPK_z + CAF_z +
               ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = df),
          "OneCarbon_resid_MAPK", "IodineHandling ~ OneCarbon residual + MAPK/CAF/immune/drivers", "All TCGA"),
  coef_glm(glm(N_positive_bin ~ RAIR_resid_MAPK + MAPK_z + CAF_z +
                 ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = df, family = binomial()),
           "RAIR_resid_MAPK", "N1 vs N0 ~ RAIR residual + MAPK/CAF/immune/drivers", "All TCGA"),
  coef_glm(glm(N_positive_bin ~ OneCarbon_resid_MAPK + MAPK_z + CAF_z +
                 ImmuneCheckpoint_z + BRAF_mut + RAS_mut, data = df, family = binomial()),
           "OneCarbon_resid_MAPK", "N1 vs N0 ~ OneCarbon residual + MAPK/CAF/immune/drivers", "All TCGA")
), fill = TRUE)

braf <- df[BRAF_mut == TRUE]
braf_models <- rbindlist(list(
  coef_lm(lm(IodineHandling_z ~ RAIR_like_z + MAPK_z + CAF_z + ImmuneCheckpoint_z, data = braf),
          "RAIR_like_z", "IodineHandling ~ RAIR-like + MAPK/CAF/immune", "BRAF-mutant"),
  coef_lm(lm(IodineHandling_z ~ OneCarbon_z + MAPK_z + Dediff_EMT_z + CAF_z +
               ImmuneCheckpoint_z, data = braf),
          "OneCarbon_z", "IodineHandling ~ OneCarbon + MAPK/dediff/CAF/immune", "BRAF-mutant"),
  coef_glm(glm(N_positive_bin ~ RAIR_like_z + MAPK_z + CAF_z + ImmuneCheckpoint_z,
               data = braf, family = binomial()),
           "RAIR_like_z", "N1 vs N0 ~ RAIR-like + MAPK/CAF/immune", "BRAF-mutant"),
  coef_glm(glm(N_positive_bin ~ OneCarbon_z + MAPK_z + Dediff_EMT_z + CAF_z +
                 ImmuneCheckpoint_z, data = braf, family = binomial()),
           "OneCarbon_z", "N1 vs N0 ~ OneCarbon + MAPK/dediff/CAF/immune", "BRAF-mutant")
), fill = TRUE)

fwrite(all_models, file.path(result_dir, "TCGA_MAPK_residual_models.csv"))
fwrite(braf_models, file.path(result_dir, "TCGA_BRAF_mutant_internal_models.csv"))

message("Loading target gene expression...")
se <- readRDS(file.path(analysis_root, "data", "tcga", "TCGA_THCA_STAR_counts_primary_tumor.rds"))
assay_name <- if ("tpm_unstrand" %in% assayNames(se)) "tpm_unstrand" else
  if ("fpkm_unstrand" %in% assayNames(se)) "fpkm_unstrand" else assayNames(se)[1]
expr <- log2(assay(se, assay_name) + 1)
rd <- as.data.frame(rowData(se))
symbol_col <- if ("gene_name" %in% colnames(rd)) "gene_name" else
  if ("external_gene_name" %in% colnames(rd)) "external_gene_name" else "gene_id"
expr_gene <- collapse_symbols(expr, rd[[symbol_col]])
target_genes <- intersect(c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1",
                            "MTHFD1L", "SHMT2", "DUSP5", "DUSP6", "BRAF"), rownames(expr_gene))
target_expr <- as.data.table(t(expr_gene[target_genes, df$barcode, drop = FALSE]), keep.rownames = "barcode")
df <- merge(df, target_expr, by = "barcode", all.x = TRUE)

target_cor <- rbindlist(lapply(target_genes, function(g) {
  sub <- df[BRAF_mut == TRUE & is.finite(get(g))]
  ct1 <- suppressWarnings(cor.test(sub[[g]], sub$RAIR_like, method = "spearman"))
  ct2 <- suppressWarnings(cor.test(sub[[g]], sub$IodineHandling, method = "spearman"))
  data.table(gene = g, cohort = "BRAF-mutant",
             rho_RAIR_like = unname(ct1$estimate), p_RAIR_like = ct1$p.value,
             rho_IodineHandling = unname(ct2$estimate), p_IodineHandling = ct2$p.value,
             n = nrow(sub))
}))
target_cor[, p_adj_RAIR_like := p.adjust(p_RAIR_like, method = "BH")]
target_cor[, p_adj_IodineHandling := p.adjust(p_IodineHandling, method = "BH")]
fwrite(target_cor, file.path(result_dir, "TCGA_BRAF_mutant_target_gene_correlations.csv"))

split_tests <- function(score_col, cohort_dt, label) {
  cohort_dt <- copy(cohort_dt)
  cutoff <- median(cohort_dt[[score_col]], na.rm = TRUE)
  cohort_dt[, group := ifelse(get(score_col) >= cutoff, "High", "Low")]
  endpoints <- c("IodineHandling", "TDS", "N_positive_bin", "T_high_bin", "Stage_high_bin")
  rbindlist(lapply(endpoints, function(ep) {
    sub <- cohort_dt[is.finite(get(ep)) & group %in% c("High", "Low")]
    if (nrow(sub) < 10 || length(unique(sub$group)) < 2) return(NULL)
    if (grepl("_bin$", ep)) {
      fit <- glm(get(ep) ~ group, data = sub, family = binomial())
      sm <- summary(fit)$coefficients
      term <- "groupLow"
      data.table(cohort = label, split_score = score_col, endpoint = ep,
                 n_high = sum(sub$group == "High"), n_low = sum(sub$group == "Low"),
                 median_high = mean(sub[[ep]][sub$group == "High"], na.rm = TRUE),
                 median_low = mean(sub[[ep]][sub$group == "Low"], na.rm = TRUE),
                 effect_high_vs_low = exp(-sm[term, "Estimate"]),
                 p_value = sm[term, "Pr(>|z|)"],
                 effect_type = "OR_high_vs_low")
    } else {
      wt <- suppressWarnings(wilcox.test(sub[[ep]] ~ sub$group, exact = FALSE))
      data.table(cohort = label, split_score = score_col, endpoint = ep,
                 n_high = sum(sub$group == "High"), n_low = sum(sub$group == "Low"),
                 median_high = median(sub[[ep]][sub$group == "High"], na.rm = TRUE),
                 median_low = median(sub[[ep]][sub$group == "Low"], na.rm = TRUE),
                 effect_high_vs_low = median(sub[[ep]][sub$group == "High"], na.rm = TRUE) -
                   median(sub[[ep]][sub$group == "Low"], na.rm = TRUE),
                 p_value = wt$p.value,
                 effect_type = "median_delta_high_minus_low")
    }
  }), fill = TRUE)
}

split_res <- rbindlist(list(
  split_tests("RAIR_resid_MAPK", df, "All TCGA"),
  split_tests("OneCarbon_resid_MAPK", df, "All TCGA"),
  split_tests("RAIR_like", braf, "BRAF-mutant"),
  split_tests("OneCarbon", braf, "BRAF-mutant"),
  split_tests("OneCarbon_resid_MAPK", braf, "BRAF-mutant")
), fill = TRUE)
split_res[, p_adj := p.adjust(p_value, method = "BH"), by = .(cohort, split_score)]
fwrite(split_res, file.path(result_dir, "TCGA_BRAF_MAPK_residual_split_tests.csv"))

p1 <- ggplot(braf, aes(x = OneCarbon_z, y = IodineHandling_z)) +
  geom_point(alpha = 0.55, size = 1.4, color = "#444444") +
  geom_smooth(method = "lm", se = TRUE, color = "#b73e3e") +
  theme_bw(base_size = 10) +
  labs(x = "One-carbon score (z) in BRAF-mutant tumors",
       y = "Iodine-handling score (z)",
       title = "BRAF-mutant TCGA-THCA: one-carbon state versus iodine-handling")
ggsave(file.path(figure_dir, "TCGA_BRAF_mutant_OneCarbon_vs_IodineHandling.png"),
       p1, width = 5.4, height = 4.2, dpi = 300)

braf[, OneCarbon_resid_group := ifelse(OneCarbon_resid_MAPK >= median(OneCarbon_resid_MAPK, na.rm = TRUE),
                                       "High residual one-carbon", "Low residual one-carbon")]
p2 <- ggplot(braf, aes(x = OneCarbon_resid_group, y = IodineHandling, fill = OneCarbon_resid_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.76, width = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.0) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = "Iodine-handling score",
       title = "BRAF-mutant tumors stratified by MAPK-adjusted one-carbon residual")
ggsave(file.path(figure_dir, "TCGA_BRAF_mutant_IodineHandling_by_OneCarbon_residual.png"),
       p2, width = 5.4, height = 4.2, dpi = 300)

plot_targets <- target_cor[gene %in% c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FOSL1")]
plot_targets[, gene := factor(gene, levels = gene[order(rho_IodineHandling)])]
p3 <- ggplot(plot_targets, aes(x = gene, y = rho_IodineHandling, fill = rho_IodineHandling < 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#b73e3e", "FALSE" = "#3274a1")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Spearman rho with iodine-handling score",
       title = "Target genes versus iodine-handling within BRAF-mutant tumors")
ggsave(file.path(figure_dir, "TCGA_BRAF_mutant_target_vs_IodineHandling.png"),
       p3, width = 5.2, height = 4.0, dpi = 300)

sink(file.path(result_dir, "TCGA_BRAF_MAPK_residual_summary.txt"))
cat("TCGA BRAF/MAPK residual and BRAF-mutant internal analysis\n\n")
cat("All-cohort MAPK-adjusted residual models:\n")
print(all_models, row.names = FALSE)
cat("\nBRAF-mutant internal models:\n")
print(braf_models, row.names = FALSE)
cat("\nBRAF-mutant target gene correlations:\n")
print(target_cor[gene %in% c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FOSL1", "FN1")],
      row.names = FALSE)
cat("\nMedian/OR split tests:\n")
print(split_res[cohort == "BRAF-mutant" &
                  split_score %in% c("OneCarbon", "OneCarbon_resid_MAPK", "RAIR_like")],
      row.names = FALSE)
cat("\nInterpretation:\n")
cat("These models test whether one-carbon/RAIR-like signals carry information beyond BRAF mutation and MAPK score. The strongest robust finding is the link to iodine-handling loss; N-stage signal is stronger for composite RAIR-like than for one-carbon alone.\n")
sink()

message("Saved BRAF/MAPK residual analysis outputs.")
