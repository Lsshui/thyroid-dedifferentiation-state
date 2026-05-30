options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(pROC)
})

root <- "D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/analysis"
data_dir <- file.path(root, "data", "tcga")
result_dir <- file.path(root, "results")
figure_dir <- file.path(root, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

gene_sets <- list(
  TDS = c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1",
          "PAX8", "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB",
          "TPO", "TSHR"),
  OneCarbon = c("SHMT2", "MTHFD2", "MTHFD1L", "PHGDH", "PSAT1", "PSPH",
                "SLC1A5", "SLC7A5", "TYMS", "DHFR", "GART"),
  MAPK = c("DUSP4", "DUSP5", "DUSP6", "ETV4", "ETV5", "SPRY1", "SPRY2",
           "SPRY4", "CCND1", "FOSL1"),
  Dediff_EMT = c("HMGA2", "AXL", "VIM", "FN1", "CD44", "ITGA6", "KRT19",
                 "EPCAM", "LAMC2", "ITGB1"),
  IodineHandling = c("SLC5A5", "TPO", "TG", "TSHR", "SLC26A4", "DUOX1", "DUOX2", "DIO1", "DIO2"),
  ImmuneCheckpoint = c("CD274", "PDCD1LG2", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  CAF = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "PDGFRA", "TAGLN")
)

target_genes <- unique(c(
  gene_sets$OneCarbon, gene_sets$MAPK, gene_sets$Dediff_EMT,
  "BRAF", "RET", "NTRK1", "NTRK3", "ALK", "MET", "VEGFA", "KDR",
  "FGFR1", "FGFR2", "PIK3CA", "AKT1", "MTOR", "TERT", "TP53"
))

module_score <- function(expr_gene, genes) {
  genes <- intersect(genes, rownames(expr_gene))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_gene)))
  z <- t(scale(t(expr_gene[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

collapse_symbols <- function(expr, symbols) {
  keep <- !is.na(symbols) & symbols != "" & symbols != "NA"
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  v <- apply(expr, 1, var, na.rm = TRUE)
  ord <- order(symbols, -v)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  expr <- expr[!duplicated(symbols), , drop = FALSE]
  rownames(expr) <- symbols[!duplicated(symbols)]
  expr
}

plot_box <- function(df, x, y, filename, title = NULL) {
  df[[x]] <- factor(df[[x]])
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
    geom_jitter(width = 0.12, size = 1.1, alpha = 0.55) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(x = NULL, y = y, title = title)
  ggsave(file.path(figure_dir, filename), p, width = 5.4, height = 4.1, dpi = 300)
}

load_or_download_se <- function() {
  rds <- file.path(data_dir, "TCGA_THCA_STAR_counts_primary_tumor.rds")
  if (file.exists(rds)) return(readRDS(rds))
  q <- GDCquery(
    project = "TCGA-THCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  GDCdownload(q, directory = data_dir, method = "api", files.per.chunk = 50)
  se <- GDCprepare(q, directory = data_dir, summarizedExperiment = TRUE)
  saveRDS(se, rds)
  se
}

message("Loading TCGA-THCA RNA-seq...")
se <- load_or_download_se()
assay_names <- assayNames(se)
message("Assays: ", paste(assay_names, collapse = ", "))

assay_name <- if ("tpm_unstrand" %in% assay_names) "tpm_unstrand" else
  if ("fpkm_unstrand" %in% assay_names) "fpkm_unstrand" else assay_names[1]
expr <- assay(se, assay_name)
expr <- log2(expr + 1)
rd <- as.data.frame(rowData(se))
symbol_col <- if ("gene_name" %in% colnames(rd)) "gene_name" else
  if ("external_gene_name" %in% colnames(rd)) "external_gene_name" else
    if ("gene_id" %in% colnames(rd)) "gene_id" else stop("No gene symbol column in rowData")
expr_gene <- collapse_symbols(expr, rd[[symbol_col]])

barcode <- colnames(expr_gene)
patient <- substr(barcode, 1, 12)
sample_short <- substr(barcode, 1, 16)

scores <- data.frame(
  barcode = barcode,
  sample = sample_short,
  patient = patient,
  TDS = module_score(expr_gene, gene_sets$TDS),
  OneCarbon = module_score(expr_gene, gene_sets$OneCarbon),
  MAPK = module_score(expr_gene, gene_sets$MAPK),
  Dediff_EMT = module_score(expr_gene, gene_sets$Dediff_EMT),
  IodineHandling = module_score(expr_gene, gene_sets$IodineHandling),
  ImmuneCheckpoint = module_score(expr_gene, gene_sets$ImmuneCheckpoint),
  CAF = module_score(expr_gene, gene_sets$CAF),
  stringsAsFactors = FALSE
)
scores$RAIR_like <- as.numeric(scale(scores$OneCarbon)) +
  as.numeric(scale(scores$MAPK)) +
  as.numeric(scale(scores$Dediff_EMT)) -
  as.numeric(scale(scores$TDS))
scores$RAIR_like_group <- ifelse(scores$RAIR_like >= median(scores$RAIR_like, na.rm = TRUE),
                                 "RAIR_like_high", "RAIR_like_low")
scores$TDS_group <- ifelse(scores$TDS <= median(scores$TDS, na.rm = TRUE),
                           "TDS_low", "TDS_high")

coverage <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  data.frame(module = nm,
             n_available = length(intersect(gene_sets[[nm]], rownames(expr_gene))),
             available_genes = paste(intersect(gene_sets[[nm]], rownames(expr_gene)), collapse = ";"),
             stringsAsFactors = FALSE)
}))

message("Loading clinical data...")
clinical <- GDCquery_clinic(project = "TCGA-THCA", type = "clinical")
if (!"submitter_id" %in% colnames(clinical)) clinical$submitter_id <- clinical$bcr_patient_barcode
clinical$patient <- clinical$submitter_id
clin_keep <- intersect(c(
  "patient", "submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up",
  "age_at_diagnosis", "gender", "ajcc_pathologic_stage", "ajcc_pathologic_t",
  "ajcc_pathologic_n", "ajcc_pathologic_m", "tumor_stage", "disease_type",
  "primary_diagnosis", "morphology", "prior_malignancy", "synchronous_malignancy"
), colnames(clinical))
clinical_slim <- clinical[, clin_keep, drop = FALSE]
df <- merge(scores, clinical_slim, by = "patient", all.x = TRUE)

## Optional progression-like endpoint from clinical follow-up if available.
df$os_time <- suppressWarnings(as.numeric(ifelse(!is.na(df$days_to_death), df$days_to_death, df$days_to_last_follow_up)))
df$os_event <- ifelse(df$vital_status == "Dead", 1, 0)
df$N_positive <- ifelse(grepl("N1", df$ajcc_pathologic_n), "N1",
                        ifelse(grepl("N0", df$ajcc_pathologic_n), "N0", NA))
df$T_high <- ifelse(grepl("T3|T4", df$ajcc_pathologic_t), "T3_T4",
                    ifelse(grepl("T1|T2", df$ajcc_pathologic_t), "T1_T2", NA))
df$stage_high <- ifelse(grepl("III|IV", df$ajcc_pathologic_stage), "Stage_III_IV",
                        ifelse(grepl("I|II", df$ajcc_pathologic_stage), "Stage_I_II", NA))

mut_genes <- c("BRAF", "NRAS", "HRAS", "KRAS", "TERT", "TP53", "EIF1AX", "PIK3CA", "PTEN")
message("Loading mutation data from cBioPortal API...")
mut_file <- file.path(data_dir, "TCGA_THCA_cbioportal_mutations_core_genes.rds")
if (file.exists(mut_file)) {
  maf <- readRDS(mut_file)
} else {
  gene_entrez <- c(BRAF = 673, NRAS = 4893, HRAS = 3265, KRAS = 3845, TERT = 7015,
                   TP53 = 7157, EIF1AX = 1964, PIK3CA = 5290, PTEN = 5728)
  body <- jsonlite::toJSON(list(sampleListId = "thca_tcga_sequenced",
                                entrezGeneIds = unname(gene_entrez)),
                           auto_unbox = TRUE)
  ## Base R cannot POST with body through url(); use curl.exe if httr is absent.
  tmp_body <- tempfile(fileext = ".json")
  tmp_out <- tempfile(fileext = ".json")
  writeLines(body, tmp_body, useBytes = TRUE)
  cmd <- sprintf(
    'curl.exe -s -X POST -H "Content-Type: application/json" -H "Accept: application/json" --data-binary "@%s" -o "%s" "https://www.cbioportal.org/api/molecular-profiles/thca_tcga_mutations/mutations/fetch"',
    normalizePath(tmp_body, winslash = "/"), normalizePath(tmp_out, winslash = "/")
  )
  system(cmd)
  maf <- jsonlite::fromJSON(tmp_out)
  entrez_to_gene <- data.frame(entrezGeneId = as.integer(gene_entrez),
                               Hugo_Symbol = names(gene_entrez),
                               stringsAsFactors = FALSE)
  maf <- merge(maf, entrez_to_gene, by = "entrezGeneId", all.x = TRUE)
  saveRDS(maf, mut_file)
}
maf$patient <- substr(maf$sampleId, 1, 12)
mut_matrix <- data.frame(patient = unique(df$patient), stringsAsFactors = FALSE)
for (g in mut_genes) {
  mut_matrix[[paste0(g, "_mut")]] <- mut_matrix$patient %in% unique(maf$patient[maf$Hugo_Symbol == g])
}
mut_matrix$RAS_mut <- mut_matrix$NRAS_mut | mut_matrix$HRAS_mut | mut_matrix$KRAS_mut
mut_matrix$Any_driver <- mut_matrix$BRAF_mut | mut_matrix$RAS_mut | mut_matrix$TERT_mut |
  mut_matrix$TP53_mut | mut_matrix$EIF1AX_mut | mut_matrix$PIK3CA_mut | mut_matrix$PTEN_mut
df <- merge(df, mut_matrix, by = "patient", all.x = TRUE)

## Association tests.
score_vars <- c("TDS", "OneCarbon", "MAPK", "Dediff_EMT", "IodineHandling",
                "ImmuneCheckpoint", "CAF", "RAIR_like")
cor_tests <- do.call(rbind, lapply(score_vars, function(v) {
  ct <- suppressWarnings(cor.test(df[[v]], df$RAIR_like, method = "spearman"))
  data.frame(score = v, spearman_rho_vs_RAIR_like = unname(ct$estimate),
             p_value = ct$p.value, stringsAsFactors = FALSE)
}))

group_tests <- list()
add_wilcox <- function(var, group_col, positive, negative, expected = "higher") {
  sub <- df[df[[group_col]] %in% c(positive, negative), ]
  sub[[group_col]] <- factor(sub[[group_col]], levels = c(negative, positive))
  if (nrow(sub) < 5 || length(unique(sub[[group_col]])) < 2) return(NULL)
  data.frame(
    comparison = paste(group_col, positive, "vs", negative, sep = "_"),
    score = var,
    positive = positive,
    negative = negative,
    n_positive = sum(sub[[group_col]] == positive),
    n_negative = sum(sub[[group_col]] == negative),
    median_positive = median(sub[[var]][sub[[group_col]] == positive], na.rm = TRUE),
    median_negative = median(sub[[var]][sub[[group_col]] == negative], na.rm = TRUE),
    delta_median = median(sub[[var]][sub[[group_col]] == positive], na.rm = TRUE) -
      median(sub[[var]][sub[[group_col]] == negative], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(sub[[var]] ~ sub[[group_col]])$p.value),
    expected = expected,
    stringsAsFactors = FALSE
  )
}
for (v in score_vars) {
  group_tests[[length(group_tests) + 1]] <- add_wilcox(v, "N_positive", "N1", "N0", expected = "higher")
  group_tests[[length(group_tests) + 1]] <- add_wilcox(v, "T_high", "T3_T4", "T1_T2", expected = "higher")
  group_tests[[length(group_tests) + 1]] <- add_wilcox(v, "stage_high", "Stage_III_IV", "Stage_I_II", expected = "higher")
  group_tests[[length(group_tests) + 1]] <- add_wilcox(v, "BRAF_mut", TRUE, FALSE, expected = "higher")
  group_tests[[length(group_tests) + 1]] <- add_wilcox(v, "RAS_mut", TRUE, FALSE, expected = "lower")
  group_tests[[length(group_tests) + 1]] <- add_wilcox(v, "TERT_mut", TRUE, FALSE, expected = "higher")
}
group_tests <- do.call(rbind, group_tests)

## Logistic models for clinical aggressiveness.
logit_one <- function(outcome, predictor) {
  sub <- df[!is.na(df[[outcome]]) & !is.na(df[[predictor]]), ]
  if (nrow(sub) < 50 || length(unique(sub[[outcome]])) < 2) return(NULL)
  y <- ifelse(sub[[outcome]] %in% c("N1", "T3_T4", "Stage_III_IV"), 1, 0)
  x <- scale(sub[[predictor]])
  fit <- glm(y ~ x, family = binomial())
  ci <- suppressWarnings(confint.default(fit))
  data.frame(
    outcome = outcome,
    predictor = predictor,
    n = nrow(sub),
    events = sum(y),
    OR_per_SD = exp(coef(fit)[2]),
    CI_low = exp(ci[2, 1]),
    CI_high = exp(ci[2, 2]),
    p_value = summary(fit)$coefficients[2, 4],
    stringsAsFactors = FALSE
  )
}
logit_results <- do.call(rbind, lapply(c("N_positive", "T_high", "stage_high"), function(out) {
  do.call(rbind, lapply(score_vars, function(v) logit_one(out, v)))
}))

## Survival is expected to be underpowered in THCA; report but do not overuse.
surv_results <- data.frame()
if (sum(!is.na(df$os_time) & df$os_time > 0 & !is.na(df$os_event)) > 50 &&
    sum(df$os_event == 1, na.rm = TRUE) >= 5) {
  for (v in c("TDS", "RAIR_like", "OneCarbon")) {
    sub <- df[!is.na(df[[v]]) & !is.na(df$os_time) & df$os_time > 0 & !is.na(df$os_event), ]
    fit <- coxph(Surv(os_time, os_event) ~ scale(sub[[v]]), data = sub)
    s <- summary(fit)
    surv_results <- rbind(surv_results, data.frame(
      endpoint = "OS",
      score = v,
      n = nrow(sub),
      events = sum(sub$os_event),
      HR_per_SD = s$coefficients[1, "exp(coef)"],
      CI_low = s$conf.int[1, "lower .95"],
      CI_high = s$conf.int[1, "upper .95"],
      p_value = s$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    ))
  }
}

## Target prioritization: differential abundance in RAIR-like high vs low,
## anti-correlation with TDS, and consistency with one-carbon/dediff signal.
target_available <- intersect(target_genes, rownames(expr_gene))
target_stats <- do.call(rbind, lapply(target_available, function(g) {
  e <- as.numeric(expr_gene[g, match(df$barcode, colnames(expr_gene))])
  high <- df$RAIR_like_group == "RAIR_like_high"
  ct <- suppressWarnings(cor.test(e, df$RAIR_like, method = "spearman"))
  ct_tds <- suppressWarnings(cor.test(e, df$TDS, method = "spearman"))
  data.frame(
    gene = g,
    median_high = median(e[high], na.rm = TRUE),
    median_low = median(e[!high], na.rm = TRUE),
    delta_high_low = median(e[high], na.rm = TRUE) - median(e[!high], na.rm = TRUE),
    p_wilcox = suppressWarnings(wilcox.test(e ~ df$RAIR_like_group)$p.value),
    rho_RAIR_like = unname(ct$estimate),
    p_rho_RAIR_like = ct$p.value,
    rho_TDS = unname(ct_tds$estimate),
    p_rho_TDS = ct_tds$p.value,
    stringsAsFactors = FALSE
  )
}))
target_stats$priority_score <- rank(-target_stats$delta_high_low, ties.method = "average") +
  rank(-target_stats$rho_RAIR_like, ties.method = "average") +
  rank(target_stats$rho_TDS, ties.method = "average")
target_stats <- target_stats[order(target_stats$priority_score), ]

write.csv(df, file.path(result_dir, "TCGA_THCA_scores_clinical_mutation.csv"), row.names = FALSE)
write.csv(coverage, file.path(result_dir, "TCGA_THCA_gene_set_coverage.csv"), row.names = FALSE)
write.csv(cor_tests, file.path(result_dir, "TCGA_THCA_score_correlations.csv"), row.names = FALSE)
write.csv(group_tests, file.path(result_dir, "TCGA_THCA_group_tests.csv"), row.names = FALSE)
write.csv(logit_results, file.path(result_dir, "TCGA_THCA_logistic_clinical_associations.csv"), row.names = FALSE)
write.csv(surv_results, file.path(result_dir, "TCGA_THCA_survival_exploratory.csv"), row.names = FALSE)
write.csv(target_stats, file.path(result_dir, "TCGA_THCA_target_priority_table.csv"), row.names = FALSE)

plot_box(df[!is.na(df$N_positive), ], "N_positive", "RAIR_like",
         "TCGA_RAIR_like_by_N_stage.png", "TCGA-THCA RAIR-like by nodal status")
plot_box(df[!is.na(df$T_high), ], "T_high", "RAIR_like",
         "TCGA_RAIR_like_by_T_stage.png", "TCGA-THCA RAIR-like by T category")
plot_box(df[!is.na(df$stage_high), ], "stage_high", "RAIR_like",
         "TCGA_RAIR_like_by_pathologic_stage.png", "TCGA-THCA RAIR-like by stage")
plot_box(df, "BRAF_mut", "RAIR_like", "TCGA_RAIR_like_by_BRAF.png", "TCGA-THCA RAIR-like by BRAF mutation")
plot_box(df, "RAS_mut", "RAIR_like", "TCGA_RAIR_like_by_RAS.png", "TCGA-THCA RAIR-like by RAS mutation")

p_scatter <- ggplot(df, aes(x = TDS, y = OneCarbon, color = RAIR_like_group)) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.4) +
  theme_bw(base_size = 11) +
  labs(title = "TCGA-THCA TDS vs one-carbon metabolism", x = "TDS", y = "One-carbon score", color = NULL)
ggsave(file.path(figure_dir, "TCGA_TDS_vs_OneCarbon.png"), p_scatter, width = 5.6, height = 4.3, dpi = 300)

top_targets <- head(target_stats, 15)
p_target <- ggplot(top_targets, aes(x = reorder(gene, delta_high_low), y = delta_high_low, fill = rho_RAIR_like)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = "Candidate target genes enriched in RAIR-like-high TCGA tumors",
       x = NULL, y = "Median expression difference: high - low")
ggsave(file.path(figure_dir, "TCGA_top_target_priority.png"), p_target, width = 6.2, height = 4.8, dpi = 300)

sink(file.path(result_dir, "TCGA_THCA_discovery_summary.txt"))
cat("TCGA-THCA discovery analysis\n\n")
cat("Expression assay:", assay_name, "\n")
cat("Samples:", nrow(df), "primary tumors\n\n")
cat("Gene-set coverage:\n")
print(coverage)
cat("\nScore correlations with RAIR-like:\n")
print(cor_tests)
cat("\nClinical/mutation group tests, selected RAIR-like rows:\n")
print(group_tests[group_tests$score == "RAIR_like", ])
cat("\nLogistic clinical associations, selected RAIR-like/TDS/OneCarbon rows:\n")
print(logit_results[logit_results$predictor %in% c("RAIR_like", "TDS", "OneCarbon"), ])
cat("\nExploratory survival:\n")
print(surv_results)
cat("\nTop target priority table:\n")
print(head(target_stats, 20))
sink()

message("TCGA discovery analysis complete: ", result_dir)
