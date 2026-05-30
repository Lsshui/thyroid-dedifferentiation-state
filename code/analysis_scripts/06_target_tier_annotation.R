options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath("analysis", winslash = "/", mustWork = TRUE)
result_dir <- file.path(root, "results")
figure_dir <- file.path(root, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

priority <- fread(file.path(result_dir, "TCGA_THCA_target_priority_table.csv"))
external <- fread(file.path(result_dir, "external_target_validation_long.csv"))

external_summary <- external[, .(
  n_external_up_validations = sum(delta_positive_negative > 0 & p_adj < 0.05, na.rm = TRUE),
  best_external_padj = min(p_adj, na.rm = TRUE),
  GSE33630_delta = delta_positive_negative[dataset == "GSE33630"][1],
  GSE33630_padj = p_adj[dataset == "GSE33630"][1],
  GSE76039_delta = delta_positive_negative[dataset == "GSE76039"][1],
  GSE76039_padj = p_adj[dataset == "GSE76039"][1]
), by = gene]

annotation <- data.table(
  gene = c("MTHFD2", "SLC1A5", "SLC7A5", "AXL", "FN1", "FOSL1", "ITGB1",
           "MTHFD1L", "SHMT2", "DUSP5", "KRT19", "LAMC2",
           "MET", "DUSP6", "DUSP4", "HMGA2", "CCND1", "ETV4", "PIK3CA", "ALK", "ETV5"),
  evidence_tier = c(rep("Tier 1A_core_actionable_axis", 4),
                    rep("Tier 1B_ecology_state_axis", 3),
                    rep("Tier 2_supportive_marker_axis", 5),
                    rep("Tier 3_deprioritize_as_primary_marker", 9)),
  biological_axis = c("mitochondrial one-carbon", "amino-acid uptake", "amino-acid uptake", "RTK/EMT survival",
                      "ECM/CAF-epithelial interface", "AP-1/MAPK transcription", "ECM adhesion",
                      "mitochondrial one-carbon", "mitochondrial one-carbon", "MAPK feedback",
                      "dedifferentiated epithelial", "basement membrane/EMT",
                      "RTK", "MAPK feedback", "MAPK feedback", "EMT/chromatin",
                      "cell cycle/MAPK", "ETS/MAPK", "PI3K", "RTK", "ETS/MAPK"),
  actionability_note = c(
    "enzyme vulnerability; strongest one-carbon gene across TCGA and GEO",
    "transporter vulnerability; robust across TCGA and GEO",
    "transporter vulnerability; robust across TCGA and GEO",
    "clinically druggable kinase class; use as actionable ecology-linked node",
    "robust ECM/state anchor; better as ecology marker than direct tumor-cell target",
    "transcriptional MAPK/AP-1 readout; marker or upstream-pathway rationale",
    "adhesion/ECM interface; ecology marker and combination-rationale node",
    "one-carbon support gene; validated mainly in PDTC-to-ATC contrast",
    "one-carbon support gene; TCGA signal, external GEO borderline after FDR",
    "MAPK feedback readout; validated in PTC-to-ATC but not PDTC-to-ATC",
    "epithelial dedifferentiation marker; one external validation",
    "EMT/basement-membrane marker; one external validation",
    "TCGA-correlated but reverses/fails in external ATC validation",
    "TCGA-correlated but not externally stable",
    "TCGA-correlated but not externally stable",
    "TCGA-correlated but not externally stable",
    "TCGA-correlated but not externally stable",
    "TCGA-correlated but not externally stable",
    "TCGA-correlated but not externally stable as expression marker",
    "TCGA-correlated but not externally stable as expression marker",
    "TCGA-correlated but not externally stable"
  )
)

tier_table <- merge(annotation, priority, by = "gene", all.x = TRUE)
tier_table <- merge(tier_table, external_summary, by = "gene", all.x = TRUE)
tier_order <- c("Tier 1A_core_actionable_axis", "Tier 1B_ecology_state_axis",
                "Tier 2_supportive_marker_axis", "Tier 3_deprioritize_as_primary_marker")
tier_table$evidence_tier <- factor(tier_table$evidence_tier, levels = tier_order)
tier_table <- tier_table[order(evidence_tier, -n_external_up_validations, -rho_RAIR_like)]
fwrite(tier_table, file.path(result_dir, "candidate_target_tiers.csv"))

plot_table <- tier_table[!is.na(rho_RAIR_like)]
plot_table$validated <- ifelse(plot_table$n_external_up_validations >= 2, "2 GEO cohorts",
                               ifelse(plot_table$n_external_up_validations == 1, "1 GEO cohort", "Not stable"))
plot_table$label <- ifelse(plot_table$evidence_tier %in%
                             c("Tier 1A_core_actionable_axis", "Tier 1B_ecology_state_axis"),
                           plot_table$gene, "")

p <- ggplot(plot_table, aes(x = rho_RAIR_like, y = delta_high_low,
                            color = validated, shape = evidence_tier)) +
  geom_hline(yintercept = 0, color = "grey78", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey78", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(aes(label = label), nudge_y = 0.18, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = c("2 GEO cohorts" = "#b73e3e",
                                "1 GEO cohort" = "#c58a2b",
                                "Not stable" = "#6f7a85")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right") +
  labs(x = "Spearman rho with RAIR-like score in TCGA-THCA",
       y = "Median expression delta: RAIR-like high vs low",
       color = "External validation",
       shape = "Target tier",
       title = "Prioritized targetable and ecology-linked nodes")
ggsave(file.path(figure_dir, "candidate_target_tier_evidence_map.png"),
       p, width = 7.0, height = 5.0, dpi = 300)

sink(file.path(result_dir, "candidate_target_tier_summary.txt"))
cat("Candidate target tiering for the one-carbon RAIR-like state\n\n")
print(tier_table[, .(gene, evidence_tier, biological_axis, n_external_up_validations,
                     delta_high_low, rho_RAIR_like, rho_TDS,
                     GSE33630_delta, GSE33630_padj, GSE76039_delta, GSE76039_padj,
                     actionability_note)], row.names = FALSE)
sink()

message("Saved target tiering outputs.")
