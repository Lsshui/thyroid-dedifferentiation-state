options(stringsAsFactors = FALSE, timeout = 100000)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(limma)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

root <- normalizePath("analysis", winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data", "tcga")
result_dir <- file.path(root, "results")
figure_dir <- file.path(root, "figures")
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

flatten_fgsea <- function(x) {
  x <- as.data.frame(x)
  if ("leadingEdge" %in% colnames(x)) {
    x$leadingEdge <- vapply(x$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  }
  x[order(x$padj, -abs(x$NES)), ]
}

message("Loading cached TCGA-THCA expression...")
se <- readRDS(file.path(data_dir, "TCGA_THCA_STAR_counts_primary_tumor.rds"))
assay_name <- if ("tpm_unstrand" %in% assayNames(se)) "tpm_unstrand" else
  if ("fpkm_unstrand" %in% assayNames(se)) "fpkm_unstrand" else assayNames(se)[1]
expr <- log2(assay(se, assay_name) + 1)
rd <- as.data.frame(rowData(se))
symbol_col <- if ("gene_name" %in% colnames(rd)) "gene_name" else
  if ("external_gene_name" %in% colnames(rd)) "external_gene_name" else "gene_id"
expr_gene <- collapse_symbols(expr, rd[[symbol_col]])

score_file <- file.path(result_dir, "TCGA_THCA_scores_clinical_mutation.csv")
scores <- fread(score_file)
common <- intersect(colnames(expr_gene), scores$barcode)
expr_gene <- expr_gene[, common, drop = FALSE]
scores <- as.data.frame(scores[match(common, scores$barcode), ])

keep_gene <- rowMeans(expr_gene > 0.5, na.rm = TRUE) >= 0.10 & apply(expr_gene, 1, var, na.rm = TRUE) > 0
expr_gene <- expr_gene[keep_gene, , drop = FALSE]

q <- quantile(scores$RAIR_like, probs = c(0.25, 0.75), na.rm = TRUE)
scores$RAIR_extreme <- ifelse(scores$RAIR_like <= q[[1]], "Q1_low",
                              ifelse(scores$RAIR_like >= q[[2]], "Q4_high", NA))
extreme <- !is.na(scores$RAIR_extreme)
group <- factor(scores$RAIR_extreme[extreme], levels = c("Q1_low", "Q4_high"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_gene[, extreme, drop = FALSE], design)
cont <- makeContrasts(RAIR_Q4_vs_Q1 = Q4_high - Q1_low, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont), trend = TRUE)
de <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
de$gene <- rownames(de)
de <- de[, c("gene", setdiff(colnames(de), "gene"))]
fwrite(de, file.path(result_dir, "TCGA_RAIR_high_low_limma.csv"))

message("Computing genome-wide Spearman correlations...")
rair <- scores$RAIR_like
rho <- suppressWarnings(cor(t(expr_gene), rair, method = "spearman", use = "pairwise.complete.obs"))
n <- sum(!is.na(rair))
tstat <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
pval <- 2 * pt(-abs(tstat), df = n - 2)
cor_df <- data.frame(
  gene = rownames(expr_gene),
  spearman_rho_RAIR_like = as.numeric(rho),
  p_value = as.numeric(pval),
  p_adj = p.adjust(as.numeric(pval), method = "BH"),
  stringsAsFactors = FALSE
)
cor_df <- cor_df[order(cor_df$p_adj, -abs(cor_df$spearman_rho_RAIR_like)), ]
fwrite(cor_df, file.path(result_dir, "TCGA_RAIR_gene_correlations.csv"))

signature_genes <- intersect(unique(unlist(gene_sets)), de$gene)
ranks <- de$t
names(ranks) <- de$gene
ranks <- ranks[is.finite(ranks) & !names(ranks) %in% signature_genes]
ranks <- sort(ranks, decreasing = TRUE)

message("Running MSigDB Hallmark GSEA...")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
hallmark_fgsea <- fgsea(pathways = hallmark_sets, stats = ranks, minSize = 10,
                        maxSize = 500, nPermSimple = 10000)
hallmark_fgsea <- flatten_fgsea(hallmark_fgsea)
fwrite(hallmark_fgsea, file.path(result_dir, "TCGA_RAIR_hallmark_fgsea.csv"))

message("Running MSigDB Reactome GSEA...")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome_sets <- split(reactome$gene_symbol, reactome$gs_name)
reactome_fgsea <- fgsea(pathways = reactome_sets, stats = ranks, minSize = 10,
                        maxSize = 500, nPermSimple = 10000)
reactome_fgsea <- flatten_fgsea(reactome_fgsea)
fwrite(reactome_fgsea, file.path(result_dir, "TCGA_RAIR_reactome_fgsea.csv"))

plot_fgsea <- function(tab, filename, title, n_each = 10) {
  sig <- tab[!is.na(tab$padj) & tab$padj < 0.05, , drop = FALSE]
  pos <- head(sig[order(-sig$NES), ], n_each)
  neg <- head(sig[order(sig$NES), ], n_each)
  plt <- rbind(pos, neg)
  if (nrow(plt) == 0) return(invisible(NULL))
  plt$pathway_clean <- gsub("^HALLMARK_|^REACTOME_", "", plt$pathway)
  plt$pathway_clean <- gsub("_", " ", plt$pathway_clean)
  plt$pathway_clean <- factor(plt$pathway_clean, levels = plt$pathway_clean[order(plt$NES)])
  p <- ggplot(plt, aes(x = pathway_clean, y = NES, fill = NES > 0)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#b73e3e", "FALSE" = "#3274a1")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Normalized enrichment score", title = title)
  ggsave(file.path(figure_dir, filename), p, width = 7.2, height = 5.4, dpi = 300)
}

plot_fgsea(hallmark_fgsea, "TCGA_RAIR_hallmark_fgsea.png",
           "Hallmark enrichment in RAIR-like high tumors")
plot_fgsea(reactome_fgsea, "TCGA_RAIR_reactome_fgsea.png",
           "Reactome enrichment in RAIR-like high tumors")

candidate_genes <- intersect(
  c("MTHFD2", "MTHFD1L", "SLC1A5", "SLC7A5", "SHMT2", "AXL", "FN1",
    "FOSL1", "ITGB1", "LAMC2", "KRT19", "DUSP5", "DUSP6", "MET"),
  rownames(expr_gene)
)
scores$RAIR_quartile <- cut(scores$RAIR_like,
                            breaks = quantile(scores$RAIR_like, probs = seq(0, 1, 0.25), na.rm = TRUE),
                            include.lowest = TRUE, labels = paste0("Q", 1:4))
gene_group <- do.call(rbind, lapply(candidate_genes, function(g) {
  data.frame(gene = g, quartile = levels(scores$RAIR_quartile),
             median_log2_tpm = as.numeric(tapply(expr_gene[g, ], scores$RAIR_quartile, median, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))
gene_group$z_median <- ave(gene_group$median_log2_tpm, gene_group$gene, FUN = function(x) as.numeric(scale(x)))
fwrite(gene_group, file.path(result_dir, "TCGA_candidate_gene_expression_by_RAIR_quartile.csv"))

p_heat <- ggplot(gene_group, aes(x = quartile, y = gene, fill = z_median)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#3274a1", mid = "white", high = "#b73e3e", midpoint = 0) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank()) +
  labs(x = "RAIR-like quartile", y = NULL, fill = "Row z", title = "Candidate target expression across RAIR-like quartiles")
ggsave(file.path(figure_dir, "TCGA_candidate_gene_expression_by_RAIR_quartile.png"),
       p_heat, width = 5.0, height = 4.6, dpi = 300)

top_hallmark <- head(hallmark_fgsea[order(hallmark_fgsea$padj, -abs(hallmark_fgsea$NES)), ], 15)
top_reactome <- head(reactome_fgsea[order(reactome_fgsea$padj, -abs(reactome_fgsea$NES)), ], 15)
summary_file <- file.path(result_dir, "TCGA_RAIR_pathway_summary.txt")
sink(summary_file)
cat("TCGA-THCA RAIR-like high-vs-low pathway analysis\n\n")
cat("Expression assay:", assay_name, "\n")
cat("Genes retained after expression/variance filter:", nrow(expr_gene), "\n")
cat("RAIR-like extreme comparison: Q4 high n =", sum(group == "Q4_high"),
    "; Q1 low n =", sum(group == "Q1_low"), "\n")
cat("Signature genes excluded from GSEA rank:", length(signature_genes), "\n\n")
cat("Top Hallmark pathways:\n")
print(top_hallmark[, c("pathway", "NES", "padj", "size")], row.names = FALSE)
cat("\nTop Reactome pathways:\n")
print(top_reactome[, c("pathway", "NES", "padj", "size")], row.names = FALSE)
cat("\nTop RAIR-like correlated genes:\n")
print(head(cor_df, 20), row.names = FALSE)
cat("\nTop RAIR-like high-vs-low up genes:\n")
print(head(de[order(-de$logFC), c("gene", "logFC", "P.Value", "adj.P.Val")], 20), row.names = FALSE)
sink()

message("Saved pathway analysis to ", summary_file)
