scores <- readRDS("D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/preanalysis/results/single_cell_per_cell_scores.rds")
out <- do.call(rbind, lapply(split(scores, scores$dataset), function(df) {
  cutoff <- quantile(df$RAIR_like, 0.9, na.rm = TRUE)
  top <- df[df$RAIR_like >= cutoff, ]
  tab <- as.data.frame(table(top$dataset, top$sample, top$tissue_group, top$cell_type),
                       stringsAsFactors = FALSE)
  names(tab) <- c("dataset", "sample", "tissue_group", "cell_type", "n_top_cells")
  tab <- tab[tab$n_top_cells > 0, ]
  tab$fraction_dataset_top10 <- tab$n_top_cells / sum(tab$n_top_cells)
  tab
}))
write.csv(out, "D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/preanalysis/results/single_cell_top_RAIR_like_by_sample.csv",
          row.names = FALSE)

cat("Top RAIR-like cells by dataset/sample, all cell types:\n")
print(out[order(out$dataset, -out$fraction_dataset_top10), ], row.names = FALSE)
cat("\nTop RAIR-like epithelial/thyroid cells by dataset/sample:\n")
epi <- out[out$cell_type == "Epithelial_Thyroid", ]
print(epi[order(epi$dataset, -epi$fraction_dataset_top10), ], row.names = FALSE)
