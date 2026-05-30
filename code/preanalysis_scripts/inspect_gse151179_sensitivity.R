df <- read.csv("D:/OneDrive/桌面/文章撰写/01_2026/5月/JTM选题/JTM1/preanalysis/results/GSE151179_scores_metadata.csv")
cat("RAI response by disease\n")
print(table(df$rai_response, df$disease, useNA = "ifany"))
cat("\nRAI response by metastatic uptake\n")
print(table(df$rai_response, df$rai_uptake, useNA = "ifany"))
cat("\nRAI response by tissue type\n")
print(table(df$rai_response, df$tissue_type, useNA = "ifany"))

tumor <- df[df$source == "papillary thyroid carcinoma", ]
cat("\nTumor samples\n")
print(tumor[, c("sample", "title", "tissue_type", "disease", "rai_response",
                "rai_uptake", "TDS", "RAIR_like", "OneCarbon", "MAPK", "Dediff_EMT")])

primary <- tumor[grepl("Primary", tumor$tissue_type), ]
cat("\nPrimary tumor only\n")
print(primary[, c("sample", "title", "tissue_type", "disease", "rai_response", "TDS", "RAIR_like")])
cat("\nPrimary-only Wilcoxon tests\n")
print(wilcox.test(TDS ~ rai_response, data = primary))
print(wilcox.test(RAIR_like ~ rai_response, data = primary))
