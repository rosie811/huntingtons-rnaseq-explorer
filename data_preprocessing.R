install.packages("BiocManager")
BiocManager::install("GEOquery")

library(readr)
library(GEOquery)

# Load count and DE files into R
counts <- read_tsv("raw_data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt")
de <- read_tsv("raw_data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt")

# Title first column of DE
colnames(de)[1] <- "gene_id"

# Load metadata using GEOquery
gse <- getGEO("GSE64810", GSEMatrix = TRUE)
meta <- pData(gse[[1]])

# Inspect each dataset
dim(counts)
head(counts)

dim(de)
head(de)

head(meta)
dim(meta)

# Create condition column
meta$condition <- ifelse(
  grepl("Huntington", meta$`diagnosis:ch1`, ignore.case = TRUE),
  "HD",
  "Control"
)

table(meta$condition)

# Create clean sample metadata table
meta_clean <- meta[, c("title", "geo_accession", "condition")]
colnames(meta_clean) <- c("sample_id", "gsm_id", "condition")

# Match metadata to counts
all(colnames(counts)[-1] %in% meta_clean$sample_id)
# > should return TRUE

# Fix counts format
colnames(counts)[1] <- "gene_id"

counts_mat <- as.data.frame(counts)
rownames(counts_mat) <- counts_mat$gene_id
counts_mat <- counts_mat[, -1]

# SANITY CHECKS BEFORE SAVING
# Confirm sample order matches
all(colnames(counts_mat) == meta_clean$sample_id)
# Condition check
table(meta_clean$condition)

# Save files
write.csv(meta_clean, "clean_data/sample_info.csv", row.names = FALSE)
write.csv(counts_mat, "clean_data/normalized_counts.csv")
write.csv(de, "clean_data/deseq_results.csv", row.names = FALSE)