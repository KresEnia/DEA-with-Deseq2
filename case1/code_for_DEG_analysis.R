# Load necessary libraries
library(tidyverse)
library(readxl)

# Preprocess step for summing technical replicates in counts file
# Load counts
counts <- read_csv("GSE283950_Raw_Read_Counts_Per_Gene-0136.csv")

# Select the columns that contain count data
# These are all columns except 'Gene_ID' and 'Gene_Name'.
count_cols_only <- colnames(counts)[!colnames(counts) %in% c("Gene_ID", "Gene_Name")]
count_cols_only

# Ensure count columns are numeric
counts_numeric_only <- counts %>%
  mutate(across(all_of(count_cols_only), as.numeric)) %>%
  # Select only the numeric count columns to be used in rowSums
  dplyr::select(all_of(count_cols_only)) %>%
  as.matrix() 

# Get Gene_ID as a separate vector to use as row names later
gene_ids <- counts$Gene_ID
rownames(counts_numeric_only) <- gene_ids

# Extract biological replicate names (e.g., Control_1 from Control_1-1)
bio_names <- gsub("-[12]_Raw_Read_Counts", "", count_cols_only)

# Collapse technical replicates by summing them
collapsed_counts <- sapply(unique(bio_names), function(biorep) {
  replicate_cols <- which(bio_names == biorep)
  rowSums(counts_numeric_only[, replicate_cols, drop = FALSE], na.rm = TRUE) # Corrected variable here
})

# Create a dataframe with new counts
collapsed_counts_data <- as.data.frame(collapsed_counts)

# Load metadata
meta <- read_excel("meta.xlsx")
metadata <- as.data.frame(meta) #convert to dataframe
rownames(metadata) <- metadata[[1]] #convert first column to rownames
metadata <- metadata[, -1] #remove first column used for rownames 

# Quality check
# making sure the row names matches column names
all(colnames(collapsed_counts_data) %in% rownames(metadata))
# checking order of row names and column names
all(colnames(collapsed_counts_data) == rownames(metadata))

library(DESeq2)
# Creating design for DEA
dds <- DESeqDataSetFromMatrix(countData = collapsed_counts_data, 
                              colData = metadata,
                              design = ~Group)
dds

# Removing low counts reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Setting reference for deg analysis
dds$Group <- relevel(dds$Group, ref = "Control")
deg <- DESeq(dds)

# Results checking
resultsNames(deg)

# Install new library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(deg, coef="Group_TMEM11KO_vs_Control", type="apeglm")
resLFC

res <- results(deg)
res
summary(res)

# Getting DEG at different alpha value
res0.05 <- results(deg, alpha = 0.05) 
summary(res0.05)

# PCA plot
vsd <- vst(deg, blind = FALSE)
plotPCA(vsd, intgroup = "Group")

#ggsave("pca_plot_case_1.png", dpi = 600, width = 12, height = 6, units = "in")
