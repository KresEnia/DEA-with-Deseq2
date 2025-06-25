## 1.INSTALLATION
install.packages("ggpubr")
# in case of necessity to install Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# example of installation of specific package
BiocManager::install(c("GenomicFeatures"))
# installation of DESeq2                     
BiocManager::install("DESeq2")
# check all options
BiocManager::available()

## 2.UPLOADING count matrix and metadata (working directory is already set)
cnt <- read.csv("counts.csv")
str(cnt)
met <- read.csv("metadata2.csv", row.names = 1)
str(met)

## 3.QUALITY check
# making sure the row names matches column names
all(colnames(cnt) %in% rownames(met))
# checking order of row names and column names
all(colnames(cnt) == rownames(met))

## 4.DEA
library(DESeq2)
# creating design for DEA
dds <- DESeqDataSetFromMatrix(countData = cnt, 
                              colData = met,
                              design = ~dexamethasone)
dds
# removing low counts reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# setting reference for deg analysis
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
deg <- DESeq(dds)
# getting results
res <- results(deg)
write.csv(res, "test_udemy.csv")
# summary stat of res
summary(res)
# getting DEG at different alpha value
res0.05 <- results(deg, alpha = 0.05) #or 0.01 (stand-0.1)
summary(res0.05)
# converting IDs to gene names
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
res0.05.df <- as.data.frame(res0.05) #transforming results to the dataframe
str(res0.05.df)
res0.05.df$Symbol <- mapIds(org.Hs.eg.db, rownames(res0.05.df), keytype = "ENSEMBL", column = "SYMBOL")
res0.05.df
str(res0.05.df)
write.csv(res0.05.df, "final_test_udemy.csv")

## 5.QUALITY check of rna seq data (PCA plot, dispersion plot, inspecting size factor)
# PCA plot
vsd <- vst(deg, blind = FALSE)
plotPCA(vsd, intgroup = "dexamethasone")

library(ggplot2)
#ggsave("pca_plot_for_gw.png", dpi = 600, width = 12, height = 6, units = "in")

# size factor estimation
sizeFactors(deg)
# estimating the dispersion
plotDispEsts(deg)

## 6. Analysis of Gene expression
# building MA plot
plotMA(res0.05)

# Find best genes
library(dplyr)
best_genes <- res0.05.df %>%
  arrange(padj) %>%
  head(30)
best_genes
write.csv(best_genes, "best_genes.csv")
# volcano plot
vol <- res0.05.df %>%
  filter(!is.na(padj))

ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05 & abs(log2FoldChange)>1)) + geom_point() +
  geom_text(data = best_genes, aes(label = Symbol), hjust = -0.2, vjust = 0.5)

#ggsave("plot_volcano.png", dpi = 600, width = 12, height = 6, units = "in")

# building heatmap in 3 steps (find top genes, normalise deg data, get z values)
BiocManager::install("ComplexHeatmap")
library (ComplexHeatmap)
top_genes <- res0.05.df %>% #find top 30 genes
  arrange(padj) %>%
  head(30)
mat <- counts(deg, normalized=T)[rownames(top_genes),]
head(mat, 5)
mat.z <- t(apply(mat, 1, scale))
head(mat.z, 5)
colnames(mat.z) <- rownames(met)
head(mat.z, 5)
Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), row_labels = top_genes$Symbol)


