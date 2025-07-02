library(DESeq2)
library(ggplot2)
library(rtracklayer)
library(dplyr)

# Load counts
count <- read.csv("raw_counts_NCBI.tsv", sep = '\t', row.names = 1)

# Load metadata
meta <- read.csv("metadata.csv", sep = ';', row.names = 1)

# Check row- and col- names matching
all(colnames(count) %in% rownames(meta))
# Check order of colnames and rownames
all(colnames(count) == rownames(meta))

# DEA
dds = DESeqDataSetFromMatrix(countData = count,
                             colData = meta,
                             design = ~Sample_treatment)
dds

# Low counts removal
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Setting reference for deg
dds$Sample_treatment <- relevel(dds$Sample_treatment, ref = "Control ethanol 24h")
dds <- DESeq(dds)

# Getting results
res <- results(dds)
summary(res)

# Filtering based on padj and log2foldchange
res.df <- as.data.frame(res)
filtered_df <- res.df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

#write.csv(filtered_df, "deg_values_filtered.csv")

# Convertion of gene_ids to gene names
annotation <- read.csv("Human_GRCh38_p13_annot.tsv", sep = '\t', row.names = 1)
ann <- as.data.frame(annotation)

# Extract relevant columns from annotation file
gene_mapping <- ann %>%
  select(GeneType, Symbol, Description)

# Combine using rownames as keys
filtered_with_ann <- cbind(filtered_df, gene_mapping[rownames(filtered_df), ])

#write.csv(filtered_with_ann, "filtered_deg_with_ann.csv")

# Extraction of normalised by deseq2 counts to use them to identify top 500 most variable genes
# obtain normalized counts
norm_counts <- counts(dds, normalized = TRUE)
# calculate the variance for each gene across samples
variances <- apply(norm_counts, 1, var)

# Select the top 500 most variable genes
top_500 = order(variances, decreasing = TRUE)[1:500]
top_500

# Subset the normalized counts for the top 500 most variable genes
top_500_norm_counts <- norm_counts[top_500, ]
top_500_norm_counts

# Perform PCA on the top 500 variable genes (using normalized counts)
pca_res <- prcomp(t(top_500_norm_counts), scale. = TRUE)

# Compute percentage of variance explained by the first two principal components
percentVar <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

# Convert PCA results into a dataframe
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], 
                     Group = meta$Sample_treatment)  # Add metadata Group column

# Plot PCA with percentage variance in axis labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of Top 500 Most Variable Genes (Normalized Counts)", 
       x = paste0("PC1 - ", round(percentVar[1], 1), "%"), 
       y = paste0("PC2 - ", round(percentVar[2], 1), "%")) +
  coord_fixed(ratio = sqrt(percentVar[2] / percentVar[1])) +
  theme_minimal()

#ggsave("var_pca.png", dpi = 600, width = 12, height = 6, units = "in")

# Analysis of Gene expression
# Find best genes
best_genes <- filtered_with_ann %>%
  arrange(padj) %>%
  head(50)
best_genes

#write.csv(best_genes, "best_genes.csv")

# Volcano plot for unfiltered data
vol <- res.df %>%
  filter(!is.na(padj))

ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05 & abs(log2FoldChange)>1)) + geom_point() +
  geom_text(data = best_genes, aes(label = Symbol), hjust = -0.2, vjust = 0.5)

#ggsave("volcano.png", dpi = 600, width = 12, height = 6, units = "in")

# Building heatmap in 3 steps (find top genes, normalize dds data, get z values)
library (ComplexHeatmap)
top_genes <- filtered_with_ann %>% #find top 50 genes
arrange(padj) %>%
  head(50)
mat <- counts(dds, normalized=T)[rownames(top_genes),]
head(mat, 5)
mat.z <- t(apply(mat, 1, scale))
head(mat.z, 5)
colnames(mat.z) <- meta$Sample_treatment # or rownames(meta)
head(mat.z, 5)
Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), row_labels = top_genes$Symbol)

# GO enrichment analysis and KEGG
res_with_symb <- cbind(res.df, gene_mapping[rownames(res.df), ]) #add symbols to universe 
detected_genes <- filtered_with_ann$Symbol #set list of detected genes
background_genes <- res_with_symb$Symbol #set background/universe

library(clusterProfiler)
library(org.Hs.eg.db) 
packageVersion("org.Hs.eg.db")
library(enrichplot)

#Convert SYMBOLs (from detected and background) to ENTREZID
# DE genes from SYMBOLs
de_entrez <- bitr(detected_genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# Convert background SYMBOLs to ENTREZID
bg_entrez <- bitr(background_genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# Analysis
# GO
ego <- enrichGO(
  gene = de_entrez$ENTREZID,       
  universe = bg_entrez$ENTREZID, 
  OrgDb = org.Hs.eg.db,        
  keyType = "ENTREZID",         
  ont = "BP",                  # Biological Process (BP)
  pAdjustMethod = "BH",        # Benjamin-Hochberg correction
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Visualisation
barplot(ego, showCategory = 20) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
# Or simply
barplot(ego, showCategory = 20)  # Show top 20 GO terms
dotplot(ego, showCategory = 20)

# Check the full list of go terms
go_results <- as.data.frame(ego)
#write.csv(go_results, "go_terms_enriched.csv")

# KEGG
kegg <- enrichKEGG(
  gene = de_entrez$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez$ENTREZID,
  qvalueCutoff = 0.2
)

kegg_results <- as.data.frame(kegg)
barplot(kegg, showCategory = 6)

# Semantic similarity of GO terms
library(rrvgo)

# Prepare GO terms and scores
terms <- go_results$ID
scores <- setNames(-log10(go_results$p.adjust), go_results$ID)

# Calculate similarity matrix
simMatrix <- calculateSimMatrix(terms, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")

# Reduce similar GO terms (you can adjust the threshold!)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")

# View reduced table
head(reducedTerms)
nrow(reducedTerms)

# Save treemaplot
png("semantic_similarity_Gene_treemap.png", width = 17, height = 12, units = "in", res = 600)
treemapPlot(reducedTerms, size = "score", title = "Semantic similarity of GO terms")
dev.off()

# Scatterplot
scatterPlot(
  simMatrix,
  reducedTerms,
  algorithm = c("pca"),
  onlyParents = FALSE,
  size = "score",
  addLabel = TRUE,
  labelSize = 3
)

