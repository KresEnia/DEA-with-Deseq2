# Load library
library(DESeq2)
# Read the files
counts <- read.csv("small.csv", header = TRUE, row.names = 1)
colnames(counts)
meta <- read.csv("metadata_sm_created.csv", header = TRUE, row.names = 1)
# Some checking
str(counts)
str(meta)
# Rows and col names matching
all(colnames(counts) %in% rownames(meta))
all(colnames(counts) == rownames(meta))

# Create a new column in metadata with combination of 2 factors
meta$Group <- factor(paste0(meta$Genotype, meta$CellType))

# Now use the new column in DESeq2 design
dds <- DESeqDataSetFromMatrix(countData = counts,  
                              colData = meta, 
                              design = ~ Group)
levels(dds$Group)
dds$Group <- relevel(dds$Group, ref = "Q20ESC")

# Run DESeq2
dds <- DESeq(dds)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# PCA plot
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
plotPCA(vsd, intgroup = "CellType")

#ggsave("plot_pca.png", dpi = 600, width = 15, height = 10, units = "in")

# Check available contrasts
resultsNames(dds)

comp_wild_NPC = results(dds, contrast=c("Group", "WTNPC", "KONPC"))
comp_wild_ESC = results(dds, contrast=c("Group", "WTESC", "KOESC"))

comp_Q50_ESC = results(dds, contrast=c("Group", "Q50ESC", "Q20ESC"))
comp_Q89_ESC = results(dds, contrast=c("Group", "Q89ESC", "Q20ESC"))
comp_Q111_ESC = results(dds, contrast=c("Group", "Q111ESC", "Q20ESC"))

comp_Q50_NPC = results(dds, contrast=c("Group", "Q50NPC", "Q20NPC"))
comp_Q89_NPC = results(dds, contrast=c("Group", "Q89NPC", "Q20NPC"))
comp_Q111_NPC = results(dds, contrast=c("Group", "Q111NPC", "Q20NPC"))

#Filtering based on padj and log2FoldChange
filtered_wild_NPC <- comp_wild_NPC[which(comp_wild_NPC$padj < 0.05 & abs(comp_wild_NPC$log2FoldChange) > 0.5), ]
filtered_wild_ESC <- comp_wild_ESC[which(comp_wild_ESC$padj < 0.05 & abs(comp_wild_ESC$log2FoldChange) > 0.5), ]
filtered_wild_NPC
filtered_wild_ESC

filtered_Q50_ESC <- comp_Q50_ESC[which(comp_Q50_ESC$padj < 0.05 & abs(comp_Q50_ESC$log2FoldChange) > 0.5), ]
filtered_Q89_ESC <- comp_Q89_ESC[which(comp_Q89_ESC$padj < 0.05 & abs(comp_Q89_ESC$log2FoldChange) > 0.5), ]
filtered_Q111_ESC <- comp_Q111_ESC[which(comp_Q111_ESC$padj < 0.05 & abs(comp_Q111_ESC$log2FoldChange) > 0.5), ]
filtered_Q50_ESC
filtered_Q89_ESC
filtered_Q111_ESC

filtered_Q50_NPC <- comp_Q50_NPC[which(comp_Q50_NPC$padj < 0.05 & abs(comp_Q50_NPC$log2FoldChange) > 0.5), ]
filtered_Q89_NPC <- comp_Q89_NPC[which(comp_Q89_NPC$padj < 0.05 & abs(comp_Q89_NPC$log2FoldChange) > 0.5), ]
filtered_Q111_NPC <- comp_Q111_NPC[which(comp_Q111_NPC$padj < 0.05 & abs(comp_Q111_NPC$log2FoldChange) > 0.5), ]
filtered_Q50_NPC
filtered_Q89_NPC
filtered_Q111_NPC

# Count upregulated and downregulated genes for Q..NPC
upregulated_Q50_NPC <- sum(filtered_Q50_NPC$log2FoldChange > 0.5)
downregulated_Q50_NPC <- sum(filtered_Q50_NPC$log2FoldChange < -0.5)

upregulated_Q89_NPC <- sum(filtered_Q89_NPC$log2FoldChange > 0.5)
downregulated_Q89_NPC <- sum(filtered_Q89_NPC$log2FoldChange < -0.5)

upregulated_Q111_NPC <- sum(filtered_Q111_NPC$log2FoldChange > 0.5)
downregulated_Q111_NPC <- sum(filtered_Q111_NPC$log2FoldChange < -0.5)

cat("Q50_NPC - Upregulated:", upregulated_Q50_NPC, "Downregulated:", downregulated_Q50_NPC, "\n")
cat("Q89_NPC - Upregulated:", upregulated_Q89_NPC, "Downregulated:", downregulated_Q89_NPC, "\n")
cat("Q111_NPC - Upregulated:", upregulated_Q111_NPC, "Downregulated:", downregulated_Q111_NPC, "\n")

# Create the bar plot 
library(ggplot2)
# Create a data frame with obtained results
df <- data.frame(
  Condition = rep(c("Q50_NPC", "Q89_NPC", "Q111_NPC"), each = 2),
  Regulation = rep(c("Upregulated", "Downregulated"), 3),
  Count = c(15, 19, 4, 26, 3, 11)
)
# Set desired order
df$Condition <- factor(df$Condition, levels = c("Q50_NPC", "Q89_NPC", "Q111_NPC"))

# Plot
ggplot(df, aes(x = Condition, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Differentially Expressed Genes per Genotype in NPC",
    y = "Number of Genes",
    x = "Genotypes in NP cells"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Upregulated" = "steelblue", "Downregulated" = "firebrick"))

#ggsave("plot_deg_NPC.png", dpi = 600, width = 15, height = 10, units = "in")

# Getting filtered DEG from comparison Q89, Q50, Q111 vs Q20 in NPC
library(dplyr)
filtered_deg_89.df <- as.data.frame(filtered_Q89_NPC)
str(filtered_deg_89.df)
#write.csv(filtered_deg_89.df, "deg_Q89_NPC_values_smallRNA.csv")

filtered_deg_50.df <- as.data.frame(filtered_Q50_NPC)
str(filtered_deg_50.df)
#write.csv(filtered_deg_50.df, "deg_Q50_NPC_values_smallRNA.csv")

filtered_deg_111.df <- as.data.frame(filtered_Q111_NPC)
str(filtered_deg_111.df)
#write.csv(filtered_deg_111.df, "deg_Q111_NPC_values_smallRNA.csv")

# Convertion of gene_ids to gene names
# use a GTF file instead of querying gene information from an annotation database like ensemble
library(rtracklayer)
# Load the GTF file
# downloaded from https://www.gencodegenes.org/mouse/release_M14.html (Release M14 (GRCm38.p5))
gtf_file <- "gencode_M14.annotation.gtf"  #used Comprehensive gene annotation (all) GTF file
gtf_data <- import(gtf_file)
# Extract relevant gene information
gtf_df <- as.data.frame(gtf_data)
# Keep only relevant columns (e.g., gene_id and gene_name)
gene_mapping <- gtf_df %>%
  filter(type == "gene") %>%  # Keep only gene entries
  select(gene_id, gene_name, gene_type)  # Adjust these column names based on GTF file structure

# Map ENSEMBL IDs to gene symbols
# Merge the results with the gene mapping table
filtered_deg_89.df$Symbol <- gene_mapping$gene_name[match(rownames(filtered_deg_89.df), gene_mapping$gene_id)]
filtered_deg_89.df$Symbol

filtered_deg_50.df$Symbol <- gene_mapping$gene_name[match(rownames(filtered_deg_50.df), gene_mapping$gene_id)]
filtered_deg_50.df$Symbol

filtered_deg_111.df$Symbol <- gene_mapping$gene_name[match(rownames(filtered_deg_111.df), gene_mapping$gene_id)]
filtered_deg_111.df$Symbol

# Verify the mapping
head(filtered_deg_89.df$Symbol)
head(filtered_deg_50.df$Symbol)
head(filtered_deg_111.df$Symbol)

# Look into gene_types 
filtered_deg_89.df$gene_type <- gene_mapping$gene_type[match(rownames(filtered_deg_89.df), gene_mapping$gene_id)]
filtered_deg_50.df$gene_type <- gene_mapping$gene_type[match(rownames(filtered_deg_50.df), gene_mapping$gene_id)]
filtered_deg_111.df$gene_type <- gene_mapping$gene_type[match(rownames(filtered_deg_111.df), gene_mapping$gene_id)]

# Gene type distribution
# Plotting the gene types with coloring and legend, keeping x-axis title and removing labels
library(forcats) # help to have NA values on the graph
ggplot(filtered_deg_89.df, aes(x = fct_na_value_to_level(gene_type), fill = fct_na_value_to_level(gene_type))) +
  geom_bar() +
  labs(title = "Distribution of Gene Types in NPC Q89 vs Q20",
       x = "Gene Type",
       y = "Count",
       fill = "Gene Type") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text()) +
  scale_fill_brewer(palette = "Set3", na.value = "grey50")

#ggsave("plot_distribution_NPC.png", dpi = 600, width = 15, height = 10, units = "in")
