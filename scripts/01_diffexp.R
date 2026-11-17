# Differential gene expression analysis
# Load packages
library(tidyverse)
library(DESeq2)
library(MetaVolcanoR)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)
library(biomaRt)

#Load row counts
counts_data <- readxl::read_excel("data/raw_counts/GSE135668_raw_counts.xlsx")
counts_data <- counts_data |>
  column_to_rownames("Geneid")

#load metadata
metadata <- read.csv("data/metadata/GSE135668_metadata.csv")
glimpse(metadata)

# Prepare sample information
col_data <-data.frame(condition = as.factor(metadata$Conditions),
                      row.names = colnames(counts_data))

# Making sure the row names in `col_data` 
colnames(counts_data)
rownames(col_data)
all(rownames(col_data) %in% colnames(counts_data))

# Are they in same order
all(rownames(col_data) == colnames(counts_data))

#Prepare count matrix
counts_data <- counts_data |>
  mutate_all(as.numeric) |>
  as.matrix()

# Create DESeq2 data abject
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~condition
)

# filter any counts less than 10
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Perform variance-stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# Perform PCA on the transformed data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each PC
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors for each conditions
custom_colors <- c("control" = "#1c9099",
                   "heat stress" = "#de2d26")

# Plot PCA with custom styling
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(aes(label = name),
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("GSE135668") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(title = "condition"))

# Export the plot
ggsave(filename = "figures/pca_plot/GSE135668.png",
       plot = pca_plot,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)

# Analyze the dataset and compile results
res <- results(dds)
res$Gene_ID <- sub("gene:",  "", rownames(res))
res_df <- as.data.frame(res)

# Convert to dataframe
res_df <- as.data.frame(res)
rownames(res_df) <- NULL

# 1. Load required libraries and query AnnotationHub
library(AnnotationHub)
library(AnnotationDbi)

ah <- AnnotationHub()
hits <- query(ah, "Oncorhynchus mykiss")

# Load the newest annotation database (OrgDb)
ah_id <- tail(names(hits), 1)
edb <- ah[[ah_id]]

# 2. Prepare gene IDs
ids <- unique(na.omit(res_df$Gene_ID))
# Remove version suffix (e.g., ".1") if present
ids <- sub("\\..*$", "", ids)

# 3. Look up gene symbols and descriptions using OrgDb keytypes and columns
annotation <- AnnotationDbi::select(
  edb,
  keys = ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

# Check the first few rows of the annotation table
head(annotation)

# Join res_df and annotations
anno_df <- res_df |>
  left_join(annotation, by = c("Gene_ID" = "ENSEMBL"))

#Export results to a csv file
write.csv(anno_df, "outputs/deseq2/GSE135668.csv", row.names = FALSE)

#One-line conversion for metavolcano 
deg_table <- prepare_deseq2(anno_df, gene_column = "SYMBOL")
write.csv(deg_table, "outputs/prepare_deseq2/GSE135668.csv", row.names = FALSE)

# for GSE253274
#Load row counts
counts_data <- readxl::read_excel("data/raw_counts/GSE253274_raw_counts.xlsx")
counts_data <- counts_data |>
  column_to_rownames("Geneid")

#load metadata
metadata <- read.csv("data/metadata/GSE253274_metadata.csv")
glimpse(metadata)

# Prepare sample information
col_data <-data.frame(condition = as.factor(metadata$Conditions),
                      row.names = colnames(counts_data))

# Making sure the row names in `col_data` 
colnames(counts_data)
rownames(col_data)
all(rownames(col_data) %in% colnames(counts_data))

# Are they in same order
all(rownames(col_data) == colnames(counts_data))

#Prepare count matrix
counts_data <- counts_data |>
  mutate_all(as.numeric) |>
  as.matrix()

# Create DESeq2 data abject
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~condition
)

# filter any counts less than 10
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Perform variance-stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# Perform PCA on the transformed data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each PC
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors for each conditions
custom_colors <- c("control" = "#1c9099",
                   "heat stress" = "#de2d26")

# Plot PCA with custom styling
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(aes(label = name),
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("GSE135668") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(title = "condition"))

# Export the plot
ggsave(filename = "figures/pca_plot/GSE253274.png",
       plot = pca_plot,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)

# Analyze the dataset and compile results
res <- results(dds)
res$Gene_ID <- sub("gene:",  "", rownames(res))
res_df <- as.data.frame(res)

# Convert to dataframe
res_df <- as.data.frame(res)
rownames(res_df) <- NULL

# 1. Load required libraries and query AnnotationHub
library(AnnotationHub)
library(AnnotationDbi)

ah <- AnnotationHub()
hits <- query(ah, "Oncorhynchus mykiss")

# Load the newest annotation database (OrgDb)
ah_id <- tail(names(hits), 1)
edb <- ah[[ah_id]]

# 2. Prepare gene IDs
ids <- unique(na.omit(res_df$Gene_ID))
# Remove version suffix (e.g., ".1") if present
ids <- sub("\\..*$", "", ids)

# 3. Look up gene symbols and descriptions using OrgDb keytypes and columns
annotation <- AnnotationDbi::select(
  edb,
  keys = ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

# Check the first few rows of the annotation table
head(annotation)

# Join res_df and annotations
anno_df <- res_df |>
  left_join(annotation, by = c("Gene_ID" = "ENSEMBL"))

#Export results to a csv file
write.csv(anno_df, "outputs/deseq2/GSE253274.csv", row.names = FALSE)

#One-line conversion for metavolcano 
deg_table <- prepare_deseq2(anno_df, gene_column = "SYMBOL")
write.csv(deg_table, "outputs/prepare_deseq2/GSE253274.csv", row.names = FALSE)

# for GSE239806
#Load row counts
counts_data <- readxl::read_excel("data/raw_counts/GSE239806_raw_counts.xlsx")
counts_data <- counts_data |>
  column_to_rownames("Geneid")

#load metadata
metadata <- read.csv("data/metadata/GSE239806_metadata.csv")
glimpse(metadata)

# Prepare sample information
col_data <-data.frame(condition = as.factor(metadata$Conditions),
                      row.names = colnames(counts_data))

# Making sure the row names in `col_data` 
colnames(counts_data)
rownames(col_data)
all(rownames(col_data) %in% colnames(counts_data))

# Are they in same order
all(rownames(col_data) == colnames(counts_data))

#Prepare count matrix
counts_data <- counts_data |>
  mutate_all(as.numeric) |>
  as.matrix()

# Create DESeq2 data abject
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~condition
)

# filter any counts less than 10
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Perform variance-stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# Perform PCA on the transformed data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each PC
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors for each conditions
custom_colors <- c("control" = "#1c9099",
                   "heat stress" = "#de2d26")

# Plot PCA with custom styling
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(aes(label = name),
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("GSE239806") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(title = "condition"))

# Export the plot
ggsave(filename = "figures/pca_plot/GSE239806.png",
       plot = pca_plot,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)

# Analyze the dataset and compile results
res <- results(dds)
res$Gene_ID <- sub("gene:",  "", rownames(res))
res_df <- as.data.frame(res)

# Convert to dataframe
res_df <- as.data.frame(res)
rownames(res_df) <- NULL

# 1. Load required libraries and query AnnotationHub
library(AnnotationHub)
library(AnnotationDbi)

ah <- AnnotationHub()
hits <- query(ah, "Oncorhynchus mykiss")

# Load the newest annotation database (OrgDb)
ah_id <- tail(names(hits), 1)
edb <- ah[[ah_id]]

# 2. Prepare gene IDs
ids <- unique(na.omit(res_df$Gene_ID))
# Remove version suffix (e.g., ".1") if present
ids <- sub("\\..*$", "", ids)

# 3. Look up gene symbols and descriptions using OrgDb keytypes and columns
annotation <- AnnotationDbi::select(
  edb,
  keys = ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

# Check the first few rows of the annotation table
head(annotation)

# Join res_df and annotations
anno_df <- res_df |>
  left_join(annotation, by = c("Gene_ID" = "ENSEMBL"))

#Export results to a csv file
write.csv(anno_df, "outputs/deseq2/GSE239806.csv", row.names = FALSE)

#One-line conversion for metavolcano 
deg_table <- prepare_deseq2(anno_df, gene_column = "SYMBOL")
write.csv(deg_table, "outputs/prepare_deseq2/GSE239806.csv", row.names = FALSE)

# for GSE262612
# 1. Load raw count matrix properly and handle any potential missing values
counts_data <- readxl::read_excel("data/raw_counts/GSE262612_raw_counts.xlsx")
counts_data <- as.data.frame(counts_data)

# Set Geneid as row names safely
counts_data <- counts_data |>
  column_to_rownames("Geneid")

# Clean counts: replace NA with 0, convert everything to numeric, and convert to matrix
counts_data[is.na(counts_data)] <- 0
counts_data <- counts_data |>
  mutate(across(everything(), as.numeric)) |>
  as.matrix()

# 2. Load sample metadata
metadata <- read.csv("data/metadata/GSE262612_metadata.csv")
glimpse(metadata) 

# 3. Prepare and align sample information with the count matrix columns securely
rownames(metadata) <- metadata$Samples
matched_metadata <- metadata[match(colnames(counts_data), metadata$Samples), ]

# 4. Create col_data dataframe for DESeq2 input
col_data <- data.frame(
  condition = as.factor(matched_metadata$Conditions),
  row.names = colnames(counts_data)
)

# Verify alignment
all(rownames(col_data) == colnames(counts_data))

# 5. Create DESeq2 data object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~condition
)

# Filter out low counts (less than 10 total counts across samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Perform variance-stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# Perform PCA on the transformed data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each PC
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors for each condition
custom_colors <- c("control" = "#1c9099",
                   "heat stress" = "#de2d26")

# Plot PCA with custom styling
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(aes(label = name),
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("GSE262612") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(title = "condition"))

# Export the PCA plot
ggsave(filename = "figures/pca_plot/GSE262612.png",
       plot = pca_plot,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)

# Analyze the dataset and compile results
res <- results(dds)
res$Gene_ID <- sub("gene:", "", rownames(res))
res_df <- as.data.frame(res)
rownames(res_df) <- NULL

# 1. Load required libraries and query AnnotationHub for Oncorhynchus mykiss
library(AnnotationHub)
library(AnnotationDbi)

ah <- AnnotationHub()
hits <- query(ah, "Oncorhynchus mykiss")

# Load the newest annotation database (OrgDb)
ah_id <- tail(names(hits), 1)
edb <- ah[[ah_id]]

# 2. Prepare gene IDs
ids <- unique(na.omit(res_df$Gene_ID))
ids <- sub("\\..*$", "", ids)

# 3. Look up gene symbols and descriptions using OrgDb keytypes and columns
annotation <- AnnotationDbi::select(
  edb,
  keys = ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME")
)

# Check the first few rows of the annotation table
head(annotation)

# Join res_df and annotations
anno_df <- res_df |>
  left_join(annotation, by = c("Gene_ID" = "ENSEMBL"))

# Export results to a csv file
write.csv(anno_df, "outputs/deseq2/GSE262612.csv", row.names = FALSE)

# One-line conversion for metavolcano 
deg_table <- prepare_deseq2(anno_df, gene_column = "SYMBOL")
write.csv(deg_table, "outputs/prepare_deseq2/GSE262612.csv", row.names = FALSE)

