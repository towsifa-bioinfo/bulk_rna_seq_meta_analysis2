# Load packages
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)

# Create output folders if missing
dir.create("figures/venn_diagram", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/venn",         recursive = TRUE, showWarnings = FALSE)

# Load DESeq2 results
GSE262612_diffexp <- read.csv("outputs/deseq2/GSE262612.csv")
GSE253274_diffexp <- read.csv("outputs/deseq2/GSE253274.csv")
GSE239806_diffexp <- read.csv("outputs/deseq2/GSE239806.csv")
GSE135668_diffexp <- read.csv("outputs/deseq2/GSE135668.csv")

# Keep significant genes
GSE262612_diffexp <- GSE262612_diffexp |>
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE253274_diffexp <- GSE253274_diffexp |>
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE239806_diffexp <- GSE239806_diffexp |>
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

GSE135668_diffexp <- GSE135668_diffexp |>
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Helper: remove NA / empty symbols and duplicates
clean_set <- function(x) {
  x <- as.character(x)
  unique(x[!is.na(x) & x != ""])
}

# 1. All four datasets
gene_sets_all <- list(
  GSE262612 = clean_set(GSE262612_diffexp$SYMBOL),
  GSE253274 = clean_set(GSE253274_diffexp$SYMBOL),
  GSE239806 = clean_set(GSE239806_diffexp$SYMBOL),
  GSE135668 = clean_set(GSE135668_diffexp$SYMBOL)
)
print(sapply(gene_sets_all, length))  

venn_plot_all <- venn.diagram(
  x = gene_sets_all,
  filename = NULL,
  fill = c("#99d8c9", "#addd8e", "#bcbddc", "#fec44f"),
  alpha = 0.6,
  cex = 0.8,
  cat.cex = 1,
  cat.col = c("#2ca25f", "#31a354", "#756bb1", "#e6550d"),
  cat.pos = c(-20, 20, -20, 20),
  cat.dist = c(0.22, 0.22, 0.1, 0.09),
  margin = 0.1,
  lwd = 2,
  disable.logging = TRUE
)

png("figures/venn_diagram/VennDiagram.png",
    units = "in", width = 6, height = 6, res = 600)
grid::grid.draw(venn_plot_all)
dev.off()

# Genes common to all four studies
intersect_all <- process_region_data(Venn(gene_sets_all))

intersect_genes_all <- intersect_all |>
  dplyr::filter(id == "1/2/3/4") |>
  dplyr::pull(item) |>
  as.data.frame()
colnames(intersect_genes_all) <- "SYMBOL"

write.csv(intersect_genes_all, "outputs/venn/intersect_genes_all.csv",
          row.names = FALSE)

# 2. Heat datasets (GSE262612 vs GSE135668)
gene_sets_heat <- list(
  GSE262612 = clean_set(GSE262612_diffexp$SYMBOL),
  GSE135668 = clean_set(GSE135668_diffexp$SYMBOL)
)

venn_plot_heat <- venn.diagram(
  x = gene_sets_heat,
  filename = NULL,
  fill = c("#99d8c9", "#fec44f"),
  alpha = 0.6,
  cex = 0.8,
  cat.cex = 1,
  cat.col = c("#31a354", "#e6550d"),
  cat.pos = c(0, 0),
  cat.dist = c(0.02, 0.02),
  margin = 0.1,
  lwd = 2,
  disable.logging = TRUE
)

png("figures/venn_diagram/VennDiagram_Heat.png",
    units = "in", width = 6, height = 6, res = 600)
grid::grid.draw(venn_plot_heat)
dev.off()

intersect_heat <- process_region_data(Venn(gene_sets_heat))

intersect_genes_heat <- intersect_heat |>
  dplyr::filter(id == "1/2") |>
  dplyr::pull(item) |>
  as.data.frame()
colnames(intersect_genes_heat) <- "SYMBOL"

write.csv(intersect_genes_heat, "outputs/venn/intersect_genes_heat.csv",
          row.names = FALSE)

# 3. Hypoxia datasets (GSE253274 vs GSE239806)
gene_sets_hypoxia <- list(
  GSE253274 = clean_set(GSE253274_diffexp$SYMBOL),
  GSE239806 = clean_set(GSE239806_diffexp$SYMBOL)
)

venn_plot_hypoxia <- venn.diagram(
  x = gene_sets_hypoxia,
  filename = NULL,
  fill = c("#bcbddc", "#addd8e"),
  alpha = 0.6,
  cex = 0.8,
  cat.cex = 1,
  cat.col = c("#756bb1", "#2ca25f"),
  cat.pos = c(0, 0),
  cat.dist = c(0.02, 0.02),
  margin = 0.1,
  lwd = 2,
  scaled = FALSE,
  disable.logging = TRUE
)

png("figures/venn_diagram/VennDiagram_Hypoxia.png",
    units = "in", width = 6, height = 6, res = 600)
grid::grid.draw(venn_plot_hypoxia)
dev.off()

intersect_hypoxia <- process_region_data(Venn(gene_sets_hypoxia))

intersect_genes_hypoxia <- intersect_hypoxia |>
  dplyr::filter(id == "1/2") |>
  dplyr::pull(item) |>
  as.data.frame()
colnames(intersect_genes_hypoxia) <- "SYMBOL"

write.csv(intersect_genes_hypoxia, "outputs/venn/intersect_genes_hypoxia.csv",
          row.names = FALSE)

