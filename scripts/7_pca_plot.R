#!/usr/bin/env Rscript

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Uso: 6_pca_plot.R <counts_file> <output_plot>")
}

input_file <- args[1]
output_plot <- args[2]

# Leer conteos
counts <- read.table(file = input_file, header = TRUE, row.names = 1, sep = "\t")

# Condiciones (ajústalas si cambia el orden o el número de muestras)
colData <- data.frame(
  row.names = colnames(counts),
  group = c(rep("Patología", 7), rep("Control", 7))
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
dds <- DESeq(dds)
rld <- rlog(dds, blind = TRUE)

# PCA
pca_data <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Gráfico
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = rownames(pca_data)), size = 3, box.padding = 0.5) +
  labs(
    title = "PCA of Gene Expression Data",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_bw() +
  scale_color_manual(values = c("Patología" = "red", "Control" = "blue"))

ggsave(output_plot, plot = p)
