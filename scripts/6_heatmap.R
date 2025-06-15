#!/usr/bin/env Rscript

suppressMessages({
  library(DESeq2)
  library(pheatmap)
})

# Leer los conteos
counts <- read.table("results/featurecounts/cleaned_gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Definir condiciones
colData <- data.frame(row.names = colnames(counts), group = c(rep("Patología", 7), rep("Control", 7)))

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)

# Filtrar por padj < 0.05
results_filtered <- subset(res, padj < 0.05)
counts_filtered <- counts[rownames(results_filtered), ]

# Normalización y selección
rld <- rlog(dds, blind = FALSE)
mat <- assay(rld)[rownames(results_filtered), ]

# Guardar heatmap
png("results/deseq2/heatmap_genes_significativos.png", width = 800, height = 600)
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = colData,
         main = "Heatmap de genes diferencialmente expresados")
dev.off()
