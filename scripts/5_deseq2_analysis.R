#!/usr/bin/env Rscript

# Cargar librerías
suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Uso: 5_deseq2_analysis.R <counts_file> <output_table> <output_plot>")
}

counts_file <- args[1]
output_file <- args[2]
volcano_file <- args[3]

# Leer conteos
counts <- read.table(file = counts_file, header = TRUE, row.names = 1, sep = "\t")

# Condiciones experimentales (ajústalas si cambia el orden)
colData <- data.frame(
  row.names = colnames(counts),
  group = c(rep("Patología", 7), rep("Control", 7))
)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)

# Guardar resultados
write.table(res, file = output_file, sep = "\t", quote = FALSE, col.names = NA)

# Volcano plot
res$log10padj <- -log10(res$padj)
res$category <- "Not Significant"
res$category[res$log2FoldChange > 1 & res$padj < 0.05] <- "Upregulated"
res$category[res$log2FoldChange < -1 & res$padj < 0.05] <- "Downregulated"
res$gene <- rownames(res)
res$label <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, res$gene, NA)

p <- ggplot(res, aes(x = log2FoldChange, y = log10padj, color = category)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = label), size = 3, box.padding = 0.3, max.overlaps = 15) +
  labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(p-adj)") +
  theme_bw()

ggsave(volcano_file, plot = p)
