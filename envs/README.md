# Entornos Conda para el Pipeline RNA-Seq

Esta carpeta contiene los archivos de definición de entornos Conda utilizados en el pipeline de análisis RNA-Seq automatizado con Snakemake. Cada archivo YAML define un entorno específico con las herramientas necesarias para ejecutar distintos pasos del flujo de trabajo.

## Archivos

- **rnaseq-pipeline.yaml**  
  Este entorno es el principal del pipeline e incluye herramientas clave como `snakemake`, `hisat2`, `samtools`, `featureCounts`, y paquetes de R como `DESeq2`, `ggplot2`, y `pheatmap` para el análisis de expresión diferencial y visualizaciones.

- **fastqc.yaml**  
  Entorno más ligero, diseñado específicamente para ejecutar FastQC, una herramienta de control de calidad de datos de secuenciación en formato FASTQ.

## Uso

Los entornos se activan automáticamente cuando se ejecuta Snakemake con la opción `--use-conda`, siempre que las reglas del Snakefile incluyan la directiva `conda:` apuntando a uno de estos archivos.

```bash
snakemake --cores 4 --use-conda
