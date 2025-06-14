# GRUPO_2-3
Repositorio colaborativo del grupo 2-3, para la resolución de la actividad grupal de la asignatura: Introducción a la Programación Científica.

## Descripción del proyecto

Se presenta un un proyecto colaborativo por los integrantes del grupo 2-3 (alumnos de la UNIR en el **Master en Bioinformática**) para el análisis de RNA-seq.

El objetivo de este proyecto es la construcción de un **pipeline automatizado mediante Snakemake**, compuesto por siete scripts que cubren todas las etapas necesarias para el procesamiento y análisis de datos transcriptómicos. A través de este *workflow*, se ejecutará el análisis completo desde cero, estructurando de forma clara y reproducible el orden y la lógica de ejecución de cada uno de los scripts mediante un archivo *snakefile*.

Todos los integrantes del grupo y creadores del repositorio colaborativo quedan reflejados en el archivo *Integrantes.txt*.


## Estructura de carpetas del proyecto: 
├── config/
│ └── samples.txt # Lista de URLs de archivos FASTQ
├── data/ # Archivos descargados (FASTQ)
├── envs/ # (opcional) Entornos Conda
├── results/
│ ├── fastqc/ # Informes FastQC
│ ├── hisat2/ # Archivos BAM alineados
│ ├── featurecounts/ # Matriz de conteos
│ └── deseq2/ # Resultados DE, volcano plot, PCA, heatmap
├── scripts/
│ ├── 1_download.py
│ ├── 2_quality_control.py
│ ├── 3_mapping.py
│ ├── 4_quantification.py
│ ├── 5_deseq2_analysis.R
│ ├── 6_heatmap.R
│ └── 7_pca_plot.R
├── Snakefile
└── README.md

## Requisitos

Snakemake
Python 3.x
R y las librerías: DESeq2, ggplot2, ggrepel, pheatmap
Herramientas externas: wget, fastqc, hisat2, samtools, subread (featureCounts)
