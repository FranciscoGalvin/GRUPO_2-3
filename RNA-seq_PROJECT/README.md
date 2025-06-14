#RNA-seq Analysis pipeline. 
Este repositorio contiene un pipeline completo para el análisis de datos de RNA-seq basado en snakemake. El flujo incluye la descarga de datos, control de calidad, alineamiento, cuantificación, análisis de expresión diferencial y visualizaciones.


#Estructura de carpetas del proyecto:
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


##Requisitos 
- [Snakemake](https://snakemake.readthedocs.io/)
- Python 3.x
- R y las librerías: `DESeq2`, `ggplot2`, `ggrepel`, `pheatmap`
- Herramientas externas: `wget`, `fastqc`, `hisat2`, `samtools`, `subread` (featureCounts)
