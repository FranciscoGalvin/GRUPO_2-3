# Diagramas del Pipeline

En esta carpeta se encuentra almacenada la imagen del diagrama del pipeline de análisis RNA-Seq. Este diagrama ilustra de manera esquemática las etapas principales del flujo de trabajo, desde la descarga de datos hasta el análisis de expresión diferencial y la visualización de resultados.

El archivo `pipeline_diagram.svg` fue generado utilizando Mermaid y puede integrarse en presentaciones o documentos para facilitar la comprensión del proceso automatizado.

El codigo de mermaid para generar el diagrama es el siguiente:
---
config:
  layout: dagre
---
flowchart TD
    A["Inicio: samples.txt con URLs"] --> B["Descarga de datos 1_download.py"]
    B --> C["Control de calidad 2_quality_control.py"]
    C --> D["Mapeo con HISAT2 3_mapping.py"]
    D --> E["Cuantificación de genes 4_quantification.py"]
    E --> F["Análisis de expresión diferencial 5_deseq2_analysis.R"] & H["Heatmap 6_heatmap.R"] & I["Análisis PCA 7_pca_plot.R"]
    F --> G["Volcano Plot"]
    style A fill:#f9f,stroke:#333,stroke-width:1px
    style B fill:#bbf,stroke:#333,stroke-width:1px
    style C fill:#bbf,stroke:#333,stroke-width:1px
    style D fill:#bbf,stroke:#333,stroke-width:1px
    style E fill:#bbf,stroke:#333,stroke-width:1px
    style F fill:#bbf,stroke:#333,stroke-width:1px
    style H fill:#bfb,stroke:#333,stroke-width:1px
    style I fill:#bfb,stroke:#333,stroke-width:1px
    style G fill:#bfb,stroke:#333,stroke-width:1px
