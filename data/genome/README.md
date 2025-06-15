# genome_index

Esta carpeta contiene los archivos del índice del genoma de referencia, generados con `hisat2-build`. Son necesarios para alinear las lecturas FASTQ al genoma durante el paso de mapeo del pipeline.

Los archivos `.ht2` que se generan forman el índice completo de HISAT2 y deben estar todos en esta carpeta. Por ejemplo:

genome_index.1.ht2
genome_index.2.ht2
...
genome_index.8.ht2

Debido a su tamaño, estos archivos no se incluyen directamente en el repositorio. Para reproducir el pipeline, deben generarse manualmente a partir del archivo FASTA del genoma, por ejemplo:

```bash
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa genome_index

