# -------------------------------
# CONFIGURACIÓN DE MUESTRAS
# -------------------------------

SAMPLE_FILES = [l.strip().split("/")[-1] for l in open("config/samples.txt") if l.strip()]
FASTQ_FILES = [f"data/{f}" for f in SAMPLE_FILES]

def get_sample_base(filename):
    return filename.replace("_f.fastq.gz", "").replace("_r.fastq.gz", "")

SAMPLES = sorted(set(get_sample_base(f) for f in SAMPLE_FILES if f.endswith("_f.fastq.gz")))

FASTQC_HTMLS = [f"results/fastqc/{f.replace('.fastq.gz', '')}_fastqc.html" for f in SAMPLE_FILES]
MAPPED_BAMS = [f"results/hisat2/{sample}_sorted.bam" for sample in SAMPLES]
MAPPED_BAIS = [f"{bam}.bai" for bam in MAPPED_BAMS]

FEATURECOUNTS_CLEANED = "results/featurecounts/cleaned_gene_counts.txt"
DESEQ2_TABLE = "results/deseq2/differential_expression_results.txt"
DESEQ2_PLOT = "results/deseq2/volcano_plot_with_labels.png"
PCA_PLOT = "results/deseq2/pca_plot.png"
HEATMAP = "results/deseq2/heatmap_genes_significativos.png"

# -------------------------------
# OBJETIVO FINAL
# -------------------------------

rule all:
    input:
        FASTQC_HTMLS,
        MAPPED_BAMS,
        MAPPED_BAIS,
        FEATURECOUNTS_CLEANED,
        DESEQ2_TABLE,
        DESEQ2_PLOT,
        PCA_PLOT,
        HEATMAP

# -------------------------------
# 1. DESCARGA DE DATOS
# -------------------------------

rule download:
    input:
        "config/samples.txt"
    output:
        FASTQ_FILES
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        "python3 scripts/1_download.py -l {input} -o data"

# -------------------------------
# 2. CONTROL DE CALIDAD (FASTQC)
# -------------------------------

rule quality_control:
    input:
        FASTQ_FILES
    output:
        FASTQC_HTMLS
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        """
        echo "{input}" | tr ' ' '\\n' > temp_fastq_list.txt
        python3 scripts/2_quality_control.py -i temp_fastq_list.txt -o results/fastqc
        rm temp_fastq_list.txt
        """

# -------------------------------
# 3. MAPEADO CON HISAT2
# -------------------------------

rule mapping:
    input:
        r1 = "data/{sample}_f.fastq.gz",
        r2 = "data/{sample}_r.fastq.gz"
    output:
        bam = "results/hisat2/{sample}_sorted.bam",
        bai = "results/hisat2/{sample}_sorted.bam.bai"
    params:
        genome_index = "data/genome/genome_index/genome_index",
        out_prefix = "{sample}",
        output_dir = "results/hisat2"
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        """
        mkdir -p {params.output_dir}
        python3 scripts/3_mapping.py \
            --genome_index {params.genome_index} \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --out_prefix {params.out_prefix} \
            --output_dir {params.output_dir}
        """

# -------------------------------
# 4. CUANTIFICACIÓN CON featureCounts
# -------------------------------

rule quantification:
    input:
        bams = MAPPED_BAMS
    output:
        cleaned = FEATURECOUNTS_CLEANED
    params:
        annotation = "data/annotation_files/Homo_sapiens.GRCh38.110.gtf",
        output_dir = "results/featurecounts"
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        """
        mkdir -p {params.output_dir}
        echo "{input.bams}" | tr ' ' '\\n' > bam_list.txt
        python3 scripts/4_quantification.py \
            -i bam_list.txt \
            -a {params.annotation} \
            -o {params.output_dir}
        rm bam_list.txt
        """

# -------------------------------
# 5. ANÁLISIS DE EXPRESIÓN (DESeq2)
# -------------------------------

rule deseq2:
    input:
        counts = FEATURECOUNTS_CLEANED
    output:
        table = DESEQ2_TABLE,
        volcano = DESEQ2_PLOT
    params:
        script = "scripts/5_deseq2_analysis.R"
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        """
        mkdir -p results/deseq2
        Rscript {params.script} {input.counts} {output.table} {output.volcano}
        """

# -------------------------------
# 6. PCA PLOT
# -------------------------------

rule pca_plot:
    input:
        counts = FEATURECOUNTS_CLEANED
    output:
        pca = PCA_PLOT
    params:
        script = "scripts/7_pca_plot.R"
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        """
        mkdir -p results/deseq2
        Rscript {params.script} {input.counts} {output.pca}
        """

# -------------------------------
# 7. HEATMAP
# -------------------------------

rule heatmap:
    input:
        counts = FEATURECOUNTS_CLEANED
    output:
        HEATMAP
    params:
        script = "scripts/6_heatmap.R"
    conda:
        "envs/rnaseq-pipeline.yaml"
    shell:
        """
        mkdir -p results/deseq2
        Rscript {params.script}
        """
