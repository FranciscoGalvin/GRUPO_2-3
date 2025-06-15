#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

def check_status(returncode, paso):
    if returncode != 0:
        print(f"Error en '{paso}'. Abortando.")
        sys.exit(1)
    else:
        print(f"{paso} ejecutado con éxito.")

def main():
    parser = argparse.ArgumentParser(description="Ejecuta featureCounts sobre varios BAM y limpia el output.")
    parser.add_argument("-i", "--input_list", required=True, help="Archivo .txt con rutas a archivos BAM (uno por línea)")
    parser.add_argument("-a", "--annotation", required=True, help="Archivo GTF con anotaciones")
    parser.add_argument("-o", "--output_dir", required=True, help="Directorio de salida")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Leer lista de BAMs
    with open(args.input_list, "r") as f:
        bam_files = [line.strip() for line in f if line.strip()]

    bam_args = " ".join(bam_files)
    output_file = os.path.join(args.output_dir, "combined_gene_counts.txt")
    cleaned_output_file = os.path.join(args.output_dir, "cleaned_gene_counts.txt")

    # Ejecutar featureCounts
    featurecounts_cmd = f"""
    featureCounts -T 8 -p -B -a {args.annotation} -o {output_file} {bam_args}
    """
    check_status(subprocess.call(featurecounts_cmd, shell=True), "FeatureCounts")

    # Limpiar y renombrar columnas del output
    limpiar_cmd = f"""
    tail -n +2 {output_file} | cut --complement -f 2-6 > {cleaned_output_file}
    """
    check_status(subprocess.call(limpiar_cmd, shell=True), "Limpieza columnas")

    # Cambiar encabezado
    num_muestras = len(bam_files)
    nombres_col = ["GeneID"] + [f"sample_{i+1}" for i in range(num_muestras)]
    header = "\\t".join(nombres_col)
    sed_cmd = f"sed -i '1s/.*/{header}/' {cleaned_output_file}"
    check_status(subprocess.call(sed_cmd, shell=True), "Renombrar columnas")

if __name__ == "__main__":
    main()
