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
        print(f" {paso} ejecutado con éxito.")

def main():
    parser = argparse.ArgumentParser(description="Mapeo con HISAT2 + samtools sort + index.")
    parser.add_argument("--genome_index", required=True, help="Ruta base al índice de HISAT2")
    parser.add_argument("--r1", required=True, help="FASTQ de lectura forward")
    parser.add_argument("--r2", required=True, help="FASTQ de lectura reverse")
    parser.add_argument("--out_prefix", required=True, help="Prefijo para los archivos de salida (.bam)")
    parser.add_argument("--output_dir", required=True, help="Directorio de salida")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Paso 1: alineamiento
    print(f"📌 Mapeando: {args.r1} + {args.r2}")
    bam_tmp = os.path.join(args.output_dir, f"{args.out_prefix}_mapped_reads.bam")
    hisat2_cmd = [
        "hisat2", "-x", args.genome_index, "-1", args.r1, "-2", args.r2
    ]
    with open(bam_tmp, "wb") as out_bam:
        p1 = subprocess.Popen(hisat2_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=p1.stdout, stdout=out_bam)
        p1.stdout.close()
        p2.communicate()
    check_status(p2.returncode, f"HISAT2 -> BAM: {args.out_prefix}")

    # Paso 2: ordenar
    bam_sorted = os.path.join(args.output_dir, f"{args.out_prefix}_sorted.bam")
    sort_cmd = ["samtools", "sort", "-o", bam_sorted, bam_tmp]
    check_status(subprocess.call(sort_cmd), f"samtools sort: {args.out_prefix}")

    # Paso 3: indexar
    index_cmd = ["samtools", "index", bam_sorted]
    check_status(subprocess.call(index_cmd), f"samtools index: {args.out_prefix}")

    # Limpiar intermedio
    os.remove(bam_tmp)

if __name__ == "__main__":
    main()
