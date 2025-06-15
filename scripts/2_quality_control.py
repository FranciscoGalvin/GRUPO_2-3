#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

def check_status(return_code, paso):
    if return_code == 0:
        print(f"Paso '{paso}' ejecutado con éxito.")
    else:
        print(f"Error en el paso '{paso}'. Abortando.")
        sys.exit(1)

def ejecutar_fastqc(archivo_listado, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    check_status(0, "Crear directorio de salida")

    with open(archivo_listado, "r") as f:
        archivos = [line.strip() for line in f if line.strip()]

    print("🧪 Ejecutando FastQC...")
    for archivo in archivos:
        print(f" -> Analizando {archivo}")
        result = subprocess.run(["fastqc", archivo, "-o", output_dir])
        check_status(result.returncode, f"FastQC: {os.path.basename(archivo)}")

def main():
    parser = argparse.ArgumentParser(description="Ejecuta FastQC sobre una lista de archivos FASTQ.")
    parser.add_argument("-i", "--input", required=True, help="Archivo .txt con rutas de archivos FASTQ (una por línea)")
    parser.add_argument("-o", "--output", required=True, help="Ruta del directorio de salida")

    args = parser.parse_args()
    ejecutar_fastqc(args.input, args.output)

if __name__ == "__main__":
    main()
