#!/usr/bin/env python3

import argparse
import os
import subprocess

def descargar_datos(urls, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for url in urls:
        url = url.strip()
        if not url:
            continue
        print(f"Descargando: {url}")
        try:
            subprocess.run(["wget", "-P", output_dir, url], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error al descargar {url}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Descargar múltiples archivos RNA-Seq desde una lista de URLs.")
    parser.add_argument("-o", "--output", required=True, help="Directorio de salida")
    parser.add_argument("-l", "--list", required=True, help="Archivo de texto con URLs (una por línea)")

    args = parser.parse_args()

    with open(args.list, "r") as f:
        urls = f.readlines()

    descargar_datos(urls, args.output)

if __name__ == "__main__":
    main()
