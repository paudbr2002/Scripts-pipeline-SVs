#!/bin/bash
# author: Paula de Blas
# transformar en bed los .vcf de los BED
carpeta_entrada="$1"
# Procesar cada subcarpeta (delly, GRIDSS, Manta)
for subcarpeta in "$carpeta_entrada"/*; do
    if [ -d "$subcarpeta" ]; then
        pass_dir="${subcarpeta}/${subcarpeta##*/}_PASS"
        if [ -d "$pass_dir" ]; then
            cd "$pass_dir"
            # Convertir archivos VCF a BED
            for archivo in *.vcf; do
                [ -e "$archivo" ] || continue  # saltar si no hay archivos
                nombre_archivo=$(basename "$archivo" .vcf)
                vcf2bed < "$archivo" > "${nombre_archivo}.bed"
            done
        fi
    fi
done
