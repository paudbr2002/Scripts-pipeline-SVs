#!/bin/bash
# author: Paula de Blas
# Transformar los BCF en VCF con bcftools y quitar las variantes en el cromosoma mitocondrial solo en la carpeta de delly, y separar variantes en PASS y No_PASS en todas las subcarpetas.

# Procesar solo la carpeta "delly" para convertir BCF a VCF y quitar variantes en el cromosoma mitocondrial
delly_dir="$carpeta_entrada/delly"
if [ -d "$delly_dir" ]; then
    cd "$delly_dir"

    # Transformar los archivos BCF en VCF
    for archivo in *.bcf; do
        nombre_archivo=$(basename "$archivo" .bcf)
        bcftools view -o "${nombre_archivo}.vcf" "$archivo"
    done

    # Quitar las variantes en el cromosoma mitocondrial
    for archivo in *.vcf; do
        nombre_archivo=$(basename "$archivo" .vcf)
        grep -v "chrM" "$archivo" > "${nombre_archivo}_sin_chrM.vcf"
    done
fi

# Procesar todas las subcarpetas (delly, GRIDSS, Manta) para separar variantes en PASS y No_PASS
for subcarpeta in "$carpeta_entrada"/*; do
    if [ -d "$subcarpeta" ]; then
        cd "$subcarpeta"

        # Crear carpetas de salida
        mkdir -p "${subcarpeta}/${subcarpeta##*/}_PASS"
        mkdir -p "${subcarpeta}/No_PASS"

        # Usar los archivos _sin_chrM.vcf si están presentes (solo en delly), de lo contrario usar los archivos .vcf
        for archivo in *_sin_chrM.vcf *.vcf; do
            [ -e "$archivo" ] || continue  # saltar si no hay archivos

            # Evitar procesar dos veces los archivos .vcf en la carpeta delly
            if [[ "$archivo" == *_sin_chrM.vcf ]] || [[ "$subcarpeta" != *"delly" ]]; then
                nombre_archivo=$(basename "$archivo" .vcf)
                nombre_archivo=$(basename "$nombre_archivo" _sin_chrM)

                # Guardar el encabezado del archivo VCF en una variable
                header=$(awk '/^#CHROM/ { print }' "$archivo")

                # Escribir el encabezado en los archivos de salida
                echo "$header" > "${subcarpeta}/${subcarpeta##*/}_PASS/${nombre_archivo}_PASS.vcf"
                echo "$header" > "${subcarpeta}/No_PASS/${nombre_archivo}_noPASS.vcf"

                # Filtrar las líneas del archivo según "PASS" y "No_PASS"
                awk -F"\t" -v carpeta_pass="${subcarpeta##*/}_PASS" -v nombre_archivo="$nombre_archivo" -v header="$header" '{ if ($7 == "PASS") { print $0 >> carpeta_pass "/" nombre_archivo"_PASS.vcf" } else { print $0 >> "No_PASS/" nombre_archivo"_noPASS.vcf" } }' "$archivo"
            fi
        done
    fi
done

