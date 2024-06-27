#!/bin/bash
# author: Paula de Blas
# Transformar los BCF en VCF con bcftools, quitar las variantes en el cromosoma mitocondrial solo en la carpeta de delly, y separar variantes en PASS y No_PASS en todas las subcarpetas.
carpeta_entrada="$1"
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
        mkdir -p "${subcarpeta}/CONTIGS"

        # Usar los archivos _sin_chrM.vcf si están presentes (solo en delly), de lo contrario usar los archivos .vcf
        for archivo in *_sin_chrM.vcf *.vcf; do
            [ -e "$archivo" ] || continue  # saltar si no hay archivos

            # Evitar procesar dos veces los archivos .vcf en la carpeta delly
            if [[ "$archivo" == *_sin_chrM.vcf ]] || [[ "$subcarpeta" != *"delly" ]]; then
                nombre_archivo=$(basename "$archivo" .vcf)
                nombre_archivo=$(basename "$nombre_archivo" _sin_chrM)

                # Guardar el encabezado del archivo VCF en una variable
                header1=$(awk '/^##/ { print }' "$archivo")
                header2=$(awk '/^#CHROM/ { print }' "$archivo")

                # Escribir el encabezado en los archivos de salida
                echo "$header1" > "${subcarpeta}/${subcarpeta##*/}_PASS/${nombre_archivo}_PASS.vcf"
                echo "$header2" >> "${subcarpeta}/${subcarpeta##*/}_PASS/${nombre_archivo}_PASS.vcf"
                echo "$header1" > "${subcarpeta}/No_PASS/${nombre_archivo}_noPASS.vcf"
                echo "$header2" >> "${subcarpeta}/No_PASS/${nombre_archivo}_noPASS.vcf"

                # Filtrar las líneas del archivo según "PASS" y "No_PASS"
                awk -F"\t" -v carpeta_pass="${subcarpeta}/${subcarpeta##*/}_PASS" -v carpeta_no_pass="${subcarpeta}/No_PASS" -v nombre_archivo="$nombre_archivo" '{ if ($7 == "PASS") { print $0 >> carpeta_pass "/" nombre_archivo"_PASS.vcf" } else { print $0 >> carpeta_no_pass "/" nombre_archivo"_noPASS.vcf" } }' "$archivo"

                # Filtrar las líneas que contienen "chr" seguido de números o letras y un guion bajo "_", y moverlas a CONTIGS
                echo "$header1" > "${subcarpeta}/CONTIGS/${nombre_archivo}_CONTIGS.vcf"
                echo "$header2" >> "${subcarpeta}/CONTIGS/${nombre_archivo}_CONTIGS.vcf"
                awk -F"\t" -v carpeta_contigs="${subcarpeta}/CONTIGS" -v nombre_archivo="$nombre_archivo" '{ if ($7 == "PASS" && ($1 ~ /^chr[0-9A-Za-z]+_/)) { print $0 >> carpeta_contigs "/" nombre_archivo"_CONTIGS.vcf" } }' "${subcarpeta}/${subcarpeta##*/}_PASS/${nombre_archivo}_PASS.vcf"

                # Remover las líneas que contienen "chr" seguido de números o letras y un guion bajo "_" del archivo PASS original y añadir los encabezados de nuevo
                echo "$header1" > "${nombre_archivo}_temp.vcf"
                echo "$header2" >> "${nombre_archivo}_temp.vcf"
                awk -F"\t" '{ if ($7 == "PASS" && !($1 ~ /^chr[0-9A-Za-z]+_/)) { print $0 } }' "${subcarpeta}/${subcarpeta##*/}_PASS/${nombre_archivo}_PASS.vcf" >> "${nombre_archivo}_temp.vcf"
                mv "${nombre_archivo}_temp.vcf" "${subcarpeta}/${subcarpeta##*/}_PASS/${nombre_archivo}_PASS.vcf"
            fi
        done
    fi
done

