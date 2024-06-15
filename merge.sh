#!/bin/bash
####
#author:paula
####
#Archivo main para el modulo "merge" en el que se interseccionan las variantes detectadas por los tres programas
####
if [ $# -lt 1 ]; then
    echo "Uso: $0 <carpeta_de_entrada>"
    exit 1
fi

# Carpeta de entrada
carpeta_entrada="$1"

# Verificar que la carpeta de entrada existe
if [ ! -d "$carpeta_entrada" ]; then
    echo "La carpeta $carpeta_entrada no existe."
    exit 1
fi

carpeta_bedM="$carpeta_entrada/Manta/Manta_PASS"
carpeta_bedD="$carpeta_entrada/delly/delly_PASS"
carpeta_bedG="$carpeta_entrada/GRIDSS/GRIDSS_PASS"

tipos_SV=("DEL" "INS" "DUP" "INV" "TRA")

for tipo in "${tipos_SV[@]}"; do
    mkdir -p "$carpeta_entrada/$tipo"
done

carpetas=("$carpeta_bedM" "$carpeta_bedD" "$carpeta_bedG")

for carpeta in "${carpetas[@]}"; do
    for tipo in "${tipos_SV[@]}"; do
        archivos_bedtipo=$(find "$carpeta" -type f \( -name "*_${tipo}_PASS.bed" -o -name "*PASS_${tipo}.bed" \))
        if [ -n "$archivos_bedtipo" ]; then
            for archivo in $archivos_bedtipo; do
                python /mnt/tblab/paula/scripts_pipe/merge/split_gt.py -i "$archivo" -o "$carpeta_entrada/$tipo/"
            done
        else
            echo "No se encontraron archivos para $tipo en $carpeta."
        fi
    done
done

for tipo in "${tipos_SV[@]}"; do
    carpeta_tipo="$carpeta_entrada/$tipo"
    if [ "$tipo" == "DEL" ]; then
        for archivo in "$carpeta_tipo"/*; do
            if [ -f "$archivo" ]; then
                sorted_file="${archivo%.bed}.sorted.bed"
                bedtools sort -i "$archivo" > "$sorted_file"
                
                collapsed_file="${sorted_file%.sorted.bed}.collapsed.bed"
                bedtools merge -c 12 -o collapse -i "$sorted_file" > "$collapsed_file"
                
                # ruta al directorio de todos los scripts de merge
                python3 /mnt/tblab/paula/scripts_pipe/merge/del_merge.py -i "$sorted_file" -m "$collapsed_file" -o "$carpeta_tipo"
            fi 
        done
        output_file="$carpeta_tipo/DEL_merged_final.bed"
        final_files=($(find "$carpeta_tipo" -type f -name "*merged_results.bed"))
        cat "${final_files[0]}" > "$output_file"
        grep -v "^#" "${final_files[1]}" >> "$output_file"
	
        echo "Todo el merge de delecciones se encuentra en $output_file."
    elif [ "$tipo" == "DUP" ]; then
        for archivo in "$carpeta_tipo"/*; do
            if [ -f "$archivo" ]; then
                python3 /mnt/tblab/paula/scripts_pipe/merge/dup_merge.py -i "$archivo" -o "$carpeta_tipo"
                echo "Archivo de duplicaciones merged"
            fi
        done
        output_file="$carpeta_tipo/DUP_merged_final.bed"
        final_files=($(find "$carpeta_tipo" -type f -name "*merged_results.bed"))
        cat "${final_files[0]}" > "$output_file"
        grep -v "^#" "${final_files[1]}" >> "$output_file"
    elif [ "$tipo" == "INS" ]; then
        for archivo in "$carpeta_tipo"/*; do
            if [ -f "$archivo" ]; then
                python3 /mnt/tblab/paula/scripts_pipe/merge/ins_merge.py -i "$archivo" -o "$carpeta_tipo"
                echo "Archivo de inserciones merged"
            fi
        done
        output_file="$carpeta_tipo/INS_merged_final.bed"
        final_files=($(find "$carpeta_tipo" -type f -name "*merged_results.bed"))
        cat "${final_files[0]}" > "$output_file"
        grep -v "^#" "${final_files[1]}" >> "$output_file"
        
    elif [ "$tipo" == "TRA" ]; then
        for archivo in "$carpeta_tipo"/*; do
            if [ -f "$archivo" ]; then
                python3 /mnt/tblab/paula/scripts_pipe/merge/tra_merge.py -i "$archivo" -o "$carpeta_tipo"
                echo "Archivo de translocaciones merged"
            fi
        done
        output_file="$carpeta_tipo/TRA_merged_final.bed"
        final_files=($(find "$carpeta_tipo" -type f -name "*merged_results.bed"))
        cat "${final_files[0]}" > "$output_file"
        grep -v "^#" "${final_files[1]}" >> "$output_file"
        
    elif [ "$tipo" == "INV" ]; then
        for archivo in "$carpeta_tipo"/*; do
            if [ -f "$archivo" ]; then
                sorted_file="${archivo%.bed}.sorted.bed"
                bedtools sort -i "$archivo" > "$sorted_file"
                
                collapsed_file="${sorted_file%.sorted.bed}.collapsed.bed"
                bedtools merge -c 12 -o collapse -i "$sorted_file" > "$collapsed_file"
                
                python3 /mnt/tblab/paula/scripts_pipe/merge/inv_merge.py -i "$sorted_file" -m "$collapsed_file" -o "$carpeta_tipo"
            fi 
        done
        output_file="$carpeta_tipo/INV_merged_final.bed"
        final_files=($(find "$carpeta_tipo" -type f -name "*merged_results.bed"))
        cat "${final_files[0]}" > "$output_file"
        grep -v "^#" "${final_files[1]}" >> "$output_file"
	
        echo "Todo el merge de inversiones se encuentra en $output_file."
    
    fi
done


# Obtener el nombre de la carpeta de entrada
nombre_carpeta=$(basename "$carpeta_entrada")

# Crear el archivo de salida final basado en el nombre de la carpeta de entrada
output_final="$carpeta_entrada/${nombre_carpeta}_SV_merged.bed"

# Tipos de variantes estructurales
tipos_SV=("DEL" "INS" "DUP" "INV" "TRA")

# Inicializar el archivo de salida final vacío
> "$output_final"

# Variable para almacenar si el header ha sido agregado
header_added=false

# Procesar cada tipo de variante estructural
for tipo in "${tipos_SV[@]}"; do
    carpeta_tipo="$carpeta_entrada/$tipo"
    echo "Procesando carpeta: $carpeta_tipo"
    
    # Buscar archivos merged_final.bed en la carpeta del tipo actual
    resultados_merge=($(find "$carpeta_tipo" -type f -name "*merged_final.bed"))
    
    if [ ${#resultados_merge[@]} -gt 0 ]; then
        for archivo in "${resultados_merge[@]}"; do
            echo "Agregando archivo: $archivo"
            if [ "$header_added" = false ]; then
                # Agregar el header y el contenido del primer archivo
                cat "$archivo" >> "$output_final"
                header_added=true
            else
                # Agregar solo el contenido de los archivos restantes, excluyendo el header
                tail -n +2 "$archivo" >> "$output_final"
            fi
        done
    else
        echo "No se encontraron archivos merged_final.bed en $carpeta_tipo"
    fi
done

echo "Archivo final de SV mergeado se encuentra en $output_final."
