#!/bin/bash

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

# Carpeta de entrada
carpeta_entrada="$1"

/mnt/tblab/paula/scripts_pipe/procesado/procesado_inicial.sh "$carpeta_entrada"

echo "Archivos separados en PASS y no PASS"

carpeta_bedM="$carpeta_entrada/Manta/Manta_PASS"

python3 /mnt/tblab/paula/scripts_pipe/procesado/separar_portiposM.py -i "$carpeta_bedM"

echo "Archivos de Manta separados por tipos"

/mnt/tblab/paula/scripts_pipe/procesado/vcf2bed.sh "$carpeta_entrada"

echo "Archivos transformados en .bed correctamente"

carpeta_bedD="$carpeta_entrada/delly/delly_PASS"

python3 /mnt/tblab/paula/scripts_pipe/procesado/delly_processing.py -p "$carpeta_bedD"

carpeta_bedM_2="$carpeta_entrada/Manta/Manta_PASS/"
python3 /mnt/tblab/paula/scripts_pipe/procesado/separar_BND_Manta.py -p "$carpeta_bedM_2"

python3 /mnt/tblab/paula/scripts_pipe/procesado/pre_processing_Manta2.py -p "$carpeta_bedM"

carpeta_bedG="$carpeta_entrada/GRIDSS/GRIDSS_PASS/"

archivo_vcf=$(find "$carpeta_bedG" -type f -name "*.vcf")
archivo_bed_original=$(find "$carpeta_bedG" -type f -name "*.bed")

sed -i '1i##fileformat=VCFv4.2' "$archivo_vcf"

Rscript /mnt/tblab/paula/scripts_pipe/procesado/script_GridssVariantAnnotation.R "$archivo_vcf"

archivo_bed_nuevo=$(find "$carpeta_bedG" -type f -name "*_portipos.bed")

python3 /mnt/tblab/paula/scripts_pipe/procesado/pre_procesado_GRIDSS.py -o "$archivo_bed_original" -n "$archivo_bed_nuevo"

python3 /mnt/tblab/paula/scripts_pipe/procesado/separar_portiposG.py -i "$archivo_bed_original" -o "$carpeta_bedG"


/mnt/tblab/paula/scripts_pipe/merge/merge.sh "$carpeta_entrada"


