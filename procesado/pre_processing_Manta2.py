#!/usr/bin/python3
# -*- coding: utf-8 -*- 
import argparse 
import os
import pandas as pd
import re

# Definir argumentos para el script
parser = argparse.ArgumentParser(
    prog='preprocessing_manta',
    description='Ayuda a obtener el formato correcto de cada archivo .bed de Manta para los siguientes pasos del procesamiento de datos',
    epilog=''
)

parser.add_argument("-p", "--path", help="Ruta al directorio donde se encuentran todos los archivos .bed de las salidas .vcf de Manta")

def procesar_deleciones(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "DEL" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'DEL' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'DEL' en el nombre.")
        return

    # Procesa cada archivo .bed encontrado
    for bed_file in bed_files:
        print("Procesando:", bed_file)
        bed_file_path = os.path.join(path, bed_file)
        archivo_bed = pd.read_csv(bed_file_path, sep='\t', header=None)

        # Itera sobre cada fila del DataFrame
        for indice, fila in archivo_bed.iterrows():
            # Columna INFO
            valor_columna_9 = fila[8].split(';')  # Accede a la columna INFO del primer programa
                
            # Filtra solo los valores que comienzan con "END="
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]
            archivo_bed.loc[indice, 11] = "M" + str(indice) + "-DEL"
            
            if real_end:
                archivo_bed.loc[indice, 1] = fila[2]
                archivo_bed.loc[indice, 2] = int(real_end[0])

    # Sobrescribe el archivo original con las modificaciones realizadas
    print(archivo_bed.head())
    archivo_bed.to_csv(bed_file_path, sep='\t', index=False, header=False)
    print("Archivo de deleciones modificado y sobrescrito exitosamente.")

def procesar_inserciones(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "INS" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'INS' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'INS' en el nombre.")
        return

    for bed_file in bed_files:
        print("Procesando:", bed_file)
        bed_file_path = os.path.join(path, bed_file)
        archivo_bed = pd.read_csv(bed_file_path, sep='\t', header=None)

        # Itera sobre cada fila del DataFrame
        for indice, fila in archivo_bed.iterrows():
            # Columna INFO
            valor_columna_9 = fila[8].split(';')  # Accede a la columna INFO del primer programa
                
            # Filtra solo los valores que comienzan con "END="
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]
            archivo_bed.loc[indice, 11] = "M" + str(indice) + "-INS"
            
            if real_end:
                archivo_bed.loc[indice, 1] = fila[2]
                archivo_bed.loc[indice, 2] = int(real_end[0])

    # Sobrescribe el archivo original con las modificaciones realizadas
    print(archivo_bed.head())
    archivo_bed.to_csv(bed_file_path, sep='\t', index=False, header=False)
    print("Archivo de inserciones modificado y sobrescrito exitosamente.")

def procesar_inversiones(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "INV" en el nombre
    bed_files = [file for file in files if file.endswith('INV.bed')]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'INV' en el nombre.")
        return

    for bed_file in bed_files:
        print("Procesando:", bed_file)
        bed_file_path = os.path.join(path, bed_file)
        archivo_bed = pd.read_csv(bed_file_path, sep='\t', header=None)

        for indice, fila in archivo_bed.iterrows():
            valor_alt = fila[6].split(':')[1]
            valor_end = re.findall(r'\d+', valor_alt)[0]
            
            valor_columna_9 = fila[8].split(';')
            sv_len = int(fila[2]) - int(valor_end)
            if int(valor_end) < int(fila[2]):
                valor_columna_9 = [re.sub(r'^SVTYPE=BND$', 'SVTYPE=INV2', item) for item in valor_columna_9]
                        # Unir los elementos de la lista valor_columna_9 de nuevo en una cadena
                archivo_bed.loc[indice, 1] = valor_end
                archivo_bed.loc[indice, 2] = fila[2]
                archivo_bed.loc[indice, 11] = "M" + str(indice) + "-INV"
            else:
                valor_columna_9 = [re.sub(r'^SVTYPE=BND$', 'SVTYPE=INV1', item) for item in valor_columna_9]
 
                archivo_bed.loc[indice, 1] = fila[2]
                archivo_bed.loc[indice, 2] = int(valor_end)
                archivo_bed.loc[indice, 11] = "M" + str(indice) + "-INV"
            
            valor_columna_9.insert(1, f"SVLEN={sv_len}")
            valor_columna_9.insert(2, f"CHR2={fila[0]}")
            valor_columna_9.insert(3, f"END={valor_end}")
            info_actualizada = ';'.join(valor_columna_9)
            archivo_bed.iloc[indice, 8] = info_actualizada
  
            
    # Sobrescribe el archivo original con las modificaciones realizadas
    print(archivo_bed.head())
    archivo_bed.to_csv(bed_file_path, sep='\t', index=False, header=False)
    print("Archivo de inversiones modificado y sobrescrito exitosamente.")

def procesar_duplicaciones(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "DUP" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'DUP' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'DUP' en el nombre.")
        return

    for bed_file in bed_files:
        print("Procesando:", bed_file)
        bed_file_path = os.path.join(path, bed_file)
        archivo_bed = pd.read_csv(bed_file_path, sep='\t', header=None)

        for indice, fila in archivo_bed.iterrows():
            valor_columna_9 = fila[8].split(';')
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]
            archivo_bed.loc[indice, 1] = fila[2]
            archivo_bed.loc[indice, 2] = int(real_end[0])
            archivo_bed.loc[indice, 11] = "M" + str(indice) + "-DUP"

    print(archivo_bed.head())
    archivo_bed.to_csv(bed_file_path, sep='\t', index=False, header=False)
    print("Archivo de duplicaciones modificado y sobrescrito exitosamente.")

def procesar_translocaciones(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "TRA" en el nombre
    bed_files = [file for file in files if file.endswith('TRA.bed')]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'TRA' en el nombre.")
        return

    for bed_file in bed_files:
        print("Procesando:", bed_file)
        bed_file_path = os.path.join(path, bed_file)
        archivo_bed = pd.read_csv(bed_file_path, sep='\t', header=None)

        for indice, fila in archivo_bed.iterrows():
            valor_alt = fila[6].split(':')[1]
            valor_end = re.findall(r'\d+', valor_alt)[0]
            valor_split = fila[6].split(':')[0]
            match = re.search(r'chr([a-zA-Z0-9]+)', valor_split)
            cromosoma_info = match.group(1)
            
            valor_columna_9 = fila[8].split(';')
            valor_columna_9.insert(1, f"CHR2=chr{cromosoma_info}")
            valor_columna_9.insert(2, f"END={valor_end}")
            info_actualizada = ';'.join(valor_columna_9)
            archivo_bed.iloc[indice, 8] = info_actualizada
            archivo_bed.loc[indice, 11] = "M" + str(indice) + "-TRA"

    print(archivo_bed.head())
    archivo_bed.to_csv(bed_file_path, sep='\t', index=False, header=False)
    print("Archivo de translocaciones modificado y sobrescrito exitosamente.")

def main():
    args = None
    try:
        args = parser.parse_args()
        if not args.path:
            raise argparse.ArgumentError(None, "Por favor, proporciona la ruta al directorio usando el argumento -p o --path.")
    except argparse.ArgumentError as e:
        print(e)
        parser.print_help()
        return

    # Si se proporciona el argumento 'path', llama a las funciones de procesamiento
    procesar_deleciones(args.path)
    procesar_inserciones(args.path)
    procesar_inversiones(args.path)
    procesar_duplicaciones(args.path)
    procesar_translocaciones(args.path)

if __name__ == '__main__':
    main()
