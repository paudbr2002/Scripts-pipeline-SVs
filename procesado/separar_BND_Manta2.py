#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import re
import argparse
import os

def procesar_archivo_bed(path):
    print('Leyendo el archivo:', path)
    try:
        # Listar los archivos en el directorio especificado
        files = os.listdir(path)
    
        # Filtra solo los archivos .bed que contienen "BND" en el nombre
        bed_files = [file for file in files if file.endswith('.bed') and 'BND' in file]
        
        if not bed_files:
            print("No se encontraron archivos .bed que contengan 'BND' en el nombre.")
            return

        # Procesa cada archivo .bed encontrado
        for bed_file in bed_files:
            bed_file_path = os.path.join(path, bed_file)
            print("Processing:", bed_file_path)
            
            archivo_bed = pd.read_csv(bed_file_path, sep='\t', header=None)
            
            # Crear DataFrames vacíos con las mismas columnas que el archivo original
            df_cromosomasigual = pd.DataFrame(columns=archivo_bed.columns)
            df_cromosomasdif = pd.DataFrame(columns=archivo_bed.columns)
            df_unk= pd.DataFrame(columns=archivo_bed.columns)

            # Itera sobre cada fila del DataFrame
            for indice, fila in archivo_bed.iterrows():
                try:
                    # Columnas INFO
                    valor_columna_9 = fila[6]  # Accede a la columna ALT
                    valor_chr1 = fila[0].split("chr")[1]  # Accede a la primera columna (chromosome 1)

                    if "decoy" in valor_columna_9 or "alt" in valor_columna_9:
                        df_contigs = pd.concat([df_unk, fila.to_frame().T])

                    else: 

                        # Buscar la primera coincidencia en la cadena de la columna INFO
                        match = re.search(r'chr([a-zA-Z0-9]+)', valor_columna_9)
                        if match:
                            cromosoma_info = match.group(1)
                        

                            if cromosoma_info == valor_chr1:
                                # Si coinciden, agregar la fila al DataFrame de cromosomas iguales
                                df_cromosomasigual = pd.concat([df_cromosomasigual, fila.to_frame().T])
                            elif cromosoma_info == "Un":
                                df_unk = pd.concat([df_unk, fila.to_frame().T])
                            else:
                                # Si no coinciden, agregar la fila al DataFrame de cromosomas diferentes
                                df_cromosomasdif = pd.concat([df_cromosomasdif, fila.to_frame().T])
                        else:
                            print(f"Error: No se encontró una coincidencia en la columna INFO en la fila {indice}")
                except Exception as fila_error:
                    print(f"Error procesando la fila {indice}: {fila_error}")

            # Generar nuevas rutas de archivo cambiando "BND" por "INV" y "TRA"
            nuevo_path_inv = bed_file_path.replace("BND", "INV")
            nuevo_path_tra = bed_file_path.replace("BND", "TRA")
            nuevo_path_unk = bed_file_path.replace("BND", "BND_CONTIGS")

            # Guardar los DataFrames resultantes en nuevos archivos
            df_cromosomasigual.to_csv(nuevo_path_inv, sep='\t', index=False, header=False)
            df_cromosomasdif.to_csv(nuevo_path_tra, sep='\t', index=False, header=False)
            df_unk.to_csv(nuevo_path_unk, sep='\t', index=False, header=False)
            print("Archivo separado en inversiones y translocaciones exitosamente")

    except Exception as e:
        print("Error:", e)

def main():
    parser = argparse.ArgumentParser(description="Procesa archivos BED y los separa según si son INV o TRA.")
    parser.add_argument("-p", "--path", required=True, help="Ruta al directorio con archivos BED de los BNDs")
    
    args = parser.parse_args()
    
    procesar_archivo_bed(args.path)

if __name__ == "__main__":
    main()
