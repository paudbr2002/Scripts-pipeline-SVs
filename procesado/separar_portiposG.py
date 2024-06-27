#!/usr/bin/python3
import pandas as pd
import argparse
import os

def separar_archivo_bed(path, output_dir):
    try:
        # Lee el archivo BED en un DataFrame de pandas
        archivo_bed = pd.read_csv(path, sep='\t', header=None)

        for indice, fila in archivo_bed.iterrows():
            if len(fila) <= 8:
                print(f"Fila {indice} no tiene suficientes columnas: {fila}")
                continue

            info = fila[8]
            if not isinstance(info, str):
                print(f"Fila {indice} tiene un tipo incorrecto en la columna 9: {type(info)}")
                continue
            
            svtype_list = [item.split("=")[1] for item in info.split(";") if item.startswith("SVTYPE=")]
            
            if svtype_list:
                svtype = svtype_list[0]
                
                if svtype in ["INV1", "INV2"]:
                    svtype = "INV"


                base_name = path.split('/')[-1].split('.')[0]
                # Crea el nombre de archivo en el directorio de salida especificado
                filename = os.path.join(output_dir, f"{base_name}_{svtype}.bed")

                with open(filename, "a") as f:
                    fila = list(fila)
                    num_linea = sum(1 for _ in open(filename)) + 1
                    fila.append(f"G{num_linea}-{svtype}")
                    f.write("\t".join(map(str, fila)) + "\n")

            else:
                print(f"No se encontró SVTYPE en la fila {indice}.")

        print("Archivos separados creados exitosamente.")
    
    except Exception as e:
        print("Error:", e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Separar archivo BED por tipo de SV y agregar ID en la última columna")
    parser.add_argument("-i", "--input", type=str, required=True, help="Ruta al archivo BED original obtenido tras el formateo con pre_procesado_GRIDSS.py")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Directorio donde se guardarán los archivos de salida")
    args = parser.parse_args()

    # Verifica si el directorio de salida existe, si no, lo crea
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    separar_archivo_bed(args.input, args.output_dir)

