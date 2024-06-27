#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

def procesar_archivo_vcf(directorio_entrada, directorio_salida):
    # Iterar sobre cada archivo en el directorio de entrada
    for archivo_entrada in os.listdir(directorio_entrada):
        # Comprobar si es un archivo .vcf
        if archivo_entrada.endswith(".vcf"):
            print(f"Procesando archivo: {archivo_entrada}")
            
            # Diccionario para almacenar los archivos de salida según SVtype
            archivos_salida = {}

            # Abrir el archivo de entrada
            with open(os.path.join(directorio_entrada, archivo_entrada), "r") as archivo:
                # Iterar sobre cada línea del archivo
                for linea in archivo:
                    # Si la línea comienza con "#" (encabezado), escribirla en todos los archivos de salida
                    if linea.startswith("#"):
                        for archivo_salida in archivos_salida.values():
                            archivo_salida.write(linea)
                    else:
                        # Extraer el valor de SVtype de la columna INFO
                        campos = linea.split("\t")
                        if len(campos) < 8:
                            print(f"Error: línea no tiene suficientes campos: {linea}")
                            continue
                        info = campos[7]
                        svtype = [item.split("=")[1] for item in info.split(";") if item.startswith("SVTYPE=")]
                        if not svtype:
                            print(f"Error: no se pudo encontrar SVTYPE en la línea: {linea}")
                            continue
                        svtype = svtype[0]

                        # Verificar si ya tenemos un archivo de salida para este tipo de SVtype
                        if svtype not in archivos_salida:
                            # Si no existe, crear un nuevo archivo de salida
                            nombre_archivo_salida = os.path.join(directorio_salida, f"{os.path.splitext(archivo_entrada)[0]}_{svtype}.vcf")
                            archivo_salida_nuevo = open(nombre_archivo_salida, "w")
                            # Escribir el encabezado en el nuevo archivo de salida
                            with open(os.path.join(directorio_entrada, archivo_entrada), "r") as archivo_entrada_lectura:
                                for linea_encabezado in archivo_entrada_lectura:
                                    if linea_encabezado.startswith("#"):
                                        archivo_salida_nuevo.write(linea_encabezado)
                            archivos_salida[svtype] = archivo_salida_nuevo

                        # Escribir la línea en el archivo de salida correspondiente al SVtype
                        archivos_salida[svtype].write(linea)

            # Cerrar todos los archivos de salida para este archivo de entrada
            for archivo_salida in archivos_salida.values():
                archivo_salida.close()

def main():
    parser = argparse.ArgumentParser(description="Procesa archivos VCF y los separa según el tipo de variante estructural (SVTYPE).")
    parser.add_argument("-i", "--entrada", help="Directorio de entrada que contiene los archivos .vcf")
    parser.add_argument("-o", "--salida", help="Directorio de salida para los archivos procesados. Si no se especifica, se usará el mismo directorio de entrada.", default=None)
    
    args = parser.parse_args()
    
    directorio_entrada = args.entrada
    directorio_salida = args.salida if args.salida else directorio_entrada
    
    procesar_archivo_vcf(directorio_entrada, directorio_salida)

if __name__ == "__main__":
    main()


