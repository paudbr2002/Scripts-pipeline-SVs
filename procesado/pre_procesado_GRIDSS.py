#!/usr/bin/python3
import pandas as pd
import argparse

def actualizar_archivo_bed(path_original, path_nuevo):
    try:
        print('Leyendo los archivos:', path_original, 'y', path_nuevo)
        
        # Lee los archivos BED en DataFrames de pandas
        archivo_original = pd.read_csv(path_original, sep='\t', header=None)
        archivo_nuevo = pd.read_csv(path_nuevo, sep='\t', header=None)

        # Itera sobre las filas del archivo original
        for indice, fila_original in archivo_original.iterrows():
            # Obtiene el ID del archivo original
            id_original = fila_original[3]

            # Verifica si el ID termina en "o"
            if id_original.endswith('o'):
                # Busca el ID correspondiente en el archivo nuevo
                fila_nuevo = archivo_nuevo[archivo_nuevo[9] == id_original]
                valor_columna_9 = fila_original[8].split(';') #valores columna info
                valor_col_format = fila_original[10].split(":")

                # Si se encuentra el ID en el archivo nuevo, actualiza las columnas en el archivo original
                if not fila_nuevo.empty:
                    if fila_nuevo.iloc[0, 6] != "TRA": 
                        if fila_nuevo.iloc[0, 6] == "INV":
                            valor_columna2 = fila_nuevo.iloc[0, 1]  # Columna 2 del archivo nuevo
                            valor_end = fila_nuevo.iloc[0, 5]  # Columna 6 del archivo nuevo
                            valor_chr2= fila_nuevo.iloc[0, 3]
                            tipoSV_nuevo = fila_nuevo.iloc[0, 6]
                            longitud= int(valor_columna2) - int(valor_end)
                            VF= int(valor_col_format [-1])
                            REF= int(valor_col_format [-7])
                            REFPAIR= int(valor_col_format [-6])
                            VAF= round(VF/(VF + REF + REFPAIR), 3)
                            if int(valor_end) < int(valor_columna2):
                                archivo_original.iloc[indice, 1] = valor_end #Columna 2 del archivo original
                                archivo_original.iloc[indice, 2] = valor_columna2  # Columna 11 del archivo original
                                valor_columna_9.insert(0, f"SVTYPE={tipoSV_nuevo}2")
                                valor_columna_9.insert(1, f"SVLEN={longitud}")
                                valor_columna_9.insert(2, f"VAF={VAF}")
                                valor_columna_9.insert(3, f"CHR2={valor_chr2}")
                                valor_columna_9.insert(4, f"END={valor_end}")
                                info_actualizada = ';'.join(valor_columna_9)

                            else:
                                valor_columna_9.insert(0, f"SVTYPE={tipoSV_nuevo}1")
                                valor_columna_9.insert(1, f"SVLEN={longitud}")
                                valor_columna_9.insert(2, f"VAF={VAF}")
                                valor_columna_9.insert(3, f"CHR2={valor_chr2}")
                                valor_columna_9.insert(4, f"END={valor_end}")
                                info_actualizada = ';'.join(valor_columna_9)
                                archivo_original.iloc[indice, 8] = info_actualizada
                                archivo_original.iloc[indice, 1] = valor_columna2  # Columna 2 del archivo original
                                archivo_original.iloc[indice, 2] = valor_end  # Columna 11 del archivo original
                      


                        else: 
                            valor_columna2 = fila_nuevo.iloc[0, 1]  # Columna 2 del archivo nuevo
                            valor_end = fila_nuevo.iloc[0, 5]  # Columna 6 del archivo nuevo
                            valor_chr2= fila_nuevo.iloc[0, 3]
                            tipoSV_nuevo = fila_nuevo.iloc[0, 6]
                            longitud= int(valor_columna2) - int(valor_end)
                            VF= int(valor_col_format [-1])
                            REF= int(valor_col_format [-7])
                            REFPAIR= int(valor_col_format [-6])
                            VAF= round(VF/(VF + REF + REFPAIR), 3)
                            valor_columna_9.insert(0, f"SVTYPE={tipoSV_nuevo}")
                            valor_columna_9.insert(1, f"SVLEN={longitud}")
                            valor_columna_9.insert(2, f"VAF={VAF}")
                            valor_columna_9.insert(3, f"CHR2={valor_chr2}")
                            valor_columna_9.insert(4, f"END={valor_end}")
                            info_actualizada = ';'.join(valor_columna_9)
                            archivo_original.iloc[indice, 8] = info_actualizada
                            archivo_original.iloc[indice, 1] = valor_columna2  # Columna 2 del archivo original
                            archivo_original.iloc[indice, 2] = valor_end  # Columna 11 del archivo original
                       
                    else: 
                        tipoSV_nuevo = fila_nuevo.iloc[0, 6]
                        longitud= fila_nuevo.iloc[0, 7]
                        VF= int(valor_col_format [-1])
                        REF= int(valor_col_format [-7])
                        REFPAIR= int(valor_col_format [-6])
                        VAF= round(VF/(VF + REF + REFPAIR), 3)
                        valor_columna_9.insert(0, f"SVTYPE={tipoSV_nuevo}")
                        valor_columna_9.insert(1, f"SVLEN={longitud}")
                        valor_columna_9.insert(2, f"VAF={VAF}")
                        valor_columna_9.insert(3, f"CHR2={valor_chr2}")
                        valor_columna_9.insert(4, f"END={valor_end}")
                        info_actualizada = ';'.join(valor_columna_9)
                        archivo_original.iloc[indice, 8] = info_actualizada


        # Sobrescribe el archivo original con las columnas actualizadas
        archivo_original.to_csv(path_original, sep='\t', index=False, header=False)

        print("Actualización completada.")

    except Exception as e:
        print("Error:", e)


def filtrar_archivo_bed(path):
    try:
        print('Leyendo el archivo:', path)
        
        # Lee el archivo BED en un DataFrame de pandas
        archivo_bed = pd.read_csv(path, sep='\t', header=None)
        
        # Filtra las filas que terminan en 'o' en la columna 4
        archivo_filtrado = archivo_bed[archivo_bed[3].str.endswith('o')]

        # Sobrescribe el archivo original con las filas filtradas
        archivo_filtrado.to_csv(path, sep='\t', index=False, header=False)

        print("Eliminados id que acaban en h y archivo sobrescrito exitosamente.")

    except Exception as e:
        print("Error:", e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Actualizar archivo BED")
    parser.add_argument("-o", "--original", type=str, help="Ruta al archivo BED original")
    parser.add_argument("-n", "--nuevo", type=str, help="Ruta al archivo BED nuevo generado a partir del script script_GridssVariantAnnotation.R")
    args = parser.parse_args()

    actualizar_archivo_bed(args.original, args.nuevo)
    filtrar_archivo_bed(args.original)
