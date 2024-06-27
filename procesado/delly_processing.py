#!/usr/bin/python3 
#Importar 
import argparse 
import os
import pandas as pd

## definimos argumentos script
parser = argparse.ArgumentParser(
                    prog='preprocessing_delly',
                    description='Helps obtaining the correct format of each delly.bed file for next steps of data processing',
                    epilog='')

parser.add_argument("-p", "--path", help="path to directory where there are located all the .bed from delly output .vcf files")
##

def deletion_processing (path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "DEL" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'DEL' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'DEL' en el nombre.")
        return
    # Procesa cada archivo .bed encontrado
    for bed_file in bed_files:
        print("Processing:", bed_file)
        
        # Carga el archivo .bed en un DataFrame
        archivo_bed = pd.read_csv(os.path.join(path, bed_file), sep='\t', header=None)
        
        for indice, fila in archivo_bed.iterrows():
            
            # Columnas INFO
            valor_columna_9 = fila[8].split(';')  # Accede a col INFO del primer programa
                
                # Filtra solo los valores que comienzan con "END="
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]       
            if real_end:
                indice_svtype = next((i for i, item in enumerate(valor_columna_9) if item.startswith('SVTYPE')), None)
                if indice_svtype is not None:
                    sv_length = int(fila[2])- int(real_end[0])    
                    # Insertar "SVLEN" después de "SVTYPE" en la lista valor_columna_9
                    valor_columna_9.insert(indice_svtype + 1, f"SVLEN={sv_length}")  # Insertar un espacio para SVLEN
                    # Unir los elementos de la lista valor_columna_9 de nuevo en una cadena
                    info_actualizada = ';'.join(valor_columna_9)
                    archivo_bed.loc[indice, 8] = info_actualizada
                    archivo_bed.loc[indice, 1] = fila[2]
                    archivo_bed.loc[indice, 2] = real_end[0]
                    archivo_bed.loc[indice, 11] = "D" + str(indice) +"-DEL"                       
    # Sobrescribe el archivo original con las modificaciones realizadas
    print(archivo_bed.head())
    archivo_bed.to_csv(os.path.join(path, bed_file), sep='\t', index=False, header=False)
    print("Archivo de delecciones modificado y sobrescrito exitosamente.")

def insertion_processing(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "DEL" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'INS' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'INS' en el nombre.")
        return
    for bed_file in bed_files:
        print("Processing:", bed_file)        
        # Carga el archivo .bed en un DataFrame
        archivo_bed = pd.read_csv(os.path.join(path, bed_file), sep='\t', header=None)
            # Itera sobre cada fila del DataFrame
        for indice, fila in archivo_bed.iterrows():

            # Columnas INFO
            valor_columna_9 = fila[8].split(';')  # Accede a col INFO del primer programa
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]        
            if real_end:
                indice_svtype = next((i for i, item in enumerate(valor_columna_9) if item.startswith('SVTYPE')), None)
                if indice_svtype is not None:
                    sv_length = (int(real_end[0]) - int(fila[2]))    
                    # Insertar "SVLEN" después de "SVTYPE" en la lista valor_columna_9
                    valor_columna_9.insert(indice_svtype + 1, f"SVLEN={sv_length}")  # Insertar un espacio para SVLEN
                    # Unir los elementos de la lista valor_columna_9 de nuevo en una cadena
                    info_actualizada = ';'.join(valor_columna_9)
                    archivo_bed.loc[indice, 8] = info_actualizada
                    archivo_bed.loc[indice, 1] = fila[2]
                    archivo_bed.loc[indice, 2] = real_end[0]
                    archivo_bed.loc[indice, 11] = "D" + str(indice) + "-INS"
    print(archivo_bed.head())
    archivo_bed.to_csv(os.path.join(path, bed_file), sep='\t', index=False, header=False)
    print("Archivo de inserciones modificado y sobrescrito exitosamente.")

def inversion_processing(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "DEL" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'INV' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'INV' en el nombre.")
        return

    # Procesa cada archivo .bed encontrado
    for bed_file in bed_files:
        print("Processing: ", bed_file)
        archivo_bed = pd.read_csv(os.path.join(path, bed_file), sep='\t', header=None)

            # Itera sobre cada fila del DataFrame
        for indice, fila in archivo_bed.iterrows():

            # Columnas INFO
            valor_columna_9 = fila[8].split(';')  # Accede a col INFO del primer programa
                
                # Filtra solo los valores que comienzan con "END="
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]
            
            
            if real_end:
                if int(real_end[0]) < int(fila[2]):
                    indice_svtype = next((i for i, item in enumerate(valor_columna_9) if item.startswith('SVTYPE')), None)
                    if indice_svtype is not None:
                            sv_length = int(fila[2])- int(real_end[0])    
                            # Insertar "SVLEN" después de "SVTYPE" en la lista valor_columna_9
                            valor_columna_9.insert(indice_svtype + 1, f"SVLEN={sv_length}")  # Insertar un espacio para SVLEN
                            # Unir los elementos de la lista valor_columna_9 de nuevo en una cadena
                            valor_columna_9[indice_svtype] += '2'
                            info_actualizada = ';'.join(valor_columna_9)
                            archivo_bed.loc[indice, 8] = info_actualizada
                            archivo_bed.loc[indice, 1] = real_end
                            archivo_bed.loc[indice, 2] = fila[2]
                            archivo_bed.loc[indice, 11] = "D" + str(indice) + "-INV2"
                else:
                    indice_svtype = next((i for i, item in enumerate(valor_columna_9) if item.startswith('SVTYPE')), None)
                    if indice_svtype is not None:
                        sv_length = int(fila[2])- int(real_end[0])    
                        # Insertar "SVLEN" después de "SVTYPE" en la lista valor_columna_9
                        valor_columna_9.insert(indice_svtype + 1, f"SVLEN={sv_length}")  # Insertar un espacio para SVLEN
                        # Unir los elementos de la lista valor_columna_9 de nuevo en una cadena
                        valor_columna_9[indice_svtype] += '1'
                        info_actualizada = ';'.join(valor_columna_9)
                        archivo_bed.loc[indice, 8] = info_actualizada
                        archivo_bed.loc[indice, 1] = fila[2]
                        archivo_bed.loc[indice, 2] = real_end[0]
                        archivo_bed.loc[indice, 11] = "D" + str(indice) + "-INV1"
    print(archivo_bed.head())
    archivo_bed.to_csv(os.path.join(path, bed_file), sep='\t', index=False, header=False)
    print("Archivo de inversiones modificado y sobrescrito exitosamente.")                   

def duplication_processing(path):
    files = os.listdir(path)
    
    # Filtra solo los archivos .bed que contienen "DEL" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'DUP' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'DUP' en el nombre.")
        return
    for bed_file in bed_files:
        print("Processing:", bed_file)        
        # Carga el archivo .bed en un DataFrame
        archivo_bed = pd.read_csv(os.path.join(path, bed_file), sep='\t', header=None)
        for indice, fila in archivo_bed.iterrows():
        # Columnas INFO
            valor_columna_9 = fila[8].split(';')  # Accede a col INFO del primer programa           
                # Filtra solo los valores que comienzan con "END="
            real_end = [item.split('=')[1] for item in valor_columna_9 if item.startswith('END=')]        
            if real_end:
                indice_svtype = next((i for i, item in enumerate(valor_columna_9) if item.startswith('SVTYPE')), None)
                if indice_svtype is not None:
                    sv_length = (int(real_end[0]) - int(fila[2]))    
                    # Insertar "SVLEN" después de "SVTYPE" en la lista valor_columna_9
                    valor_columna_9.insert(indice_svtype + 1, f"SVLEN={sv_length}")  # Insertar un espacio para SVLEN
                    # Unir los elementos de la lista valor_columna_9 de nuevo en una cadena
                    info_actualizada = ';'.join(valor_columna_9)
                    archivo_bed.loc[indice, 8] = info_actualizada
                archivo_bed.loc[indice, 1] = fila[2]
                archivo_bed.loc[indice, 2] = real_end[0]
                archivo_bed.loc[indice, 11] = "D" + str(indice) + "-DUP"
    print(archivo_bed.head())
    archivo_bed.to_csv(os.path.join(path, bed_file), sep='\t', index=False, header=False)
    print("Archivo de duplicaciones modificado y sobrescrito exitosamente.")      

def translocation_processing(path):
    files = os.listdir(path)    
    # Filtra solo los archivos .bed que contienen "DEL" en el nombre
    bed_files = [file for file in files if file.endswith('.bed') and 'TRA' in file]
    if not bed_files:
        print("No se encontraron archivos .bed que contengan 'TRA' en el nombre.")
        return
    for bed_file in bed_files:
        print("Processing:", bed_file)         
        archivo_bed = pd.read_csv(os.path.join(path, bed_file), sep='\t', header=None)
        # Itera sobre cada fila del DataFrame
        for indice, fila in archivo_bed.iterrows():
            archivo_bed.loc[indice, 11] = "D" + str(indice) + "-TRA"
    print(archivo_bed.head())
    archivo_bed.to_csv(os.path.join(path, bed_file), sep='\t', index=False, header=False)
    print("Archivo de translocaciones modificado y sobrescrito exitosamente.")                


def main():
    args = None
    try:
        args = parser.parse_args()
        # Si no se proporciona el argumento 'path', muestra la ayuda y sale del script
        if not args.path:
            raise argparse.ArgumentError(None, "Please provide the path to the directory using -p or --path argument.")
    except argparse.ArgumentError as e:
        print(e)
        parser.print_help()
        return

    # Si se proporciona el argumento 'path', llama a las funciones de procesamiento
    deletion_processing(args.path)
    insertion_processing(args.path)
    duplication_processing(args.path)
    inversion_processing(args.path)
    translocation_processing(args.path)


if __name__ == '__main__':
    main()
