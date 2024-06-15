#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
@author: paula/yolanda

"""
#Importar librerías 
import argparse #modulo para incluir comandos en el script
import pandas as pd #modulo para manipular datos 
import os #modulo para interactuar con el sistema operativo y obtener directorios

#Función para asegurar que el argumento proporcionado es un directorio válido
def directory(raw_path):
    if not os.path.isdir(raw_path):
        raise argparse.ArgumentTypeError('"{}" is not an existing directory'.format(raw_path))
    return os.path.abspath(raw_path)

parser = argparse.ArgumentParser(
                    prog='merge_DEL',
                    description='Merge deletions from multiple structural variants callers./ Fusionar las deleciones detectadas por varios programas de detección de SVs.',
                    epilog='')

parser.add_argument("-i", "--SVs_file", help="input file containing SVs from all variant callers used/archivo que contiene las SVs ordenadas de todos los programas")
parser.add_argument("-m", "--merged_file", help="input file containing merged SVs with bedtools merge/archivo que contiene las variantes tras utilizar bedtools merge")
parser.add_argument('-o', '--output_dir', help="output directory for merged results/directorio para dar el archivo fusionado", default=os.path.curdir)


def ParseMergedDels(path_merged, path_SVs, output_dir):
    # Leer archivos de entrada como DataFrames de Pandas
    bed_SVs = pd.read_csv(path_SVs, sep='\t', header=None)
    bed_merged = pd.read_csv(path_merged, sep='\t', header=None)
    
    # Asignar nombres a las columnas del DataFrame bed_merged
    bed_merged.columns = ["#chr", "start", "end", "ID_PROGRAM"]
    
    # Añadir columnas adicionales para bed_merged
    bed_merged["SVTYPE"] = "DEL" #el tipo de SV
    bed_merged["DEL_ID"] = bed_merged['#chr'].astype(str) +"_" + bed_merged["start"].astype(str) +"_" + bed_merged['end'].astype(str) #ID basado en las coordenadas
    bed_merged["merged_len"] = bed_merged["end"] - bed_merged["start"]#longitud delecion
    
    
    # Calcular la longitud de las deleciones en bed_SVs
    bed_SVs["del_len"] = bed_SVs[2] - bed_SVs[1]
    bed_SVs["Individual_coord"] =  bed_SVs[11].astype(str) +":"  +bed_SVs[0].astype(str) +":" + bed_SVs[1].astype(str) +"-" + bed_SVs[2].astype(str)
    bed_merged = bed_merged.set_index('DEL_ID')
    df_deletions = pd.DataFrame(columns = ['DEL_ID', 'Program', 'ProgramID', 'Individual_coord','Merged_DEL_LEN', 'Individual_DEL_LEN'])

    for ind in bed_merged.index:
        ID_PROGRAM = bed_merged.loc[ind, "ID_PROGRAM"].split(",")
        for program_id in ID_PROGRAM:
            new_row = [ind, program_id[0], program_id, list(bed_SVs[bed_SVs[11] == program_id]["Individual_coord"])[0], bed_merged["merged_len"][ind], list(bed_SVs[bed_SVs[11] == program_id]["del_len"])[0]]
            df_deletions.loc[len(df_deletions.index)] = new_row
    
    # Calcular la longitud total de deleciones por programa
    df_deletions["Program_DEL_LEN_TOTAL"]=df_deletions.groupby(["DEL_ID","Program"],as_index = False)["Individual_DEL_LEN"].transform('sum')
    
    # Calcular el solapamiento individual y de programa   
    df_deletions["Individual_Overlap"] = round(df_deletions["Individual_DEL_LEN"] /  df_deletions["Merged_DEL_LEN"] *100)   
    df_deletions.loc[df_deletions["Individual_Overlap"]<0.01, "Individual_Overlap"] = "<0.01"
    df_deletions["Program_Overlap_SUM"] = round(df_deletions["Program_DEL_LEN_TOTAL"] /  df_deletions["Merged_DEL_LEN"] *100)
    df_deletions.loc[df_deletions["Program_Overlap_SUM"]<0.01, "Program_Overlap_SUM"] = "<0.01"
    df_deletions["Individual_Overlap_ID"] = df_deletions["ProgramID"] +":"+ df_deletions["Individual_Overlap"].astype("str")
    df_deletions["GT"] = bed_SVs[10].apply(lambda x: x.split(":")[0])
 
    df_deletions["GT"] = df_deletions["GT"].apply(lambda x: x if x in ["0/1", "1/1"] else None)
    df_deletions["VAF"] = bed_SVs[10].apply(lambda x: x.split(":")[1] if df_deletions.loc[df_deletions["GT"].notnull(), "GT"].any() else None)


    n_programs = []
    programs = []
    
    for x in bed_merged["ID_PROGRAM"]:
        ids_ind = (x.split(","))
        ind_prog = [x[0] for x in ids_ind]
        n_programs.append((len(set(ind_prog))))
        programs.append(','.join(set(ind_prog)))
    
    bed_merged["N_programs"] = n_programs
    bed_merged["Programs"] = programs
    bed_merged["GT_MD"] = df_deletions.groupby("DEL_ID")["GT"].first()
    for program in list(set(df_deletions["Program"])):
        bed_merged[str(program + "_overlap")] = df_deletions[df_deletions["Program"] == program].groupby("DEL_ID")[["Individual_Overlap_ID"]].agg(lambda x: x.astype(str).str.cat(sep=","))
        bed_merged[str(program + "_Total_overlap")] = df_deletions[df_deletions["Program"] == program].groupby("DEL_ID")[["Program_Overlap_SUM"]].first().astype("str")
        bed_merged[str(program + "_overlap_INFO")] = bed_merged[str(program + "_Total_overlap")] + "=" + bed_merged[str(program + "_overlap")]
        bed_merged[str(program + "_coord")] = df_deletions[df_deletions["Program"] == program].groupby("DEL_ID")[["Individual_coord"]].agg(lambda x: x.astype(str).str.cat(sep=","))
        if program == "G":
            bed_merged[str("GRIDSS_VAF")] = df_deletions[df_deletions["Program"] == program].groupby("DEL_ID")[["VAF"]].agg(lambda x: x.astype(str).str.cat(sep=","))
    
    # Define the order of columns
    column_order = [
        "#chr", "start", "end", "ID_PROGRAM", "SVTYPE", "merged_len", "N_programs", "Programs", "GT_MD",
        "G_overlap", "G_Total_overlap", "G_overlap_INFO", "G_coord", "GRIDSS_VAF",
        "D_overlap", "D_Total_overlap", "D_overlap_INFO", "D_coord",
        "M_overlap", "M_Total_overlap", "M_overlap_INFO", "M_coord"
    ]

    # Ensure all columns are present
    for col in column_order:
        if col not in bed_merged.columns:
            bed_merged[col] = pd.NA
    
    bed_merged = bed_merged[column_order]
    
    nombre_base= os.path.splitext(os.path.basename(path_SVs))[0]
    nombre_final=f"{nombre_base}_merged_results.bed"

    output_file = os.path.join(output_dir, nombre_final)
    bed_merged.to_csv(output_file, index=False, sep="\t")
    print(f"Output written to {output_file}")


def main():
    args = parser.parse_args()
    output_dir = args.output_dir

    try:
        print("Working dir: " + output_dir)
        print("SVs file after cat all individual bed files: " + args.SVs_file)
        print("Bed file after bedtools merge all programs: " + args.merged_file)

        path_merged = args.merged_file
        path_SVs = args.SVs_file

    except TypeError:
        print("Invalid argument.")

    ParseMergedDels(path_merged, path_SVs, output_dir)


if __name__ == '__main__':
    main()
