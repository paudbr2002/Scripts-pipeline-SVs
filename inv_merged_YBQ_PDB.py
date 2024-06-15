#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 16:35:32 2024

@author: yolanda/paula
"""

import argparse 
import pandas as pd
import os


def directory(raw_path):
    if not os.path.isdir(raw_path):
        raise argparse.ArgumentTypeError('"{}" is not an existing directory'.format(raw_path))
    return os.path.abspath(raw_path)

parser = argparse.ArgumentParser(
                    prog='merge_DEL',
                    description='Merge deletions from multiple structural variants callers.',
                    epilog='')

parser.add_argument("-i", "--SVs_file", help="input file containing SVs from all variant callers used")
parser.add_argument("-m", "--merged_file", help="input file containing merged SVs with bedtools merge")
parser.add_argument('-o', '--output_dir', help="output directory for merged results", default=os.path.curdir)

#path_MGD = "/home/yolanda/Documents/paula/beds_paula/BED_MERGED_PRUEBA.sorted.bed"
#path_merged = "/home/yolanda/Documents/paula/beds_paula/BED_MERGED_PRUEBA.collapsed.bed"
def swap_values(row):
    svtype = row[8]  # Assuming column 8 contains SVTYPE information
    if 'SVTYPE=INV2' in svtype:
        return [row[2], row[1]]
    else:
        return [row[1], row[2]]

def ParseMergedDels(path_merged, path_SVs, output_dir):
    # El input son los paths al bed merged y al bed con todas las deleciones. 

    bed_SVs = pd.read_csv(path_SVs, sep='\t', header=None)
    bed_SVs[[1, 2]] = bed_SVs.apply(swap_values, axis=1, result_type='expand')
    bed_merged = pd.read_csv(path_merged, sep='\t', header=None)
    #print(bed_SVs.head())
    bed_merged.columns = ["#chr", "start", "end", "ID_PROGRAM"]
    
    bed_merged["SVTYPE"] = "INV"
    
    bed_merged["DEL_ID"] = bed_merged['#chr'].astype(str) +"_" + bed_merged["start"].astype(str) +"_" + bed_merged['end'].astype(str)
    bed_merged["merged_len"] = bed_merged["end"] - bed_merged["start"] 
    
    
    
    bed_SVs["del_len"] = bed_SVs[2] - bed_SVs[1]
    bed_SVs["Individual_coord"] =  bed_SVs[11].astype(str) +":"  +bed_SVs[0].astype(str) +":" + bed_SVs[1].astype(str) +"-" + bed_SVs[2].astype(str)
    
    # crear un diccionario con el DEL_ID con la info de cada programa. 
    
    bed_merged = bed_merged.set_index('DEL_ID')
    df_inversions = pd.DataFrame(columns = ['DEL_ID', 'Program', 'ProgramID', 'Individual_coord','Merged_DEL_LEN', 'Individual_DEL_LEN']) # creamos el df vacío
    
    
    for ind in bed_merged.index:
        
        #print(bed_merged[3][ind])
        
        #ID_PROGRAM = bed_merged["ID_PROGRAM"][ind].split(",") # IDs de los programas que la soportan (por ejemplo: M1, D1, D2; D1 y D2 solapan con M1)
        #print(bed_merged.loc[ind, "ID_PROGRAM"])
        ID_PROGRAM = bed_merged.loc[ind, "ID_PROGRAM"].split(",")
        
        for program_id in ID_PROGRAM:          
            #list(bed_SVs[bed_SVs[3] == program_id]["Individual_coord"])[0]
            new_row = [ind, program_id[0], program_id, list(bed_SVs[bed_SVs[11] == program_id]["Individual_coord"])[0], bed_merged["merged_len"][ind], list(bed_SVs[bed_SVs[11] == program_id]["del_len"])[0]]
            df_inversions.loc[len(df_inversions.index)] = new_row
            
   # con la información que nos interese. Con el program ID podemos recuperar la info. Darle el formato adecuado para pasarlo por AnnotSV. 

    df_inversions["Program_DEL_LEN_TOTAL"]=df_inversions.groupby(["DEL_ID","Program"],as_index = False)["Individual_DEL_LEN"].transform('sum')
    df_inversions["Individual_Overlap"] = round(df_inversions["Individual_DEL_LEN"] /  df_inversions["Merged_DEL_LEN"] *100) 
    df_inversions.loc[df_inversions["Individual_Overlap"]<0.01, "Individual_Overlap"] = "<0.01"
    
    df_inversions["Program_Overlap_SUM"] = round(df_inversions["Program_DEL_LEN_TOTAL"] /  df_inversions["Merged_DEL_LEN"] *100) # ahora tendríamos que coger esta tabla, y hacer una línea por deleción 
    df_inversions.loc[df_inversions["Program_Overlap_SUM"]<0.01, "Program_Overlap_SUM"] = "<0.01"
    # cuando se pasa del 100% de overlap es porque hay varias deleciones que solapan en la misma deleción. 
    df_inversions["Individual_Overlap_ID"] = df_inversions["ProgramID"] +":"+ df_inversions["Individual_Overlap"].astype("str")
    
    
    df_inversions["GT"] = bed_SVs[10].apply(lambda x: x.split(":")[0])
    df_inversions["GT"] = df_inversions["GT"].apply(lambda x: x if x in ["0/1", "1/1"] else None)
    df_inversions["VAF"] = bed_SVs[10].apply(lambda x: x.split(":")[1] if df_inversions.loc[df_inversions["GT"].notnull(), "GT"].any() else None)

    
    #list_programs = [(x.split(",")) for x in bed_merged["ID_PROGRAM"]]
    
    n_programs = []
    programs = []
    
    for x in bed_merged["ID_PROGRAM"]: 
        ids_ind = (x.split(","))
        ind_prog = [x[0] for x in ids_ind]
        n_programs.append((len(set(ind_prog))))
        programs.append(','.join(set(ind_prog)))
    
    bed_merged["N_programs"] = n_programs
    bed_merged["Programs"] = programs
    bed_merged["GT_MD"] = df_inversions.groupby("DEL_ID")["GT"].first()
    
    for program in list(set(df_inversions["Program"])): # creamos las nuevas columnas. 
        #print(program + "_coord")
        bed_merged[str(program + "_overlap")] = df_inversions[df_inversions["Program"] == program].groupby("DEL_ID")[["Individual_Overlap_ID"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
        bed_merged[str(program + "_Total_overlap")] =df_inversions[df_inversions["Program"] == program].groupby("DEL_ID")[["Program_Overlap_SUM"]].first().astype("str") 
        bed_merged[str(program + "_overlap_INFO")] = bed_merged[str(program + "_Total_overlap")] + "=" + bed_merged[str(program + "_overlap")]
        #list(set(df_inversions[df_inversions["Program"] == program]))
        
        #print(program + "_overlap")
        bed_merged[str(program + "_coord")] = df_inversions[df_inversions["Program"] == program].groupby("DEL_ID")[["Individual_coord"]].agg(lambda x: x.astype(str).str.cat(sep=",") )
        if program == "G":
            bed_merged[str("GRIDSS_VAF")] = df_inversions[df_inversions["Program"] == program].groupby("DEL_ID")[["VAF"]].agg(lambda x: x.astype(str).str.cat(sep=",") ) 
    
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

