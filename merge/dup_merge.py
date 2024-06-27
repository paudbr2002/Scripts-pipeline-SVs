#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:26:11 2024

@author: paula
"""

import os
import pandas as pd
import argparse

def directory(raw_path):
    if not os.path.isdir(raw_path):
        raise argparse.ArgumentTypeError('"{}" is not an existing directory'.format(raw_path))
    return os.path.abspath(raw_path)


def ParseMergedInsertions(path_input, output_dir):
    df = pd.read_csv(path_input, sep="\t", header=None)
    df.columns = ["chr", "start", "end", "ins_ids", "num", "REF", "ALT", "FILTER", "INFO", "FORMAT", "FORMAT_2", "ID_PROGRAM"]
    df["SV_ID"] = df['chr'].astype(str) + "_" + df["start"].astype(str) + "_" + (df['start']+1).astype(str)

    df_ins = pd.DataFrame(columns=['SV_ID', 'Program', 'ProgramID', 'Individual_coord', 'Merged_INS_LEN', 'Individual_INS_LEN'])
    df_ins["SV_ID"] = df['chr'].astype(str) + "_" + df["start"].astype(str) + "_" + (df['start']+1).astype(str)
    df_ins["Program"] = df["ID_PROGRAM"].apply(lambda x: x[0])
    df_ins["ProgramID"] = df["ID_PROGRAM"]
    df_ins["Individual_coord"] = "chr" + df["chr"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    df_ins["Individual_INS_LEN"] = df["end"] - df["start"]
    df_ins["Merged_INS_LEN"] = df_ins.groupby(["SV_ID"], as_index=False)["Individual_INS_LEN"].transform('max')
    df_ins["Program_DEL_LEN_TOTAL"] = df_ins.groupby(["SV_ID", "Program"], as_index=False)["Individual_INS_LEN"].transform('sum')
    df_ins["Individual_Overlap"] = round(df_ins["Individual_INS_LEN"] / df_ins["Merged_INS_LEN"] * 100)
    df_ins.loc[df_ins["Individual_Overlap"] < 0.01, "Individual_Overlap"] = "<0.01"
    df_ins["Program_Overlap_SUM"] = round(df_ins["Program_DEL_LEN_TOTAL"] / df_ins["Merged_INS_LEN"] * 100)
    df_ins.loc[df_ins["Program_Overlap_SUM"] < 0.01, "Program_Overlap_SUM"] = "<0.01"
    df_ins["Individual_Overlap_ID"] = df_ins["ProgramID"] + ":" + df_ins["Individual_Overlap"].astype("str")
    df_ins["GT"] = df["FORMAT_2"].apply(lambda x: x.split(":")[0])
    df_ins["GT"] = df_ins["GT"].apply(lambda x: x if x in ["0/1", "1/1"] else None)
    df_ins["VAF"] = df.apply(lambda row: f'{row["ID_PROGRAM"]}:{row["FORMAT_2"].split(":")[1]}' if row["FORMAT_2"].split(":")[0] not in ["0/1", "1/1"] else None, axis=1)
 

    merged = df.groupby("SV_ID")["ID_PROGRAM"].apply(lambda x: ','.join(x.unique())).reset_index()
    merged["SVTYPE"] = "DUP"
    merged["merged_len"] = df_ins.groupby(["SV_ID"], as_index=False)["Individual_INS_LEN"].transform('max')
    

    n_programs = []
    programs = []

    for x in merged["ID_PROGRAM"]:
        ids_ind = (x.split(","))
        inds_prog = [x[0] for x in ids_ind]
        n_programs.append((len(set(inds_prog))))
        programs.append(','.join(set(inds_prog)))

    merged["N_programs"] = n_programs
    merged["Programs"] = programs
    merged["GT_MD"] = list(df_ins.groupby("SV_ID")["GT"].first())
    if len(set(df_ins["Program"])) < 3:        
            if "G" not in set(df_ins["Program"]):
                # Si no está presente, agregamos "G" a la lista de programas únicos
                programs = list(set(df_ins["Program"])) + ["G"]
            elif "M" not in set(df_ins["Program"]):
                # Si no está presente, agregamos "G" a la lista de programas únicos
                programs = list(set(df_ins["Program"])) + ["M"]
            elif "D" not in set(df_ins["Program"]):
                    # Si no está presente, agregamos "G" a la lista de programas únicos
                    programs = list(set(df_ins["Program"])) + ["D"]
    else:
        # Si "G" ya está presente, simplemente usamos la lista existente de programas únicos
        programs = list(set(df_ins["Program"]))

    df_prog = pd.DataFrame()
    for program in programs:
        program_df = df_ins[df_ins["Program"] == program]
        df_prog_temp = program_df.groupby("SV_ID").agg({
            "Individual_Overlap_ID": lambda x: ','.join(x),
            "Program_Overlap_SUM": 'first',
            "Individual_coord": lambda x: ','.join(x),
            "VAF": lambda x: ','.join(x.dropna().astype(str)),
        })

        df_prog_temp = df_prog_temp.rename(columns={
            "Individual_Overlap_ID": str(program + "_overlap"),
            "Program_Overlap_SUM": str(program + "_Total_overlap"),
            "Individual_coord": str(program + "_coord"),
            "VAF": str(program + "_VAF")
        })

        df_prog_temp[str(program + "_overlap_INFO")] = df_prog_temp[str(program + "_Total_overlap")].astype(str) + "=" + df_prog_temp[str(program + "_overlap")]

        if df_prog.empty:
            df_prog = df_prog_temp
        else:
            df_prog = df_prog.join(df_prog_temp, how='outer')

    merged = merged.merge(df_prog, on='SV_ID', how='left')

    merged[['#chr', 'start', 'end']] = merged['SV_ID'].str.split('_', expand=True)
    merged.drop(columns=['SV_ID'], inplace=True)
  

    cols = ['#chr', 'start', 'end'] + [col for col in merged.columns if col not in ['#chr', 'start', 'end']]
    merged = merged[cols]

    if 'G_VAF' in merged.columns:
        merged['GRIDSS_VAF'] = merged['G_VAF']
    
    for program in list(set(df_ins["Program"])):
        merged.drop(columns=[str(program)+"_VAF"], inplace=True)


     
    # Define the order of columns
    column_order = [
        "#chr", "start", "end", "ID_PROGRAM", "SVTYPE", "merged_len", "N_programs", "Programs", "GT_MD",
        "G_overlap", "G_Total_overlap", "G_overlap_INFO", "G_coord", "GRIDSS_VAF",
        "D_overlap", "D_Total_overlap", "D_overlap_INFO", "D_coord",
        "M_overlap", "M_Total_overlap", "M_overlap_INFO", "M_coord"
    ]

    # Ensure all columns are present
    for col in column_order:
        if col not in merged.columns:
            merged[col] = pd.NA
    
    merged = merged[column_order]

    input_base = os.path.basename(path_input)
    output_file = os.path.join(output_dir, os.path.splitext(input_base)[0] + "_merged_results.bed")

    merged.to_csv(output_file, index=False, sep="\t")




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="input BED file containing insertion variants", required=True)
    parser.add_argument("-o", "--output_dir", help="output directory for merged insertions", type=directory)
    args = parser.parse_args()

    try:
        print("Input file: " + args.input_file)
        print("Output directory: " + args.output_dir)
        ParseMergedInsertions(args.input_file, args.output_dir)
    except TypeError:
        print("Invalid argument.")

if __name__ == '__main__':
    main()
