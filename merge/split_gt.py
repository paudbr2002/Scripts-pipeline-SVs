#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 16:04:44 2024

@author: paula/yolanda

Script para separar las variantes del bed PASS de delly y de Manta en 0/1 y 1/1, de forma que estemos comparando con bedtools merge solo los que tienen el mismo gt
"""
import argparse 
import os


def directory(raw_path):
    if not os.path.isdir(raw_path):
        raise argparse.ArgumentTypeError('"{}" is not an existing directory'.format(raw_path))
    return os.path.abspath(raw_path)

parser = argparse.ArgumentParser(
                    prog='gt_splitter',
                    description='Separate variants into homozygosity and heterozygosity.',
                    epilog='')

parser.add_argument("-i", "--SVs_file", help="path of the input file containing SVs")
#parser.add_argument('--working-dir', type=directory, default=os.path.curdir) 
parser.add_argument('-o', '--output_dir', type=directory, default=os.path.curdir) 

def Split_GT(SVs_file, outdir):
    file_name = os.path.basename(SVs_file).split(".")[0]
    SV_TYPE = os.path.basename(outdir)
    with open(SVs_file, "r") as archivo:
        for linea in archivo:            
            gt = linea.split("\t")[10].split(":")[0] #sacamos el gt
            if str(gt) == "1/1": 
                with open(str(outdir + "/" + SV_TYPE + "_homo.bed"),"a") as archivo_salida:
                    archivo_salida.write(linea)
                    archivo_salida.close()
            elif str(gt) == "0/1":
                with open(str(outdir + "/" + SV_TYPE + "_het.bed"),"a") as archivo_salida:
                    archivo_salida.write(linea)    
                    archivo_salida.close()
            elif str(gt) == ".":
                    with open(str(outdir + "/" + SV_TYPE + "_homo.bed"),"a") as archivo_homo:
                        archivo_homo.write(linea)
                        archivo_homo.close()
                    with open(str(outdir + "/" + SV_TYPE + "_het.bed"),"a") as archivo_het:
                        archivo_het.write(linea)
                        archivo_het.close()
    print("Archivo separado en gt exitosamente")
    

def main():
    args = parser.parse_args()
    
    try:
        print("Output dir: " + args.output_dir)
        print("SVs file: " + args.SVs_file)
                
        path_SVs = args.SVs_file

    except TypeError:
        print("Invalid argument.")
    
    Split_GT(path_SVs, args.output_dir)


if __name__ == '__main__':
    main()

