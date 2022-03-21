#!/usr/bin/python

import sys

def read_inversion(inputFile):
    with open(inputFile,"r") as inv:
        inversion = inv.readlines()
        inv_identify = []
        for line in inversion:
            if line[0] != "#":
                inv_info = line.split()
                chr = inv_info[0]
                pos_start = inv_info[1]
                pos_info = inv_info[7].split(";")
                pos_end = pos_info[3][4:]
                inv_identify.append(chr)
                inv_identify.append(pos_start)
                inv_identify.append(pos_end)
    inv.close()
    return inv_identify

def merge_inv_TE(inv_input,inv_TE):
    with open(inv_TE,"r") as inv_TE:
        TE_input = inv_TE.readlines()
        TE_qvchong = []
        feature = []
        for i in TE_input:
            if i not in TE_qvchong:
                TE_qvchong.append(i)
        for i in TE_qvchong:
            line = i.split()
            chr_pos = str(line[0]+line[1])
            feature.append(chr_pos)
        for i in range(0,len(inv_input)//3):
            inv_chr_pos = str(inv_input[3*i]+inv_input[3*i+1])
            if inv_chr_pos not in feature:
                print(inv_input[3*i],inv_input[3*i+1],inv_input[3*i+2],sep="\t")
        for i in TE_qvchong:
            line = i.split()
            print(*line,sep = "\t")
    inv_TE.close()



if __name__ == "__main__":
    inv = read_inversion(sys.argv[1])
    merge_inv_TE(inv,sys.argv[2])
#    inv = read_inversion("D:\\DESKTOP\\Labratory\\Quercus_a\\Delly\\DW_1_12.vcf.gz")
#    merge_inv_TE(inv,"D:\\DESKTOP\\Labratory\\Quercus_a\\Delly\\DW_TE_define.txt")

