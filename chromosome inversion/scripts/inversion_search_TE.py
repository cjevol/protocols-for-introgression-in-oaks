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
def read_TEanno(TEfile):
    with open(TEfile,"r") as TE:
        TEanno = TE.readlines()
        TE_identify = []
        for line in TEanno:
            TE_info = line.split()
            chr = TE_info[0]
            TE_class = TE_info[2]
            pos_start = TE_info[3]
            pos_end = TE_info[4]
            description = TE_info[8]
            TE_identify.append(chr)
            TE_identify.append(TE_class)
            TE_identify.append(pos_start)
            TE_identify.append(pos_end)
            TE_identify.append(description)
    TE.close()
    return TE_identify

def search_TE(inv,TE):
    for i in range(0,len(inv)//3):
        for j in range(0,len(TE)//5):
            if inv[3*i] == TE[5*j]:
                inv1_start = int(inv[3*i+1])-1000
                inv1_end = int(inv[3*i+1])+1000
                inv2_start =int(inv[3*i+2])-1000
                inv2_end = int(inv[3*i+2])+1000
                TE_start = int(TE[5*j+2])
                TE_end = int(TE[5*j+3])
                if (TE_start>inv1_start and TE_start<inv1_end) or (TE_start<inv1_start and TE_end>inv1_end) or (TE_start>inv1_start and TE_end<inv1_end) or (TE_start<inv1_start and TE_end>inv1_start):
                    info = [inv[3*i],inv[3*i+1],inv[3*i+2],"left", TE[5*j+1],TE[5*j+2],TE[5*j+3],TE[5*j+4]]
                    print(*info, sep='\t')
                if (TE_start>inv2_start and TE_start<inv2_end) or (TE_start<inv2_start and TE_end>inv2_end) or (TE_start>inv2_start and TE_end<inv2_end) or (TE_start<inv2_start and TE_end>inv2_start):
                    info = [inv[3 * i], inv[3 * i + 1], inv[3 * i + 2], "right", TE[5 * j + 1], TE[5 * j + 2], TE[5 * j + 3], TE[5 * j + 4]]
                    print(*info, sep='\t')


if __name__ == "__main__":
    #inv = read_inversion("D:\\DESKTOP\\Labratory\\Quercus_a\\Delly\\20220112补充\\DW\DW.select.vcf")
    #TE = read_TEanno("D:\\DESKTOP\\Labratory\\Quercus_a\\Delly\\Quercus_a_chr1_12_TEanno.gff")
    inv = read_inversion(sys.argv[1])
    TE = read_TEanno(sys.argv[2])
    search_TE(inv,TE)



