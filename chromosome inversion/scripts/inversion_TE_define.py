#!/usr/bin/python

import sys

def define_TE(TE_search):
    with open(TE_search,"r") as inv_TE:
        TE_identify = inv_TE.readlines()
        for i in TE_identify:
            for j in TE_identify:
                    line1 = i.split()
                    line2 = j.split()
                    chr1 = line1[0]
                    chr2 = line2[0]
                    pos1_start = line1[1]
                    pos2_start = line2[1]
                    pos1_end = line1[2]
                    pos2_end = line2[2]
                    direction1 = line1[3]
                    direction2 = line2[3]
                    pos1 = line1[5]
                    pos2 = line2[5]
                    name1 = line1[4]
                    name2 = line2[4]

                    if chr1 == chr2 and pos1_start == pos2_start and pos1_end == pos2_end and pos1 != pos2 and direction1 != direction2 and name1 == name2 :
                        print(*line1,sep = "\t")




if __name__ == "__main__":
    #define_TE("D:\\DESKTOP\\Labratory\\Quercus_a\\Delly\\test.txt")
    define_TE(sys.argv[1])


