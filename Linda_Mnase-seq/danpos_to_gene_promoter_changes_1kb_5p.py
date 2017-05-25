__author__ = 'mjohnpayne'

import glob
import subprocess
import os
import sys
import openpyxl as opxl

from itertools import izip,islice,tee
from time import sleep as sl



ing = "/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.fasta"

def fasta_to_dict(infasta):
    genome ={}
    gen = open(infasta,"r").read()
    gen = gen.split(">")
    for i in gen:
        j = i.split('\n')
        name = j[0][:j[0].find(" ")]
        genome[name] = "".join(j[1:])
    return genome

genome = fasta_to_dict(ing)

# for i in genome:
#     print i,len(genome[i])




ingff = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.sorted.gff","r")

def gene_pos_dict(gff):
    dict = {}
    for line in gff:
        if line[0] != "#":
            col = line.strip('\n').split('\t')
            try:
                if col[2] == 'gene':
                    cont = col[0]
                    st = col[3]
                    en = col[4]
                    pmaa = col[8].split(';')[1].replace("Name=","")
                    orient = col[6]
                    if cont not in dict:
                        dict[cont] = {st:[en,pmaa,orient]}
                    else:
                        dict[cont][st] = [en,pmaa,orient]
            except:
                pass
    return dict

gendict = gene_pos_dict(ingff)
# s = 0
# for i in gendict:
#     for j in gendict[i]:
#         print j,gendict[i][j]
#         sl(0.4)
# print s


def ident_env(name,cont,pos,gdict):
    contlis = gdict[cont]
    poslis = sorted(map(int,contlis.keys()))
    ups = []
    dwn = []
    ingene = []
    for j in range(len(poslis)):
        try:
            key = str(poslis[j])
            if j == len(poslis)-1:
                key_p_1 = "terminus"
            else:
                key_p_1 = str(poslis[j+1])
            gene_inf = contlis[key]
            gene1_inf = contlis[key_p_1]
            if int(poslis[j]) < int(pos) < int(gene_inf[0]):
                continue
                # print "contig",cont
                # print "start",poslis[j]
                # print "insertion",pos
                # print "end",gene_inf[0]

                # size = int(gene_inf[0]) - poslis[j]
                # lost = int(gene_inf[0]) - pos
                # perc = (float(lost)/size)*100
                # if contlis[key_p_1][2] == "+":
                #     #outf.write("%s\t%s\t%s\t%s\tYes: %s\t%s\tN/A\tN/A\n" %(name,cont,pos,site,gene_inf[1],str(perc)[:5]))
                #     #print "%s peak is in %s %s%% through the gene" %(name,gene_inf[1],str(perc)[:5])
                #     ingene.append(gene_inf[1])
                # elif contlis[key_p_1][2] == "-":
                #     perc = 100-perc
                #     #outf.write("%s\t%s\t%s\t%s\tYes: %s\t%s\tN/A\tN/A\n" %(name,cont,pos,site,gene_inf[1],str(perc)[:5]))
                #     #print "%s peak is in %s %s%% through the gene" %(name,gene_inf[1],str(perc)[:5])
                #     ingene.append(gene_inf[1])
            elif int(gene_inf[0]) < int(pos) < int(key_p_1):
                # print "contig",cont
                # print "left gene pos",str(gene_inf[0])
                # print "insertion",pos
                # print "right gene pos",str(key_p_1)
                left_d = int(pos) - int(gene_inf[0])
                right_d = int(key_p_1) - int(pos)
                if gene_inf[2] == "-":
                    if contlis[key_p_1][2] == "-" and int(gene_inf[0]) > int(pos) > int(gene_inf[0])+1000:
                        #outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp upstream of %s\t%sbp downstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                        #print "%s site is %sbp upstream of %s and %sbp downstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])
                        ups.append(gene_inf[1])
                    elif contlis[key_p_1][2] == "+" and gene_inf[0] > int(pos) > int(gene_inf[0])+1000 and int(poslis[j+1])-1000 > int(pos) > int(poslis[j+1]):
                        #outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp upstream of %s\t%sbp upstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                        #print "%s site is %sbp upstream of %s and %sbp upstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])
                        ups.append(gene_inf[1])
                        ups.append(contlis[key_p_1][1])
                elif gene_inf[2] == "+":
                    if contlis[key_p_1][2] == "+" and int(poslis[j+1])-1000 > int(pos) > int(poslis[j+1]):
                        #outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp downstream of %s\t%sbp upstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                        #print "%s site is %sbp downstream of %s and %sbp upstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])
                        ups.append(contlis[key_p_1][1])
        except:
            print key
    return ups


peaks = open("/Volumes/MP_HD/Linda_MNase_Seq/danpos_r1_out/Volumes_MP_HD_Linda_MNase_Seq_bam_11H-Volumes_MP_HD_Linda_MNase_Seq_bam_15H.peaks.integrative.txt","r").read().split('\r')

keys = peaks[0].split('\t')

#print keys[0],keys[3],keys[16]
up_5_genes = []
down_5_genes = []
for i in peaks[1:]:
    col = i.split('\t')
    if float(col[-1]) < 0.01:
        dir = ''
        if float(col[16]) < 0:
            dir = "up at 15H"
            u = ident_env(dir,col[0],float(col[3]),gendict)
            if len(u) > 0:
                up_5_genes += u
        elif float(col[16]) > 0:
            dir = "down at 15H"
            u = ident_env(dir,col[0],float(col[3]),gendict)
            if len(u) > 0:
                down_5_genes += u
up_5_genes = set(up_5_genes)
down_5_genes = set(down_5_genes)
exclusive_5_up = up_5_genes.difference(down_5_genes)
exclusive_5_down = down_5_genes.difference(up_5_genes)
all_5_change = set(list(up_5_genes)+list(down_5_genes))
lsls = [up_5_genes,down_5_genes,exclusive_5_up,exclusive_5_down,all_5_change]
lsnames = ["up_5_genes","down_5_genes","exclusive_5_up","exclusive_5_down","all_5_change"]

print "increase occupancy: 5'  decreased occupancy 5'\n"
print len(up_5_genes),len(down_5_genes)

for i in range(len(lsls)):
    outf = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/"+lsnames[i]+"_genes.txt","w")
    for j in lsls[i]:
        outf.write(j+'\n')
    outf.close()


