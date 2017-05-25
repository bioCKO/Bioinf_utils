__author__ = 'mjohnpayne'

from glob import glob
from Bio import SeqIO
from time import sleep as sl

## interproscan output files

intpros_in = glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/*.tsv")

## protein sequences

prots = SeqIO.parse("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/int_files/goodProteins.fasta","fasta")

##dict of prot id and prot seq object
prot = {}
for i in prots:
    prot[i.id] = i

## Output fasta file

outfasta = "/Users/mjohnpayne/Documents/PhD/ASPS/MEGA_analysis_CLEANUP/asp_interpro_all/IPR021109_all_Talaro.fasta"

## Target Interpro ID

id = "IPR021109"

## iterate over tsvs and collect list of ids with target

genelist = []

for i in intpros_in:
    if "rename" not in i:
        tmp = open(i,"r")
        for line in tmp:
            col = line.split('\t')
            if len(col) > 11:
                if col[11] == id:
                    if col[0][:3] in ["Tfl","Pfu","Tst","Pma"]:
                    #print col[0]
                        genelist.append(prot[col[0]])
        tmp.close()
print len(genelist)
genelist = list(set(genelist))
print len(genelist)


SeqIO.write(genelist,outfasta,"fasta")
