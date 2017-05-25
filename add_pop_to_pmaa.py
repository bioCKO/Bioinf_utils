import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq


##infile = open(sys.argv[1],"r")
##outfile = open(sys.argv[2],"w")
##inconv = open(sys.argv[3],"r")


outfile = open("/Users/mjohnpayne/Documents/Uni/phd/Asp_sequences/MEGA_analysis/500bp_5prime_alignment/1000bp_5prime_seq_asps_mod.fasta","w")
inconv = open("/Users/mjohnpayne/Documents/Uni/phd/asps_to_pops_fix.txt","r")

convert = {}

for line in inconv:
    columns = line.strip("\n").split("\t")
    adpop = str(columns[0] + "_" + columns[1] + " ")
    convert[columns[0]] = adpop

infile = {}
for record in SeqIO.parse("/Users/mjohnpayne/Documents/Uni/phd/Asp_sequences/MEGA_analysis/500bp_5prime_alignment/1000bp_5prime_seq_asps.fasta", "fasta"):
        infile[record.id] = record.seq

for asp in convert:
    for name in infile:
        if asp in str(name):
            outfile.writelines(">" + convert[asp] + "\n" + infile[name] + "\n")

outfile.close()
inconv.close()
