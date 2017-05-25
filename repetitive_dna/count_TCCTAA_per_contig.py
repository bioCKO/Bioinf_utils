__author__ = 'mjohnpayne'

from Bio import SeqIO


genome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome_mod.fasta","fasta")

cA = 0
cB = 0
for i in genome:
    c1 = str(i.seq).count("TCCTAATCCTAA")
    c2 = str(i.seq).count("TTAGGATTAGGA")
    c3 = str(i.seq).count("CCCTAACCCTAA")
    c4 = str(i.seq).count("TTAGGGTTAGGG")
    cA += c1+c2
    cB += c3+c4
print "TCCTAATCCTAA","CCCTAACCCTAA"
print cA,cB