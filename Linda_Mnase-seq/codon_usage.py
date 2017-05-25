__author__ = 'mjohnpayne'

from Bio.SeqUtils import CodonUsage as co
from Bio import SeqIO
from time import sleep as sl

infasta = SeqIO.parse("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13_orf_coding.fasta","fasta")

codon_d = co.CodonsDict
codon_index = co.SynonymousCodons


total = 0

for i in infasta:
    for j in range(0,len(i.seq),3):
        codon = i.seq[j:j+3]
        if len(codon) == 3:
            codon_d[codon] += 1
            total +=1

ncodon = {}

for i in codon_d:
    ncodon[i] = (float(codon_d[i])/total)*100



outfile = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13_codon_usage.txt","w")
outfile.write("Amino Acid\tCodon\tcount\tpercentage\tcodon percentage of AA\n")
aa_codon_p = {}

for i in codon_index:
    sum = 0
    for j in codon_index[i]:
        sum += ncodon[j]
    aa_codon_p[i] = sum


for i in codon_index:
    for j in codon_index[i]:
        outfile.write(i+'\t'+j+"\t"+str(codon_d[j])+'\t'+str(ncodon[j])+'\t'+str((ncodon[j]/aa_codon_p[i])*100)+'\n')

outfile.close()
