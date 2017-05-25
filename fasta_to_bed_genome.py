
from Bio import SeqIO

in_fasta = SeqIO.parse('/Volumes/MP_HD/BAM/S_pombe.fasta','fasta')
outfile = open('/Volumes/MP_HD/BAM/pombe_cov.bed.txt','w')

for contig in in_fasta:
    outfile.write(contig.id + '\t' + str(len(contig.seq)) + '\n')


in_fasta.close()
outfile.close()
    
