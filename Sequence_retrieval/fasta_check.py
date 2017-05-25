__author__ = 'mjohnpayne'

from Bio import SeqIO

fasta = SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pm_proteins_with_byss_2.fasta','fasta')

for i in fasta:
    if len(i.seq) == 0:
        print i.id