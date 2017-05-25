__author__ = 'mjohnpayne'


from Bio import SeqIO
import re
from time import sleep as sl


outfile = open("/Users/mjohnpayne/Documents/PhD/transposon insertion sites/all_te_sites",'w')

genome = {}

for record in SeqIO.parse("/Users/mjohnpayne/Documents/Phd/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

tot = []

outfile.write('Contig\tPosition\n')

for contig in genome:
    con = str(genome[contig])
    for m in re.finditer('TTAA', con):
        pos = m.start()
        outfile.write(contig+'\t'+str(pos)+'\n')
#print tot[:100]