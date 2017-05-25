
import re
from Bio import SeqIO
from Bio.Seq import Seq
import sys

gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3","r")
inlis = open(sys.argv[1],'r').read().split('\r')

## Make dictionary of contig numbers and sequences, contig number key

genome = {}

for record in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively

length = 2000


gene = {}
for i in inlis:
    print i + '\n'
    for line in gff:
        if 'ID=gene' in line and i in line:
            col = line.strip('\n').split('\t')
            det = col[8].split(';')
            pmaa = det[1][5:]
            if col[6] == '+':
                seq = str(Seq(str(genome[col[0]][int(col[3])-length:int(col[3])])))
                gene[pmaa] = seq
            elif col[6] == '-':
                seq = str(Seq(str(genome[col[0]][int(col[4]):int(col[3])+length])).reverse_complement())
                gene[pmaa] = seq


outfile = open(sys.argv[2],'w')

for i in gene:
    outfile.write('>' + i + '\n' + gene[i] + '\n')


# print gene['PMAA_090410']
#
# outlong.write('Gene\tLong motif occurences\tShort motif occurences\n')
#
# for prom in gene:
#     long_count = len(re.findall('ATTTGGC[CT]GG[GC]CC',gene[prom])) + len(re.findall('GG[GC]CC[GA]GCCAAAT',gene[prom]))
#     short_count = len(re.findall('TTGGC[CT]GG',gene[prom])) + len(re.findall('CC[GA]GCCAA',gene[prom]))
#     if short_count > 0 or long_count > 0:
#         outlong.write(prom + '\t' + str(long_count) + '\t' + str(short_count) + '\n')
    
                                  
                                  
gff.close()
outfile.close()


