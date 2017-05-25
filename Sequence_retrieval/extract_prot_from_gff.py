__author__ = 'mjohnpayne'



import re
from Bio import SeqIO
from Bio.Seq import Seq
import regex
from time import sleep as sl

#inacc = open('/Users/mjohnpayne/Documents/PhD/bys/bys_promoter_analysis/complete_list_bys_ids.txt','r').read()
outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Pm_Ts_augustus/Ts_augustus_gene.fasta','w')
type = 'prot'


gff = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Pm_Ts_augustus/Ts_augustus.gff","r").readlines()
#outlong = open("/Users/mjohnpayne/Documents/PhD/bys/bys_promoter_analysis/bys_motif_counts","w")
#out5primes = open("/Users/mjohnpayne/Documents/PhD/wt_genome/2161_5prime_intergenics.fa","w")

## Make dictionary of contig numbers and sequences, contig number key

class loci:
    def __init__(self,gene,mrna,cds,orient,prot):
        self.gene = gene
        self.mrna = mrna
        self.cds = cds
        self.orient = orient
        self.prot = prot


def revcomp_re(str):
    new = []
    for i in str:
        if i == 'T':
            new = ['A'] + new
        elif i == 'A':
            new = ['T'] + new
        elif i == 'G':
            new = ['C'] + new
        elif i == 'C':
            new = ['G'] + new
        elif i == '[':
            new = [']'] + new
        elif i == ']':
            new = ['['] + new
    return ''.join(new)



genome = {}

for record in SeqIO.parse("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Pm_Ts_augustus/Ts_all_genome.fasta", "fasta"):
        genome[record.id] = record.seq

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively

contigs = {}
gene = {}
for i in range(len(gff)-1):
        col = gff[i].strip('\n').split('\t')
        if '\tgene\t' in gff[i]:
            det = col[8]
            pmaa = det[3:].strip('\n')
            gene[pmaa] = loci('','','','','')
            st = int(col[3])
            en = int(col[4])
            orient = col[6]
            gene[pmaa].orient = orient
            cont = col[0]
            if cont in contigs:
                    contigs[cont] += [pmaa]
            elif cont not in contigs:
                    contigs[cont] = [pmaa]
            if orient == '+':
                gene[pmaa].gene = genome[cont][st-1:en]
            elif orient == '-':
                gene[pmaa].gene = str(Seq(str(genome[cont][st:en])).reverse_complement())
        elif '\ttranscript\t' in gff[i]:
            st = int(col[3])
            en = int(col[4])
            if orient == '+':
                gene[pmaa].mrna = genome[cont][st-1:en]
            elif orient == '-':
                gene[pmaa].mrna = str(Seq(str(genome[cont][st:en])).reverse_complement())
        elif '\tCDS\t' in gff[i]:
            if '\tCDS\t' not in gff[i-1] and '\tCDS\t' not in gff[i+1]:
                st = int(col[3])
                en = int(col[4])
                curcds = genome[cont][st-1:en]
                if orient == '+':
                    gene[pmaa].cds = curcds
                    gene[pmaa].prot = str(Seq(str(curcds)).translate())
                elif orient == '-':
                    revcds = Seq(str(curcds)).reverse_complement()
                    gene[pmaa].cds = str(revcds)
                    gene[pmaa].prot = str(Seq(str(revcds)).translate())
            elif '\tCDS\t' not in gff[i-1] and '\tCDS\t' in gff[i+1]:
                curcds = ''
                st = int(col[3])
                en = int(col[4])
                curcds = genome[cont][st-1:en]
            elif '\tCDS\t' in gff[i-1] and '\tCDS\t' in gff[i+1]:
                st = int(col[3])
                en = int(col[4])
                curcds += genome[cont][st-1:en]
            elif '\tCDS\t' in gff[i-1] and '\tCDS\t' not in gff[i+1]:
                st = int(col[3])
                en = int(col[4])
                curcds += genome[cont][st-1:en]
                if orient == '+':
                    gene[pmaa].cds = curcds
                    gene[pmaa].prot = str(Seq(str(curcds)).translate())
                elif orient == '-':
                    revcds = Seq(str(curcds)).reverse_complement()
                    gene[pmaa].cds = str(revcds)
                    gene[pmaa].prot = str(Seq(str(revcds)).translate())




for cont in contigs:
        contigs[cont] = sorted(contigs[cont])

# for i in inacc.split('\r'):
#     if type == 'prot':
#         outfile.write('>' + i + '\n' + gene[i].prot + '\n')
#     elif type == 'cds':
#         outfile.write('>' + i + '\n' + str(gene[i].cds) + '\n')
#     elif type == 'mrna':
#         outfile.write('>' + i + '\n' + str(gene[i].mrna) + '\n')
#     elif type == 'gene':
#         outfile.write('>' + i + '\n' + str(gene[i].gene) + '\n')
#
# outfile.close()


for i in gene:
    # print i
    # print gene[i].gene
    # sl(0.2)
    outfile.write('>' + str(i) + '\n' + str(gene[i].gene) + '\n')

outfile.close()


#print gene['PMAA_030490'].cds
#print gene['PMAA_030490'].prot