__author__ = 'mjohnpayne'


import re
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import regex

limit = int(raw_input('5 prime limit: '))
inacc = open(raw_input('accession number list path: '),'r')
outfile = raw_input('output file path: ')
fasta = '/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta'
gff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')



#outlong = open('/Users/mjohnpayne/Documents/PhD/bys/bys_promoter_analysis/bys_motif_counts','w')
#out5primes = open('/Users/mjohnpayne/Documents/PhD/wt_genome/2161_5prime_intergenics.fa','w')

## Make dictionary of contig numbers and sequences, contig number key
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


def get_genome(fasta_str):
    gnm = {}
    for record in SeqIO.parse(fasta_str, 'fasta'):
            gnm[record.id] = record.seq
    return gnm

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively


def get_gene_conigs(gf):
    conts = {}
    gen = {}
    for line in gf:
            col = line.strip('\n').split('\t')
            if 'ID=gene' in line:
                    det = col[8].split(';')
                    pmaa = det[1][5:]
                    st = int(col[3])
                    en = int(col[4])
                    orient = col[6]
                    cont = col[0]
                    if cont in conts:
                            conts[cont] += [pmaa]
                    elif cont not in conts:
                            conts[cont] = [pmaa]
                    gen[pmaa] = [st,en,orient]
    for cont in conts:
            conts[cont] = sorted(conts[cont])
    return conts,gen

def getproms(cont,genme,lim,genes):
    proms = {}
    for cont in cont:
            for i in range(0,len(cont[cont])):
                    g = cont[cont][i]
                    wins = 0
                    wine = 0
                    winlen = 0
                    if genes[g][2] == '+':
                            try:
                                    if genes[g][0] < genes[cont[cont][i-1]][1]:
                                            wins = genes[cont[cont][i-2]][1]
                                            wine = genes[g][0]
                                            winlen = wine - wins
                                    else:
                                            wins = genes[cont[cont][i-1]][1]
                                            wine = genes[g][0]
                                            winlen = wine - wins
                            except:
                                    wins = genes[g][0]-lim
                                    wine = genes[g][0]
                                    if wins < 0:
                                            wins = 0
                                    winlen = wine - wins
                            if winlen > lim:
                                    wins = wine-lim
                            proms[g] = str(Seq(str(genme[cont][wins:wine])))
                    elif genes[g][2] == '-':
                            try:
                                    if genes[g][1] > genes[cont[cont][i+1]][0]:
                                            wins = genes[g][1]
                                            wine = genes[cont[cont][i+2]][0]
                                            winlen = wine - wins
                                    else:
                                            wins = genes[g][1]
                                            wine = genes[cont[cont][i+1]][0]
                                            winlen = wine - wins
                            except:
                                    wins = genes[g][1]
                                    wine = genes[g][1]+lim
                                    winlen = wine - wins
                            if winlen > lim:
                                    wine = wins+lim
                            proms[g] = str(Seq(str(genme[cont][wins:wine])).reverse_complement())
                    if winlen < 0:
                            continue
    return proms

def main(fas,gf,lim,of,inacs):
    genome = get_genome(fas)
    contigs,gene = get_gene_conigs(gf)
    proms = getproms(contigs,genome,lim,gene)
    #outlong.write('Gene\t5 prime intergenic length\tLong motif occurences\n')
    of = open(of,'w')
    for i in inacs.read().split('\r'):
        print i
        of.write('>' + i + '\n' + proms[i] + '\n')
    of.close()

main(fasta,gff,limit,outfile,inacc)