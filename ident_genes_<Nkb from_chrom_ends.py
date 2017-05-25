__author__ = 'mjohnpayne'

from Bio import SeqIO

ingenome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome.fasta","fasta")

genome = {}

for i in ingenome:
    genome[i.id[:-2]] = i.seq

N = 100000

outfile = open("/Users/mjohnpayne/Documents/PhD/wt_genome/genes_closer_than_" + str(N/1000) + "kb_from_chrom_ends.txt",'w')


gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')


for line in gff:
    if 'ID=gene' in line:
        col = line.strip('\n').split('\t')
        det = col[8].split(';')
        pmaa = det[1][5:]
        st = int(col[3])
        en = int(col[4])
        cont = col[0]
        if st < N:
            outfile.write(pmaa+'\t'+"Y" + '\n')
        elif en > (len(genome[cont])-N):
            outfile.write(pmaa+'\t'+"Y" + '\n')
        else:
            outfile.write(pmaa+'\t'+"N" + '\n')

outfile.close()