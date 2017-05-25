__author__ = 'mjohnpayne'

from Bio import SeqIO

ingenome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome.fasta","fasta")

genome = {}

for i in ingenome:
    genome[i.id[:-2]] = i.seq

N = 100000

#outfile = open("/Users/mjohnpayne/Documents/PhD/wt_genome/genes_closer_than_" + str(N/1000) + "kb_from_chrom_ends.txt",'w')

genlist = open("/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids_f.txt","r")
glist = []
for i in genlist:
    i = i.strip('\n')
    glist.append(i)

gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')

y = 0
t = 0
for line in gff:
    if 'ID=gene' in line:
        col = line.strip('\n').split('\t')
        det = col[8].split(';')
        pmaa = det[1][5:]
        st = int(col[3])
        en = int(col[4])
        cont = col[0]
        dist = min([st,len(genome[cont])-en])
        # if pmaa in glist:
        if st < N:
            print pmaa,"Y",dist
            y+=1
        elif en > (len(genome[cont])-N):
            print pmaa,"Y",dist
            y+=1
        else:
            print pmaa,"N",dist
            t+=1
print y
print str((y/float(t+y))*100)+"%"
print t+y
# print len(glist)
# print str((y/float(len(glist)))*100)+"%"
