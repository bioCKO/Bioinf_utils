__author__ = 'mjohnpayne'


from Bio import SeqIO
from time import sleep as sl


species = 'Ts'

####
gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')

genomefasta = "/Volumes/MP_HD/repDNA_data/"+species+"_all_genome.fasta"

ingenome = SeqIO.parse(genomefasta,"fasta")
####

genome = {}
for i in ingenome:
    genome[i.id[:-2]] = i.seq

####
rep = "/Volumes/MP_HD/repDNA_data/repeatmasker/"+species+"_all_genome.fasta.out"

inreps = open(rep,'r').readlines()
####

repeatpos = {}

for i in inreps[3:]:
    col = i.split()
    st = col[11].replace(")","").replace("(","")
    en = col[12].replace(")","").replace("(","")
    if col[14] not in repeatpos:
        repeatpos[col[14]] = [int(st),int(en),col[9],col[4][:-2],int(col[0])]
    else:
        repeatpos[col[14]][1] = int(en)
        if col[0] < repeatpos[col[14]][4]:
            repeatpos[col[14]][4] = col[0]


repeats = {}

####
minlength = 1000
minscore = 2500
####

for i in repeatpos:
    if repeatpos[i][1]-repeatpos[i][0] > minlength and repeatpos[i][4] > minscore:
        id = repeatpos[i][2]
        if id not in repeats:
            repeats[id] = [(repeatpos[i][3],repeatpos[i][0],repeatpos[i][1])]
        else:
            repeats[id].append((repeatpos[i][3],repeatpos[i][0],repeatpos[i][1]))



for i in repeats:

    ####
    outfile = open('/Volumes/MP_HD/repDNA_data/repeatmasker/repeat_fastas_more_stringent/' + i + '_'+species+'_sequences.fasta','w')
    ####

    for j in repeats[i]:
        outfile.write(">" + i + "_" + j[0] + "_" + str(j[1]) + ":" + str(j[2]) + '\n' + str(genome[j[0]][j[1]:j[2]]) + '\n')
    outfile.close()