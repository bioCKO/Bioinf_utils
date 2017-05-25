import sys
import re
import regex
from Bio import SeqIO
from Bio.Seq import Seq
import operator
from time import sleep as sl

motif = sys.argv[1]#raw_input('motif in regular expression format (ie CYG = C[CT]G): ')
limit = int(sys.argv[2])#int(raw_input('max length of intergenic region to search(ie 1000): '))
mismat = sys.argv[3]#raw_input('number of mismatches allowed: ')
outlong = open(sys.argv[4],'w')#open(raw_input('path to output file: '),'w')
#gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3","r")
gen = sys.argv[5]
gff = open(sys.argv[6],'r')#open("/Users/mjohnpayne/Documents/PhD/wt_genome/tsta1_working_models.gff3","r")




## Make dictionary of contig numbers and sequences, contig number key

genome = {}

# for record in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
#         genome[record.id] = record.seq

for record in SeqIO.parse(gen, "fasta"):
        genome[record.id] = record.seq

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively

# def motif_find(mot,seq,mismat):
#     st = 0
#     en = st + len(mot)
#     while en < len(seq):
#         region = seq[st:en]
#         count = 0
#         score = 0
#         mot =
#         for char in mot:
#             if region[count]in char:
#                 score += 1
#             count += 1
#         if score >= (len(pat) - mis):
#             hits.append(start)



contigs = {}
gene = {}
for line in gff:
    if line[0] != '#':
        col = line.strip('\n').split('\t')
        # for pm and ts
        # if col[2] =='gene':
        #         det = col[8].split(';')
        #         pmaa = det[1][5:]
        #         st = int(col[3])
        #         en = int(col[4])
        #         orient = col[6]
        #         cont = col[0]
        #         if cont in contigs:
        #                 contigs[cont] += [pmaa]
        #         elif cont not in contigs:
        #                 contigs[cont] = [pmaa]
        #         gene[pmaa] = [st,en,orient]
        # for Pf and Tf
        if col[2] =='gene':
                det = col[8].split(';')
                pmaa = det[1][5:]
                st = int(col[3])
                en = int(col[4])
                orient = col[6]
                cont = col[0]
                if cont in contigs:
                        contigs[cont] += [pmaa]
                elif cont not in contigs:
                        contigs[cont] = [pmaa]
                gene[pmaa] = [st,en,orient]



for cont in contigs:
        contigs[cont] = sorted(contigs[cont])


proms = {}
for cont in contigs:
        for i in range(0,len(contigs[cont])):
                g = contigs[cont][i]
                wins = 0
                wine = 0
                winlen = 0
                if gene[g][2] == '+':
                        try:
                                if gene[g][0] < gene[contigs[cont][i-1]][1]:
                                        wins = gene[contigs[cont][i-2]][1]
                                        wine = gene[g][0]
                                        winlen = wine - wins                                        
                                else:
                                        wins = gene[contigs[cont][i-1]][1]
                                        wine = gene[g][0]
                                        winlen = wine - wins
                        except:
                                wins = gene[g][0]-2000
                                wine = gene[g][0]
                                if wins < 0:
                                        wins = 0
                                winlen = wine - wins
                        if winlen > limit:
                                wins = wine-limit
                        proms[g] = str(Seq(str(genome[cont][wins:wine])))
                elif gene[g][2] == '-':
                        try:
                                if gene[g][1] > gene[contigs[cont][i+1]][0]:
                                        wins = gene[g][1]
                                        wine = gene[contigs[cont][i+2]][0]
                                        winlen = wine - wins
                                else:
                                        wins = gene[g][1]
                                        wine = gene[contigs[cont][i+1]][0]
                                        winlen = wine - wins
                        except:
                                wins = gene[g][1]
                                wine = gene[g][1]+2000
                                winlen = wine - wins
                        if winlen > limit:
                                wine = wins+limit
                        proms[g] = str(Seq(str(genome[cont][wins:wine])).reverse_complement())



outlong.write('Gene\t5 prime intergenic length\tMotif count\tMotif positions(plus strand)\tMotif sequences(plus strand)\tMotif positions(minus strand)\tMotif sequences(minus strand)\tpromoter score\n')

allmot = {}
for prom in proms:
    if int(mismat) > 0:
        p = regex.compile(r'(?:'+motif+'){s<' + mismat + '}')##r''motif{s<=2})
    else:
        p = regex.compile(r'(?:'+motif+')')
    lsp = []
    dicp = {}
    lsn = []
    dicn = {}
    for m in p.finditer(proms[prom]):
        pos = int(m.start()-len(proms[prom]))
        mot = proms[prom][m.start():m.end()]
        dicp[pos] = mot
        if mot not in allmot:
            allmot[mot] = 1
        else:
            allmot[mot] += 1
    for m in p.finditer(str(Seq(str(proms[prom])).reverse_complement())):
        pos = 0 - m.start()
        mot = str(Seq(str(proms[prom])).reverse_complement())[m.start():m.end()]
        dicn[pos] = mot
        if mot not in allmot:
            allmot[mot] = 1
        else:
            allmot[mot] += 1
    lsp = sorted(dicp.keys(), key=int,reverse=True)
    lsn = sorted(dicn.keys(), key=int,reverse=True)
    #lsp = map(str,lsp)
    #lsn = map(str,lsn)
    motlisp = []
    for i in lsp:
        motlisp += [dicp[i]]
    motlisn = []
    for j in lsn:
        motlisn += [dicn[j]]

    mots = lsp + lsn



    score = 0
    score += len(mots)
    # for m in mots:
    #     score -= 2
    #     for n in mots:
    #         if m-40 < n < m+40:
    #             score +=4
    for m in mots:
        if m > -500:
            score +=2
            score -= 4
            for n in mots:
                if m-40 < n < m+40:
                    score +=4
        elif -1000 < m < -500:
            score +=1
            score -= 2
            for n in mots:
                if m-40 < n < m+40:
                    score +=2

    #print prom, lsp, lsn, mots, score
    #sl(0.5)
    if len(lsn) > 0 or len(lsp) > 0:
         outlong.write(prom + '\t' + str(len(proms[prom])) + '\t' + str(len(lsp)+len(lsn)) + '\t' + ','.join(map(str,lsp)) + '\t' + ','.join(motlisp) + '\t' + ','.join(map(str,lsn)) + '\t' + ','.join(motlisn)+ '\t'+ str(score) + '\n')
outlong.write('\n\n\n\nNumber of occurences of individual motif sequences\n')

sorted_x = sorted(allmot.items(), key=operator.itemgetter(1))

for i in sorted_x:
    outlong.write(str(i[0]) + ': ' + str(i[1]) + '\n')
                                  
gff.close()
outlong.close()



# each motif +1
# two motifs within 40bp +2 each motif
# 0-500 +2
# 500-1000 +1

# m2 -- m4 -- m2   +8


## ATTTGGC[CT]GG[GC]CC


######/Volumes/MP_HD/for Kylie/new_trial_8-9-14/test_long3
#
# inp = open(raw_input('inpath: '),'r').readlines()
#
# pat = inp[0].strip('\n')
# seq = inp[1].strip('\n')
# mis = int(inp[2])
#
# start = 0
#
# end = start + len(pat)
#
# hits = []
#
# while end < len(seq) + 1:
#     hit = seq[start:end]
#     count = 0
#     score = 0
#     for char in pat:
#         if hit[count] == char:
#             score += 1
#         count += 1
#     if score >= (len(pat) - mis):
#         hits.append(start)
#     start += 1
#     end += 1
#
# print ' '.join(map(str,hits))