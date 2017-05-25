
import re
from Bio import SeqIO
from Bio.Seq import Seq
import regex

limit = int(raw_input('Length of 5 prime region to search(int): '))
acc = raw_input('Regular expression of motif: ')
mismat = int(raw_input('Number of mismatches to allow(int): '))


gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3","r")
#outlong = open("/Users/mjohnpayne/Documents/PhD/bys/bys_promoter_analysis/bys_motif_counts","w")
#out5primes = open("/Users/mjohnpayne/Documents/PhD/wt_genome/2161_5prime_intergenics.fa","w")

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



genome = {}

for record in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively

contigs = {}
gene = {}
for line in gff:
        col = line.strip('\n').split('\t')
        if 'ID=gene' in line:
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
                                wins = gene[g][0]-limit
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
                                wine = gene[g][1]+limit
                                winlen = wine - wins
                        if winlen > limit:
                                wine = wins+limit
                        proms[g] = str(Seq(str(genome[cont][wins:wine])).reverse_complement())
                if winlen < 0:
                        continue


print 'Gene\t5 prime intergenic length\tLong motif occurences'
#outlong.write('Gene\t5 prime intergenic length\tLong motif occurences\n')




def repos(reg,string,num):
    pos = []
    for i in list(set(regex.findall(r'(?=('+reg+'){s<='+str(num)+'})',string))):
        for j in re.finditer(i,string):
            pos.append(j.start())
    return pos



for prom in proms:
    fexp = regex.findall(r'(?=('+acc+'){s<='+str(mismat)+'})',proms[prom])
    rexp = regex.findall(r'(?=('+revcomp_re(acc)+'){s<='+str(mismat)+'})',proms[prom])
    #out5primes.write('>'+prom+'_' + str(len(proms[prom])) +'\n'+proms[prom]+'\n')
    posit = repos(acc,proms[prom],mismat) + repos(revcomp_re(acc),proms[prom],mismat)
    posit[:] = [x - len(proms[prom]) for x in posit]
    posit = ','.join(map(str,sorted(posit)))
    long_count = len(fexp) + len(rexp)
    if long_count > 0:
        print prom + '\t' + str(len(proms[prom])) + '\t' + str(long_count) + '\t' + posit
        if len(fexp) > 0:
            print fexp
        if len(rexp) > 0:
            print rexp
        #outlong.write(prom + '\t' + str(len(proms[prom])) + '\t' + str(long_count) + '\n')



gff.close()
#outlong.close()
#out5primes.close()

