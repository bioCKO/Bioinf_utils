__author__ = 'mjohnpayne'


from Bio import SeqIO
import subprocess
from time import sleep as sl


ingenome = SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/repDNA_data/Pm_all_genome.fasta","fasta")

outfile = open("/Users/mjohnpayne/Documents/PhD/wt_genome/repDNA_data/Pm_all_genome_rep_counts_20mer.txt",'w')

## Make concatenated genome for db to search
concat = ''

for i in ingenome:
    concat += i.seq

ingenome.close()
ingenome = SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/repDNA_data/Pm_all_genome.fasta","fasta")

print "Genome length: " + str(len(concat))

## use repeat-match -n 50 /Users/mjohnpayne/Documents/PhD/wt_genome/repDNA_data/Pm_all_genome.fasta instead of own algorithm

merlen = 20

merdict = {}
ls = []

repargs = '/usr/bin/repeat-match -n 50' +



# outfile.write(str(merlen) + "mer start position" + '\t' + "number of occurences\n")
#
# for i in ingenome:
#     print i.id
#     st = 0
#     while st < len(i.seq)-merlen:
#         window = i.seq[st:st+merlen]
#         if window not in merdict:
#             count = concat.count(window)
#             merdict[window] = count
#             outfile.write(str(st) + '\t' + str(count) + '\n')
#         else:
#             outfile.write(str(st) + '\t' + str(merdict[window]) + '\n')
#         st += merlen
#
# outfile.close()

