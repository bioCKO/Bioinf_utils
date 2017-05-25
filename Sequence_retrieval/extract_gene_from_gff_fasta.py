from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
import sys
from time import sleep as sl

in_gff = open(sys.argv[1],'r')#open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_denovo.gff','r')

in_fas = SeqIO.parse(sys.argv[2],'fasta')#SeqIO.parse('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/PF_SSPACE.final.scaffolds.fasta','fasta')

outfile = open(sys.argv[3],'w')#open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_cds.fa','w')





#
# tf_fas = {}
# for i in in_fas:
#     tf_fas[i.id] = str(i.seq)
# cds = ''
# name = ''
# orient = '+'
# scaf = ''
# for line in in_gff:
#     if line[0] != "#":
#         #print line
#         #sl(0.5)
#         col = line.strip('\n').split('\t')
#         if col[2] == 'gene':
#             inf = col[8].split(';')
#             name = inf[1].strip('Name=')
#             for i in tf_fas:
#                 if i == col[0]:
#                     scaf = tf_fas[i]
#             orient = col[6]
#             st = int(col[3])
#             en = int(col[4])
#             if orient == "+":
#                 gene = scaf[st-1:en]
#                 outfile.write('>' + name + '\n' + gene + '\n')
#             elif orient == "-":
#                 gene = str(Seq(scaf[st-1:en],generic_nucleotide).reverse_complement())
#                 outfile.write('>' + name + '\n' + gene + '\n')
#
# outfile.close()

tf_fas = {}
for i in in_fas:
    tf_fas[i.id] = str(i.seq)
cds = ''
name = ''
orient = '+'
scaf = ''
for line in in_gff:
    if line[0] != "#":
        #print line
        #sl(0.5)
        col = line.strip('\n').split('\t')
        if col[2] == 'gene':
            inf = col[8].split(';')
            name = inf[1].strip('Name=')
        if col[2] == 'CDS':
            for i in tf_fas:
                if i == col[0]:
                    scaf = tf_fas[i]
            orient = col[6]
            st = int(col[3])
            en = int(col[4])
            if orient == "+":
                gene = scaf[st-1:en]
                outfile.write('>' + name + '\n' + gene + '\n')
            elif orient == "-":
                gene = str(Seq(scaf[st-1:en],generic_nucleotide).reverse_complement())
                outfile.write('>' + name + '\n' + gene + '\n')

outfile.close()