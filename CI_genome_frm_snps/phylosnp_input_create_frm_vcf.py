__author__ = 'mjohnpayne'


from glob import glob

from Bio import SeqIO

# filelis = glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_snps/*.vcf")
# #outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_phylosnp/CI_phylosnp_input.csv","w")
# print filelis
#
#
# #outfile.write("Chrom,Position,change,frequency,coverage,quality,Name\n")
#
# for file in filelis:
#     tmp = open(file,'r')
#     name = file.split('/')[-1].split('_')[0]
#     outfile= open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_phylosnp/"+name+"_phylosnp_input.csv","w")
#     outfile.write("Chrom,Position,change,frequency,coverage,quality,Name\n")
#     for line in tmp:
#         if line[0] != '#':
#             col = line.split('\t')
#             if col[0] != '104':
#                 inf = col[7].split(';')
#                 freq = '1'#inf[1].replace('AF=',"")
#                 cov = ''
#                 for i in inf:
#                     if i.startswith('DP='):
#                         cov = i.replace('DP=','')
#                 outfile.write(col[0]+','+col[1] + ',' + col[3] + '|' + col[4] + ',' + freq + ',' + cov + ',' + col[5] + ',' + name + '\n')
#     outfile.close()
#     tmp.close()

##remove excess fastas from genome-shrink.pl
ls = ['2161','3841','4059','3840','F4','BR2SD2','BR2','BR2SD','043','702','203SD4','203SD3','203','203SD','027SD2','027SD1','027','3871','3482','HR2','012']


seqs = SeqIO.parse("/Volumes/MP_HD/CI_GENOME_SEQ/CI_phylosnp/output_file.fasta","fasta")
new = []
for i in seqs:
    if i.id in ls and i not in new:
        new.append(i)

SeqIO.write(new,"/Volumes/MP_HD/CI_GENOME_SEQ/CI_phylosnp/output_file_fixed.fasta","fasta")