import sys
import re
import subprocess
import shlex
import sys
import time
import glob

inp = sys.argv[1]
out = inp[:-27] + 'provean_input.txt'
print out

infile = open(inp,'r')

outfile = open(out,'w')

#infile = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_mpileup_VCFs/rep_qual_filtered_SNPS_only/012_functional_consequences.txt','r')
#outfile = open('/Volumes/MP_HD/provean_analysis/012_out','w')
aa = open('/Volumes/MP_HD/panther_scoring_lib/amino_acid_abbrev.txt','r')

abrev = {}

for line in aa:
    col = line.strip('\n').split('     ')
    abrev[col[0]] = col[1]
    
for line in infile:
    pos = ''
    change = ''
    snpname = ''
    prot = ''
    if line[0] == '"':
        continue
    else:
        col = line.strip('\n').split('\t')
        if col[12] == 'Yes' and '*' not in col[11]:
            try:
                snpname = col[0] + '_' + col[1]
                prot = col[10][:11]
                if '[' in col[11]:
                    inf = col[11].split(';')
                    pos = inf[1][4:-4]
                    change = abrev[inf[1][1:4]] + pos + abrev[inf[1][-4:-1]]
                else:
                    pos = col[11][17:-3]
                    change = abrev[col[11][14:17]] + pos + abrev[col[11][-3:]]
            except:
                pass
            outfile.writelines(prot + '\t' + change + '\n')
outfile.close()
        
        
    
