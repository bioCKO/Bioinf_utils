__author__ = 'mjohnpayne'

import sys
import subprocess
import glob


inmapped = open("/Volumes/MP_HD/RNASeq_data/RNAseq_mapped_read_counts.txt",'r').read().split('\r')


coverage = {}
for i in inmapped:
    col = i.replace(',','').split('\t')
    coverage[col[0]] = int(col[1])

print coverage

filelis = glob.glob('/Volumes/MP_HD/RNASeq_data/*sort.bam')

print filelis

samplis = [x.split('/')[-1].replace('_sort.bam','') for x in filelis]

print samplis



poslis = ['103:3315219-3315872']


genlis = ['PMAA_102265']
outfile = open('/Volumes/MP_HD/RNASeq_data/RPKM_generation/PMAA_102265.txt','w')


RPKMS = {}
outs = {}
for i in range(len(poslis)):
    for j in range(len(samplis)):
        print poslis[i]
        count_args = '/opt/local/bin/samtools view -c ' + filelis[j] + ' ' + poslis[i]
        count_out = subprocess.check_output(count_args, shell=True)
        count_out = int(count_out.strip('\n'))
        size = poslis[i].replace('-',':').split(':')
        size = int(size[2])-int(size[1])
        print samplis[j]
        print count_out
        if count_out == 0:
            rpkm = 0
        else:
            rpkm = float(count_out)/((float(size)/1000)*(float((coverage[samplis[j]]))/1000000))
        RPKMS[samplis[j]] = rpkm

sortr =  sorted(samplis)

outfile.write('\t'.join(sortr) + '\n')
outfile.write('\t'.join(map(str,[RPKMS[x] for x in sortr])) + '\n')

outfile.close()
# # samtools view -c $file 93:71937-72046 >> /Users/mjohnpayne/Documents/PhD/bys/PMAA_076860\(missannot_fused_bys\)/PMAA_000915