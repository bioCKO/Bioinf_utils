__author__ = 'mjohnpayne'

import sys
import subprocess
import glob
import time



inmedian = open('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/depthofcoverage_stats/CI_median_readdepths.txt','r')

outfile = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/unnannot_bys_gene_read_nos/test.txt','w')

filelis = glob.glob('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/*_GATK_processed.bam')

#geneinf = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/unnannot_bys_gene_read_nos/unannot_byss_pos.txt",'r')

ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')

poslis = []
genlis = []

for i in ingff:
    if '##' not in i:
        col = i.strip('\n').split('\t')
        if col[2] == 'gene':
            pos = col[0] + ':' + col[3] + '-' + col[4]
            poslis.append(pos)
            name = col[8].split(';')[1][5:]
            genlis.append(name)


median = {}
for i in inmedian:
    col = i.strip('\n').split('\t')
    median[col[0]] = float(col[1])

print median
files = {}

for i in filelis:
    ID = i.split('/')[-1].replace("_GATK_processed.bam","")
    print ID
    files[ID] = i


print filelis

samplis = [x.split('/')[-1].replace('_GATK_processed.bam','') for x in filelis]

#print samplis

sortr =  sorted(samplis)

print sortr

outfile.write("GeneID\t" + '\t'.join(sortr) + '\n')

## example pos format in: PMAA_000775\t93:25661-26155\n


# for i in geneinf:
#     col = i.strip('\n').split('\t')
#     poslis.append(col[1])
#     genlis.append(col[0])

#poslis = ['103:618251-619362']
#genlis = ['PMAA_092600']

start = time.time()
print start

RPKMS = {}
outs = {}
for i in range(len(poslis)):
    if i%100 == 0:
        print i
        e = time.time()
        t = e - start
        print t
    outfile.write(genlis[i] + '\t')
    for j in range(len(sortr)):
        f = files[sortr[j]]
        f = f.replace("(","\(")
        f = f.replace(")","\)")
        count_args = '/opt/local/bin/samtools view -c ' + f + ' ' + poslis[i]
        count_out = subprocess.check_output(count_args, shell=True)
        count_out = int(count_out.strip('\n'))
        size = poslis[i].replace('-',':').split(':')
        # print size
        size = int(size[2].replace(',',''))-int(size[1].replace(',',''))
        # print sortr[j]
        # print count_out
        # print size
        rpkm = 0
        if count_out == 0:
            rpkm = 0
        else:
            rpkm = ((float(count_out)*100)/(float(size)))/median[sortr[j]]*100
            #rpkm = float(count_out)/(float(size))*(float((mapped[sortr[j]])))
        RPKMS[sortr[j]] = rpkm
    outfile.write('\t'.join(map(str,[RPKMS[x] for x in sortr])) + '\n')

outfile.close()
# samtools view -c $file 93:71937-72046 >> /Users/mjohnpayne/Documents/PhD/bys/PMAA_076860\(missannot_fused_bys\)/PMAA_000915