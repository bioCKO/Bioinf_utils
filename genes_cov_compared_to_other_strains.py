import subprocess
import shlex
import sys
import re
import time
import glob

inp = '/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage/gene coverage/'# sys.argv[1]
inmapped = open('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam/mapped_read_counts.txt','r')

outfile = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage/gene coverage/cov_comparison_between_strains_uncorrected_for_tot_reads.txt','w')

mapped = {}
for i in inmapped:
    col = i.strip('\n').split('\t')
    mapped[col[0]] = float(col[1])/28234821


file_lis = glob.glob(inp + '/*')
file_lis.sort()
genes = {}

#create class made up of 2 lists

class gencov:
    strain = []
    cov = []

#for each cov file add a gencov object for each gene then extend the lists of strain and coverage for each strains file
for j in file_lis:
    if 'gene_cov.txt' in j:
        strain = j.split('/')[-1].replace('_gene_cov.txt','')
        inf = open(j,'r').readlines()
        for i in range(1,len(inf)):
            col = inf[i].strip('\n').split('\t')
            #normalise cov by total mapped read coverage per strain
            cov = float(col[13])/mapped[strain]
            gene = col[2]
            if gene not in genes:
                newclass = gencov()
                newclass.strain += [strain]
                newclass.cov += [cov]
                genes[gene] = newclass
                gencov.strain = []
                gencov.cov = []
            else:
                genes[gene].strain.append(strain)
                genes[gene].cov.append(cov)

outfile.write('Gene\t' + '\t'.join(genes['PMAA_090410'].strain) + '\n')

for i in genes:
    outfile.write(i + '\t' + '\t'.join(map(str,genes[i].cov)) + '\n')

outfile.close()
