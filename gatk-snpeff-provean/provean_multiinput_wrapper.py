import subprocess
import shlex
import sys
import re
import time
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import os
from datetime import datetime
from time import sleep as sl

print time.strftime('%X')

my_env = os.environ
my_env["PATH"] = "/opt/local/bin:/opt/local/sbin:" + my_env["PATH"]

file_lis = glob.glob('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_indels/*indel_provean_inp.txt')
file_lis.sort()
print file_lis

#inp = sys.argv[1]

#insnp = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/CLC_cSNPs/012_provean_input.txt','r')


infastas = SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pm_proteins_with_byss.fasta','fasta')

prots = {}
snps = {}
genes = {}

for seq in infastas:
    prots[seq.id] = seq
    snps[seq.id] = []
    genes[seq.id] = []

strains = []
for i in file_lis:
    end = i.split('/')[-1]
    strain = end[:end.find('_')]
    print strain
    strains.append(strain)
    insnp = open(i,'r')
    for line in insnp:
        col = line.strip('\n').split('\t')
        pmaa = col[0]
        if pmaa not in snps:
            continue
        else:
            snp = col[1]
            snps[pmaa].append(snp)
    for i in genes:
        if len(genes[i]) == 0:
            genes[i].append(len(snps[i]))
        else:
            genes[i].append(len(snps[i])-sum(genes[i]))
print time.strftime('%X')

## provean run

# print '\n'.join(snps['PMAA_058850'])
# print prots['PMAA_058850'].seq
# SeqIO.write(prots['PMAA_058850'],'/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/tmp/PMAA_058850.fasta', "fasta")
# outsnp = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/tmp/PMAA_058850.var','w')
# outsnp.write('\n'.join(snps['PMAA_058850']))
# outsnp.close()

# subprocess.Popen('cd /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/', shell=True).wait()
# tot = len(prots)
# count = 1
#
# for pmaa in snps:
#        isdone = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/done_ls.txt','r').readlines()
#        if pmaa + '\n' in isdone:
#            continue
#        else:
#            if len(snps[pmaa]) > 0:
#                print str(float(count)/tot*100) + '%'
#                SeqIO.write(prots[pmaa],'/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.fasta', "fasta")
#                outsnp = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.var','w')
#                outsnp.write('\n'.join(snps[pmaa]))
#                outsnp.close()
#                protean_args = '/usr/local/bin/provean.sh -q /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.fasta -v /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.var --save_supporting_set /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.sss'
#                outp = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeloutput/'+pmaa+'.out','w')
#                subprocess.Popen(protean_args, shell=True,stdout=outp,env=my_env).wait()
#                outp.close()
#                subprocess.Popen('rm /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.fasta', shell=True).wait()
#                subprocess.Popen('rm /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.var', shell=True).wait()
#                subprocess.Popen('rm /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.sss', shell=True).wait()
#                subprocess.Popen('rm /Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/'+pmaa+'.sss.fasta', shell=True).wait()
#                print time.strftime('%X')
#                donels = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeltmp/done_ls.txt','a')
#                donels.write(pmaa + '\n')
#                donels.close()
#        count += 1


### provean output processing

file_lis2 = glob.glob('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeloutput/*.out')
file_lis2.sort()
outfile = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeloutput/combined_output.txt','w')
outfile2 = open('/Volumes/MP_HD/provean_analysis_derivative_SNP_analysis/indeloutput/combined_neg_scores_output.txt','w')

outfile.write('Accession\t' + '\t'.join(strains) + '\n')
outfile2.write('Accession\t' + '\t'.join(strains) + '\n')

for f in file_lis2:
    name = f.split('/')[-1][:11]
    outfile.write(name + '\t')
    outfile2.write(name + '\t')
    ifile = open(f,'r')
    outsnps = []
    for line in ifile:
        if '[' in line or '#' in line or 'No variations'in line:
            continue
        else:
            outsnps.append(line.strip('\n').replace('\t',':'))
    pos = 0
    for i in range(len(strains)):
        num = genes[name][i]
        if num == 0:
            outfile.write('-\t')
            outfile2.write('-\t')
        else:
            sig = outsnps[pos:pos+num]
            nsig = []
            comb_score = 0
            for j in sig:
                score = float(j.split(':')[1])
                if score <= -2.5:
                    nsig.append(j)
                    comb_score += score
            if len(nsig) > 0:
                outfile.write(','.join(nsig) + '\t')
                outfile2.write(str(comb_score) + '\t')
            else:
                outfile.write('-\t')
                outfile2.write('-\t')
            pos = pos + num
    outfile.write('\n')
    outfile2.write('\n')

outfile.close()
outfile2.close()

    
    
    
    
    
    
    
                 
    
    



