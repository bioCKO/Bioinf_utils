import subprocess
import shlex
import sys
import re
import time
import glob


 
def aln_run(in_fastq):
    
    sai_path = in_fastq

    sai_path = sai_path.strip('.fastq')

    sai_path = sai_path + '.sai'

    #outname = sam_path[len(sam_path)-1].strip('.fastq')

    #sam_path = sam_path.pop(len(sam_path)-1])

    #sam_path = sam_path.join('/')

#    aln_args = '/usr/bin/bwa aln -t 6 -q 20 /Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta ' + in_fastq + ' > ' + sai_path
    aln_args = '/usr/local/bin/bwa aln -t 6 -q 20 /Volumes/MP_HD/Linda_MNase_Seq/A_nidulans_FGSC_A4_version_s10-m03-r13.fasta ' + in_fastq + ' > ' + sai_path

    aln_out = subprocess.Popen(aln_args, shell=True).wait()

t0= time.time()

def time_out():
    t= time.time()
    diff = t - t0
    timer = float(diff)
    timer = timer/60
    print '%1.2f'%timer + ' minutes'


inp = sys.argv[1]

file_lis = glob.glob(inp + '/*.fastq')
file_lis.sort()

for f in file_lis:
    print f.split('/')[-1].strip('.fastq')
a = 0
b = 1

while b < len(file_lis):
    f1 = file_lis[a]
    f2 = file_lis[b]

    print '\nRunning ' + f1.split('/')[-1].strip('.fastq')

    print '\n\nRunning bwa aln on 1st pair\n'
    aln_run(f1)
    time_out()


    print '\n\nRunning bwa aln on 2nd pair\n'
    aln_run(f2)
    time_out()

    print '\n\nRunning bwa sampe on aln outputs\n'

    fastq1 = f1
    fastq2 = f2

    sai1 = fastq1.strip('.fastq') + '.sai'
    sai2 = fastq2.strip('.fastq') + '.sai'


    sam_path = fastq1.strip('_1.fastq')

    sam_path = sam_path + '.sam'

    sampe_args = '/usr/local/bin/bwa sampe /Volumes/MP_HD/Linda_MNase_Seq/A_nidulans_FGSC_A4_version_s10-m03-r13.fasta ' + sai1 + ' ' + sai2 + ' ' + fastq1 + ' ' + fastq2 + ' > ' + sam_path

    sampe_out = subprocess.Popen(sampe_args, shell=True).wait()

    print '\n\n'

    time_out()

    bam_path = sam_path.strip('sam') + 'bam'

    to_bam_args = '/usr/local/bin/samtools view -S -b -o ' + bam_path + ' ' + sam_path

    print '\n\nConverting sam to bam file\n'

    bam_out = subprocess.Popen(to_bam_args, shell=True).wait()

    print '\n\n'

    time_out()

    sort_bam_path = bam_path.strip('.bam') + '_sort'

    sort_bam_args = '/usr/local/bin/samtools sort ' + bam_path + ' ' + sort_bam_path

    print '\n\nSorting bam file\n'

    bam_sort = subprocess.Popen(sort_bam_args, shell=True).wait()

    print '\n\n'

    time_out()
    a+=2
    b+=2


