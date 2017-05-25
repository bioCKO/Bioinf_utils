from Bio import SeqIO
from Bio.Seq import Seq


infile = open('/Volumes/MP HD/PhD/CLC_output/CI-rerun/3842_fixed_snps.fa','r')
outfile = open('/Volumes/MP HD/PhD/CLC_output/CI-rerun/3842_fixed_snps_old_contigs.fasta','w')

conversions = open('/Users/mjohnpayne/Documents/PhD/wt_genome/contig name conversion pm_fix.txt','r')

conv_list = {}
new_list = []
for line in conversions:
    columns = line.strip('\n').split('\t')
    new_list.append(columns[1])
    conv_list[columns[1]] = columns[0]

seqs = {}

for seq in SeqIO.parse(infile, 'fasta'):
    seqs[seq.id] = seq.seq

for name in new_list:
    for seq in seqs:
        if seq == conv_list[name]:
            outfile.write('>' + str(name) + '_3842' + '\n' + str(seqs[seq]) + '\n')

infile.close()
outfile.close()
conversions.close()
