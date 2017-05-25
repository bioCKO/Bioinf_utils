import re
from Bio import SeqIO
from Bio.Seq import Seq
import regex

motif = 'TTAGGG'

outfile = open("/Users/mjohnpayne/Documents/PhD/wt_genome/Poss_telomere_window/" + motif + "_exact_tab_file_10kbwindow.txt","w")

## Make dictionary of contig numbers and sequences, contig number key

genome = {}

for record in SeqIO.parse("/Users/mjohnpayne/Documents/Phd/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq


def revcomp_re(str):
    new = []
    for i in str:
        if i == 'T':
            new = ['A'] + new
        elif i == 'A':
            new = ['T'] + new
        elif i == 'G':
            new = ['C'] + new
        elif i == 'C':
            new = ['G'] + new
        elif i == '[':
            new = [']'] + new
        elif i == ']':
            new = ['['] + new
    return ''.join(new)


for contig in genome:
    winstart = 0
    winend = 10000
    step = 1000
    window_no = 0
    if len(genome[contig]) < winend:
        sequence = genome[contig]
        count = len(regex.findall(r'(?=('+motif+'){s<=0})',str(sequence)))
        count += len(regex.findall(r'(?=('+revcomp_re(motif)+'){s<=1})',str(sequence)))
        outfile.writelines(str(contig) + "\t" + "1" + "\t" + str(len(genome[contig])) + "\t" + str(count) + "\n")
    else:
        while winend < len(genome[contig]):
            window = genome[contig][winstart:winend]
            wincount = len(regex.findall(r'(?=('+motif+'){s<=0})',str(window)))
            wincount += len(regex.findall(r'(?=('+revcomp_re(motif)+'){s<=1})',str(window)))
            outfile.writelines(contig + "\t" + str(winstart) + "\t" + str(winend) + "\t" + str(wincount) + "\n")
            winstart = winstart + step
            winend = winend + step
            window_no = window_no + 1

outfile.close()