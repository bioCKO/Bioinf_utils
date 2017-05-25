import sys
from Bio import SeqIO
import re

## python extract_subset_seqs_from_fasta.py path/to/larger/fasta path/to/acc/list path/to/output/fasta

fs = sys.argv[1]
ac = open(sys.argv[2],'r').read()
of = open(sys.argv[3],'w')

def fasta_extract(fasta,accs,outfasta):
    inf = SeqIO.parse(fasta, "fasta")
    acclist = []
    for line in re.split('[\n\r]',accs):
        acclist.append(line)

    #acclist.append(sys.argv[2])
    #print acclist

    for i in inf:
        for a in acclist:
            if i.id == a:
                SeqIO.write(i,outfasta, "fasta")
    outfasta.close()
            
fasta_extract(fs,ac,of)