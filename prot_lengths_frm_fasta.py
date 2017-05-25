import sys
from Bio import SeqIO
from Bio import GenBank
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC


inp = sys.argv[1]
out = inp.replace('.fasta','_lengths.txt')

inf = SeqIO.parse(inp, "fasta")
outfile = open(out,'w')
ls = []
for i in inf:
    if i.id not in ls:
        outfile.write(i.id + '\t' + str(len(i.seq)) + '\n')
        ls.append(i.id)

outfile.close()
