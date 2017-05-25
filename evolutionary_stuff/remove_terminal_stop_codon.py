__author__ = 'mjohnpayne'


import sys
from Bio import SeqIO

infile = sys.argv[1]
outfile = sys.argv[1][:-6] + '_notstop.fasta'

fasta = SeqIO.parse(infile,'fasta')
out = []

for i in fasta:
    new = SeqIO.SeqRecord(i.seq[:-3],id=i.id,name=i.name,description=i.description)
    out.append(new)

SeqIO.write(out,outfile,'fasta')