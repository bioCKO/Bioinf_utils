__author__ = 'mjohnpayne'

from Bio import SeqIO
import sys

in_gb = SeqIO.parse(sys.argv[1],'genbank')

for i in in_gb:
    print i.id
    out_embl = '/'.join(sys.argv[1].split('/')[:-1]) + i.id + '.embl'
    SeqIO.write(i,out_embl,'embl')