__author__ = 'mjohnpayne'

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import subprocess


ingenome = SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta','fasta')


gen = {}

for i in ingenome:
    gen[i.id] = str(i.seq)

def findOccurences(s, ch):
    pos = []
    win = len(ch)
    st = 0
    en = st + win
    while en < len(s):
        if s[st:en] == ch:
            pos.append(str(st))
        st += 1
        en += 1
    return pos

for i in gen:
    pos = findOccurences(gen[i],'CCCGGG')
    print i + '\t' + '\t'.join(pos)
