import subprocess
import shlex
import sys
import re
import time
import glob
import os
from cogent import LoadSeqs, DNA
from cogent.parse.fasta import MinimalFastaParser
from Bio import SeqIO
import os



inp = sys.argv[1]


def check3div(fasta):
    seqs = SeqIO.parse(fasta,'fasta')
    outseq = []
    for i in seqs:
        if len(i.seq)%3 != 0:
            i.seq +=('N')*(3-(len(i.seq)%3))
            outseq += [i]
        else:
            outseq += [i]
    return outseq



def runPrank(input):
    out = input.split('/')
    out[-2] = 'Prank_fasta_alns'
    out = '/'.join(out).replace('.fasta','_prank')
    outlis = check3div(inp)
    outfile = open(input,'w')
    SeqIO.write(outlis,outfile,'fasta')
    outfile.close()
    prank_inp = '/Users/mjohnpayne/Documents/PhD/bioinformatics_tools/prank/bin/prank -codon -quiet -showtree -f=fasta -o=' + out + ' -d=' + input
    subprocess.Popen(prank_inp, shell=True).wait()
    # newp = inp.split('/')
    # newp[-2] = 'Prank_finished_fastas'
    # os.rename(inp,'/'.join(newp))


runPrank(inp)

