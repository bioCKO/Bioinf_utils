__author__ = 'mjohnpayne'

from Bio import SeqIO
import sys


fs = sys.argv[1]
insignalp = open(sys.argv[2],'r')
of = open(sys.argv[3],'w')

def get_ids_frm_signalp_out(sigp):
    accs = []
    for i in sigp:
        if i[0] != '#':
            col = i.strip('\n').split()
            if col[9] == 'Y':
                accs.append(col[0])
    return accs




def fasta_extract(fasta,accs,outfasta):
    inf = SeqIO.parse(fasta, "fasta")
    acclist = accs
    for i in inf:
        for a in acclist:
            if i.id == a:
                SeqIO.write(i,outfasta, "fasta")
    outfasta.close()

ac = get_ids_frm_signalp_out(insignalp)

fasta_extract(fs,ac,of)