__author__ = 'mjohnpayne'

import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq

acc = sys.argv[1]
out = open(sys.argv[2],'a')
intpro_files = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan'
fasta_files = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/compliant'


intpro_ls = glob.glob(intpro_files + '/*.tsv')

fasta_ls = glob.glob(fasta_files + '/*')

fastas = {}
for f in fasta_ls:
    for record in SeqIO.parse(f, "fasta"):
        fastas[record.id] = record


iprs = {}
for f in intpro_ls:
    intpr = open(f,'r')
    for line in intpr:
        col = line.split('\t')
        if len(col) > 11:
            ipr = col[11]
            if ipr not in iprs:
                iprs[ipr] = [col[0]]
            else:
                iprs[ipr] += [col[0]]
    intpr.close()
iprs2 = {}
for i in iprs:
    iprs2[i] = set(iprs[i])

seqs = []
for i in iprs2[acc]:
    seqs += [fastas[i]]

SeqIO.write(seqs,out,"fasta")

out.close()