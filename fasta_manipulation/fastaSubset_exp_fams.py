from Bio import SeqIO
import sys
from time import sleep as sl

fasta_file = sys.argv[1]#"/Users/lukemn/Documents/prc/pABH1_13Illumina/velvet/velk40/contigs.fa"
#subset_file = sys.argv[2]#"/Users/lukemn/Documents/prc/pABH1_13Illumina/velvet/velk40/contigRetrieve.txt"
output = sys.argv[3]#"/Users/lukemn/Documents/prc/pABH1_13Illumina/velvet/velk40/Retrieved.fa"
newnames = sys.argv[2]


wanted = set()
with open(newnames) as f:
    for line in f:
        line = line.strip()
        col = line.split('\t')
        #print line
        if line != "":
            #print line
            wanted.add(col[0])

newnames = open(sys.argv[2],'r')
newname = {}
for i in newnames.readlines():
    col = i.strip('\n').split('\t')
    if len(col) > 2:
        newname[col[0]] = col[0] + "_" + col[1] + "_(" + col[2] + ")"
    else:
        newname[col[0]] = col[0] + "_" + col[1]

fasta_sequences = SeqIO.parse(open(fasta_file),"fasta")
with open(output, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            seq.id = newname[seq.id]
            seq.name = ""
            seq.description = ""
            SeqIO.write(seq, f, "fasta")

fasta_sequences.close()
