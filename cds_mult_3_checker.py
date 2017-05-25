__author__ = 'mjohnpayne'

from Bio import SeqIO
from Bio.Seq import Seq


incds = '/Users/mjohnpayne/Documents/PhD/bys/bys_orthologue_stuff/2nd_round_with_unannot/complete_set_bys_cds_with_4spec/complete_set_bys_cds_with_4spec.fasta'

cds = {}

for record in SeqIO.parse(incds, "fasta"):
        cds[record.id] = record.seq


for i in cds:
    print i
    print float(len(cds[i]))/3