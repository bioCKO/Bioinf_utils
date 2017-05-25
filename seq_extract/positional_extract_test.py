import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

gb_path = '/Users/mjohnpayne/Documents/PhD/random_python_scripts/seq_extract/Pm_from_gb_old_acc.gb'
in_gb = open(gb_path,'rU')

sequences = SeqIO.parse(in_gb, "genbank")

contig = ''
for record in sequences:
    contig = record.id
    for i,feature in enumerate(record.features):
#        if feature.type == 'source':
#            print feature.id
        if feature.type == 'gene':
           geneseq = feature.extract(record.seq)
           pmaa = feature.qualifiers['locus_tag'][0]
#           print contig
           pos = str(feature.location)
           print pos
           pos = pos.translate(None, '<>[]()')
           print pos
           pos = pos.split(':')
           print pos
#           gene[pmaa] = geneseq
