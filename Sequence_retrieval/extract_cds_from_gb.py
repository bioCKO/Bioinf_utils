import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq


gb_path = sys.argv[1]
outfile = open(sys.argv[2],'w')

in_gb = open(gb_path,'rU')

sequences = SeqIO.parse(in_gb, "genbank")

gene = ''
count = 0
for record in sequences:
    for i,feature in enumerate(record.features):
        try:
            if feature.type == 'gene':
               pmaa = feature.qualifiers['locus_tag'][0]
               count += 1
               #print count
               #print pmaa
            elif feature.type == 'CDS':
                cdsseq = feature.extract(record.seq)
                pmaa = feature.qualifiers['protein_id'][0]
                #cdsseq = str(Seq(str(cdsseq)).translate())
                outfile.write('>' + str(pmaa) + '\n' + str(cdsseq) + '\n')
        except:
            pass

outfile.close()
