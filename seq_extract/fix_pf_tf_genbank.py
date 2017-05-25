__author__ = 'mjohnpayne'

from Bio import SeqIO

in_gb = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf_reorder.gb'

sequences = SeqIO.parse(in_gb, "genbank")

pmaa = ''

for record in sequences:
    geneno = 0
    for i,feature in enumerate(record.features):
        if feature.type == 'gene':

           geneno = 0
        elif feature.type == 'transcript':
            print feature
        elif feature.type == 'CDS':
            if geneno == 0:
                #print feature.location
                geneno = 1
