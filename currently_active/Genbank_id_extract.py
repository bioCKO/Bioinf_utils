__author__ = 'mjohnpayne'

from Bio import SeqIO
from Bio.Alphabet import IUPAC
import subprocess
ids = raw_input('enter sequence id: ')#PMAA_101030'
#spec = raw_input('enter species identifier (Pm,Ts,Tf,Pf): ')#Pm'

left = 4000
right = 4000
in_gb = ''

if ids[:2] =='PM':
    in_gb = SeqIO.parse('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pm_genbank_gene_containing_only.gb','genbank')
elif ids[:2] == 'TS':
    in_gb = SeqIO.parse('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Ts_reorder.gb','genbank')
elif ids[:2] == 'TF':
    in_gb = SeqIO.parse('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Tf_reorder.gb','genbank')
elif ids[:2] == 'PF':
    in_gb = SeqIO.parse('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf_reorder.gb','genbank')


outseq = open('/Users/mjohnpayne/Desktop/HW_sod_loci/' + ids + '.fasta','w')

for i in in_gb:
    for j in i.features:
        if j.type == 'gene':
            if ids == j.qualifiers['locus_tag'][0]:
                print i.id
                loc = str(j.location).replace('>','')
                loc = loc.replace('<','')
                st = int(loc[1:loc.find(':')])
                en = int(loc[loc.find(':')+1:loc.find(']')])
                outseq.write('>' + str(ids) + '\n' + str(i.seq[st-left:en+right]) + '\n')
