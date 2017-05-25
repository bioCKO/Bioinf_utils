__author__ = 'mjohnpayne'


import re

inf = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf_reorder.gb'

infile = open(inf,'r').readlines()
outfile = open(inf[:-3] + '_add_tag.gb','w')

print infile[0:100]
out = []

for i in infile:
     if '                     /ID' in i:
         out.append(i)
         ids = i[26:37]
         tag = "                     /locus_tag=\"" + ids + "\"\n"
         out.append(tag)
     else:
         out.append(i)

outfile.write(''.join(out))

outfile.close()