__author__ = 'mjohnpayne'

import re

inf = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf_reorder.gb'

infile = open(inf,'r')
outfile = open(inf[:-3] + '_fix.gb','w')

for line in infile:
    if line.count("..") == 1 and "join(" in line:
        line = line.replace("join(","")
        if "))" in line:
            line = line.replace("))",")")
        else:
            line = line.replace(")","")
        outfile.write(line)
    else:
        outfile.write(line)

outfile.close()
