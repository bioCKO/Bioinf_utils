__author__ = 'mjohnpayne'

import re

inf = '/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gb'

infile = open(inf,'r').read()
outfile = open(inf[:-3] + '_fix.gb','w')

join = [m.start() for m in re.finditer('join\(join\(', infile)]



infile = infile.split('join(join(')
new = []
for i in infile:
    i = i.replace('))',')',1)
    new.append(i)

infile = 'join('.join(new)

outfile.write(infile)

outfile.close()