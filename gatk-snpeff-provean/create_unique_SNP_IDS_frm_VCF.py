__author__ = 'mjohnpayne'

import sys

inv = sys.argv[1]
invcf = open(inv,'r')
outfile = open(sys.argv[2],'w')

for line in invcf:
    if '#' not in line:
        col = line.split('\t')
        a_freq = col[-1].split(';')[1].replace('AF=','')
        name = col[0] + "_" + col[1] + '_' + col[3] + col[4]
        outfile.write(name + '\n')

outfile.close()