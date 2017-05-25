__author__ = 'mjohnpayne'

import sys
import re

infile = open(sys.argv[1],'r').readlines()

outfile = open(sys.argv[2],'w')


outfile.write('\n'.join(infile[:3]))

types = ["Gypsy","Mariner","Tad1","Copia","AFUT1","I-5","Harbinger","I-1","I-2","I-3","I-4","I-6","PiggyBac"]

for i in infile:
    for j in types:
        if j in i:
            #print i[:re.search(j, i,).start()+len(j)]+"\""
            outfile.write(i[:re.search(j, i,).start()+len(j)]+"\"" + '\n')

outfile.close()
