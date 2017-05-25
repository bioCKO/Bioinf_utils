__author__ = 'mjohnpayne'

import sys

infile = open(sys.argv[1],'r')

column = int(sys.argv[2])-1

type_to_remove = sys.argv[3]

outfile = open(sys.argv[4],'w')

if len(sys.argv) < 5:
    print "script.py infile column threshold outfile"

for i in infile:
    if i[0] != "#":
        col = i.strip('\n').split('\t')
        if col[column] != type_to_remove:
            outfile.write(i)
    else:
        outfile.write(i)
outfile.close()



