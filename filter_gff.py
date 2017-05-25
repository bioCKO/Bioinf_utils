__author__ = 'mjohnpayne'

import sys

infile = open(sys.argv[1],'r')

column = int(sys.argv[2])-1

threshold = float(sys.argv[3])

outfile = open(sys.argv[4],'w')

if len(sys.argv) < 5:
    print "script.py infile column threshold outfile"

for i in infile:
    col = i.strip('\n').split('\t')
    if float(col[column]) > threshold:
        outfile.write(i)
outfile.close()



