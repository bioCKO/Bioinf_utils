__author__ = 'mjohnpayne'


import sys

infile = open(sys.argv[1],'r')

outfile = open(sys.argv[2],'w')

outfile.write("!gaf-version: 2.0\n")

for i in infile:
    cols = i.split("\t")
    if len(cols) > 4:
        if cols[12] == "taxon:162425":
            outfile.write(i)

outfile.close()
