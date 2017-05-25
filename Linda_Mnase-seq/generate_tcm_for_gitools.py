__author__ = 'mjohnpayne'


import sys

infile = open(sys.argv[1],'r')

genecol = int(sys.argv[2])
annotcol = int(sys.argv[3])

outfile = open(sys.argv[4],'w')

for i in infile:
    cols = i.split("\t")
    if len(cols) > 4:
        if cols[12] == "taxon:162425":
            gene = cols[10].split("|")[0]
            annot = cols[annotcol]
            outfile.write(gene+'\t'+annot+'\n')

outfile.close()
