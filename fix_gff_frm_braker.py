__author__ = 'mjohnpayne'

import sys

infile = open(sys.argv[1],'r')

outfile = open(sys.argv[2],'w')

gene=""
mrna=""
for i in infile:
    col = i.strip('\n').split('\t')
    if i[0] != "#":
        if col[2] == "gene":
            gene = col[8]
            col[8] = "ID=" + col[8]
            outfile.write("\t".join(col)+'\n')
        elif col[2] == "transcript":
            col[2] = 'mRNA'
            mrna = col[8]
            col[8] = "ID=" + col[8] + ";Parent=" + gene
            outfile.write("\t".join(col)+'\n')
        else:
            col[8] = col[8].split(";")[0].replace("transcript_id ","Parent=").replace("\"","")
            outfile.write("\t".join(col) + '\n')
    else:
        outfile.write(i)

outfile.close()



