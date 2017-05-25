__author__ = 'mjohnpayne'
import sys

infile = open(sys.argv[1],"r")
outfile = open(sys.argv[2],"w")

for i in infile:
    n = i.replace("\\","")
    outfile.write(n)

outfile.close()