__author__ = 'mjohnpayne'

import xlrd
import csv
from os import sys
from time import sleep as sl
import sys


inf = sys.argv[1]
infile = open(inf,"r").readlines()
outf = open(inf.replace(".xls",".txt"),"w")


for i in infile:
    col = i.split('\t')
    outf.write(col[0] + '\t' + "\t".join(col[4:]))

outf.close()