__author__ = 'mjohnpayne'

### fix multiple annotations per PMAA so that each pmaa has only one output


import sys
from time import sleep as sl
import re
cazyout = open(sys.argv[1],'r').readlines()

cazy_class = open(sys.argv[2],'r').read().split('\r')

class_dict = {}

out = open(sys.argv[3],'w')
out.write("GeneID\tdbCAN ID\tOther Clustering\tdbCAN description\n")

for i in cazy_class[1:]:
    col = i.split('\t')
    class_dict[col[0]] = [col[1],col[4]]

inpdict = {}

for j in cazyout:
    col = j.split('\t')
    eval = float(col[2])
    if eval < float(1e-10):
        if col[1][:-4] not in class_dict:
            out.write(col[0] + '\t' + col[1][:-4] + '\tNo Pfam\tNo Description\n')
        else:
            out.write(col[0] + '\t' + col[1][:-4] + '\t' + class_dict[col[1][:-4]][0] + '\t' + class_dict[col[1][:-4]][1] + '\n')


out.close()