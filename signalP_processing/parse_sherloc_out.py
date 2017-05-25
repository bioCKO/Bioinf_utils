__author__ = 'mjohnpayne'

import glob
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec
from pylab import *
from numpy import array
import numpy as np


inlist = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/*_result.txt")

outf = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/"

def parse_sheloc2out(inls,outpref):
    outd = {}
    for i in inls:
        species = i.split("/")[-1][:-20]
        inf = open(i,'r').readlines()[4:]
        types = {}
        outfile = open(outpref + species + "_sherloc2_parsed_out.txt","w")
        outfile.write("GeneID\ttop_prediction\tscore")
        names = []
        for line in inf:
            col= line.split('\t')
            name = col[0]
            if name not in names:
                pos = col[1].find(":")
                top = col[1][:pos]
                prob = float(col[1][pos+2:])
                outfile.write("%s\t%s\t%s\n" %(name,top,str(prob)))
                if prob > 0.5:
                    if top not in types:
                        types[top] = [name]
                    else:
                        types[top].append(name)
                else:
                    if "Inconclusive" not in types:
                        types["Inconclusive"] = [name]
                    else:
                        types["Inconclusive"].append(name)
                names.append(name)
        outfile.close()
        typels = []
        countls = []
        outsum = open(outpref + species + "_sherloc2_summary.txt","w")
        outsum.write("Location\tCount\n")
        for j in types:
            typels.append(j)
            ct = len(types[j])
            outsum.write(j + '\t' + str(ct) + '\n')
            countls.append(ct)
        outsum.close()
        outd[species] = countls


    # for i in
    # return outd,typels


parse_sheloc2out(inlist,outf)
