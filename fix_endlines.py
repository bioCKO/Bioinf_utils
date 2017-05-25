__author__ = 'mjohnpayne'

import sys

infile = "/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/enrichment_trial/point<0.005_ext_prom_list.txt"

def fixlines(inpath):
    inf = open(inpath,'r').read()
    outf = '/'.join(inpath.split('/')[:-1])+'/' + inpath.split('/')[-1][:-4] + "_fix_newlines" + inpath.split('/')[-1][-4:]
    outfile = open(outf,'w')
    inf=inf.replace("\r","\n")
    outfile.write(inf)
    outfile.close()
    print inpath
    print outf

fixlines(infile)