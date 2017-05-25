__author__ = 'mjohnpayne'

import glob
from time import sleep as sl



inlist = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/TMHMM/*_tmhmm_web_out.txt")

outf = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/TMHMM/"

def parse_tmhmmout(inls,outpref):
    for i in inls:
        inf = open(i,"r").readlines()
        spec = i.split('/')[-1][:-18]
        out = open(outf+spec+"tmhmm_parsed.txt","w")
        out.write('GeneID\tnumber helices\n')
        for j in inf:
            col = j.split('\t')
            pmaa = col[0]
            hel_no = col[4][8:]
            out.write(pmaa+'\t'+hel_no+'\n')
        out.close()


parse_tmhmmout(inlist,outf)