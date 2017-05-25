__author__ = 'mjohnpayne'

import sharepathway as sp


filein="/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/enrichment_trial/point<0.005_ext_prom_list_fix_newlines.txt"

fileout="/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/enrichment_trial/point<0.005_ext_prom_list_enrich_out.txt"



sp.Run(fi=filein,fo=fileout,species='ani',r=0.1)