__author__ = 'mjohnpayne'


import glob

inls = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/tm_outputs/*_sherloc2_out.txt")

outfile = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/tm_sherloc2_result.txt","w")

outfile.write("SherLoc2 Prediction Result\n\nOrigin = fungal\n\n")

for i in inls:
    inp = open(i,"r").read().split('\n')[4]
    outfile.write(inp+"\n")

outfile.close()