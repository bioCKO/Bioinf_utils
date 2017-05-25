__author__ = 'mjohnpayne'

### cd /Users/mjohnpayne/Documents/PhD/bioinformatics_tools/SherLoc2-26-10-2009


from Bio import SeqIO
import os
import sys
import subprocess
from time import sleep as sl

inf = sys.argv[1]

infasta = SeqIO.parse(inf,"fasta")


name = inf.split('/')[-1][:-6]
donels = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/python_tmp/"+name+"_donels.txt"

if not os.path.exists(donels):
    file(donels, 'w').close()

c = 0
for i in infasta:
    done = open(donels,"r").read().split('\n')
    print c, i.id
    c += 1
    if i.id not in done:
        fs = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/python_tmp/"+name+"_tmp.fasta"
        SeqIO.write(i,fs,"fasta")
        sherlocargs = "/usr/local/bin/python ./src/sherloc2_prediction.py -fasta="+fs+" -origin=fungal -result=/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/sherloc2_outputs/tm_outputs/" + i.id + "_sherloc2_out.txt"
        subprocess.Popen(sherlocargs, shell=True).wait()
        done = open(donels,"a")
        done.write(i.id+"\n")
        done.close()


