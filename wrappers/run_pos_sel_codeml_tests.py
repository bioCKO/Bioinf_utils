import subprocess
import shlex
import sys
import re
import time
import glob

start = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/inparanoids/codeml_pos_sel_tests/startingpoint','r')

alnin = sys.argv[1]
treein = sys.argv[2]

aln_lis = glob.glob(alnin + '/*')
nalnlist = []
for i in aln_lis:
    if 'original1' not in i and '.phy' in i:
        nalnlist.append(i)
aln_lis = nalnlist
aln_lis.sort()


tree_lis = glob.glob(treein + '/*')
nalnlist = []
for i in tree_lis:
    if 'original1' not in i and '.dnd' in i:
        nalnlist.append(i)
tree_lis = nalnlist
tree_lis.sort()

st = start.read()
st = int(st)
start.close()
print '\n\n\n\n' + str(st)

for j in range(st,len(tree_lis)):
    inalign = aln_lis[j]
    intree = tree_lis[j]
    out = intree.split('/')[-1][:-13] + '_pos_sel_results.txt'
    print out
    codeml_inp = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/inparanoids/codeml_pos_sel_tests/sample_codeml.ctl.txt','r').readlines()
    codeml_ctl = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/inparanoids/codeml_pos_sel_tests/outputs/tempcodeml.ctl','w')
    codeml_out = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/inparanoids/codeml_pos_sel_tests/outputs/' + out
    for i in range(len(codeml_inp)):
        if i == 0:
            codeml_ctl.writelines('      seqfile = ' + inalign + '\n')
        elif i == 1:
            codeml_ctl.writelines('      treefile = ' + intree + '\n')
        elif i == 2:
            codeml_ctl.writelines('      outfile = ' + codeml_out + '\n')
        else:
            codeml_ctl.writelines(codeml_inp[i])
    codeml_ctl.close()
    codeml_args = '/Users/mjohnpayne/Documents/PhD/bioinformatics_tools/paml4.7/bin/codeml /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/inparanoids/codeml_pos_sel_tests/outputs/tempcodeml.ctl'
    subprocess.Popen(codeml_args, shell=True).wait()
    st += 1
    start = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/inparanoids/codeml_pos_sel_tests/startingpoint','w')
    start.write(str(st))
    start.close()
    
    
    
