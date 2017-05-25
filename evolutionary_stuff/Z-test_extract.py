__author__ = 'mjohnpayne'


import numpy
import glob
import matplotlib.pyplot as plt

pat = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/Z-test_output'

indata = glob.glob(pat + '/*.meg')

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/Z-test-stats.txt','w')

pfpml = []
pftfl = []
pmtfl = []
pftsl = []
pmtsl = []
tftsl = []

name = ['pfpm','pftf','pmtf','pfts','pmts','tfts']
outfile.write('ID\t' + '\t'.join(name)+'\n')


c = 0
for i in indata:
    b=0
    # if c%1000 == 0:
    #    print c
    infile = open(i,'r').read()
    pos = infile.index('#Pma|')
    pmaa = infile[pos+5:pos+17].strip('\r\n')
    if infile[pos-3] != '2':
        b=1
    infile = infile.split('\r\n')
    l1 = infile[35][infile[35][6:].index('[')+7:].strip('\r\n')
    l1 = l1.replace(']','').split('[')
    l2 = infile[36][infile[36][6:].index('[')+7:].strip('\r\n')
    l2 = l2.replace(']','').split('[')
    l3 = infile[37][infile[37][6:].index('[')+7:].strip('\r\n')
    l3 = l3.replace(']','')
    pfpm = float(l1[0])#float(infile[35][21:32].strip('\n'))
    pfpml.append(pfpm)
    pftf = float(l1[1])#float(infile[35][37:49].strip('\n'))
    pftfl.append(pftf)
    pmtf = float(l2[0])#float(infile[36][37:49].strip('\n'))
    pmtfl.append(pmtf)
    pfts = float(l1[2])#float(infile[35][52:64].strip('\n'))
    pftsl.append(pfts)
    pmts = float(l2[1])#float(infile[36][52:64].strip('\n'))
    pmtsl.append(pmts)
    tfts = float(l3)#float(infile[37][52:64].strip('\n'))
    tftsl.append(tfts)
    gene = [pfpm,pftf,pmtf,pfts,pmts,tfts]
    if b==1:
        outfile.write(pmaa + '_mixed\t' + '\t'.join(map(str,gene))+'\n')
    else:
        outfile.write(pmaa + '\t' + '\t'.join(map(str,gene))+'\n')
    # if pmtf > -2:
    #     print pmaa
    c +=1

inf = [pfpml,pftfl,pmtfl,pftsl,pmtsl,tftsl]
data = {'pfpm':pfpml,'pftf':pftfl,'pmtf':pmtfl,'pfts':pftsl,'pmts':pmtsl,'tfts':tfts}
name = ['pfpm','pftf','pmtf','pfts','pmts','tfts']
av = []
std = []
for x in data:
    # name += [x]
    av += [numpy.average(data[x])]
    std += [numpy.std(data[x])]

ax = plt.axes()
plt.boxplot(inf)
ax.set_xticklabels(name)
plt.savefig('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/boxplots_z-score.pdf',dpi=400)

# c = 1
# for i in infile:
#     print c
#     print itfts
#     c+=1
