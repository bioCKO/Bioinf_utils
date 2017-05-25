__author__ = 'mjohnpayne'

import numpy
import glob
import matplotlib.pyplot as plt

pat = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/fishers_exact_prop_out'

indata = glob.glob(pat + '/*.meg')


pfpml = []
pftfl = []
pmtfl = []
pftsl = []
pmtsl = []
tftsl = []

c = 0
for i in indata:
    if c%1000 == 0:
        print c
    infile = open(i,'r').readlines()
    pmaa = infile[25][9:].strip('\n')
    pfpm = float(infile[31][6:16].strip('\n'))
    pfpml.append(pfpm)
    pftf = float(infile[32][6:16].strip('\n'))
    pftfl.append(pftf)
    pmtf = float(infile[32][17:27].strip('\n'))
    pmtfl.append(pmtf)
    pfts = float(infile[33][6:16].strip('\n'))
    pftsl.append(pfts)
    pmts = float(infile[33][17:27].strip('\n'))
    pmtsl.append(pmts)
    tfts = float(infile[33][28:38].strip('\n'))
    tftsl.append(tfts)
    gene = [pfpm,pftf,pmtf,pfts,pmts,tfts]
    if sum(gene)<6:
        print pmaa
        print gene
    c +=1

# inf = [pfpml,pftfl,pmtfl,pftsl,pmtsl,tftsl]
# data = {'pfpm':pfpml,'pftf':pftfl,'pmtf':pmtfl,'pfts':pftsl,'pmts':pmtsl,'tfts':tfts}
name = ['pfpm','pftf','pmtf','pfts','pmts','tfts']
# av = []
# std = []
# for x in data:
#     # name += [x]
#     av += [numpy.average(data[x])]
#     std += [numpy.std(data[x])]
#
# ax = plt.axes()
# plt.boxplot(inf)
# ax.set_xticklabels(name)
# plt.show()#plt.savefig('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/boxplots_syn_mut_prop.pdf',dpi=400)

print name
# c = 1
# for i in infile:
#     print c
#     print itfts
#     c+=1
