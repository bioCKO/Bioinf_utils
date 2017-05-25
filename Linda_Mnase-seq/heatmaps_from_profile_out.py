__author__ = 'mjohnpayne'


import matplotlib.pyplot as plt
import numpy as np
from time import sleep as sl
import itertools

# Generate Data

data = open("/Volumes/MP_HD/Linda_MNase_Seq/danpos_r1_out/profile/tRNAs/tRNA_genes_diff_TSS_heatmap/tRNA_genes_tss_11H-15H_heatmap.csv",'r').read().split('\r')

datals = []
colnames = data[0].split(',')[4:]
rownames = []
for i in data[1:]:
    i=i.split(',')
    rownames.append(i[0])
    datals.append(map(float,i[4:]))
dataar = np.asarray(datals)
print rownames[:20]
print colnames[:20]
print dataar[:10]

# data=temp1.values
# test=data.reshape(7,24)
# #data = np.random.rand(7,24)
rows = rownames
columns = ' '.join(colnames).split()[::50]
print columns
pos = np.asarray([x+150 for x in map(int,columns)])
print pos


plt.pcolor(dataar,cmap=plt.cm.RdYlBu,vmin=0,vmax=50)
plt.colorbar()
plt.yticks([])
plt.xticks(np.asarray(pos)+0.5,columns)
plt.show()
plt.close()