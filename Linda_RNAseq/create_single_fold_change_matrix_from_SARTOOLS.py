__author__ = 'mjohnpayne'


import glob
from time import sleep as sl


inlis = glob.glob("/Volumes/MP_HD/Linda_MNase_Seq/RNAseq/SARTools_diff_analysis/All/tables/*.complete.txt")

outfile = open("/Volumes/MP_HD/Linda_MNase_Seq/RNAseq/SARTools_diff_analysis/All/tables/all_norm_counts_with_FC.txt",'w')
count = 0
row1 = ''
values = {}
for i in inlis:
    comp = i.split('/')[-1].strip(".complete.txt")
    inf = open(i,'r').readlines()
    if count < 1:
        print i
        head = inf[0].split('\t')
        row1 = "Gene_ID\t"+head[45]+"\t"+head[46]+"\t"+"\t".join(head[19:37])+"\t" + "\t".join(head[38:45])
        for j in inf[1:]:
            vals = j.split('\t')
            values[vals[0]] = [vals[45],vals[46]] + vals[19:37] + vals[38:45]
        count = 1
    head = inf[0].split('\t')
    row1+="\t" + comp + "_" + head[47] + "\t" + comp + "_" + head[50]
    for x in inf[1:]:
        vals = x.split('\t')
        if vals[50] == "NA":
            values[vals[0]]+= [vals[47],1]
        else:
            values[vals[0]]+= [vals[47],vals[50]]
row1 += '\n'
print row1
print values["AN0005"]

outfile.write(row1)

for i in list(sorted(values.keys())):
    outfile.write(i+'\t'+"\t".join(map(str,values[i]))+'\n')

outfile.close()


