__author__ = 'mjohnpayne'

import sys,re
from time import sleep as sl
inf = open(sys.argv[1],'r')
#inf = open('/Volumes/MP_HD/Blast2Go/Tf/Tfl_f.fasta.xml','r')
outfile = open(sys.argv[2],'w')

def rn(gene):
    if 'TF' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif 'Pf' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    return gene


def rename(file,change):
    newfile = []
    ofile = file.readlines()
    for i in ofile:
        if change in i:
            pos = i.find(change) #+ 1
            end = i.find('\t')
            acc = i[pos:end]
            nacc = rn(acc)
            nline = i[:pos] + nacc + i[end:]
            newfile.append(nline)
        else:
            newfile.append(i)
    return ''.join(newfile)

    # no = 5-len(gene[3:])
    # gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    # elif 'Pf' in gene:
    #     no = 5-len(gene[3:])
    #     gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    # return gene



outfile.write(rename(inf,'TF_'))

outfile.close()