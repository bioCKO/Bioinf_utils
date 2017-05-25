__author__ = 'mjohnpayne'

import sys



intsv = open(sys.argv[1],'r')

string = sys.argv[2]

def rename(gene):
    if "|" in gene:
        gene = gene[4:]
        if 'TF' in gene:
            no = 5-len(gene[3:])
            gene = 'TFLA_' + no*'0' + gene[3:] + '0'
        elif 'Pf' in gene:
            no = 5-len(gene[3:])
            gene = 'PFUN_' + no*'0' + gene[3:] + '0'
        elif "ANID" in gene:
            gene = gene[:10]
        elif "Pc_gi_" in gene:
            gene = gene[16:]
    else:
        if 'TF' in gene:
            no = 5-len(gene[3:])
            gene = 'TFLA_' + no*'0' + gene[3:] + '0'
        elif 'Pf' in gene:
            no = 5-len(gene[3:])
            gene = 'PFUN_' + no*'0' + gene[3:] + '0'
        elif "ANID" in gene:
            gene = gene[:10]
        elif "Pc_gi_" in gene:
            gene = gene[16:]
    return gene


for i in intsv:
    if string in i:
        col = i.strip('\n')
        col = col.split('\t')
        id = rename(col[0])
        print id
