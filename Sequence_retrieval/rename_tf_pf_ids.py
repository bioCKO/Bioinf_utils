__author__ = 'mjohnpayne'

import sys

infile = sys.argv[1]
outfile = sys.argv[2]
outf = open(outfile,"w")


def rn(gene):
    if 'TF' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif "g" in gene:
        gene = gene[1:]
        no = 5-len(gene)
        gene = 'TFLA_' + no*'0' + gene + '0'
    elif 'Pf' in gene:
        gene = gene[3:]
        no = 5-len(gene)
        gene = 'PFUN_' + no*'0' + gene + '0'
    else:
        gene = gene[4:]
    return gene

with open(infile,"r") as f:
    for line in f:
        if ">" in line:
            outf.write(">" + rn(line[1:].strip())+'\n')
        else:
            outf.write(line)

outf.close()