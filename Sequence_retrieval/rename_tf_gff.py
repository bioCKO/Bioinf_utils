__author__ = 'mjohnpayne'



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

infile = open("")