__author__ = 'mjohnpayne'
## take orthomcl output and print each gene with a list of orthologues/paralogues

from time import sleep as sl

def rn(gene):
    if 'Tfl' in gene:
        gene = gene[7:]
        no = 5-len(gene)
        gene = 'TFLA_' + no*'0' + gene + '0'
    # elif "g" in gene:
    #     gene = gene[1:]
    #     no = 5-len(gene)
    #     gene = 'TFLA_' + no*'0' + gene + '0'
    elif 'Pfu' in gene:
        gene = gene[7:]
        no = 5-len(gene)
        gene = 'PFUN_' + no*'0' + gene + '0'
    else:
        gene = gene[4:]
    return gene

infile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups.txt','r')

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/4_species_eurot_groups_by_gene.txt','w')

dict = {}

for i in infile:
    ls = i.strip('\n').split(' ')
    genes = ls[1:]
    genes = [rn(x) for x in genes]
    groupid = ls[0][:-1]
    for j in range(len(genes)):
        gene = genes[j]
        list = genes[:j] + genes[j+1:]
        acc = gene[:4]
        para = [t for t in list if acc in t]
        para = [s[s.find("|")+1:] for s in para]
        orthos = [t for t in list if acc not in t]
        orthos = [s[s.find("|")+1:] for s in orthos]
        dict[gene[gene.find("|")+1:]] = [para,orthos,groupid]

outfile.write("GeneID\tGroupID\tParalogues\tOrthologues\n")

poss = ["PMAA","TSTA","PFUN","TFLA"]





for j in dict:
    if j[:4] in poss:
        outfile.write(j + '\t' + dict[j][2] +'\t' +','.join(dict[j][0]) + '\t' + ','.join(dict[j][1]) + '\n')

outfile.close()