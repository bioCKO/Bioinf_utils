__author__ = 'mjohnpayne'

# 1 - generate dict of clusters with clusterid as key and gene ids as data
# 2 - generate dict of interpro ID:description objects with gene id as key and id-descriptions as data
# 2 - generate dict of GOTermobjects with gene id as key and id-descriptions as data
# 3 - merge by having clusterid link to collapsed list of interpro ids (set)
# 4 - print out

import sys
from time import sleep as sl

clusterinf = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups.txt','r')

intpro = ['/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/Pma_f.fasta.tsv','/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/Tfl_f.fasta.tsv','/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/Pfu_f.fasta.tsv','/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/Tst_f.fasta.tsv']


def rename(gene):
    if 'TF' in gene:
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif 'Pf' in gene:
        no = 5-len(gene[3:])
        gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    return gene


def ortho_clust_to_dict(inf):
    clust = {}
    for line in inf:
        col = line.strip('\n').split(' ')
        clust[col[0][:-1]] = [rename(x[4:]) for x in col[1:]]
    return clust



def intpro_annot_to_dict(intproin):
    ipr = {}
    go = {}
    for f in intproin:
        tmp = open(f,'r')
        for line in tmp:
            col = line.strip('\n').split('\t')
            gene = rename(col[0][4:])
            if len(col) > 11:
                if gene not in ipr:
                    ipr[gene] = [col[11] + ' = ' + col[12]]
                else:
                    ipr[gene].append(col[11] + ' = ' + col[12])
            if len(col) > 13:
                if gene not in go:
                    go[gene] = col[13].split('|')
                else:
                    go[gene] += col[13].split('|')
    nipr = {}
    ngo = {}
    for i in ipr.keys():
        nipr[i] = list(set(ipr[i]))
        if len(nipr[i]) == 0:
            nipr[i] = ['none']
    for i in go.keys():
        ngo[i] = list(set(go[i]))
        if len(ngo[i]) == 0:
            ngo[i] = ['none']

    return nipr,ngo

def merge(orthos,ints,gos):
    annots = {}
    for i in orthos:
        iprs = []
        go = []
        for j in orthos[i]:
            if j in ints:
                iprs += ints[j]
                go += gos[j]
        iprs = list(set(iprs))
        go = list(set(go))
        niprs = []
        ngo = []
        if '' in go:
            del(go[0])
        if len(iprs) == 0:
            niprs = ['none']
        else:
            niprs = iprs
        if len(go) == 0:
            ngo = ['none']
        else:
            ngo = go
        annots[i] = (niprs,ngo)
    return annots



def main(cl,int):
    clu = ortho_clust_to_dict(cl)
    interpro,goterms = intpro_annot_to_dict(int)
    annots = merge(clu,interpro,goterms)
    return annots,goterms

annotations,gos = main(clusterinf,intpro)

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/ortho_groups_by_species/4_species_ortho_group_annots.txt','w')

for i in annotations:
    outfile.write(i + '\t' + ';'.join(annotations[i][0]) + '\t' + ';'.join(annotations[i][1]) + '\n')

outfile.close()