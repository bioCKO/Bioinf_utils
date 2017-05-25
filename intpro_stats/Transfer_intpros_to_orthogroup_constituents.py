__author__ = 'mjohnpayne'


## 1 link all interpros to respective pmaa genes

## 2 link each pmaa gene to orthologues

## 3 transfer interpros from all pmaa to their orthologues

import sys
from time import sleep as sl

clusterinf = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/Tm_CI_orthogroups.txt','r')

intpro = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/CI_orthogroup_annots.txt','r')


def intpro_annot_to_dict(intproin):
    intpro_dict = {}
    for i in intproin:
        col = i.strip('\n').split('\t')
        group = col[0]
        intrpros = [x[:9] for x in col[1].split(";")]
        gos = col[2].split(';')
        intpro_dict[group] = (intrpros,gos)
    return intpro_dict


def ortho_clust_to_dict(inf):
    clust = {}
    for line in inf:
        col = line.strip('\n').split(' ')
        clust[col[0][:-1]] = [x.split('|')[:2] for x in col[1:]]
    return clust


def annot_genes(clusters,diction):
    gene_ints = {}
    gene_gos = {}
    for i in clusters:
        for j in clusters[i]:
            gene_ints[j] = diction[i][0]
            gene_gos[j] = diction[i][1]
    return gene_ints,gene_gos

def counts_per_strain()

def main(cl,int):
    clu = ortho_clust_to_dict(cl)
    cludict = intpro_annot_to_dict(int)
    genint,gengo = annot_genes(clu,cludict)



annotations,gos = main(clusterinf,intpro)