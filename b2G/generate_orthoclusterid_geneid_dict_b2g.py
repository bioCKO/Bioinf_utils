__author__ = 'mjohnpayne'

# 1 - generate dict of clusters with clusterid as key and gene ids as data
# 2 - generate dict of interpro ID:description objects with gene id as key and id-descriptions as data
# 2 - generate dict of GOTermobjects with gene id as key and id-descriptions as data
# 3 - merge by having clusterid link to collapsed list of interpro ids (set)
# 4 - print out

import sys
from time import sleep as sl

clusterinf = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups.txt','r')

descriptions = ['/Volumes/MP_HD/Blast2Go/Pm/exported_data/Pm_b2g_descriptions.txt','/Volumes/MP_HD/Blast2Go/Tf/exported data/Tf_b2g_descriptions.txt','/Volumes/MP_HD/Blast2Go/Pf/exported_data/Pf_b2g_descriptions.txt','/Volumes/MP_HD/Blast2Go/Ts/exported data/Ts_b2g_descriptions.txt']


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



def desc_to_dict(intproin):
    sd = {}
    gt = {}
    enz = {}
    kgm = {}
    for f in intproin:
        tmp = open(f,'r').readlines()
        for line in tmp[1:]:
            col = line.strip('\r\n').split('\t')
            gene = col[0]#rename(col[0][4:])
            sd[gene] = col[1]
            gt[gene] = col[2]
            enz[gene] = col[3]
            kgm[gene] = col[4]
    return sd,gt,enz,kgm

def merge(orthos,sd,gt,enz,kgm):
    annots = {}
    for i in orthos:
        sdd = []
        gtd = []
        enzd = []
        kgmd = []
        for j in orthos[i]:
            if j in sd:
                sdd.append(sd[j])
            if j in gt:
                gtd.append(gt[j])
            if j in enz:
                enzd.append(enz[j])
            if j in kgm:
                kgmd.append(kgm[j])
        sdl = list(set(sdd))
        gtl = list(set(gtd))
        enzl = list(set(enzd))
        kgml = list(set(kgmd))
        lst = [sdl,gtl,enzl,kgml]
        for k in lst:
            if len(k) == 1 and k[0] == '':
                k[0] = 'no annot'
        annots[i] = lst
    return annots



def main(cl,int):
    clu = ortho_clust_to_dict(cl)
    sd,gt,enz,kgm = desc_to_dict(int)
    annots = merge(clu,sd,gt,enz,kgm)
    return annots

annotations = main(clusterinf,descriptions)

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/ortho_groups_by_species/4_species_ortho_group_b2g_annots.txt','w')

outfile.write('ID\tDescription\tGOterms\tEnzymes\tKeggmap\n')

for i in annotations:
    outfile.write(i + '\t' + ';'.join(annotations[i][0]) + '\t' + ';'.join(annotations[i][1]) + '\t' + ';'.join(annotations[i][2]) + '\t' + ';'.join(annotations[i][3]) + '\n')

outfile.close()