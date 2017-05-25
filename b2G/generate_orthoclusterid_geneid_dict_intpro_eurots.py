__author__ = 'mjohnpayne'

# 1 - generate dict of clusters with clusterid as key and gene ids as data
# 2 - generate dict of interpro ID:description objects with gene id as key and id-descriptions as data
# 2 - generate dict of GOTermobjects with gene id as key and id-descriptions as data
# 3 - merge by having clusterid link to collapsed list of interpro ids (set)
# 4 - print out

import sys
from time import sleep as sl
import glob

clusterinf = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups.txt','r')

descriptions = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/*.tsv")  #['/Volumes/MP_HD/Blast2Go/Pm/exported_data/Pm_b2g_descriptions.txt','/Volumes/MP_HD/Blast2Go/Tf/exported data/Tf_b2g_descriptions.txt','/Volumes/MP_HD/Blast2Go/Pf/exported_data/Pf_b2g_descriptions.txt','/Volumes/MP_HD/Blast2Go/Ts/exported data/Ts_b2g_descriptions.txt']


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
        clust[col[0][:-1]] = [rename(x[x.find("|")+1:]) for x in col[1:]]
    return clust



def desc_to_dict(intproin):
    ipro = {}
    for f in intproin:
        tmp = open(f,'r').readlines()
        for line in tmp[1:]:
            if "\tIPR" in line:
                col = line.strip('\r\n').split('\t')
                gene = rename(col[0][col[0].find("|")+1:])
                ipr = ""
                desc = ""
                for i in range(len(col)):
                    if col[i][:3] == "IPR":
                        ipr = col[i]
                        desc = col[i+1]
                if gene not in ipro:
                    ipro[gene] = [ipr+":"+desc]
                else:
                    if ipr+":"+desc not in ipro[gene]:
                        ipro[gene].append(ipr+":"+desc)
    return ipro

def merge(orthos,ip):
    annots = {}
    for i in orthos:
        ipd = []
        for j in orthos[i]:
            if j in ip:
                ipd+=ip[j]
        ipl = list(set(ipd))
        if len(ipl) == 0:
            annots[i] = ["None:No Description"]
        else:
            annots[i] = ipl
    return annots



def main(cl,int):
    clu = ortho_clust_to_dict(cl)
    # for i in clu:
    #     print i,clu[i]
    #     sl(0.5)
    ipro = desc_to_dict(int)
    # for i in ipro:
    #     print i,ipro[i]
    #     sl(0.5)
    annots = merge(clu,ipro)
    # for i in annots:
    #     print i,annots[i]
    #     sl(0.5)
    return annots

annotations = main(clusterinf,descriptions)

# for i in annotations:
#     print i,annotations[i]
#     sl(0.5)
outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_ortho_group_intpros.txt','w')

outfile.write('ID\tInterpro IDs and descriptions\n')

for i in annotations:
    outfile.write(i + '\t' + ' , '.join(annotations[i]) + '\n')

outfile.close()