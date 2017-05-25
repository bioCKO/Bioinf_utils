__author__ = 'mjohnpayne'

from time import sleep as sl

clusterinf = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups.txt','r')

aftab = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/gene description files/A_fumigatus_Af293_current_chromosomal_feature.tab",'r')
antab = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/gene description files/A_nidulans_FGSC_A4_current_chromosomal_feature.tab",'r')


pmgff = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/gene description files/pmfa1_working_models_fix.gff",'r')
tsgff = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/gene description files/tsta1_working_models.gff3",'r')

t_gffs=[pmgff,tsgff]

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups_aggregated_descriptions.txt','w')

def rename(gene):
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

def ortho_clust_to_dict(inf,genlst):
    clust = {}
    for line in inf:
        col = line.strip('\n').split(' ')
#        good = ["ANID","Pma","CIM","Tst","Tfl",'Pfu',"Hca","Pbr","Pch","Afu","Nfi","Acl","Ate","Afl","Aor","Ani"]
#        asper = ["ANID","Afu","Nfi","Acl","Ate","Afl","Aor","Ani"]
        tested = ["ANID","Afu","Pma","Tst"]
        genes = [rename(x[x.find("|")+1:]) for x in col[1:] if x[:x.find("|")] in tested]
        clust[col[0][:-1]] = genes
#        aspgenes = [rename(x[x.find("|")+1:]) for x in col[1:] if x[:x.find("|")] in asper]
        genlst += genes
    return clust,genlst

def gff_to_desc_dict(gfflst):
    genes = {}
    count = 0
    spec = ["Pm","Ts"]
    for i in gfflst:
        for line in i:
            if '\tgene\t' in line:
                desc = line.split('\t')[-1].strip('\n').split(';')
                id = desc[1][5:]
                inf = desc[2][6:].replace('\r','')
                genes[id] = [spec[count],"",inf]
        count +=1
    return genes

def tab_to_desc_dict(af,an,genes):
    for line in an:
        if 'ANID_' in line:
            col = line.split('\t')
            nid = col[2].split('|')[-1]
            name = col[1]
            desc = col[10]
            genes[nid] = ["An",name,desc]
    for line in af:
        if "!" not in line:
            col = line.split('\t')
            nid = col[0]
            name = col[1]
            desc = col[10]
            genes[nid] = ["Af",name,desc]
    return genes

def ortho_inf_dict(cluster,inf_dict):
    clusterdict = {}
    for i in cluster:
        for j in cluster[i]:
            if j in inf_dict:
                if i not in clusterdict:
                    clusterdict[i] = [inf_dict[j]]
                else:
                    clusterdict[i].append(inf_dict[j])
    return clusterdict


def main(aft,ant,gfflst,clust):
    clustdict,glist = ortho_clust_to_dict(clust,[])
    infdict = gff_to_desc_dict(gfflst)
    infdict = tab_to_desc_dict(aft,ant,infdict)
    cluster_inf = ortho_inf_dict(clustdict,infdict)
    # for i in cluster_inf:
    #     print i,cluster_inf[i],"\n\n"
    #     sl(0.2)
    return cluster_inf

out_dict = main(aftab,antab,t_gffs,clusterinf)

for i in out_dict:
    outfile.write(i + '\t' + '\t'.join([' '.join(x) for x in out_dict[i]])+'\n')

outfile.close()