__author__ = 'mjohnpayne'


# generate list of all genes and dictionary of which genes are in which orthogroup - done
# generate dict of genes with GO IDS from GAFs (may need accession tweaking) - done
# merge gene dictionaries - done
# Use orthogroup gene list to generate orthogroup GO ids




import sys, traceback
from time import sleep as sl
import glob


clusterinf = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups.txt','r')

ingaf = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/GAF_files/*.gaf")


intsv = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/GAF_files/*.tsv")

aspgd = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/GAF_files/aspgd_annot.gf"

outfile = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/GAF_files/eurot_gafs_aggregated_GO_terms.gaf",'w')
out_simple_assoc = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/GAF_files/eurot_gafs_aggregated_GO_terms.assoc",'w')


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
        good = ["ANID","Pma","CIM","Tst","Tfl",'Pfu',"Hca","Pbr","Pch","Afu","Nfi","Acl","Ate","Afl","Aor","Ani"]
        asper = ["ANID","Afu","Nfi","Acl","Ate","Afl","Aor","Ani"]
        genes = [rename(x[x.find("|")+1:]) for x in col[1:] if x[:x.find("|")] in good]
        clust[col[0][:-1]] = genes
        aspgenes = [rename(x[x.find("|")+1:]) for x in col[1:] if x[:x.find("|")] in asper]
        genlst += aspgenes
    return clust,genlst

def correct_taxa(string):
    col = string.split('\t')
    taxals = ["162425",'746128','5061','5062','341663','331117','344612','332952']
    taxa = col[12][6:]
    if taxa in taxals:
        return (True,taxa)
    else:
        return (False,"")

def make_go_dict_frm_aspgd(infiles,genes,gos):
    file = open(infiles,'r').readlines()
    for line in file:
        if line[0] != "!" and correct_taxa(line)[0]:
            col = line.split('\t')
            gene = ''
            inf = (col[4],col[8])
            if correct_taxa(line)[1] == "162425":
                if len(col[10].split("|")) > 1:
                    gene = col[10].split("|")[1]
                else:
                    pass
            elif correct_taxa(line)[1] == "746128":
                gene = col[10].split('|')[0]
            elif correct_taxa(line)[1] == "5061":
                if len(col[10].split("|")) > 1:
                    gene = col[10].split("|")[1]
                else:
                    pass
            elif correct_taxa(line)[1] == "5062":
                gene = col[10].split('|')[0]
            elif correct_taxa(line)[1] == "341663":
                gene = col[2]
            elif correct_taxa(line)[1] == "331117":
                gene = col[2]
            elif correct_taxa(line)[1] == "344612":
                gene = col[2]
            elif correct_taxa(line)[1] == "332952":
                gene = col[2]
            if gene not in gos:
                gos[gene] = [inf]
            else:
                if inf not in gos[gene]:
                    gos[gene].append(inf)
    return gos


def parse_gafs(gafls):
    gos = {}
    for i in gafls:
        tmp = open(i,'r')
        for j in tmp:
            if j[0] != "!":
                col = j.split('\t')
                gene = col[2]
                inf = (col[4],col[8])
                if gene not in gos:
                    gos[gene] = [inf]
                else:
                    if inf not in gos[gene]:
                        gos[gene].append(inf)
        tmp.close()
    return gos

def fixname(string):
    if string[:3] == "PCH":
        string = string[4:]
    return string


def parse_tsvs(tsvls,gos):
    for i in tsvls:
        tmp = open(i,'r')
        for j in tmp:
            if j[0] != "D":
                col = j.split('\t')
                gene = fixname(col[-1].split('|')[-1].strip('\n'))
                inf = (col[4],col[6][0])
                if gene not in gos:
                    gos[gene] = [inf]
                else:
                    if inf not in gos[gene]:
                        gos[gene].append(inf)
        tmp.close()
    return gos

def ortho_go_dict(cluster,godict):
    clusterdict = {}
    for i in cluster:
        for j in cluster[i]:
            if j in godict:
                if i not in clusterdict:
                    clusterdict[i] = godict[j]
                else:
                    for x in godict[j]:
                        if x not in clusterdict[i]:
                            clusterdict[i].append(x)
    return clusterdict




def main(cl,int,asp,tsv):
    clu,glist = ortho_clust_to_dict(cl,[])
    go = parse_gafs(int)
    #go = parse_tsvs(tsv,go)
    # go = {}
    # go = make_go_dict_frm_aspgd(asp,glist,go)
    cludict = ortho_go_dict(clu,go)
    return cludict

annotations = main(clusterinf,ingaf, aspgd,intsv)

outfile.write("!Generated by aggregate_go_terms_for_eurots.py\n!Date: 20150416\n!gaf-version: 2.0\n")

for i in annotations:
    out_simple_assoc.write(i + '\t' + ';'.join([x[0] for x in annotations[i]]) + '\n')
    for j in annotations[i]:
        outfile.write('MJP\t' + i + '\t' + i + '\t\t' + j[0] + '\t' + i + '\tISS\t\t' + j[1] + '\t\t\tprotein\tTalaromyces marneffei\t20150416\tmpayne\t\t\n')
out_simple_assoc.close()
outfile.close()