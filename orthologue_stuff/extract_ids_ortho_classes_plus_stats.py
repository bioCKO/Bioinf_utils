__author__ = 'mjohnpayne'

from numpy import average as av
from numpy import std
from numpy import sqrt
from time import sleep as sl

orthogroups = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups.txt'

of = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/cluster_class_data.txt','w')

gffs = {'Pm':'/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pmfa1_working_models_fix.gff','Ts':'/Users/mjohnpayne/Documents/PhD/wt_genome/ts_wt_dbs/tsta1_working_models.gff3',"Tf":"/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/T_flavus/vel_denovo_genome/TF_vel_pfams_para_genome_fix_rename.gff","Pf":"/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_denovo_rename.gff"}

group_cats = ['ab_conserved','conserved','PM_only','TF_only','PF_only','TS_only','other']


def rn(gene):
    if 'TF' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif 'Pf' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    else:
        gene = gene[4:]
    return gene


def get_orthoids(og,type):
    orthos = []
    allgenes = []
    og = open(og,'r')
    for line in og:
        id = ''
        genes = []

        col = line.strip('\n').split(' ')
        id = col[0].replace(':','')
        genes = col[1:]
        ngenes = []
        for i in genes:
            i = rn(i)
            ngenes.append(i)
        if count_ortho_nos(line) == type:
            orthos += ngenes
        allgenes += ngenes
    og.close()
    return orthos, allgenes


def count_ortho_nos(ln):
    if ln.count('Pma|') == 1 and ln.count('Tst|') == 1 and ln.count('Tfl|') == 1 and ln.count('Pfu|') == 1:
        return 'ab_conserved'
    elif ln.count('Pma|') >= 1 and ln.count('Tst|') >= 1 and ln.count('Tfl|') >= 1 and ln.count('Pfu|') >= 1:
        return 'conserved'
    elif ln.count('Pma|') >= 1 and ln.count('Tst|') == 0 and ln.count('Tfl|') == 0 and ln.count('Pfu|') == 0:
        return 'PM_only'
    elif ln.count('Pma|') == 0 and ln.count('Tst|') >= 1 and ln.count('Tfl|') == 0 and ln.count('Pfu|') == 0:
        return 'TS_only'
    elif ln.count('Pma|') == 0 and ln.count('Tst|') == 0 and ln.count('Tfl|') >= 1 and ln.count('Pfu|') == 0:
        return 'TF_only'
    elif ln.count('Pma|') == 0 and ln.count('Tst|') == 0 and ln.count('Tfl|') == 0 and ln.count('Pfu|') >= 1:
        return 'PF_only'
    else:
        return 'other'


def get_gene_conigs(gfin):
    gf = open(gfin,'r')
    conts = {}
    gen = {}
    for line in gf:
        if not line.startswith('#'):
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                det = col[8].split(';')
                pmaa=''
                if col[8].startswith('ID=gene'):
                    pmaa = det[1][5:]
                else:
                    pmaa = det[0][3:]
                st = int(col[3])
                en = int(col[4])
                orient = col[6]
                cont = ''
                if 'TF_vel' in gfin:
                    cont = 'Tf' + col[0]
                else:
                    cont = col[0]
                if cont in conts:
                        conts[cont] += [pmaa]
                elif cont not in conts:
                        conts[cont] = [pmaa]
                gen[pmaa] = [st,en,orient]
    for cont in conts:
        conts[cont] = sorted(conts[cont])
    gf.close()
    gf = open(gfin,'r').read().split('\tgene\t')
    for i in gf:
        if 'CDS' in i:
            cdss = i.count('\tCDS\t')
            pmaa = ''
            i = i.replace('\n','\t')
            col = i.split('\t')
            det = col[5].split(';')
            if col[5].startswith('ID=gene'):
                pmaa = det[1][5:]
            else:
                pmaa = det[0][3:]
            gen[pmaa].append(cdss)

    return conts,gen


def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

def get_stats(gendict,genlis,spec,type,algen):
    genlen = []
    cdsno = []
    for i in genlis:
        if spec in i:
            if i in gendict:
                gene = gendict[i]
                genlen.append(gene[1]-gene[0])
                cdsno.append(gene[3])
    if 'only' in type:
        for i in gendict:
            if spec in i and type[:2] in i:
                if i not in algen:
                    # if spec == 'PM':
                    #     print i
                    gene = gendict[i]
                    genlen.append(gene[1]-gene[0])
                    cdsno.append(gene[3])
    if len(genlen) > 0:
        return len(genlen), av(genlen), std(genlen), std(genlen)/sqrt(len(genlen)), len(cdsno), av(cdsno), std(cdsno), std(cdsno)/sqrt(len(cdsno))
    else:
        return 'N/A', 'N/A', 'N/A', 'N/A','N/A', 'N/A', 'N/A', 'N/A'


def stats(og,gfls,type,spec):
    contigs = {}
    genes = {}

    for gff in gfls:
        c2,g2 = get_gene_conigs(gffs[gff])
        contigs = merge_two_dicts(contigs,c2)
        genes = merge_two_dicts(genes,g2)

    geneids,allgen = get_orthoids(og,type)
    l,a,sd,s,cl,ca,csd,cs = get_stats(genes,geneids,spec,type,allgen)
    return l,a,sd,s,cl,ca,csd,cs

species = ['PM','TF','TS','PF']

of.write('Group type\tspecies\tgene number\taverage length\tstdev\tsem\texon number average\texon number stdev\texon number sem\n')

for i in group_cats:
    for j in species:
        length,avg,sd,sem,clen,cavg,csd,csem = stats(orthogroups,gffs,i,j)
        if length == 'N/A':
            continue
        else:
            #print i +'_'+ j + '\t' + str(length) + '\t' + str(avg) + '\t' + str(sd) + '\t' + str(sem)
            of.write(i +'_'+ j + '\t' + str(length) + '\t' + str(avg) + '\t' + str(sd) + '\t' + str(sem) + '\t' + str(cavg) + '\t' + str(csd) + '\t' + str(csem) +'\n')

of.close()