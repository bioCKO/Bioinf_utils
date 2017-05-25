__author__ = 'mjohnpayne'

from time import sleep as sl
from Bio import SeqIO

orthos_by_gene = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/eurot_groups_by_gene.txt','r')

all_talaro_cds = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/all_talaro_cds.fasta"

outfolder = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_and_all_talaro_orthos/"

allc = SeqIO.parse(all_talaro_cds,"fasta")
all_cds_seqs = {}

for i in allc:
    all_cds_seqs[i.id] = i

def rn(gene):
    if 'TF' in gene:
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif 'Pf' in gene:
        no = 5-len(gene[3:])
        gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    else:
        gene = gene[4:]
    return gene



for line in orthos_by_gene:
    col = line.strip('\n').split('\t')
    if "PMAA" in col[0]:
        outgenes = []
        outgenes.append(col[0])
        paras = col[1].split(',')
        if paras[0] != "":
            outgenes += paras
        if len(col) > 2:
            orthos = col[2].split(',')
            for i in orthos:
                if "TSTA" in i:
                    outgenes.append(i)
                elif "TF_" in i:
                    outgenes.append(rn(i))
                elif "Pf_" in i:
                    outgenes.append(rn(i))
        write_genes = [all_cds_seqs[x] for x in outgenes]
        SeqIO.write(write_genes,outfolder + outgenes[0] + '_talaro_orthos.fasta',"fasta")