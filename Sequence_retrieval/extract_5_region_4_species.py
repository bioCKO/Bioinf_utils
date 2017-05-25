__author__ = 'mjohnpayne'


import re
from Bio import SeqIO
from Bio.Seq import Seq
import sys
from time import sleep as sl
#import regex

limit = 5000#int(raw_input('5 prime limit: '))
inacc = open('/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids_f.txt','r')#open(raw_input('accession number list path: '),'r')
outfile = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/4_spec_conserved_orthos_5000bp_proms/'#raw_input('output file path: ')
fastas = {'Pm':'/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta','Ts':'/Users/mjohnpayne/Documents/PhD/wt_genome/tsta1_annot_scaf.fasta','Tf':'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/T_flavus/vel_denovo_genome/TF_vel_genome.fasta','Pf':'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Pf_databases/PF_scaff_genome.fasta'}
gffs = {'Pm':'/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff','Ts':'/Users/mjohnpayne/Documents/PhD/wt_genome/ts_wt_dbs/tsta1_working_models.gff3',"Tf":"/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/T_flavus/vel_denovo_genome/TF_vel_pfams_para_genome_fix_rename.gff","Pf":"/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_denovo_rename.gff"}

orthogroups = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups.txt','r')


#outlong = open('/Users/mjohnpayne/Documents/PhD/bys/bys_promoter_analysis/bys_motif_counts','w')
#out5primes = open('/Users/mjohnpayne/Documents/PhD/wt_genome/2161_5prime_intergenics.fa','w')

## Make dictionary of contig numbers and sequences, contig number key
def revcomp_re(str):
    new = []
    for i in str:
        if i == 'T':
            new = ['A'] + new
        elif i == 'A':
            new = ['T'] + new
        elif i == 'G':
            new = ['C'] + new
        elif i == 'C':
            new = ['G'] + new
        elif i == '[':
            new = [']'] + new
        elif i == ']':
            new = ['['] + new
    return ''.join(new)


def get_genome(fasta_str):
    gnm = {}
    for record in SeqIO.parse(fasta_str, 'fasta'):
        if 'TF_vel_genome' in fasta_str:
            id = str(record.id)[:str(record.id).find('|')]
            gnm['Tf' + id] = record.seq
        else:
            gnm[record.id] = record.seq
    return gnm

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively


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
    return conts,gen

def getproms(cont,genme,lim,genes):
    proms = {}
    for conti in cont:
            for i in range(0,len(cont[conti])):
                    g = cont[conti][i]
                    wins = 0
                    wine = 0
                    winlen = 0
                    if genes[g][2] == '+':
                            try:
                                    if genes[g][0] < genes[cont[conti][i-1]][1]:
                                            wins = genes[cont[conti][i-2]][1]
                                            wine = genes[g][0]
                                            winlen = wine - wins
                                    else:
                                            wins = genes[cont[conti][i-1]][1]
                                            wine = genes[g][0]
                                            winlen = wine - wins
                            except:
                                    wins = genes[g][0]-lim
                                    wine = genes[g][0]
                                    if wins < 0:
                                            wins = 0
                                    winlen = wine - wins
                            if winlen > lim:
                                    wins = wine-lim
                            proms[g] = str(Seq(str(genme[conti][wins:wine])))
                    elif genes[g][2] == '-':
                            try:
                                    if genes[g][1] > genes[cont[conti][i+1]][0]:
                                            wins = genes[g][1]
                                            wine = genes[cont[conti][i+2]][0]
                                            winlen = wine - wins
                                    else:
                                            wins = genes[g][1]
                                            wine = genes[cont[conti][i+1]][0]
                                            winlen = wine - wins
                            except:
                                    wins = genes[g][1]
                                    wine = genes[g][1]+lim
                                    winlen = wine - wins
                            if winlen > lim:
                                    wine = wins+lim
                            #print conti, wins, wine, g
                            proms[g] = str(Seq(str(genme[conti][wins:wine])).reverse_complement())
                    if winlen < 0:
                            continue
    return proms

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

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

def count_ortho_nos(ln):
    if ln.count('Pma|') == 1 and ln.count('Tst|') == 1 and ln.count('Tfl|') == 1 and ln.count('Pfu|') == 1:
        return True
    else:
        return False


def get_1to1_orthoids(og):
    cons_orthos = {}
    for line in og:
        id = ''
        genes = []
        if count_ortho_nos(line):
            col = line.strip('\n').split(' ')
            id = col[0].replace(':','')
            genes = col[1:]
            ngenes = []
            for i in genes:
                i = rn(i)
                ngenes.append(i)
            cons_orthos[id] = ngenes
    print len(cons_orthos)
    return cons_orthos






def main(fasls,gfls,lim,of):
    proms = {}
    genome = {}
    contigs = {}
    gene = {}
    for i in fasls:
        print i
        genome = get_genome(fasls[i])
        #genome = genome.update(get_genome(fasls[i]))
        c,g = get_gene_conigs(gfls[i])
        #contigs = contigs.update(c)
        #gene = gene.update(c)
        proms = merge_two_dicts(proms,getproms(c,genome,lim,g))
    print len(proms)
    groups = get_1to1_orthoids(orthogroups)
    c = 0
    for i in groups:
        inacs = groups[i]
        out = open(of+'conserved_orthos_prom_' + i,'w')
        for j in inacs:
            out.write('>' + j + '\n' + proms[j] + '\n')
        out.close()
        c += 1
        if c%100 == 0:
            print c,'Done'


main(fastas,gffs,limit,outfile)