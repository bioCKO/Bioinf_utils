__author__ = 'mjohnpayne'


## counts of genes with x characteristic that are gained/lost in tm branch

## characteristics of lost gene are determined by majority type for ortho_group

## make dict of gene: interpros

## for each ortho group add all intpros in constituent genes if % greater than 40 - add lists of interpros, do counts accept if above 40% threshold

## go through eurot groups - assign gains and losses to interpros progressively i.e. as interpro is in a group add gains losses to dict for that interpro

## go through interpros and print



import glob
from time import sleep as sl
from collections import Counter
import seaborn as sns
from matplotlib import pyplot as plt


def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]


annot_lists = glob.glob("/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_functional_annot/*_descriptions.txt")

bdrate_changes = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_spec_gene_gain_loss/eurot_gene_gain_loss/eurot_gene_gain_loss_from_badirate.txt","r").readlines()
num_change = {}

outfile = open("/Users/mjohnpayne/Documents/PhD/wt_genome/4_species_functional_annot/orthogroup_interpros_gain_loss_functional.txt","w")
outfile.write("Interpro Accession\tGains\tLosses\tTotal\n")

for i in bdrate_changes[1:]:
    col = i.split('\t')
    num_change[col[0]] = int(col[3])


iprdict = {}
genenos = {}
for i in annot_lists:
    spec = i.split('/')[-1][:-17]
    print spec
    info = open(i,"r").read().split('\r')
    c = 1
    for line in info[1:]:
        col = line.split('\t')
        interpro = col[2].split(",")
        interpros = []
        for i in interpro:
            if i.startswith("IPR"):
                interpros.append(i[:9])
        if col[13] == "-":
            # continue
            iprdict[spec+"_singleton_" + str(c)] = interpros
            genenos[spec+"_singleton_" + str(c)] = 1
            c +=1
        else:
            if col[13] not in iprdict:
                iprdict[col[13]] = interpros
                genenos[col[13]] = 1
            else:
                iprdict[col[13]]+=interpros
                genenos[col[13]] += 1


annots = {}
for i in iprdict:
    annots[i] = []
    for j in set(iprdict[i]):
        count = iprdict[i].count(j)
        frac = float(count)/genenos[i]
        if frac >= 0.4:
            annots[i].append(j)

ipr_changes = {}

for i in annots:
    group = annots[i]
    for ipr in group:
        change = 0
        if "T_marneffei_singleton" in i:
            change = 1
        elif "singleton" in i:
            change = 0
        else:
            change = num_change[i]
        if ipr not in ipr_changes:
            ipr_changes[ipr] = [0,0]
            if change > 0:
                ipr_changes[ipr][0] += change
            elif change < 0:
                ipr_changes[ipr][1] += change
        else:
            if change > 0:
                ipr_changes[ipr][0] += change
            elif change < 0:
                ipr_changes[ipr][1] += change

for i in ipr_changes:
    if ipr_changes[i][0] != 0 or ipr_changes[i][1] != 0:
        g = ipr_changes[i][0]
        l = ipr_changes[i][1]
        t=g+l
        outfile.write(i+'\t'+str(g)+'\t'+str(l)+'\t'+str(t)+'\n')
outfile.close()



