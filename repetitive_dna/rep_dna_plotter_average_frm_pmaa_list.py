__author__ = 'mjohnpayne'
from matplotlib import *
from pylab import *
from numpy import *
from Bio import SeqIO
from Bio import Seq
import subprocess
from time import sleep as sl
import matplotlib.pyplot as plt
import math
import matplotlib.gridspec as gridspec
import scipy

flank = 10000

ingenome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome.fasta","fasta")

#outpos = open("/Users/mjohnpayne/Documents/PhD/wt_genome/repDNA_data/Pm_all_genome_rep_counts30mer_5step.txt",'w')


#ingene = open("/Users/mjohnpayne/Documents/PhD/ASPS/pop_ids.txt",'r')
#ingene = open("/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids_f.txt",'r')
#ingene = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/galactomannoproteins/gmps_pm_f.txt",'r')
ingene = open("/Volumes/MP_HD/repDNA_data/pops/pop_ids.txt")


genome = {}

for i in ingenome:
    genome[i.id[:-2]] = i.seq


incount = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_100mer_counts.txt",'r')

count_dict = {}

for line in incount:
    col = line.strip('\n').split('\t')
    count_dict[col[0]] = int(col[1])
print len(count_dict)


def n_random_genes_from_gff(n):
    gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')
    genes = []
    for line in gff:
        if 'ID=gene' in line:
            col = line.strip('\n').split('\t')
            det = col[8].split(';')
            pmaa = det[1][5:]
            genes.append(pmaa)
    rand_smpl = []
    gff.close()
    cnt = 0
    while cnt < n:
        random.seed(cnt+295)
        rand_smpl.append(genes[random.randint(0,len(genes)-1)])
        cnt += 1
    return rand_smpl

def extract_genome_counts_across_gene(gene,genom,prime5,prime3):
    gff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff",'r')
    seq = ""
    for line in gff:
        if 'ID=gene' in line and gene in line:
            col = line.strip('\n').split('\t')
            det = col[8].split(';')
            pmaa = det[1][5:]
            if int(col[3]) > prime5:
                p5 = len(str(genom[col[0]][int(col[3])-prime5:int(col[3])]))
            else:
                p5 = len(str(genom[col[0]][0:int(col[3])]))
            if (int(col[4]) + prime3) > len(genom[col[0]]):
                p3 = len(str(genom[col[0]][int(col[4]):len(genom[col[0]])]))
            else:
                p3 = len(str(genom[col[0]][int(col[4]):int(col[4])+prime3]))
            if int(col[3]) < prime5:
                seq = str(genom[col[0]][0:int(col[4])+prime3])
            elif (int(col[4]) + prime3) > len(genom[col[0]]):
                seq = str(genom[col[0]][int(col[3])-prime5:len(genom[col[0]])])
            else:
                seq = str(genom[col[0]][int(col[3])-prime5:int(col[4])+prime3])
            glen = len(str(genom[col[0]][int(col[3]):int(col[4])]))
    counts = []
    for i in range(len(seq)-100):
        win = seq[i:i+100]
        if "N" in win or "n" in win:
            counts.append(0)
        else:
            counts.append(count_dict[win])
    gff.close()
    return counts,p5,p3,glen

### Add counts to overall list if repetitive (base of 1) for each position where count is > 1, add 1

def merge_counts_lists(flank5,flank3,countls):
    eq_counts = {}
    for i in flank5:
        eq_len_counts = (flank-flank5[i])*[0] + countls[i][:flank5[i]] + 1000*[-1] + countls[i][-flank3[i]:] + (flank-flank3[i])*[0]
        eq_counts[i] = eq_len_counts
        print len(eq_len_counts)
    sumcounts = []
    for i in range(0,21000-1):
        sumval = average([eq_counts[j][i]*100 for j in sorted(eq_counts.keys())])
        sumcounts.append(sumval)
    return sumcounts


gcounts = {}
g5primes = {}
g3primes = {}
glens = {}

#ingene = n_random_genes_from_gff(50)

for i in ingene:
    gene = i.strip('\n')
    gcounts[gene],g5primes[gene],g3primes[gene],glens[gene] = extract_genome_counts_across_gene(gene,genome,flank,flank)
    ##### converts repetitive counts into 1 no reps and 2 reps #####
    nlist = []
    for j in gcounts[gene]:
        if j > 1:
            nlist.append(1)
        else:
            nlist.append(0)
    gcounts[gene] = nlist
    #####

sumls = merge_counts_lists(g5primes,g3primes,gcounts)




fig=plt.figure()#figsize=(8.27,11.7))

ax1 = plt.axes()
ax0 = plt.axes()

genpos = [0]*10000 + [100]*1000 + [0]*10000
ax1.plot(genpos, lw=0.1, antialiased = True, color = 'lightcoral')

ys = sumls

ax0.plot(ys, lw=0.1, antialiased = True,color = 'steelblue')

d = scipy.zeros(len(genpos))
ax1.fill_between([j for j in range(len(genpos))],genpos,where=genpos>=d,interpolate=True,color = 'lightcoral')

d = scipy.zeros(len(ys))
ax0.fill_between([j for j in range(len(ys))],ys,where=ys>=d,interpolate=True,color='steelblue')

ylim(0,30)
xlim(0,21000)
# xticks([len(genpos)/2,len(genpos)],fontsize=8)
# yticks([10,20],fontsize=8)

title("% of pop genes repetitive at position",fontsize=12)

#plt.show()
plt.savefig("/Volumes/MP_HD/repDNA_data/pops/pop_10kb_flanking_avg_repetitiveness.pdf",dpi=200)