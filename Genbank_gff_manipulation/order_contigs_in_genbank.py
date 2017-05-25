__author__ = 'mjohnpayne'

from Bio import SeqIO
from Bio.Alphabet import IUPAC
import subprocess
## import order and orientation
in_gb = SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/MAuve_reorder/alignment1/PM_1.gb','genbank')

gb_dict = {}

for rec in in_gb:
    print rec.id
    gb_dict[rec.id] = rec

in_ordering = open('/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/MAuve_reorder/alignment6/Talaromyces marneffei PM1_contigs.tab')

num = 1
def revcomp_rec(seq_rec):
        new = seq_rec.reverse_complement()
        new.seq.alphabet = IUPAC.ambiguous_dna
        new.id = seq_rec.id
        return new

def genorder(infile):
    infile = infile.readlines()
    st = 0
    en = 0
    for i in range(len(infile)):
        if 'Ordered Contigs' in infile[i]:
            st = i+2
        elif 'Contigs with conflicting ordering information' in infile[i]:
            en = i-2
    orderinf = infile[st:en]
    orderlis = []
    orient = {}
    for i in orderinf:
        col = i.strip('\n').split('\t')
        orderlis += [col[1]]
        orient[col[1]] = col[3]
    return orderlis,orient


order,direct = genorder(in_ordering)
#
# print direct

num = 1

for i in order:
    if direct[i] == 'forward':
        new = gb_dict[i + '.1']
        SeqIO.write(new,'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf/temp_gb_'+str(num)+'.gb','genbank')
        num += 1
    elif direct[i] == 'complement':
        new = gb_dict[i + '.1']
        new = revcomp_rec(new)
        SeqIO.write(new,'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf/temp_gb_'+str(num)+'.gb','genbank')
        num += 1

comb_out = open('/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/PM_1_reorder.gb','a')

for i in range(1,num):
    tempin = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf/temp_gb_'+str(i)+'.gb','r').read()
    comb_out.write(tempin)# + '\n//\n')

args = 'rm /Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pf/temp_gb_*'

subprocess.Popen(args, shell=True).wait()

comb_out.close()
