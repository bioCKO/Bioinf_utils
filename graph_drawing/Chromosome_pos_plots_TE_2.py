__author__ = 'mjohnpayne'

from Bio import SeqIO
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
import sys
from time import sleep as sl

ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')

#inTE = open('/Users/mjohnpayne/Documents/PhD/Chromosome_plots/TE/pmfa1_annot_scaf_repeats.gff','r')

inTE = open('/Volumes/MP_HD/repDNA_data/repeatmasker/Pm_all_genome.fasta.txt','r').read()


print len(sys.argv)

def remove_dupes(ls):
    ls = list(set(ls))
    return ls


# if len(sys.argv) == 3:
#     acclis1 = remove_dupes(open(sys.argv[1],'r').read().split())
# elif len(sys.argv) == 4:
#     acclis1 = remove_dupes(open(sys.argv[1],'r').read().split())
#     acclis2 = remove_dupes(open(sys.argv[2],'r').read().split())
# elif len(sys.argv) == 5:
#     acclis1 = remove_dupes(open(sys.argv[1],'r').read().split())
#     acclis2 = remove_dupes(open(sys.argv[2],'r').read().split())
#     acclis3 = remove_dupes(open(sys.argv[3],'r').read().split())
# out = sys.argv[-1]

##Extract gene position, contig and orientation info from gff
# contigs = {}
# gene = {}
# for line in ingff:
#         col = line.strip('\n').split('\t')
#         if 'ID=gene' in line:
#             det = col[8].split(';')
#             pmaa = det[1][5:]
#             st = int(col[3])
#             en = int(col[4])
#             orient = col[6]
#             cont = col[0]
#             if pmaa not in gene:
#                 if cont in contigs:
#                         contigs[cont] += [pmaa]
#                 elif cont not in contigs:
#                         contigs[cont] = [pmaa]
#                 gene[pmaa] = [cont,st,en,orient]
#
# for cont in contigs:
#         contigs[cont] = sorted(contigs[cont])

#tetype = 'TTCGGG'

tetype = "Copia"

def any_of_list(string,list):
    out = False
    for i in list:
        if i in string:
            out = True
    return out

#te_types = {'dna_te':["hobo","Tc1","IS630","Pogo","PiggyBac","Tourist","Harbinger"],'retro_te' = ["Penelope",]


## make all repmasker annots gene features

##also make lists of "gene" ids for different TE categories
dnalist = []
retro_list = []
contigs = {}
gene = {}
for line in inTE.split("\r")[3:]:
    col = line.strip('\n').split('\t')
    pmaa = str(col[4]) + '_' + str(col[5])
    st = int(col[5])
    en = int(col[6])
    if col[8] == "C":
        orient = "-"
    else:
        orient = col[8]
    cont = col[4]
    name = col[9]
    if "LTR/" in col[10]or "LINE/" in col[10]:
        retro_list.append(pmaa)
    elif "DNA/" in col[10]:
        dnalist.append(pmaa)
    if pmaa not in gene:
        if cont in contigs:
                contigs[cont] += [pmaa]
        elif cont not in contigs:
                contigs[cont] = [pmaa]
        gene[pmaa] = [cont,st,en,orient,name]


# 100 {'100_2': ['2110000', '2140000']}
# 101 {'101_1': ['0', '10000']}
# 102 {'102_3': ['2610000', '2640000']}
# 103 {'103_6': ['1520000', '1550000']}
# 61 {'61_5': ['1', '4275']}
# 93 {'93_7': ['1350000', '1370000']}
# 95 {'95_8': ['760000', '790000']}
# 96 {'96_10': ['3520000', '3550000']}
# 97 {'97_9': ['2050000', '2070000']}
# 99 {'99_4': ['2850000', '2880000']}



centromeres = {"c-100":["100.1",2110000, 2140000,None,"c100"],'c-102': ["102.1",2610000, 2640000,None,"c102"], 'c-103': ["103.1",1520000, 1550000,None,"c103"],'c-93': ["93.1",1350000, 1370000,None,"c93"],'c-95': ["95.1",760000, 790000,None,"c95"],'c-96': ["96.1",3520000, 3550000,None,"c96"],'c-99': ['99.1',2850000, 2880000,None,"c99"],'c-97': ['97.1',2050000, 2070000,None,"c97"]}

cent_list = ["c-100","c-99","c-102","c-103","c-93","c-95","c-96","c-97"]

for cont in contigs:
        contigs[cont] = sorted(contigs[cont])

## Generate list of contig,length tuples then sort in decending order
genome = {}
entries = []
lengths = []
for i in SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome.fasta", "fasta"):
    genome[i.id] = i
    entries.append((i.id,len(i.seq)))
    lengths.append(len(i.seq))

entries = sorted(entries,key=lambda x: x[1],reverse=True)
## Create gene SeqFeatures as subsets of contig seqRecord objects
for i in genome:
    for j in gene:
        if gene[j][0] == genome[i].id:
            # print j,gene[j][0],gene[j][3],genome[i].id
            # sl(0.5)
            direc = int(gene[j][3] + '1')
            genome[i].features.append(SeqFeature(FeatureLocation(gene[j][1], gene[j][2], strand = direc),type='gene',id=j,qualifiers={'locus_tag':[gene[j][4]]}))
    for j in centromeres:
        if centromeres[j][0] == genome[i].id:
            # print j,gene[j][0],gene[j][3],genome[i].id
            # sl(0.5)
            direc = None#int(gene[j][3] + '1')
            genome[i].features.append(SeqFeature(FeatureLocation(centromeres[j][1], centromeres[j][2], strand = direc),type='gene',id=j,qualifiers={'locus_tag':[centromeres[j][4]]}))

## telomere length - rounded ends of chromosome size
max_len = max(lengths)
telomere_length = 40000

chr_diagram = BasicChromosome.Organism()
chr_diagram.page_size = (60*cm, 21*cm)


for index, (name, length) in enumerate(entries):
    if length > 80000:
        features = []
        for i in cent_list:
            for f in genome[name].features:
                if f.id==i:
                    f.qualifiers['color'] = ["grey"]
                    features += [f]
        for i in dnalist:
            for f in genome[name].features:
                if f.id==i:
                    f.qualifiers['color'] = ["red"]
                    features += [f]
        for i in retro_list:
            for f in genome[name].features:
                if f.id==i:
                    f.qualifiers['color'] = ["blue"]
                    features += [f]


        #for f in features: f.qualifiers["color"] = [index+2]
        cur_chromosome = BasicChromosome.Chromosome(name)
        cur_chromosome.scale_num = max_len + 2 * telomere_length
        cur_chromosome.label_size = 0
        cur_chromosome.chr_percent = 0
        # cur_chromosome.label_sep_percent = 0



        #Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.scale = telomere_length
        cur_chromosome.add(start)

        #Add a body - using bp as the scale length here.
        body = BasicChromosome.AnnotatedChromosomeSegment(length,features)
        body.scale = length
        cur_chromosome.add(body)

        #Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = telomere_length
        cur_chromosome.add(end)

        #This chromosome is done
        chr_diagram.add(cur_chromosome)

scale = BasicChromosome.Chromosome('Legend')
scale.scale_num = max_len + 2 * telomere_length
scale.label_size = 8
scale.chr_percent = 0.1
scale.label_sep_percent = 0.12
#
#
acc1id = "DNA_transposons"
feats = [SeqFeature(FeatureLocation(100000, 100500, strand = 1),type='gene',id=acc1id,qualifiers={'locus_tag':[acc1id],'color':["red"]})]
#
acc2id = "Retroelements"
feats += [SeqFeature(FeatureLocation(400000, 400500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':["blue"]})]

acc2id = "Putative Centromeres"
feats += [SeqFeature(FeatureLocation(700000, 700500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':["grey"]})]
#
#Add a body - using bp as the scale length here.
body = BasicChromosome.AnnotatedChromosomeSegment(1000000,feats)
body.scale = 1000000
scale.add(body)


# This chromosome is done
chr_diagram.add(scale)
#chr_diagram.demo()

chr_diagram.draw('/Users/mjohnpayne/Documents/PhD/Chromosome_plots/test/2_types_TE_chromosome_positions.pdf', 'DNA and retroelement transposon positions')
