__author__ = 'mjohnpayne'

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from time import sleep as sl
import re


ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')


#ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/tsta1_working_models.gff3','r')

#record = SeqIO.read("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/4_species_mauve/reordering/Pm_genbank_gene_containing_only.gb", "genbank")

#helicases = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/recQ_helicase_ts_gene_ids.txt",'r').read().split('\r')

helicases = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/recQ_helicase_gene_ids.txt",'r').read().split('\r')
# pops = open("/Users/mjohnpayne/Documents/PhD/ASPS/pop_ids.txt",'r').read().split('\n')
# byss = open('/Users/mjohnpayne/Documents/PhD/bys/complete_list_bys_ids_f.txt','r').read().split('\n')

print helicases
#print pops

genome = {}
for i in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
    genome[i.id] = i.seq
size = 40000
ends5 = []
ends3 = []
for i in genome:
    l = len(genome[i])
    if l > 500000:
        ends5 += [(i,0,size)]
        ends3 += [(i,l-size,l)]


contigs = {}
gene = {}
for line in ingff:
        col = line.strip('\n').split('\t')
        if 'ID=gene' in line:
            det = col[8].split(';')
            pmaa = det[1][5:]
            st = int(col[3])
            en = int(col[4])
            orient = col[6]
            cont = col[0]
            if pmaa not in gene:
                if cont in contigs:
                        contigs[cont] += [pmaa]
                elif cont not in contigs:
                        contigs[cont] = [pmaa]
                gene[pmaa] = [cont,st,en,orient]


#pos = ('93',0,40000)




lab = False
lsize = 5
athick = 1


gdd = GenomeDiagram.Diagram("5_ends")
c = 1
for pos in ends5:
    gd_track_for_features = gdd.new_track(c, name=pos[0], greytrack=True)
    gd_feature_set = gd_track_for_features.new_set()
    for i in gene:
        if gene[i][0] == pos[0] and gene[i][2] < pos[2] and i in helicases:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=i, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=i, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        # elif gene[i][0] == pos[0] and gene[i][2] < pos[2] and i in pops:
        #     orient = int(gene[i][3] + '1')
        #     if orient > 0:
        #         feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
        #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=30)
        #     elif orient < 0:
        #         feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
        #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        # elif gene[i][0] == pos[0] and gene[i][2] < pos[2] and i in byss:
        #     orient = int(gene[i][3] + '1')
        #     if orient > 0:
        #         feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
        #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="green", arrowshaft_height=athick, label_size = lsize, label_angle=30)
        #     elif orient < 0:
        #         feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
        #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="green", arrowshaft_height=athick, label_size = lsize, label_angle=150)
        elif gene[i][0] == pos[0] and gene[i][2] < pos[2]:
            orient = int(gene[i][3] + '1')
            if orient > 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
            elif orient < 0:
                feature = SeqFeature(FeatureLocation(gene[i][1], gene[i][2]), strand=orient)
                gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=150)

    if pos[1] == 0:
        feature = SeqFeature(FeatureLocation(0, 200), strand=0)
        gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
    else:
        feature = SeqFeature(FeatureLocation(pos[2]-200, pos[2]), strand=0)
        gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)

    # tel_rep = 'TAAGGG'## or in nidulans 'TAAGGG'
    # telomere_sites = [m.start() for m in re.finditer(tel_rep,str(genome[pos[0]][pos[1]:pos[2]]))]
    # for i in telomere_sites:
    #     feature = SeqFeature(FeatureLocation(i, i+5))
    #     gd_feature_set.add_feature(feature, color="blue", name='telomeric_rep',label=True, label_size = lsize)

    c +=1


gdd.draw(format='linear', pagesize='A4', fragments=1, start=0, end=size,track_size=0.5)

gdd.write("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/5_prime_40kb_labels.pdf", "pdf")



#
#
# gdd = GenomeDiagram.Diagram("3_ends")
# c = 1
# for pos in ends5:
#     gd_track_for_features = gdd.new_track(c, name=pos[0], greytrack=True)
#     gd_feature_set = gd_track_for_features.new_set()
#     pst = 0
#     pen = size
#     for i in gene:
#         contig = gene[i][0]
#         st = gene[i][1]-len(genome[contig]) + size
#         en = gene[i][2]-len(genome[contig]) + size
#         if contig == pos[0] and st > 0 and i in helicases:
#             orient = int(gene[i][3] + '1')
#             if orient > 0:
#                 feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#                 gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=30)
#             elif orient < 0:
#                 feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#                 gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="teal", arrowshaft_height=athick, label_size = lsize, label_angle=150)
#         # elif contig == pos[0] and st > 0 and i in pops:
#         #     orient = int(gene[i][3] + '1')
#         #     if orient > 0:
#         #         feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#         #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=30)
#         #     elif orient < 0:
#         #         feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#         #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="red", arrowshaft_height=athick, label_size = lsize, label_angle=150)
#         # elif contig == pos[0] and st > 0 and i in byss:
#         #     orient = int(gene[i][3] + '1')
#         #     if orient > 0:
#         #         feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#         #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="green", arrowshaft_height=athick, label_size = lsize, label_angle=30)
#         #     elif orient < 0:
#         #         feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#         #         gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="green", arrowshaft_height=athick, label_size = lsize, label_angle=150)
#         elif contig == pos[0] and st > 0:
#             orient = int(gene[i][3] + '1')
#             if orient > 0:
#                 feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#                 gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
#             elif orient < 0:
#                 feature = SeqFeature(FeatureLocation(st, en), strand=orient)
#                 gd_feature_set.add_feature(feature, name=i, label=lab, sigil="ARROW", color="grey", arrowshaft_height=athick,  label_size = lsize, label_angle=150)
#     feature = SeqFeature(FeatureLocation(pen-200, pen), strand=0)
#     gd_feature_set.add_feature(feature, name="Scaffold end", label=lab, sigil="BOX", color="red", arrowshaft_height=athick,  label_size = lsize, label_angle=30)
#     c +=1
#


# gdd.draw(format='linear', pagesize='A4', fragments=1, start=0, end=size, track_size=0.5)
#
# gdd.write("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/recQ_helicases/ts_5_primes_40kb.pdf", "pdf")

tel_rep = 'TTAGGG' ##

telomere_sites = [m.start() for m in re.finditer(tel_rep,str(genome[pos[0]][pos[1]:pos[2]]))]

print telomere_sites
