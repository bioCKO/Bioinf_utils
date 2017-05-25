__author__ = 'mjohnpayne'

from Bio import SeqIO
from reportlab.lib.units import cm
import reportlab.lib.colors as colours
from Bio.Graphics import BasicChromosome
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
import sys

ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')

print len(sys.argv)

def remove_dupes(ls):
    ls = list(set(ls))
    return ls


if len(sys.argv) == 3:
    vcf1 = sys.argv[1]
elif len(sys.argv) == 4:
    vcf1 = sys.argv[1]
    vcf2 = sys.argv[2]
elif len(sys.argv) == 5:
    vcf1 = sys.argv[1]
    vcf2 = sys.argv[2]
    vcf3 = sys.argv[3]

elif len(sys.argv) == 6:
    vcf1 = sys.argv[1]
    vcf2 = sys.argv[2]
    vcf3 = sys.argv[3]
    vcf4 = sys.argv[4]
out = sys.argv[-1]

centromeres = {"c-100":["100",2110000, 2140000,None,""],'c-102': ["102",2610000, 2640000,None,""], 'c-103': ["103",1520000, 1550000,None,""],'c-93': ["93",1350000, 1370000,None,""],'c-95': ["95",760000, 790000,None,""],'c-96': ["96",3520000, 3550000,None,""],'c-99': ['99',2850000, 2880000,None,""],'c-97': ['97',2050000, 2070000,None,""]}

cent_list = ["c-100","c-99","c-102","c-103","c-93","c-95","c-96","c-97"]


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

## Generate list of contig,length tuples then sort in decending order
genome = {}
entries = []
lengths = []
for i in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
    #genome[i.id] = i
    entries.append((i.id,len(i.seq)))
    lengths.append(len(i.seq))

entries = sorted(entries,key=lambda x: x[1],reverse=True)
## Create gene SeqFeatures as subsets of contig seqRecord objects



def make_featuredicts(invcf):
    geno = {}
    for i in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        geno[i.id] = i
    inv = open(invcf,"r")
    for i in inv:
        if i[0] != "#":
            col = i.split('\t')
            cont = col[0]
            st = int(col[1])
            details = col[7].split(";")
            end = int(details[6][4:])
            ident = col[2]
            tp = invcf.split('/')[-1][-7:-4]
            endcont = details[5][5:]
            if tp == "DEL":
                geno[cont].features.append(SeqFeature(FeatureLocation(st, end, strand = None),type=tp,id=ident,qualifiers={"color":[]}))
            elif tp == "INV":
                geno[cont].features.append(SeqFeature(FeatureLocation(st, end, strand = -1),type=tp,id=ident,qualifiers={"color":[]}))
            elif tp == "DUP":
                geno[cont].features.append(SeqFeature(FeatureLocation(st, end, strand = 1),type=tp,id=ident,qualifiers={"color":[]}))
            elif tp == "TRA":
                geno[cont].features.append(SeqFeature(FeatureLocation(st, st+10, strand = 1),type=tp,id=ident,qualifiers={"color":[]}))
                geno[endcont].features.append(SeqFeature(FeatureLocation(end, end+10, strand = 1),type=tp,id=ident,qualifiers={"color":[]}))
    return geno


# for i in genome:
#     # for j in gene:
#     #     if gene[j][0] == genome[i].id:
#     #         direc = int(gene[j][3] + '1')
#     #         genome[i].features.append(SeqFeature(FeatureLocation(gene[j][1], gene[j][2], strand = direc),type='gene',id=j,qualifiers={'locus_tag':[j]}))
#     for j in centromeres:
#         if centromeres[j][0] == genome[i].id:
#             genome[i].features.append(SeqFeature(FeatureLocation(centromeres[j][1], centromeres[j][2], strand = None),type='gene',id=j,qualifiers={'locus_tag':[centromeres[j][4]]}))

#genome = make_featuredicts(genome)

## telomere length - rounded ends of chromosome size
max_len = max(lengths)
telomere_length = 40000

chr_diagram = BasicChromosome.Organism()
#chr_diagram.page_size = (60*cm, 21*cm)
chr_diagram.page_size = (40*cm, 21*cm)

fill = colours.CMYKColorSep(0.4, 0, 0, 0.3,density=0.4, spotName='PMS_7496')

for index, (name, length) in enumerate(entries):
    if length > 80000:
        features = []
        for j in centromeres:
                    if centromeres[j][0] == name:
                        features += [SeqFeature(FeatureLocation(centromeres[j][1], centromeres[j][2], strand = None),type='gene',id=j,qualifiers={'locus_tag':[centromeres[j][4]],'color':["grey"]})]
        for j in make_featuredicts(vcf1)[name].features:
                    j.qualifiers['color'] = ["red"]
                    features += [j]
        if len(sys.argv) == 4:
            for i in make_featuredicts(vcf2)[name].features:
                        i.qualifiers['color'] = ["blue"]
                        features += [i]
        elif len(sys.argv) == 5:
            for i in make_featuredicts(vcf2)[name].features:
                        i.qualifiers['color'] = ["blue"]
                        features += [i]
            for i in make_featuredicts(vcf3)[name].features:
                        i.qualifiers['color'] = ["green"]
                        features += [i]
        elif len(sys.argv) == 6:
            for i in make_featuredicts(vcf2)[name].features:
                        i.qualifiers['color'] = ["blue"]
                        features += [i]
            for i in make_featuredicts(vcf3)[name].features:
                        i.qualifiers['color'] = ["green"]
                        features += [i]
            for i in make_featuredicts(vcf4)[name].features:
                        i.qualifiers['color'] = ["yellow"]
                        features += [i]
        # for i in features[:10]:
        #     print i,i.qualifiers
        for j in centromeres:
            if centromeres[j][0] == name:
                features += [SeqFeature(FeatureLocation(centromeres[j][1], centromeres[j][2], strand = None),type='gene',id=j,qualifiers={'locus_tag':[centromeres[j][4]],'color':["grey"]})]

        # for j in features:
        #     print j.qualifiers
        #for f in features: f.qualifiers["color"] = [index+2]
        cur_chromosome = BasicChromosome.Chromosome(name)
        cur_chromosome.scale_num = max_len + 2 * telomere_length
        cur_chromosome.label_size = 8
        cur_chromosome.chr_percent = 0.2
        cur_chromosome.label_sep_percent = 0.1


        #Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.fill_color = fill
        start.scale = telomere_length
        cur_chromosome.add(start)

        #Add a body - using bp as the scale length here.
        body = BasicChromosome.AnnotatedChromosomeSegment(length,features)
        body.scale = length
        body.fill_color = fill
        cur_chromosome.add(body)

        #Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.fill_color = fill
        end.scale = telomere_length
        cur_chromosome.add(end)

        #This chromosome is done
        chr_diagram.add(cur_chromosome)

scale = BasicChromosome.Chromosome('Legend')
scale.scale_num = max_len + 2 * telomere_length
scale.label_size = 8
scale.chr_percent = 0.3
scale.label_sep_percent = 0.12


acc1id = sys.argv[1].split('/')[-1][-7:-4]
feats = [SeqFeature(FeatureLocation(100000, 100500, strand = 1),type='gene',id=acc1id,qualifiers={'locus_tag':[acc1id],'color':["red"]})]

accid = "Putative Centromeres"
feats += [SeqFeature(FeatureLocation(900000, 900500, strand = 1),type='gene',id=accid,qualifiers={'locus_tag':[accid],'color':["grey"]})]

if len(sys.argv) == 4:
    acc2id = sys.argv[2].split('/')[-1][-7:-4]
    feats += [SeqFeature(FeatureLocation(300000, 300500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':["blue"]})]

if len(sys.argv) == 5:
    acc2id = sys.argv[2].split('/')[-1][-7:-4]
    feats += [SeqFeature(FeatureLocation(300000, 300500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':["blue"]})]
    acc3id = sys.argv[3].split('/')[-1][-7:-4]
    feats += [SeqFeature(FeatureLocation(500000, 500500, strand = 1),type='gene',id=acc3id,qualifiers={'locus_tag':[acc3id],'color':["green"]})]

if len(sys.argv) == 6:
    acc2id = sys.argv[2].split('/')[-1][-7:-4]
    feats += [SeqFeature(FeatureLocation(300000, 300500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':["blue"]})]
    acc3id = sys.argv[3].split('/')[-1][-7:-4]
    feats += [SeqFeature(FeatureLocation(500000, 500500, strand = 1),type='gene',id=acc3id,qualifiers={'locus_tag':[acc3id],'color':["green"]})]
    acc4id = sys.argv[4].split('/')[-1][-7:-4]
    feats += [SeqFeature(FeatureLocation(700000, 700500, strand = 1),type='gene',id=acc4id,qualifiers={'locus_tag':[acc4id],'color':["yellow"]})]

#Add a body - using bp as the scale length here.
body = BasicChromosome.AnnotatedChromosomeSegment(1000000,feats)
body.scale = 1000000
body.fill_color = fill
scale.add(body)



# This chromosome is done
chr_diagram.add(scale)


chr_diagram.draw(out, sys.argv[-1].split('/')[-1].replace('.pdf',''))
