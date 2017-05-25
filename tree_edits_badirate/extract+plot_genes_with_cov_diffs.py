__author__ = 'mjohnpayne'

from Bio import SeqIO
from reportlab.lib.units import cm
import reportlab.lib.colors as colours
from Bio.Graphics import BasicChromosome
from Bio.SeqFeature import SeqFeature, FeatureLocation

indiffs = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/gene_counts_CIs.tsv",'r')



def find_genes(ind):
    up = []
    down = []
    neither = []
    both = []
    for i in indiffs:
        if not i.startswith('ID'):
            col = i.strip('\n').split('\t')
            id = col[0]
            levels = map(int,col[2:])
            if 0 in levels and max(levels) > 1:
                both.append(id)
            elif 0 in levels:
                down.append(id)
            elif max(levels) > 1:
                up.append(id)
            else:
                neither.append(id)
    return up,down,neither,both

def find_strain_genes(strain,ind):
    up = []
    down = []
    neither = []
    pos = 0
    for i in indiffs:
        col = i.strip('\n').split('\t')
        if col[0] == 'ID':
            strains = col[1:]
            pos = strains.index(strain)
            print pos
        else:
            id = col[0]
            level = int(col[pos+1])
            if level == 0:
                down.append(id)
            elif level > 1:
                up.append(id)
            else:
                neither.append(id)
    return up,down,neither

#up,down,neither,both = find_genes(indiffs)

strain = "HR2"

up,down,neither = find_strain_genes(strain,indiffs)

inlist = [down,up]
names = ["Deleted","Duplicated"]

print len(down), len(up)

# inlist = [down,both,up]
# names = ["Deleted","Both","Duplicated"]

# inlist = [down]
# names = ["Deleted"]
#
# inlist = [up]
# names = ["Duplicated"]
#
ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')

print len(inlist)

def remove_dupes(ls):
    ls = list(set(ls))
    return ls


if len(inlist) == 1:
    acclis1 = inlist[0]
elif len(inlist) == 2:
    acclis1 = inlist[0]
    acclis2 = inlist[1]
elif len(inlist) == 3:
    acclis1 = inlist[0]
    acclis2 = inlist[1]
    acclis3 = inlist[2]

elif len(inlist) == 4:
    acclis1 = inlist[0]
    acclis2 = inlist[1]
    acclis3 = inlist[2]
    acclis4 = inlist[3]
out = "/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/"+strain+"_change_chromosome_pos.pdf"

##Extract gene position, contig and orientation info from gff
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

for cont in contigs:
        contigs[cont] = sorted(contigs[cont])

## Generate list of contig,length tuples then sort in decending order
genome = {}
entries = []
lengths = []
for i in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
    genome[i.id] = i
    entries.append((i.id,len(i.seq)))
    lengths.append(len(i.seq))

entries = sorted(entries,key=lambda x: x[1],reverse=True)
## Create gene SeqFeatures as subsets of contig seqRecord objects
for i in genome:
    for j in gene:
        if gene[j][0] == genome[i].id:
            direc = int(gene[j][3] + '1')
            genome[i].features.append(SeqFeature(FeatureLocation(gene[j][1], gene[j][2], strand = direc),type='gene',id=j,qualifiers={'locus_tag':[j]}))

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
        for i in acclis1:
            for f in genome[name].features:
                if f.id==i:
                    f.qualifiers['color'] = [2]
                    features += [f]
        if len(inlist) == 2:
            for i in acclis2:
                for f in genome[name].features:
                    if f.id==i:
                        f.qualifiers['color'] = [4]
                        features += [f]
        elif len(inlist) == 3:
            for i in acclis2:
                for f in genome[name].features:
                    if f.id==i:
                        f.qualifiers['color'] = [4]
                        features += [f]
            for i in acclis3:
                for f in genome[name].features:
                    if f.id==i:
                        f.qualifiers['color'] = [8]
                        features += [f]
        elif len(inlist) == 4:
            for i in acclis2:
                for f in genome[name].features:
                    if f.id==i:
                        f.qualifiers['color'] = [4]
                        features += [f]
            for i in acclis3:
                for f in genome[name].features:
                    if f.id==i:
                        f.qualifiers['color'] = [8]
                        features += [f]
            for i in acclis4:
                for f in genome[name].features:
                    if f.id==i:
                        f.qualifiers['color'] = [10]
                        features += [f]
        #for f in features: f.qualifiers["color"] = [index+2]
        cur_chromosome = BasicChromosome.Chromosome(name)
        cur_chromosome.scale_num = max_len + 2 * telomere_length
        cur_chromosome.label_size = 4
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


acc1id = names[0]
feats = [SeqFeature(FeatureLocation(100000, 100500, strand = 1),type='gene',id=acc1id,qualifiers={'locus_tag':[acc1id],'color':[2]})]

if len(inlist) == 2:
    acc2id = names[1]
    feats += [SeqFeature(FeatureLocation(300000, 300500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':[4]})]

if len(inlist) == 3:
    acc2id = names[1]
    feats += [SeqFeature(FeatureLocation(300000, 300500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':[4]})]
    acc3id = names[2]
    feats += [SeqFeature(FeatureLocation(500000, 500500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc3id],'color':[8]})]

if len(inlist) == 4:
    acc2id = names[1]
    feats += [SeqFeature(FeatureLocation(300000, 300500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc2id],'color':[4]})]
    acc3id = names[2]
    feats += [SeqFeature(FeatureLocation(500000, 500500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc3id],'color':[8]})]
    acc4id = names[3]
    feats += [SeqFeature(FeatureLocation(700000, 700500, strand = 1),type='gene',id=acc2id,qualifiers={'locus_tag':[acc4id],'color':[10]})]

#Add a body - using bp as the scale length here.
body = BasicChromosome.AnnotatedChromosomeSegment(1000000,feats)
body.scale = 1000000
body.fill_color = fill
scale.add(body)



# This chromosome is done
chr_diagram.add(scale)


chr_diagram.draw(out, "Genes gained and lost in " + strain)
