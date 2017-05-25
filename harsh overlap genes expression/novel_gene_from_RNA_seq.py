import random
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC

##in_wig = open("/Volumes/Storage/PhD/phage_bams/NoPhage_1_sort.bam.wig")
##
##outfile = open("/Volumes/Storage/PhD/phage_bams/NoPhage_1_expression_genes_25bp_window.txt",'w')
##
##in_gb = SeqIO.parse('/Users/mjohnpayne/Documents/Uni/phd/wt_genome/Pm_from_gb_old_acc.gb', "genbank")

######################

# in terminal: python novel_gene_from_seq.py [wig path] [Pm_from_gb_old_acc.gb path] [outfile.txt] [pmfa1_annot_scaf.fasta path]

######################

in_wig = open(sys.argv[1],'r')

in_gb = SeqIO.parse(sys.argv[2], "genbank")

outfile = open(sys.argv[3],'w')

genome = {}
for record in SeqIO.parse(sys.argv[4], "fasta"):
        genome[record.id] = record.seq

# parse wig file
chrom = {}
contig = ""
for line in in_wig:
    if line[0:12] == "variableStep":
        columns = line.split(" ")
        contig = str(columns[1])
        contig = contig.replace("chrom=","")
        chrom[contig] = {}
    elif line[0] == "#":
        continue
    else:
        line = line.strip("\n").split("\t")
        chrom[contig][line[0]] = [float(line[1]),0]
        
genome2 = {}
counter = 1
step = 25

# remove genome contigs if not in wig

for contig in genome:
        if contig not in chrom:
                continue
        else:
                info = genome[contig]
                genome2[contig] = info

# adds [0,0] data for window not in wig file

for contig in genome2:
    counter = 1
#    outfile.writelines(contig + "\n")
    while counter < len(genome[contig]):
        if str(counter) in chrom[contig]:
            name = counter
            coverage = chrom[contig][str(counter)]
            outline = str(name) + "\t" + str(coverage) + "\n"
#            outfile.writelines(outline)           
        else:
            chrom[contig][str(counter)] = [0,0]
            name = str(counter)
            coverage = chrom[contig][str(counter)]
            outline = str(name) + "\t" + str(coverage) + "\n"
#            outfile.writelines(outline)
        counter = int(counter) + int(step)

#make dictionary of genes, contig and position

genes = {}
        
for record in in_gb:
    contig = record.name
#    print contig
    contig_len = len(record)
    index = 1
    for i,feature in enumerate(record.features):
        if feature.type == 'gene':
            loc = str(feature.location)
            loc = loc.strip("\[""\<""(-)""(+)""\]")
            loc = loc.replace(">","").split(":")
            loc = map(int, loc)
            pmaa = feature.qualifiers['locus_tag'][0]
            genes[pmaa] = [contig] + loc

# Step through window and check if gene is present and add to dict

winstart = 1
winend = 25
step = 25
windowcount = 0
index = 0
counter = {}


for contig in genome2:
    winstart = 1
    winend = 25
    step = 25
    windowcount = 0
    while winend < len(genome2[contig]):
        for gene in genes:
            if genes[gene][0] == contig and winstart > genes[gene][1] and winend < genes[gene][2]:
                chrom[contig][str(winstart)][1] = 1
            elif genes[gene][0] == contig and winstart < genes[gene][1] or winend > genes[gene][2]:
                if chrom[contig][str(winstart)][1] == 1:
                        continue
                else:
                        chrom[contig][str(winstart)][1] = 0
            else:
                    continue
        outfile.writelines(str(contig) + " " + str(winstart) + " " + str(chrom[contig][str(winstart)]) + "\n")
        windowcount = 0
        winstart = winstart + step
        winend = winend + step
    print str(contig) + "done"
        
outfile.close()




