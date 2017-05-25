__author__ = 'mjohnpayne'

import sys

ingff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/tsta1_working_models.gff3",'r')

genes = 0
exons = 0

for line in ingff:
    if "\tgene\t" in line:
        genes += 1
    elif "\tCDS\t" in line:
        exons += 1
print "Genes:", genes
print "Exons:", exons
print "Average exons per gene:",float(exons)/genes


#vary to find gene set info

ingff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/tsta1_working_models.gff3",'r') #sys.argv[1]

inls = sys.argv[2]

genes = 0
exons = 0

for line in ingff:
    if "\tgene\t" in line:
        genes += 1
    elif "\tCDS\t" in line:
        exons += 1
print "Genes:", genes
print "Exons:", exons
print "Average exons per gene:",float(exons)/genes