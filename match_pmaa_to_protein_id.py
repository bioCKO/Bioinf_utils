
import re

infile = open('/Users/mjohnpayne/Desktop/merge_pm_embl/Pm_annot.embl',"r")
outfile = open('/Users/mjohnpayne/Desktop/merge_pm_embl/pmaa_prot_ids_link.txt',"w")


count = 1
link = {}

for line in infile:
    if line.find('/locus_tag=') > 0:
        if count == 1:
            pmaa = line[33:44]
            count = 0

    if line.find('protein_id') > 0:
        prot_id = line[34:44]
        link[pmaa] = prot_id
        count = 1

for info in link:
#    print info
#    print link[info]
    outfile.writelines(info + "\t" + link[info] + "\n")


outfile.close()
infile.close()
