__author__ = 'mjohnpayne'

import sys
from time import sleep as sl

# if len(sys.argv) < 2:
#     print "\n\nTakes tab delimited detailed gene associated processed danpos outputs (made by danpos_to_promoter_changes_v3_3_cats_orig.py)\n and makes list of genes in file with less than 'cutoff' FDR in point changes"
#     print "\nUsage: python sig_changes_for_go_enrichment.py 'input.txt file' cutoff\n"
# else:
inp = sys.argv[1]

cutoff = float(sys.argv[2])



out1 = inp.replace(".txt","_sig_point_changes_increase_from_1st.txt")
outf1 = open(out1,"w")
out2 = inp.replace(".txt","_sig_point_changes_decrease_from_1st.txt")
outf2 = open(out2,"w")

for i in open(inp,"r").readlines()[1:]:
    col = i.strip().split('\t')
    if col[17] != "":
        pvals = map(float,col[7].strip(",").split(","))
        if min(pvals) < cutoff:
            fold_changes = map(float,col[6].strip(",").split(","))
            if min(fold_changes) < 0:
                outf1.write(col[1]+"\n")
            elif max(fold_changes) > 0:
                outf2.write(col[1]+"\n")