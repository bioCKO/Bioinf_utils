__author__ = 'mjohnpayne'

import sys
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace
import math

# infile = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_node_assignments_tree.nwk','r')
# outfile = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_node_assignments_tree_nos_only.nwk','w')
# infile = infile.read()

t = PhyloTree('/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_node_assignments_tree.nwk', format=1)

# ts = TreeStyle()
#
# t.show(tree_style=ts)

for node in t:#.iter_search_nodes():
    # name = node.name
    # name = name[name.find("_")+1:]
    # node.name = name
    print node.name
    if node.name == "41":
        node.dist = 5e-05

# t.write(outfile='/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_node_assignments_tree_nos.nwk',format=1,)

ts = TreeStyle()

t.show(tree_style=ts)

#
# print infile
#
# def rem_names(s):
#     if "_" in s:
#         p = s.find("_")
#         nstr = s[:p-3]+s[p+1:]
#         return rem_names(nstr)
#     else:
#         return s
#
# def rem_supports(s):
#     if ":" in s:
#         p = s.find(":")
#         nstr = s[:p]+s[p+11:]
#         return rem_supports(nstr)
#     else:
#         return s
#
# def swap_no_for_name(s):
#     if "_" in s:
#         p = s.find("_")
#         nstr = s[p+1:]
#         n = nstr.find(':')
#         id = nstr[:n]
#         print id
#         # id = id[4:] + ':' + id[:3]
#         # print id
#         # nstr = s[:p-3]+s[p+11:]
#         return swap_no_for_name(nstr)
#     else:
#         return s
#
# swap_no_for_name(infile)
#
#
# out1 = rem_names(infile)
#
# print out1
# #
# # out2 = rem_supports(out1)
# #
# # print out2
# #
# # out3 = rem_supports(infile)
# #
# # print out3
# #
# outfile.write(out1)
# #
# # outfile.close()