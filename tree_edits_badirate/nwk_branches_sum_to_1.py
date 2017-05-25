__author__ = 'mjohnpayne'


from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace
import re



def count_tot(s,c):
    if ":" in s:
        p = s.find(":")
        num = float(s[p+1:p+11])
        c += num
        nstr = s[:p]+s[p+1:]
        return count_tot(nstr,c)
    else:
        return c

def double_branches(s,c,ls):
    c += 1
    if c <= s.count(":"):
        p = ls[c-1]
        num = s[p+1:p+11]
        num = float(num)*100
        nstr = s[:p+1] + str("%0.8f" % num) + s[p+11:]
        return double_branches(nstr,c,ls)
    else:
        return s

#tree = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_node_assignments_tree_nos.nwk",'r').read()
tree = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/gain_loss_tree_frm_orthogroups/CI_badirate_branch_no_tree.nwk",'r').read()


print count_tot(tree.read(),0)

#pos = [m.start() for m in re.finditer(':',tree)]

#print double_branches(tree,0,pos)