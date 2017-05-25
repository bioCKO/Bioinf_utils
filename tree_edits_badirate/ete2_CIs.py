__author__ = 'mjohnpayne'

from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, PhyloTree, PieChartFace
import math


#
def my_layout(node):
    if node.name != "43":
        print node.name
        if node.is_leaf():
             # If terminal node, draws its name
            name_face = AttrFace("species", fsize=16)
            pie = PieChartFace([changes[node.name][0],changes[node.name][1]],changes[node.name][2],changes[node.name][2],["Green","Red"])
            pie.opacity = 0.5
            #name_face = AttrFace("name", fsize=16)
            faces.add_face_to_node(name_face, node, column=0, position="branch-right")
            faces.add_face_to_node(pie, node, column=0, position="float")
        else:
             # If internal node, draws label with smaller font size
            #name_face = AttrFace("name", fsize=10)
            pie = PieChartFace([changes[node.name][0],changes[node.name][1]],changes[node.name][2],changes[node.name][2],["Green","Red"])
            pie.opacity = 0.5
            #faces.add_face_to_node(name_face, node, column=0, position="branch-right")
            faces.add_face_to_node(pie, node, column=0, position="float")


ts = TreeStyle()
# Do not add leaf names automatically ts.show_leaf_name = False
# Use my custom layout
ts.show_leaf_name = False
ts.layout_fn = my_layout





#t = PhyloTree('/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_node_assignments_tree_nos.nwk', format=1)

t = PhyloTree('/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/gain_loss_tree_frm_orthogroups/CI_badirate_branch_no_tree_no_names.nwk', format=1)

#dataorder = ['FRR2161','FRR3841','FRR4059','FRR3840','F4','BR2SD2','BR2','BR2SD1','G09043','G11702','G11203SD4','G11203SD3','G11203','G11203SD1','G09027SD2','G09027SD1','G09027','FRR3871','FRR3482','HR2','G11012']
#nos = ["1",'2','4','5','7','10','11','13','15','17','19','20','22','24','26','27','29','32','33','35','37']

#dataorder = ['FRR2161','FRR3841','FRR3840','FRR4059','F4','BR2SD2','BR2','BR2SD1','G09043','G11702','G11203SD4','G11203SD3','G11203','G11203SD1','G09027SD2','G09027SD1','G09027','FRR3871','FRR3482','HR2','G11012']
#nos = ["1",'2','4','5','7','10','11','13','15','17','19','20','22','24','26','27','29','32','33','35','37']

#branch_to_node = {24:25,39:38,41:43,33:35,5:5,31:0,18:22,14:16,40:42,42:41,35:36,36:32,27:30,15:2,26:29,12:15,29:26,21:19,11:4,32:34,6:11,17:20,22:17,16:18,13:3,34:33,43:37,3:6,7:13,8:14,37:39,10:10,44:31,30:24,20:21,2:8,1:7,38:40,28:28,4:9,25:27,19:23,23:1,9:12}
#print branch_to_node

inchanges = open('/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/gain_loss_tree_frm_orthogroups/CI_denovo_gene_gain_loss_table.txt','r')

changes = {}
for i in inchanges:
    #print i
    if '#' not in i:
        inf = i[i.find(">")+1:].strip('\n').split('\t')
        gain = int(inf[1])
        loss = int(inf[2])
        tot = gain+loss
        pgain = gain/float(tot)*100
        ploss = 100-pgain
        size = tot
        size = math.sqrt((size/math.pi))*3
        changes[inf[0]] = (pgain,ploss,size)

#print changes




#dt = {}

dt = {'1':'FRR2161','2':'FRR3841','4':'FRR3840','5':'FRR4059','7':'F4','34':'BR2SD2','33':'BR2','36':'BR2SD1','38':'G09043','40':'G11702','13':'G11203SD4','15':'G11203SD3','10':'G11203','11':'G11203SD1','18':'G09027SD2','17':'G09027SD1','20':'G09027','23':'FRR3871','26':'FRR3482','28':'HR2','30':'G11012','24':'Pm1'}

# for i in range(len(nos)):
#     dt[nos[i]] = dataorder[i]

def get_species_name(node_name_string):
    # Species code is the first part of leaf name (separated by an # underscore character)
    spcode = node_name_string
    # We could even translate the code to complete names
    code2name = dt
    return code2name[spcode]

t.set_species_naming_function(get_species_name)

for node in t.iter_search_nodes():
    if node.name == "43":
        node.dist = 5e-05

#t.show(tree_style=ts)
# t.render("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_copy_no_tree/CI_gain_loss_tree.pdf",tree_style=ts,w=3200,h=4800,dpi=200)
t.render("/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/gain_loss_tree_frm_orthogroups/CI_denovo_gene_gain_loss_tree.pdf",tree_style=ts,w=3200,h=4800,dpi=200)