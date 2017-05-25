__author__ = 'mjohnpayne'


from time import sleep as sl

infile = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/orthologue_stuff/eurot_groups_by_gene.txt","r").readlines()
outfile = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/orthologue_stuff/Pm_An_groups_by_gene.txt","w")

genes = {}
outfile.write(infile[0])
for i in infile[1:]:
    col = i.strip('\n').split('\t')
    if "PMAA" in col[0]:
        outfile.write(col[0]+"\t"+col[1]+'\t'+",".join([x for x in col[2].split(',') if "ANID" in x])+"\n")

outfile.close()