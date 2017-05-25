__author__ = 'mjohnpayne'

from time import sleep as sl


infile = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/An_aspgd.gaf","r")
outfile = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/An_GO.map","w")

dict = {}
for line in infile:
    if line[0] != "!":
        col=line.split('\t')
        an = col[10].split('|')[0]
        if an not in dict:
            dict[an] = [col[4]]
        else:
            dict[an].append(col[4])

for i in list(sorted(dict.keys())):
    outfile.write(i+'\t'+", ".join(dict[i])+'\r')

outfile.close()

## gene \t GO term one GO term at a time
#
# outfile = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/An_GO_single.map","w")
# infile2 = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/An_aspgd.gaf","r")
#
# outfile.write('gene_id\tgo_id\n')
# for line in infile2:
#     if line[0] != "!":
#         col=line.split('\t')
#         outfile.write(col[10].split('|')[0]+'\t'+col[4]+'\n')
# outfile.close()
#



## GO \t genes all genes for each go term

# outfile3 = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/GO_An_single.map","w")
# infile3 = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/An_aspgd.gaf","r")
#
# outfile3.write('go_id\tgene_id\n')
# outdict = {}
# for line in infile3:
#     if line[0] != "!":
#         col=line.split('\t')
#         an = col[10].split('|')[0]
#         go = col[4]
#         if go not in outdict:
#             outdict[go] = [an]
#         else:
#             outdict[go].append(an)
# for i in sorted(outdict.keys()):
#     outfile3.write(i + '\t' + ",".join(outdict[i]) + '\n')
#
# outfile3.close()