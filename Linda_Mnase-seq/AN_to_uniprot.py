__author__ = 'mjohnpayne'

outfile = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/Linda_Mnase-seq/An_uniprot_links.txt","w")

outfile2 = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/Linda_Mnase-seq/An_uniprot_links_uniprots.txt","w")

ingaf = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/Linda_Mnase-seq/An_aspgd.gaf","r").readlines()

inuniprots = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/Linda_Mnase-seq/gp2protein.aspgd","r").readlines()

print ingaf[0:5],inuniprots[0]

uniprots = {}

for i in inuniprots:
    col = i.split('\t')
    uniprots[col[0][6:]] = col[1][10:16]

genelinks = {}
for i in ingaf[1:]:
    col = i.strip('\n').split('\t')
    geneid = col[10].split("|")[0]
    if col[1] in uniprots:
        genelinks[geneid] = uniprots[col[1]]
    else:
        genelinks[geneid] = "none"


for i in genelinks:
    outfile.write(i+'\t'+genelinks[i]+'\n')
    if genelinks[i] != "none":
        outfile2.write(i+'\t'+genelinks[i]+'\n')

outfile2.close()
outfile.close()
