__author__ = 'mjohnpayne'


from time import sleep as sl

ingff = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.sorted.gff","r").read().split('\n')

outfile = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13_danpos2_genefile.txt","w")
gene = ''
chrom = ''
strand = ''
txstart = ''
txend = ''
exonst = []
exonend = []
exonc = 0
cdsStart = 0
cdsEnd = 0
for i in ingff[9:]:
    print i
    col = i.split('\t')
    if col[2] == "gene":
        outfile.write(str(gene)+'\t'+str(chrom)+'\t'+str(strand)+'\t'+str(txstart)+'\t'+str(txend)+'\t'+str(cdsStart)+'\t'+str(cdsEnd)+'\t'+str(exonc)+'\t'+','.join(exonst)+'\t'+','.join(exonend)+'\n')
        exonst = []
        exonend = []
        exonc = 0
        gene = col[8].split(';')[1].replace("Name=","")
        chrom = col[0]
        strand = col[6]
        txstart = col[3]
        txend = col[4]
    elif col[2] == "CDS":
        cdsStart = col[3]
        cdsEnd = col[4]
    elif col[2] == "exon":
        exonc +=1
        exonst.append(str(col[3]))
        exonend.append(str(col[4]))

outfile.write(str(gene)+'\t'+str(chrom)+'\t'+str(strand)+'\t'+str(txstart)+'\t'+str(txend)+'\t'+str(cdsStart)+'\t'+str(cdsEnd)+'\t'+str(exonc)+'\t'+','.join(exonst)+'\t'+','.join(exonend)+'\n')

outfile.close()