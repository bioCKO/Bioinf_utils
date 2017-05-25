__author__ = 'mjohnpayne'

hrslist = ["9H","11H","15H","17H"]

#infiles = ["/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/profile/tRNA/4way_trna_TSS_heatmap/4way_trna.tss.9H.heatmap.txt","/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/profile/tRNA/4way_trna_TSS_heatmap/4way_trna.tss.11H.heatmap.txt","/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/profile/tRNA/4way_trna_TSS_heatmap/4way_trna.tss.15H.heatmap.txt","/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/profile/tRNA/4way_trna_TSS_heatmap/4way_trna.tss.17H.heatmap.txt"]

combined = open("/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/profile/all_an/all_an_genes_TSS_heatmap/all_genes.tss.9H,11H,15H,17H.heatmap_shortened.txt","w")

outdict = {}

for i in hrslist:
    con = "\t"+i+"_"
    tmp = open("/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/profile/all_an/all_an_genes_TSS_heatmap/all_genes.tss."+i+".heatmap.txt","r")
    for j in tmp:
        col = j.strip().split('\t')
        if col[0] == "name":
            if col[0] not in outdict:
                outdict[col[0]] = i+"_" + con.join(col[81:-80])
            else:
                outdict[col[0]] += "\t"+i+"G1\t"+i+"G2\t"+i+"G3\t"+i+"G4\t"+i+"G5\t"+i+"G6\t"+i+"G8\t"+i+"G8\t"+i+"G9\t"+i+"G10\t" + i+"_" +con.join(col[81:-80])
        else:
            if col[0] not in outdict:
                outdict[col[0]] = "\t".join(col[81:-80])
            else:
                outdict[col[0]] += "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" + "\t".join(col[81:-80])

combined.write("name\t" + outdict["name"] + '\n')

for i in outdict:
    if i != "name":
        combined.write(i +"\t"+outdict[i] + "\n")

combined.close()
