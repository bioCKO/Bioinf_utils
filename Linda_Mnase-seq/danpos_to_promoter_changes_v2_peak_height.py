__author__ = 'mjohnpayne'

from time import sleep as sl

cutoff = 0.05
dist = 500

pks = open("/Volumes/MP_HD/Linda_MNase_Seq/danpos_r1_out/Volumes_MP_HD_Linda_MNase_Seq_bam_11H-Volumes_MP_HD_Linda_MNase_Seq_bam_15H.peaks.integrative.txt","r").read().split('\r')


def get_sig_peaks(infile,cutoff):
    peaks = {}
    for i in infile[1:]:
        col = i.split('\t')
        if float(col[8]) < cutoff:
            dir = ''
            if float(col[6]) < 0:
                if col[0] not in peaks:
                    peaks[col[0]] = [(int(col[3]),"u",col[6])]
                else:
                    peaks[col[0]].append((int(col[3]),"u",col[6]))
            elif float(col[6]) > 0:
                if col[0] not in peaks:
                    peaks[col[0]] = [(int(col[3]),"d",col[6])]
                else:
                    peaks[col[0]].append((int(col[3]),"d",col[6]))
    return peaks


peaks = get_sig_peaks(pks,cutoff)

ingff = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.sorted.gff","r")

def gene_prom_dict(gff,pkdict,distance):
    out = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/"+str(distance)+"_upstream_"+str(cutoff)+"_cutoff_peak_height_changes.txt","w")
    out.write("Chromosome\tGene_ID\tdistance from gene start\tlog2FC(-ve = decreased occupancy at 11H relative to 15H)\n")
    genes = []
    for line in gff:
        if line[0] != "#":
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                cont = col[0]
                st = col[3]
                en = col[4]
                pmaa = col[8].split(';')[1].replace("Name=","")
                orient = col[6]
                if orient == "+":
                    pst = int(st)-distance
                    pen = int(st)
                    for i in pkdict[cont]:
                        if pst < i[0] < pen:
                            genes.append(pmaa)
                            dist = pen - i[0]
                            out.write(cont+'\t'+pmaa+'\t'+str(dist)+'\t'+i[2]+'\n')
                elif orient == "-":
                    pst = int(en)
                    pen = int(en)+distance
                    for i in pkdict[cont]:
                        if pst < i[0] < pen:
                            genes.append(pmaa)
                            dist = i[0] - pst
                            out.write(cont+'\t'+pmaa+'\t'+str(dist)+'\t'+i[2]+'\n')
    out.close()
    outgenes = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/"+str(distance)+"_upstream_"+str(cutoff)+"_cutoff_peak_height_changes_genes_only.txt","w")
    outgenes.write('\n'.join(genes))
    outgenes.close


gendict = gene_prom_dict(ingff,peaks,dist)
