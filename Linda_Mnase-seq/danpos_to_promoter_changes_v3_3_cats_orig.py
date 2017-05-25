__author__ = 'mjohnpayne'

from time import sleep as sl
import urllib
import sys

cutoff = 0.01



## 4 cats to get overall counts for, also do gene by gene counts of intragenic, promoter and pstream region - count and position relative to TSS
## 1 - intragenic (TSS+50 to end of gene)
## 2 - promoter (TSS-500 to TSS+50)
## 3 - intergenic (all other areas)
## 4 - entire upstream region of all genes (will contain overlap within set and overlap promoter and intergenic groups

## 4 stats
## 1 - point_diff_FDR  < 0.01 for general dynamic nucleosomes
## 2 - treat2control_dis between 50 and 90 and point_diff_FDR < 0.01 for position changes
## 3 - fuzziness_diff_FDR < 0.01 and point_diff_FDR < 0.01 for fuzziness changes
## 4 - smt_diff_FDR < 0.01 and point_diff_FDR < 0.01 for occupancy changes
##               Some time we may want to remove the subset that shows both low smt_diff_FDR and low fuzziness_diff_FDR

## output dict should be gene by gene with subdict of positions with subsub dict of change types gene[promoter][position_changes] = [(position relative to TSS, change statistic, FDR value),(...,...,...)]
input = sys.argv[1]

pks = open(input,"r").read().split('\r')

sample = input.split('/')[-1].replace(".positions.integrative.txt","")



#pks = open("/Volumes/MP_HD/Linda_MNase_Seq/danpos_r1_out/Volumes_MP_HD_Linda_MNase_Seq_bam_11H-Volumes_MP_HD_Linda_MNase_Seq_bam_15H.positions.integrative.txt","r").read().split('\r')

def get_all_nucs(infile):
    odict = {}
    for i in infile[1:]:
        col = i.split('\t')
        if col[0] not in odict:
            odict[col[0]] = [col[0]+"_"+str(col[3])]
    return odict

all_peaks = {}

all_nucs = get_all_nucs(pks)

def get_sig_nuc(infile,cutoff,type):
    odict = {}
    sites = []
    for i in infile[1:]:
        col = i.split('\t')
        odict[col[0]] = []
    for i in infile[1:]:
        col = i.split('\t')

        if type == "points":
            if float(col[17]) < cutoff:
                    dir = ''
                    odict[col[0]].append((int(col[3]),col[15],col[17]))
                    sites.append(col[0]+'_'+col[3])
        elif type == "positions":
            if float(col[17]) < cutoff and 50 < int(col[7]) < 90:
                    dir = ''
                    odict[col[0]].append((int(col[3]),col[7],col[17]))
                    sites.append(col[0]+'_'+col[3])
        elif type == "fuzzes":
            if float(col[22]) < cutoff and float(col[17]) < cutoff:
                    dir = ''
                    odict[col[0]].append((int(col[3]),col[20],col[22]))
                    sites.append(col[0]+'_'+col[3])
        elif type == "occupancys":
            if float(col[12]) < cutoff and float(col[17]) < cutoff:
                    dir = ''
                    odict[col[0]].append((int(col[3]),col[10],col[12]))
                    sites.append(col[0]+'_'+col[3])
    return odict,sites

ingff = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.sorted.gff","r")

def gene_prom_dict(gff,peaks_file,distance):
    #out = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/"+sample+"_"+str(distance)+"_upstream_"+str(cutoff)+"point_changes.txt","w")
    #out.write("Chromosome\tGene_ID\tNumber of point changes\tPositions of point changes\tpoint change stats\tpoint change FDR\tNumber of position changes\tPositions of position changes\tposition change stats\tposition change FDR\tNumber of fuzziness changes\tPositions of fuzziness changes\tfuzziness change stats\tfuzziness change FDR\tNumber of occupancy changes\tPositions of occupancy changes\toccupancy change stats\toccupancy change FDR\t\n")
    genes = []
    sites = {}
    points,sites["points"] = get_sig_nuc(peaks_file,0.01,"points")
    positions,sites["positions"] = get_sig_nuc(peaks_file,0.01,"positions")
    fuzzes,sites["fuzzes"] = get_sig_nuc(peaks_file,0.01,"fuzzes")
    occupancys,sites["occupancys"] = get_sig_nuc(peaks_file,0.01,"occupancys")
    p_size = distance
    genedict = {}
    ends = {}
    final_dict = {}
    gene_chrom = {}
    gene_inf = {}
    for line in gff:
        if line[0] != "#":
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                cont = col[0]
                st = col[3]
                en = col[4]
                pmaa = col[8].split(';')[1].replace("Name=","")
                notes = col[8].split(';')[2]
                gene = ''
                if notes[:5] == "Gene=":
                    gene = notes.replace("Gene=","")
                    notes = col[8].split(';')[3]
                    notes = notes.replace("Note=","")
                    notes = gene+": "+urllib.unquote_plus(notes)
                else:
                    notes = col[8].split(';')[2]
                    notes = notes.replace("Note=","")
                    notes = urllib.unquote_plus(notes)
                gene_inf[pmaa] = notes
                orient = col[6]
                gene_chrom[pmaa] = col[0]
                if cont not in genedict:
                    genedict[cont] = {}
                    genedict[cont][st] = (pmaa,en,orient)
                else:
                    genedict[cont][st] = (pmaa,en,orient)
            elif col[2] == 'chromosome':
                ends[col[0]] = col[4]
    changetypes = [points,positions,fuzzes,occupancys]
    cnames = ["points","positions","fuzzes","occupancys"]
    for i in genedict:

        starts = map(str,sorted(map(int,genedict[i].keys())))

        for j in range(len(starts)):

            st = starts[j]
            en = genedict[i][st][1]
            pmaa = genedict[i][st][0]
            orient = genedict[i][st][2]

            prevend = 0
            nextstart = 0
            if j == len(starts)-1:
                nextstart = int(ends[i])
            else:
                nextstart = int(starts[j+1])

            if j != 0:
                prevend = genedict[i][starts[j-1]][1]

            prom = []
            gene = []
            prom_extended = []

            if orient == "+":
                if int(prevend) > int(st)-p_size:
                    prom = [int(prevend),int(st)+50]
                else:
                    prom = [int(st)-p_size,int(st)+50]
                gene = [int(st)+50,en]
                prom_extended = [int(prevend),int(st)+50]

            elif orient == "-":
                if nextstart < int(en)+p_size:
                    prom = [int(en)-50,nextstart]
                else:
                    prom = [int(en)-50,int(en)+p_size]
                gene = [int(st),int(en)-50]
                prom_extended = [int(en)-50,nextstart]

            postype = [gene,prom,prom_extended]
            pnames = ["gene","prom","prom_extended"]
            final_dict[pmaa] = {"gene":{},"prom":{},"prom_extended":{}}
            for t in range(len(postype)):
                z = postype[t]
                t_name = pnames[t]
                for changep in range(len(changetypes)):
                    change_name = cnames[changep]
                    final_dict[pmaa][t_name][change_name] = []
                    for y in changetypes[changep][i]:
                        pos = int(y[0])
                        score = y[1]
                        fdr = y[2]
                        if int(z[0]) < pos < int(z[1]):
                            #print pmaa,t_name,change_name,pos,score,fdr,z[0],z[1]
                            if orient == "+":
                                final_dict[pmaa][t_name][change_name].append((pos-int(st),score,fdr,pos))
                            elif orient == "-":
                                final_dict[pmaa][t_name][change_name].append((int(en)-pos,score,fdr,pos))
    return final_dict,gene_chrom,sites,gene_inf





    # out.close()
    # outgenes = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/"+str(distance)+"_upstream_"+str(cutoff)+"_cutoff_point_changes_genes_only.txt","w")
    # outgenes.write('\n'.join(genes))
    # outgenes.close


gdict,chrom,sites,inf = gene_prom_dict(ingff,pks,500)

## Gets counts of changes (named as chromosome_position makes list of lists of change type and gene region ###

def get_counts(genedict,chroms,out,site):
    outfile = open(out,"w")
    pnames = ["gene","prom","prom_extended"]
    cnames = ["points","fuzzes","occupancys","positions"]
    regions_list = {"gene":[],"prom":[],"prom_extended":[]}
    stat_type_list = {"points":[],"fuzzes":[],"occupancys":[],"positions":[]}
    change_dict = {}
    for i in pnames:
        change_dict[i] = {}
        for j in cnames:
            change_dict[i][j] = []
    print change_dict
    for gene in genedict:
        for pos in pnames:
            for stat in cnames:
                if len(genedict[gene][pos][stat]) > 0:
                    for c in genedict[gene][pos][stat]:
                        name = chroms[gene].replace("_A_nidulans_FGSC_A4","") + "_" + str(c[3])
                        change_dict[pos][stat].append(name)
                        regions_list[pos].append(name)
                        stat_type_list[stat].append(name)
                        if name in site[stat]:
                            site[stat].remove(name)
    outfile.write('\t'+'\t'.join(cnames)+'\n')
    for pos in pnames:
        outfile.write(pos)
        for stat in cnames:
            outfile.write('\t'+str(len(set(change_dict[pos][stat]))))
        outfile.write('\n')
    outfile.write("intergenic")
    for stat in cnames:
        outfile.write('\t' + str(len(set(site[stat]))))
    outfile.write('\n')
    outfile.close()

    nregions_list = {}
    for l in regions_list:
        nregions_list[l] = list(set(regions_list[l]))
    regions_list = nregions_list
    nstat_type_list = {}
    for l in stat_type_list:
        nstat_type_list[l] = list(set(stat_type_list[l]))
    stat_type_list = nstat_type_list

    for p in regions_list:
        outf = open("/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/4_stats_4_cats/"+sample+"_all_changed_"+p+"position_nucleosomes.txt",'w')
        outf.write('\n'.join(regions_list[p]))
        outf.close()

    for p in stat_type_list:
        outf = open("/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/4_stats_4_cats/"+sample+"_"+p+"_changed_all_position_nucleosomes.txt",'w')
        outf.write('\n'.join(stat_type_list[p]))
        outf.close()

    fuz_occ_ls = set(stat_type_list["fuzzes"]).intersection(stat_type_list["occupancys"])
    outf = open("/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/4_stats_4_cats/"+sample+"_fuzzys_and_occupancys_changed_all_position_nucleosomes.txt",'w')
    outf.write('\n'.join(fuz_occ_ls))
    outf.close()



get_counts(gdict,chrom,"/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/4_stats_4_cats/"+sample+"_counts_cats_out.txt",sites)



## per gene output of counts, positions, scores and FDRs of all types of changes ##

## output dict  dict[gene][area][stat_type] = [(position relative to TSS, change statistic, FDR value),(...,...,...)]

def write_per_gene(gene_dict,cdict,outpath,type):
    of = open(outpath,"w")
    of.write("Chromosome\tGene_ID\tGene Description\tNumber of point changes\tPositions of point changes\tnucleosome_id_of_point_changes\tpoint change stats\tpoint change FDR\tNumber of fuzziness changes\tPositions of fuzziness changes\tnucleosome_id_of_fuzziness_changes\tfuzziness change stats\tfuzziness change FDR\tNumber of occupancy changes\tPositions of occupancy changes\tnucleosome_id_of_occupancy_changes\toccupancy change stats\toccupancy change FDR\tNumber of position changes\tPositions of position changes\tnucleosome_id_of_position_changes\tposition change stats\tposition change FDR\n")
    cnames = ["points","fuzzes","occupancys","positions"]


    for gene in gene_dict:
        of.write(cdict[gene]+'\t'+gene+'\t'+inf[gene])
        for stat in cnames:
            of.write("\t"+str(len(gdict[gene][type][stat])))
            if len(gdict[gene][type][stat]) > 1:
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(str(i[0])+",")
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(cdict[gene].replace("_A_nidulans_FGSC_A4","")+"_"+str(i[3])+",")
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(str(i[1])+",")
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(str(i[2])+",")
            elif len(gdict[gene][type][stat]) <= 1:
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(str(i[0]))
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(cdict[gene].replace("_A_nidulans_FGSC_A4","")+"_"+str(i[3])+",")
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(str(i[1]))
                of.write("\t")
                for i in gdict[gene][type][stat]:
                    of.write(str(i[2]))
        of.write('\n')
    of.close()

pnames = ["gene","prom","prom_extended"]

for i in pnames:
    write_per_gene(gdict,chrom,"/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/4_stats_4_cats/"+sample+"_"+i+"_out.txt",i)