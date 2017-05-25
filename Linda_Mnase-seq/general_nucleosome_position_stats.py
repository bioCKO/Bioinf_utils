__author__ = 'mjohnpayne'

import urllib
from time import sleep as sl

pks = open("/Volumes/MP_HD/Linda_MNase_Seq/danpos_r1_out/Volumes_MP_HD_Linda_MNase_Seq_bam_11H-Volumes_MP_HD_Linda_MNase_Seq_bam_15H.positions.integrative.txt","r").read().split('\r')

ingff = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13.sorted.gff","r")


def get_all_nucs(infile):
    tot_nuc_count = 0
    odict = {}
    for i in infile[1:]:
        tot_nuc_count +=1
        col = i.split('\t')
        if col[0] not in odict:
            odict[col[0]] = [col[3]]
        else:
            odict[col[0]].append(col[3])
    return odict,tot_nuc_count

all_peaks = {}

all_nucs,tot_nuc_count = get_all_nucs(pks)


def gene_prom_dict(gff,peaks_file):
    # out = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/"+str(distance)+"_upstream_"+str(cutoff)+"point_changes.txt","w")
    # out.write("Chromosome\tGene_ID\tNumber of point changes\tPositions of point changes\tpoint change stats\tpoint change FDR\tNumber of position changes\tPositions of position changes\tposition change stats\tposition change FDR\tNumber of fuzziness changes\tPositions of fuzziness changes\tfuzziness change stats\tfuzziness change FDR\tNumber of occupancy changes\tPositions of occupancy changes\toccupancy change stats\toccupancy change FDR\t\n")
    # genes = []
    p_size = 500
    sites,nccount = get_all_nucs(peaks_file)
    genedict = {}
    ends = {}
    final_dict = {}
    gene_chrom = {}
    gene_inf = {}
    tot_size = 0
    nuc_cats = {"gene":[],"prom":[],"prom_extended":[]}
    region_sizes = {"gene":0,"prom":0,"prom_extended":0}
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
                tot_size += int(col[4])
    c = 0
    for i in genedict:
        starts = map(str,sorted(map(int,genedict[i].keys())))

        for j in range(min([20000000,len(starts)])):
            if c%100 == 0:
                print c
            c+=1

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

            nextorient = ""

            if j == len(starts)-1:
                nextorient = "-"
            else:
                nextorient = genedict[i][starts[j+1]][2]

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
            final_dict[pmaa] = {"gene":[],"prom":[],"prom_extended":[]}
            for t in range(len(postype)):
                z = postype[t]
                t_name = pnames[t]
                region_sizes[t_name] += int(z[1])-int(z[0])
                if "prom" in t_name:
                    if orient == "-" and nextorient == "+":
                        if p_size < (nextstart - int(en)) < 2*p_size:
                            overlap = 2*p_size - (nextstart - int(en))
                            region_sizes[t_name] -= overlap
                        elif (nextstart - int(en)) <= p_size:
                            region_sizes[t_name] -= (nextstart - int(en))
                for n in sites[i]:
                    n = int(n)
                    if int(z[0]) < n < int(z[1]):
                        #print n,t_name,z[0],z[1]
                        #print pmaa,t_name,change_name,pos,score,fdr,z[0],z[1]
                        if orient == "+":
                            final_dict[pmaa][t_name].append(n-int(st))
                            nuc_cats[t_name].append(i.replace("_A_nidulans_FGSC_A4","")+"_"+str(n))
                        elif orient == "-":
                            final_dict[pmaa][t_name].append(int(en)-n)
                            nuc_cats[t_name].append(i.replace("_A_nidulans_FGSC_A4","")+"_"+str(n))
    nnuc_cats = {}
    for i in nuc_cats:
        nnuc_cats[i] = list(set(nuc_cats[i]))
    nuc_cats = nnuc_cats
    #### make so nucleosomes are counted once and other regions are only counted once
    return final_dict,gene_chrom,sites,gene_inf,nuc_cats,region_sizes,tot_size


dict,chrom_assign,site_pos,gene_desc,nucleosomes,reg_sizes,tot_size = gene_prom_dict(ingff,pks)

outfile = open("/Volumes/MP_HD/Linda_MNase_Seq/python_analysis/R2_4_stats_4_categories/general_nucleosome_position_stats.txt","w")
outfile.write("Location\tnucleosome count\tcumulative region length\tnucleosomes per kb\n")


ncount = 0
size_count = 0
for i in nucleosomes:
    type = i
    count = len(nucleosomes[i])
    print nucleosomes[i][:100]
    region_sum_size = reg_sizes[i]
    nucleosomes_per_kb = str(float(len(nucleosomes[i]))/(float(reg_sizes[i])/1000))
    if type != "prom_extended":
        ncount += count
        size_count += region_sum_size
    outfile.write(type+'\t'+str(count)+'\t'+str(region_sum_size)+'\t'+nucleosomes_per_kb+'\n')

extracellular = tot_nuc_count-ncount
extra_size = tot_size-size_count

outfile.write("extracellular"+'\t'+str(extracellular)+'\t'+str(extra_size)+'\t'+str(float(extracellular)/(float(extra_size)/1000))+'\n')
outfile.write("overall"+'\t'+str(tot_nuc_count)+'\t'+str(tot_size)+'\t'+str(float(tot_nuc_count)/(float(tot_size)/1000))+'\n')

outfile.close()