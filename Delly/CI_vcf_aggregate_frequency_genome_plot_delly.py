__author__ = 'mjohnpayne'


import glob
from Bio import SeqIO
from matplotlib import *
from pylab import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties


## snpeff muts with coding consequence
#inls = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_snps/*_GATK_filtered_snps_pass.ann.vcf")
#indels = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_indels/*_GATK_filtered_indels_pass.ann.vcf")

## all muts
#inls = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/Pass_filter_snps/*_GATK_filtered_snps_pass.vcf")
#indels = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/Pass_filter_indels/*_GATK_filtered_indels_pass.vcf")

## from platypus
#inls = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/Pass_filter_snps/*_GATK_filtered_snps_pass.vcf")## snps for compare
#indels = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/CI_platypus_SV/AllVariants.vcf")

##delly
indels = ["/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_dedup_all_DEL.vcf"]
ininvs = ["/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_dedup_all_INV.vcf"]
indupe = ["/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_dedup_all_DUP.vcf"]
intran = ["/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_dedup_all_TRA.vcf"]

centromeres = {"100":2125000,"102":2625000,"103":1535000,"93":1360000,"95":775000,"96":3535000,'99':2865000,'97':2060000}


genome = {}
for record in SeqIO.parse("/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
       genome[record.id] = record.seq
       print record.id + "\t" + str(len(genome[record.id]))


def make_snp_dict(vcfls):
    snps = {"93":[],"95":[],"96":[],"97":[],"99":[],"100":[],"102":[],"103":[]}
    for vcf in vcfls:
        inf = open(vcf,"r")
        for line in inf:
            if '#' not in line:
                col = line.split('\t')
                if col[0] not in snps:
                    snps[col[0]] = [int(col[1])]
                else:
                    snps[col[0]].append(int(col[1]))

        inf.close()
    nsnps = {}
    for i in snps:
        # print i
        # print len(snps[i])
        nsnps[i] = list(set(snps[i]))
        # print len(nsnps[i])
    snps = nsnps


    ### generate list of tuples with contig, winstart, winend, snp_count
    winstart = 0
    winend = 50000
    step = 1000
    hwind = winend/2
    windowcount = 0
    index = 0
    outls = []

    #contlis = ["93","95","96","97","99","100","102","103"]

    for contig in snps:
       winstart = 0
       winend = 5000
       step = 500
       if len(genome[contig]) < winend:
           count = len(snps[contig])
       else:
           while winend < len(genome[contig]):
               while index < len(snps[contig]):
                   if snps[contig][index] > winstart and snps[contig][index] < winend:
                       windowcount = windowcount + 1
                       index = index + 1
                   else:
                       index = index + 1
               outls.append((contig,str(winstart),str(winend),str(windowcount)))
               windowcount = 0
               winstart = winstart + step
               winend = winend + step
               index = 0



    ######



    ## generate corresponding dictionaries of window positions and window values

    contigs = []#"93","95","96","97","99","100","102","103"]
    values = {}#{"93":[],"95":[],"96":[],"97":[],"99":[],"100":[],"102":[],"103":[]}
    position = {}#{"93":[],"95":[],"96":[],"97":[],"99":[],"100":[],"102":[],"103":[]}

    for columns in outls:
        if columns[0] in contigs:
            values[columns[0]] = values[columns[0]] + [columns[3]]
            position[columns[0]] =  position[columns[0]] + [(int(columns[1])+hwind)/1000]
        else:
            contigs = contigs + [columns[0]]
            values[columns[0]] = [columns[3]]
            position[columns[0]] = [(int(columns[1])+hwind)/1000]

    # newpos = {}
    # newval = {}
    #
    # for contig in position:
    #     if len(genome[contig]) > 200000:
    #         x_start = position[contig][0]
    #         pos_len = len(position[contig]) - 1
    #         x_end = position[contig][pos_len]
    #         if x_end > 10:
    #             newpos[contig] = position[contig]
    #             newval[contig] = values[contig]
    # for contig in position:
    #     if int(position[contig][-1]) > 200000:
    #         newpos[contig] = position[contig]
    #         newval[contig] = values[contig]

    # position = newpos
    # values = newval
    return position,values

position,values = make_snp_dict(indels)
invposition,invvalues = make_snp_dict(ininvs)
dupposition,dupvalues = make_snp_dict(indupe)
for i in dupvalues:
    print i, dupvalues[i]
for i in dupposition:
    print i, dupposition[i]
traposition,travalues = make_snp_dict(intran)


rel_size = {}

largest = ''
size = 0
for contig in position:
    if len(position[contig]) > size:
           size = len(position[contig])
           largest = contig




for contig in position:
    rel_size[contig] = float((float(len(position[contig]))/float(len(position[largest]))))

contig_no = (len(position))



fig = plt.figure


counter = 1

clis = ["93","95","96","97","99","100","102","103"]

for contig in clis:
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    X = np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(values[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    ax1 = ax0.twinx()
    ax0.plot(X,Y, lw=1.0, antialiased = True, alpha=0.5, color = 'red', label = "Dels")
    #ax0.plot([centromeres[contig], centromeres[contig]], [0, 500], 'k-', lw=2)
    ax0.axvline(float(centromeres[contig])/1000, color='black',lw=1)
    ax0.set_ylim([0,6])
#    ax0.set_yticks([250,500,750,1000])
    Y2 = numpy.asarray(invvalues[contig])
    ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = "Inversions",)
    Y3 = numpy.asarray(dupvalues[contig])
    ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = "Duplication",)
    Y4 = numpy.asarray(travalues[contig])
    ax0.plot(X,Y4, lw=1, antialiased = True, alpha=0.5, color = 'orange', label = "Transpositions",)
    # ax1.set_ylim([0,120])
    # ax1.set_yticks([40,80,120])
    xlim(0,x_end+1)
    # ylim(0,50)
    #y_max = max([max(values[contig]),max(idvalues[contig])])
    for label in ax0.yaxis.get_majorticklabels():
        label.set_fontsize(6)
    for label in ax1.yaxis.get_majorticklabels():
        label.set_fontsize(6)
    ax0.set_ylabel(contig, rotation = 'horizontal', fontsize = 12)
    ax0.set_xticks([x_end])
    for label in ax0.xaxis.get_majorticklabels():
        label.set_fontsize(6)
    counter = counter + 1
#
# fontP = FontProperties()
# fontP.set_size('small')

#plt.figlegend([input1,input2], ("SNPs","Indels"), 'upper right',prop=fontP)
plt.show()
#outpdf = "/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/platypus_test_indel_snp_plot.pdf"
#plt.savefig(outpdf,dpi=400)
