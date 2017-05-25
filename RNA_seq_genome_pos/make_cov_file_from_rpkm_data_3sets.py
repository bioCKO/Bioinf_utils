__author__ = 'mjohnpayne'

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import numpy as np
from matplotlib import *
from pylab import *
from numpy import *

RPKM_index = 51
RPKM_index2 = 63
RPKM_index3 = 57
morm = 'median'

window = 100000

genome = {}
for record in SeqIO.parse("/Users/mjohnpayne/Documents/Phd/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

## get gene positions - gff


rpkm_path = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/RNA_Seq_data_mapped_to_genome/RNA-seq_summary_HW_Jan14.txt'



def make_cov_dict(idx):
    pos = {}
    ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')
    rpkms = open(rpkm_path,'r')
    in_rpkms = rpkms.read().split('\r')
    for line in ingff:
        if '#' not in line:
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                pmaa = col[8].split(';')[1].replace('Name=','')
                if pmaa[-1] == '0':
                    pos[pmaa] = (float(col[3]),float(col[4]),col[0])

    count = 0
    for line in in_rpkms:
        col = line.strip('\n').split('\t')
        if count == 0:
            for i in col: print i + ' ' + str(col.index(i))
            count = 1
        else:
            if col[2] in pos:
                pos[col[2]] += (col[idx],)

    out_dict = {}

    for i in genome:
        out_dict[i] = len(genome[i])*[float(-1)]

    for i in pos:
        contig = pos[i][2]
        if len(pos[i]) == 4:
            st = pos[i][0]
            en = pos[i][1]
            exp = pos[i][3]
            for j in range(int(st),int(en)+1):
                out_dict[contig][j] = float(exp)
    ingff.close()
    rpkms.close()
    return out_dict

## get RPKM values

##generate output with per position RPKM value

def make_cov(coverage):
    cov_win = {}
    pos_win = {}
    for contig in coverage:
        step = window/10
        start = 0
        end = start + window
        cov_win[contig] = []
        pos_win[contig] = []
        if len(coverage[contig]) < end:
            cov_win[contig] = [sum(coverage[contig])/len(coverage[contig])]
            pos_win[contig] = [len(coverage[contig])/2]
        else:
            while end < len(coverage[contig]):
                win = coverage[contig][start:end]
                newin = [0]
                for i in win:
                    if i >= 0:
                        newin += [i]
                if morm == 'mean':
                    win = np.average(newin)
                elif morm == 'median':
                    win = np.median(newin)
                cov_win[contig].append(win)
                pos = end - (window/2)
                pos_win[contig].append(pos)
    #            print contig
    #            print cov_win[contig]
    #            print pos_win[contig]
                start = start + step
                end = start + window
    values = cov_win
    position = pos_win
    ##Only allows contigs > 200
    newpos = {}
    newval = {}
    for contig in position:
        if int(position[contig][-1]) > 200000:
            newpos[contig] = position[contig]
            newval[contig] = values[contig]
    for contig in position:
        avg_depth = sum(values[contig])/len(values[contig])
        if avg_depth == 0:
            avg_depth = 400
        else:
            avg_depth = avg_depth
    return newpos,newval

def getavg(val_list):
    avgd = sum(val_list)/len(val_list)
    return avgd

def normalize_val(vals):
    newvals = {}
    for contig in vals:
        newvals[contig] = []
        avg = getavg(vals[contig])
        for x in vals[contig]:
            y = (float(x)/avg)*100
            newvals[contig].append(y)
    return newvals




cov = make_cov_dict(RPKM_index)
cov2 = make_cov_dict(RPKM_index2)
cov3 = make_cov_dict(RPKM_index3)



in_rpkms = open(rpkm_path,'r').read().split('\r')

name1 = in_rpkms[0].split('\t')[RPKM_index]
name2 = in_rpkms[0].split('\t')[RPKM_index2]
name3 = in_rpkms[0].split('\t')[RPKM_index3]

position,values = make_cov(cov)
pos2,val2 = make_cov(cov2)
pos3,val3 = make_cov(cov3)
#values = normalize_val(values)






# elif len(sys.argv) == 6:
#     in_cov2 = open(sys.argv[3],'r')
#     name2 = sys.argv[3].split('/')[-1][:-8]
#     in_cov3 = open(sys.argv[4],'r')
#     name3 = sys.argv[4].split('/')[-1][:-8]
#     pos2,val2 = make_cov(in_cov2)
#     pos3,val3 = make_cov(in_cov3)
#     val2 = normalize_val(val2)
#     val3 = normalize_val(val3)
# elif len(sys.argv) == 7:
#     in_cov2 = open(sys.argv[3],'r')
#     name2 = sys.argv[3].split('/')[-1][:-8]
#     in_cov3 = open(sys.argv[4],'r')
#     name3 = sys.argv[4].split('/')[-1][:-8]
#     in_cov4 = open(sys.argv[5],'r')
#     name4 = sys.argv[5].split('/')[-1][:-8]
#     pos2,val2 = make_cov(in_cov2)
#     pos3,val3 = make_cov(in_cov3)
#     pos4,val4 = make_cov(in_cov4)
#     val2 = normalize_val(val2)
#     val3 = normalize_val(val3)
#     val4 = normalize_val(val4)



rel_size = {}


### Finds largest contig

largest = ''
size = 0
for contig in position:
    if len(position[contig]) > size:
           size = len(position[contig])
           largest = contig

##finds relative sizes of all contigs

for contig in position:
    rel_size[contig] = float((float(len(position[contig]))/float(len(position[largest]))))

contig_no = (len(position))
print rel_size

##makes plot




fig = plt.figure


counter = 1

for contig in position:
    avg_depth = getavg(values[contig])
    if avg_depth == 0:
        avg_depth = 400
    else:
        avg_depth = avg_depth
    y_min = 0
    y_max = max([max(values[contig]),max(val2[contig]),max(val3[contig])])
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    X = np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(values[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    input1, = ax0.plot(X,Y, lw=1, antialiased = True, alpha=0.5, color = 'red', label = name1)
    #if RPKM_index2 > 0:
    Y2 = numpy.asarray(val2[contig])
    input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
    Y3 = numpy.asarray(val3[contig])
    input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
    # elif len(sys.argv) == 6:
    #     Y2 = numpy.asarray(val2[contig])
    #     input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
    #     Y3 = numpy.asarray(val3[contig])
    #     input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
    # elif len(sys.argv) == 7:
    #     Y2 = numpy.asarray(val2[contig])
    #     input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
    #     Y3 = numpy.asarray(val3[contig])
    #     input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
    #     Y4 = numpy.asarray(val4[contig])
    #     input4, = ax0.plot(X,Y4, lw=1, antialiased = True, alpha=0.5, color = 'orange', label = name4)
    ylabel(contig, rotation = 'horizontal', fontsize = 12)
    xlim(0,x_end+1)
    ylim(y_min,y_max)
    xticks(range(0,x_end,100000),fontsize = 1), yticks([y_min,y_max],fontsize = 7)
    counter = counter + 1

fontP = FontProperties()
fontP.set_size('small')

#if RPKM_index2 > 0:
#     plt.figlegend([input1,input2], (name1,name2), 'upper right',prop=fontP)
# elif len(sys.argv) == 6:
plt.figlegend([input1,input2,input3], (name1,name2,name3), 'upper right',prop=fontP)
# elif len(sys.argv) == 7:
#plt.figlegend([input1,input2,input3,input4], (name1,name2,name3,name4), 'upper right',prop=fontP)
# elif len(sys.argv) == 4:
#plt.figlegend([input1], (name1,), 'upper right',prop=fontP)


#plt.show()
plt.savefig('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/RNA_Seq_data_mapped_to_genome/all_invitro_RPKM_median_100k.pdf',dpi=600)


