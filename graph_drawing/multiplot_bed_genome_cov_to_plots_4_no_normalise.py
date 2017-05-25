from matplotlib import *
from pylab import *
from numpy import *
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import sys
import numpy as np



def make_cov(inc):
    contigs = []
    coverage = {}
    ## builds dictionaries for each contig of a list of cov values
    for line in inc:
        columns = line.strip('\n').split('\t')
        if columns[0] in contigs:
            coverage[columns[0]].append(float(columns[2]))
        else:
            coverage[columns[0]] = [float(columns[2])]
            contigs.append(columns[0])
    ##builds dictionaries of windows of coverage and another dictionary of contigs for position of those windows
    cov_win = {}
    pos_win = {}
    for contig in coverage:
        step = window/10
        start = 0
        end = start + window
        cov_win[contig] = []
        pos_win[contig] = []
        if len(coverage[contig]) < end:
            cov_win[contig] = np.average(coverage[contig])
            pos_win[contig] = [len(coverage[contig])/2]
        else:
            while end < len(coverage[contig]):
                win = coverage[contig][start:end]
                win = sum(win)/len(win)
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
        avg_depth = np.average(values[contig])
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


window = int(sys.argv[1])

in_cov = open(sys.argv[2],'r')
name1 = sys.argv[2].split('/')[-1][:-8]

position,values = make_cov(in_cov)

# values = normalize_val(values)

if len(sys.argv) == 5:
    in_cov2 = open(sys.argv[3],'r')
    name2 = sys.argv[3].split('/')[-1][:-8]
    pos2,val2 = make_cov(in_cov2)
    # val2 = normalize_val(val2)
elif len(sys.argv) == 6:
    in_cov2 = open(sys.argv[3],'r')
    name2 = sys.argv[3].split('/')[-1][:-8]    
    in_cov3 = open(sys.argv[4],'r')
    name3 = sys.argv[4].split('/')[-1][:-8]
    pos2,val2 = make_cov(in_cov2)
    pos3,val3 = make_cov(in_cov3)
    # val2 = normalize_val(val2)
    # val3 = normalize_val(val3)
elif len(sys.argv) == 7:
    in_cov2 = open(sys.argv[3],'r')
    name2 = sys.argv[3].split('/')[-1][:-8] 
    in_cov3 = open(sys.argv[4],'r')
    name3 = sys.argv[4].split('/')[-1][:-8]
    in_cov4 = open(sys.argv[5],'r')
    name4 = sys.argv[5].split('/')[-1][:-8]
    pos2,val2 = make_cov(in_cov2)
    pos3,val3 = make_cov(in_cov3)
    pos4,val4 = make_cov(in_cov4)
    # val2 = normalize_val(val2)
    # val3 = normalize_val(val3)
    # val4 = normalize_val(val4)



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
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    X = np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(values[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    input1, = ax0.plot(X,Y, lw=1, antialiased = True, alpha=0.5, color = 'red', label = name1)
    y_max = max([max([values[contig]])])
    if len(sys.argv) == 5:
        Y2 = numpy.asarray(val2[contig])
        input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
        y_max = max([max(values[contig]),max(val2[contig])])
    elif len(sys.argv) == 6:
        Y2 = numpy.asarray(val2[contig])
        input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
        Y3 = numpy.asarray(val3[contig])
        input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
        y_max = max([max(values[contig]),max(val2[contig]),max(val3[contig])])
    elif len(sys.argv) == 7:
        Y2 = numpy.asarray(val2[contig])
        input2, = ax0.plot(X,Y2, lw=1, antialiased = True, alpha=0.5, color = 'green', label = name2)
        Y3 = numpy.asarray(val3[contig])
        input3, = ax0.plot(X,Y3, lw=1, antialiased = True, alpha=0.5, color = 'blue', label = name3)
        Y4 = numpy.asarray(val4[contig])
        input4, = ax0.plot(X,Y4, lw=1, antialiased = True, alpha=0.5, color = 'orange', label = name4)
        y_max = max([max(values[contig]),max(val2[contig]),max(val3[contig]),max(val4[contig])])
    ylabel(contig, rotation = 'horizontal', fontsize = 12)
    xlim(0,x_end+1)
    ylim(0,min([50,y_max]))
    xticks(range(0,x_end,100000),fontsize = 1), yticks([0,min([50,y_max])],fontsize = 7)
    counter = counter + 1

fontP = FontProperties()
fontP.set_size('small')

if len(sys.argv) == 5:
    plt.figlegend([input1,input2], (name1,name2), 'upper right',prop=fontP)   
elif len(sys.argv) == 6:
    plt.figlegend([input1,input2,input3], (name1,name2,name3), 'upper right',prop=fontP)
elif len(sys.argv) == 7:
    plt.figlegend([input1,input2,input3,input4], (name1,name2,name3,name4), 'upper right',prop=fontP)
elif len(sys.argv) == 4:
    plt.figlegend([input1], (name1), 'upper right',prop=fontP)


#plt.show()
plt.savefig(sys.argv[-1],dpi=300)

close()
in_cov.close()
