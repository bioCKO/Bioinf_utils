from matplotlib import *
from pylab import *
from numpy import *
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import sys


in_cov = open(sys.argv[1],'r')

contigs = []
coverage = {}
## builds dictionaries for each contig of a list of cov values
for line in in_cov:
    columns = line.strip('\n').split('\t')
    if columns[0] in contigs:
        coverage[columns[0]].append(int(columns[2]))
    else:
        coverage[columns[0]] = [int(columns[2])]
        contigs.append(columns[0])

##builds dictionaries of windows of coverage and another dictionary of contigs for position of those windows

cov_win = {}
pos_win = {}

for contig in coverage:
    window = 1000
    step = 100
    start = 0
    end = 999
    cov_win[contig] = []
    pos_win[contig] = []
    if len(coverage[contig]) < end:
        cov_win[contig] = [sum(coverage[contig])/len(coverage[contig])]
        pos_win[contig] = [len(coverage[contig])/2]
    else:
        while end < len(coverage[contig]):
            win = coverage[contig][start:end]
            win = sum(win)/len(win)
            cov_win[contig].append(win)
            pos = end - 500
            pos_win[contig].append(pos)
#            print contig
#            print cov_win[contig]
#            print pos_win[contig]
            start = start + step
            end = start + window
            
        
        
##
##    
##
##

values = cov_win
position = pos_win

##Only allows contigs > 200

newpos = {}
newval = {}

for contig in position:
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    if x_end > 50000:
        newpos[contig] = position[contig]
        newval[contig] = values[contig]
        
position = newpos
values = newval

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
    avg_depth = sum(values[contig])/len(values[contig])
    if avg_depth == 0:
        avg_depth = 400
    else:
        avg_depth = avg_depth
    y_max = avg_depth*3
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    X = np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(values[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    ax0.plot(X,Y, lw=0.2, antialiased = True, alpha=0.9, color = 'red')
    ylabel(contig, rotation = 'horizontal', fontsize = 12)
    xlim(0,x_end+1)
    ylim(0,y_max)
    xticks([x_end],fontsize = 7), yticks([0,avg_depth],fontsize = 7)
    counter = counter + 1


#plt.show()
plt.savefig(sys.argv[2],dpi=500)

close()
in_cov.close()
