from matplotlib import *
from pylab import *
from numpy import *
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import sys


in_data = open(sys.argv[1],'r')#

contigs = []
values = {}
position = {}

for line in in_data:
    columns = line.strip('\n').split('\t')
    if columns[0] in contigs:
        values[columns[0]] = values[columns[0]] + [columns[3]]
        position[columns[0]] =  position[columns[0]] + [(int(columns[1])+5000)/1000]
    else:
        contigs = contigs + [columns[0]]
        values[columns[0]] = [columns[3]]
        position[columns[0]] = [(int(columns[1])+5000)/1000]

newpos = {}
newval = {}

for contig in position:
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    if x_end > 10:
        newpos[contig] = position[contig]
        newval[contig] = values[contig]
        
position = newpos
values = newval
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



for contig in position:
    x_start = position[contig][0]
    pos_len = len(position[contig]) - 1
    x_end = position[contig][pos_len]
    X = np.linspace(x_start, x_end, pos_len + 1)
    Y = numpy.asarray(values[contig])
    new_right = (10 + (85*rel_size[contig]))/100
    gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
    ax0 = plt.subplot(gs[counter-1])
    ax0.plot(X,Y, lw=1.0, antialiased = True, alpha=0.9, color = 'red')
    ylabel(contig, rotation = 'horizontal', fontsize = 12)
    xlim(0,x_end+1)
    ylim(0.25,0.75)
    xticks([x_end],fontsize = 7), yticks([0.5],fontsize = 2)
    counter = counter + 1


#plt.show()
outpdf = sys.argv[2]
plt.savefig(outpdf,dpi=800)

close()
in_data.close()
