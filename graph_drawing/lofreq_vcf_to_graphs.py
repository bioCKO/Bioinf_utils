__author__ = 'mjohnpayne'

from matplotlib import *
from pylab import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import sys

values = 100*[0]

inv = sys.argv[1]
invcf = open(inv,'r')


for line in invcf:
    if '#' not in line and 'INDEL' not in line:
        col = line.strip().split()
        a_freq = float(col[7].split(';')[1].replace('AF=',''))
        a_freq = round(a_freq,2)
        values[int((a_freq*100)-1)] += 1

values = values
maxval = max(values[:-3])
fig = plt.figure

title('Counts of allele frequency in ' + str(inv.split('/')[-1][:-19]),fontsize = 15)
X = np.linspace(0,len(values),len(values))
Y = numpy.asarray(values)
ax0 = plt.plot(X,Y)
ylim(0,maxval)
ylabel('Count of SNPS', rotation = 'vertical', fontsize = 12)
xlabel('Allele Frequency of SNPs', rotation = 'horizontal', fontsize = 12)
#plt.show()

outpdf = inv.replace('.vcf','.pdf')
plt.savefig(outpdf,dpi=800)

# for contig in position:
#     x_start = position[contig][0]
#     pos_len = len(position[contig]) - 1
#     x_end = position[contig][pos_len]
#     X = np.linspace(x_start, x_end, pos_len + 1)
#     Y = numpy.asarray(values[contig])
#     new_right = (10 + (85*rel_size[contig]))/100
#     gs = gridspec.GridSpec(contig_no,1, hspace = 0.5,left = 0.10, right = new_right, bottom = 0.05, top = 0.95)
#     ax0 = plt.subplot(gs[counter-1])
#     ax0.plot(X,Y, lw=1.0, antialiased = True, alpha=0.9, color = 'red')
#     ylabel(contig, rotation = 'horizontal', fontsize = 12)
#     xlim(0,x_end+1)
# #    ylim(0,1.00)
#     xticks([x_end],fontsize = 7), yticks(fontsize = 2)
#     counter = counter + 1
#
#
# #plt.show()