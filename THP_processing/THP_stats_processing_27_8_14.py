import numpy
from scipy import stats
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import re

outliers = ['11G','13D','20F','16G','13H','20N','12I']

neg_cont = ['GAPDH','OTP-NT']

indata = open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-8-14/Plate_140121-HR-OPT-001B.txt','r')#open(sys.argv[1],'r')
indata2 = []
for line in indata:
    line = line.replace('\r\r','\n')
    line = line.replace('\r\n','\n')
    line = line.replace('\r','\n')
    indata2.append(line)
indata2 = indata2[0].split('\n')

##
##for i in indata2:
##    print i

in_nuc = indata2[19:36]#open(sys.argv[1],'r') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/nuclei_count.txt','r')

#print in_nuc
infect_cells = indata2[37:54]#open(sys.argv[2],'r')#open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/infected_cell_count.txt','r')
pm_area = indata2[55:72]#open(sys.argv[3],'r') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/Pm_area.txt','r')
tot_cells = indata2[1:18]#open(sys.argv[4],'r') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/tot_cell_count.txt','r')
platemap = open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/PLate_map.txt','r')#open(sys.argv[2],'r') #

class well:
    def __init__(self, nuc_count,infected,pm_area,tot_cells,ID):
        self.nuc_count = nuc_count
        self.infected = infected
        self.pm_area = pm_area
        self.tot_cells = tot_cells
        self.ID = ID
        self.Pm_per_Infected = (nuc_count*tot_cells)/infected
        self.pmsize = pm_area*nuc_count
        self.pct_inf = (infected/tot_cells)*100
    def __str__(self):
        return "nuc_count\tinfected\tpm_area\ttot_cells\tID\tPm_per_Infected\tpmsize\tpct_inf\n" + str(self.nuc_count) + '\t' + str(self.infected) + '\t' + str(self.pm_area) + '\t' + str(self.tot_cells) + '\t' + str(self.ID) + '\t' + str(self.Pm_per_Infected) + '\t' + str(self.pmsize) + '\t' + str(self.pct_inf)

data = {}
def parse_inp(infile):
    newlines = []
    lines = infile
    del(lines[0])
    for i in lines:
        newlines += [i.split('\t')[1:]]
    return newlines

def parse_plate(infile):
    newlines = []
    lines = infile.read().split('\r')
    del(lines[0])
    for i in lines:
        newlines += [i.split('\t')[1:]]
    return newlines

inf = parse_inp(infect_cells)
pm_ar = parse_inp(pm_area)
totcell = parse_inp(tot_cells)
ids = parse_plate(platemap)
inuc = parse_inp(in_nuc)
rowcount = 0
colcount = 0
order = []
for row in inuc:
    colcount = 0
    for col in row:
        pos = str(colcount+1) + chr(rowcount + ord('A'))
        if inuc[rowcount][colcount] == "":
            pass#print pos
        else:
            data[pos] = well(float(inuc[rowcount][colcount]),float(inf[rowcount][colcount]),float(pm_ar[rowcount][colcount]),float(totcell[rowcount][colcount]),ids[rowcount][colcount])
        order.append(pos)
        colcount +=1
    rowcount += 1

for i in outliers:
   del(data[i])


outstats =open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-8-14/1st_set_with_CV_without_outliers.txt','w')# open(sys.argv[3],'w') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/stats_by_sample_type.txt','w')

outstats.write('Sample ID\tTotal cells average\tCV\tSEM\tInfected Cells average\tCV\tSEM\tPm Per infected cell average\tCV\tSEM\tPm Size average\tCV\tSEM\tNo of replicates\n')

results = {}
for i in data:
    if data[i].ID not in results:
        results[data[i].ID] = [i]
    else:
        results[data[i].ID].append(i)

#cont = data['GAPDH']
#cont2 = data['OTP-NT']

outorder = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','GAPDH','OTP-NT','PLK','M','']



totcellsav = ()
totcellssem = ()
totcellscv = ()
totcellsstd = ()
infectcellsav = ()
infectcellssem = ()
infectcellscv = ()
infectcellsstd = ()
pmperinfav = ()
pmperinfsem = ()
pmperinfcv = ()
pmperinfstd = ()
pmsizeav = ()
pmsizesem = ()
pmsizecv = ()
pmsizestd = ()
idsav = ()
# out_per_type = open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-8-14/1st_set_with_CV.txt','w')#open(sys.argv[3].replace('.txt','_per_sample_type.txt'),'w')



inftot = []
pmSize = []
Totcells = []
Infectedcells = []
outpos = []
for pos in results['GAPDH']:
    inftot.append(data[pos].Pm_per_Infected)
    pmSize.append(data[pos].pmsize)
    Totcells.append(data[pos].tot_cells)
    Infectedcells.append(data[pos].pct_inf)
    outpos.append(pos)
infcontrol = numpy.average(numpy.array(inftot))
pmSizecontrol = numpy.average(numpy.array(pmSize))
Totcellscontrol = numpy.average(numpy.array(Totcells))
Infectedcellscontrol = numpy.average(numpy.array(Infectedcells))
print Infectedcellscontrol


for i in outorder:
    print i
    inftot = []
    pmSize = []
    Totcells = []
    Infectedcells = []
    outpos = []
    outstats.write(i + '\t')
    for pos in results[i]:
        inftot.append(data[pos].Pm_per_Infected)
        pmSize.append(data[pos].pmsize)
        Totcells.append(data[pos].tot_cells)
        Infectedcells.append(data[pos].pct_inf)
        outpos.append(pos)
    inftotarr = numpy.array(inftot)
  #  inftotarr = inftotarr/infcontrol
    pmSizearr = numpy.array(pmSize)
 #   pmSizearr = pmSizearr/pmSizecontrol
    Totcellsarr = numpy.array(Totcells)
 #   Totcellsarr = Totcellsarr/Totcellscontrol
    Infectedcellsarr = numpy.array(Infectedcells)
 #   Infectedcellsarr = Infectedcellsarr/Infectedcellscontrol


    totcellsav += (numpy.average(Totcellsarr),)
    totcellsstd += (numpy.std(Totcellsarr),)
    totcellssem += (stats.sem(Totcellsarr),)
    totcellscv += (numpy.std(Totcellsarr)/totcellsav,)
    outstats.write(str(numpy.average(Totcellsarr)) + '\t')
    outstats.write(str(numpy.std(Totcellsarr)/numpy.average(totcellsav)) + '\t')
    outstats.write(str(stats.sem(Totcellsarr)) + '\t')

    infectcellsav += (numpy.average(Infectedcellsarr),)
    infectcellsstd += (numpy.std(Infectedcellsarr),)
    infectcellssem += (stats.sem(Infectedcellsarr),)
    infectcellscv += (numpy.std(Infectedcellsarr)/infectcellsav,)
    outstats.write(str(numpy.average(Infectedcellsarr)) + '\t')
    outstats.write(str(numpy.std(Infectedcellsarr)/numpy.average(infectcellsav)) + '\t')
    outstats.write(str(stats.sem(Infectedcellsarr)) + '\t')

    pmperinfav += (numpy.average(inftotarr),)
    pmperinfstd += (numpy.std(inftotarr),)
    pmperinfsem += (stats.sem(inftotarr),)
    pmperinfcv += (numpy.std(inftotarr)/pmperinfav,)
    outstats.write(str(numpy.average(inftotarr)) + '\t')
    outstats.write(str(numpy.std(inftotarr)/numpy.average(pmperinfav)) + '\t')
    outstats.write(str(stats.sem(inftotarr)) + '\t')

    pmsizeav += (numpy.average(pmSizearr),)
    pmsizestd += (numpy.std(pmSizearr),)
    pmsizesem += (stats.sem(pmSizearr),)
    pmsizecv += (numpy.std(pmSizearr)/pmsizeav,)
    outstats.write(str(numpy.average(pmSizearr)) + '\t')
    outstats.write(str(numpy.std(pmSizearr)/numpy.average(pmsizeav)) + '\t')
    outstats.write(str(stats.sem(pmSizearr)) + '\t')


    outstats.write(str(len(inftot)) + '\n')
    idsav += (i,)
    # out_per_type.write(i + '\t' + '\t'.join(outpos) + '\taverage' + '\tSEM' + '\n')
    # out_per_type.write('Total Cells' + '\t' + '\t'.join(map(str,Totcells)) + '\t' + str(numpy.average(Totcellsarr)) + '\t' + str(stats.sem(Totcellsarr)) + '\n')
    # out_per_type.write('Percent Infected' + '\t' + '\t'.join(map(str,Infectedcells)) + '\t' + str(numpy.average(Infectedcellsarr)) + '\t' + str(stats.sem(Infectedcellsarr)) + '\n')
    # out_per_type.write('Pm per Infected Cell' + '\t' + '\t'.join(map(str,inftot)) + '\t' + str(numpy.average(inftotarr)) + '\t' + str(stats.sem(inftotarr)) + '\n')
    # out_per_type.write('Pm Cell size' + '\t' + '\t'.join(map(str,pmSize)) + '\t' + str(numpy.average(pmSizearr)) + '\t' + str(stats.sem(pmSizearr)) + '\n\n')



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 7}
matplotlib.rc('font', **font)



ind = numpy.array(range(1,len(idsav)+1))
ind += 0.3
width = 0.6
fonts = 10
xlim = len(idsav)+1.5
gs = gridspec.GridSpec(2,2,hspace = 0.5)
ax1 = plt.subplot(gs[0,0])
rects1 = ax1.bar(ind, totcellsav, width, color='b', yerr=totcellsstd, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Total Cells', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)


ax1 = plt.subplot(gs[0,1])
rects1 = ax1.bar(ind, infectcellsav, width, color='b', yerr=infectcellsstd, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Percent Infected Cells', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.xlim = len(idsav)
ax1.set_xlim(0,xlim)
#ax1.set_ylim(0,100)

ax1 = plt.subplot(gs[1,0])
rects1 = ax1.bar(ind, pmperinfav, width, color='b', yerr=pmperinfstd, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Pm per Infected Cell', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)


ax1 = plt.subplot(gs[1,1])
rects1 = ax1.bar(ind, pmsizeav, width, color='b', yerr=pmsizestd, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Pm cell size', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

plt.savefig('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-8-14/1st_set_with_2std_without_outliers.pdf',dpi=500)

#plt.show()

outstats.close()
#out_per_type.close()
