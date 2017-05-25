import numpy
from scipy import stats
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import re

outliers = []#['3B','5B','3D']

#neg_cont = ['GAPDH','OTP-NT']

indata = open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/HW_test_stats/HarshMutsForMP-B plateNH.csv','r')#open(sys.argv[1],'r')

#data = indata.

def chunks(l, n):
    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]


indata2 = []
for line in indata:
    line = line.replace('\r\r','\n')
    line = line.replace('\r\n','\n')
    line = line.replace('\r','\n')
    indata2.append(line)
indata2 = indata2[0].split('\n')

#indata2 = chunks(indata2,20)

##
##for i in indata2:
##    print i

# ValidObjectCount  =  total macrophage cells
# MEAN_ROI_B_Target_I_ObjectTotalArea  =  Pm Area
# MEAN_ROI_A_Target_II_ObjectCount = ??? (nuclei count?)
# %HIGH_ROI_B_Target_I_ObjectCount  =  Percent infected
# MEAN_ROI_B_Target_I_ObjectCount  = pm count
# Plate Map
# %LOW_ROI_B_Target_I_ObjectCount = ???
# MEAN_ROI_B_Target_I_ObjectTotalInten = ???
tot_cells = []
perc_inf = []
pmcells = []
pm_area = []
in_nuc = []
perc_LOW_ROI_B_Target_I_ObjectCount = []
MEAN_ROI_B_Target_I_ObjectTotalInten = []
platemap = []
for i in range(len(indata2)):
    if 'Feature: ValidObjectCount' in indata2[i]:
        tot_cells = indata2[i:i+19]
    if '%HIGH_ROI_B_Target_I_ObjectCount' in indata2[i]:
        perc_inf = indata2[i:i+19]
    if 'MEAN_ROI_B_Target_I_ObjectCount' in indata2[i]:
        pmcells = indata2[i:i+19]
    if 'MEAN_ROI_B_Target_I_ObjectTotalArea' in indata2[i]:
        pm_area = indata2[i:i+19]
    if 'MEAN_ROI_A_Target_II_ObjectCount' in indata2[i]:
        in_nuc = indata2[i:i+19]
    if '%LOW_ROI_B_Target_I_ObjectCount' in indata2[i]:
        perc_LOW_ROI_B_Target_I_ObjectCount = indata2[i:i+19]
    if 'MEAN_ROI_B_Target_I_ObjectTotalInten' in indata2[i]:
        MEAN_ROI_B_Target_I_ObjectTotalInten = indata2[i:i+19]
    if 'Plate Map' in indata2[i]:
        platemap = indata2[i:i+19]


#
# in_nuc = indata2[3]#open(sys.argv[1],'r') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/nuclei_count.txt','r')
#
# #print in_nuc
# infect_cells = indata2[55:72]#open(sys.argv[2],'r')#open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/infected_cell_count.txt','r')
# pm_area = indata2[1]#open(sys.argv[3],'r') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/Pm_area.txt','r')
# tot_cells = indata2[0]#open(sys.argv[4],'r') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/tot_cell_count.txt','r')
# platemap = indata2[5]#open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/PLate_map.txt','r')#open(sys.argv[2],'r') #

class well:
    def __init__(self, nuc_count,perc_infected,pm_area,tot_cells,pmcells,perc_LOW_ROI_B_Target_I_ObjectCount,MEAN_ROI_B_Target_I_ObjectTotalInten,ID):
        self.MEAN_ROI_B_Target_I_ObjectTotalInten = MEAN_ROI_B_Target_I_ObjectTotalInten
        self.perc_LOW_ROI_B_Target_I_ObjectCount = perc_LOW_ROI_B_Target_I_ObjectCount
        self.nuc_count = nuc_count
        self.pmcells = pmcells
        self.perc_infected = perc_infected
        self.pm_area = pm_area
        self.tot_cells = tot_cells
        self.ID = ID
        self.Pm_per_Infected = pmcells*(float(perc_infected)/100)
        self.pmsize = pm_area/pmcells
#        self.pct_inf = (infected/tot_cells)*100
    def __str__(self):
        return "nuc_count\tpercent_infected\tpm_area\ttot_cells\tID\tPm_per_Infected\tpmsize\n" + str(self.nuc_count) + '\t' + str(self.infected) + '\t' + str(self.pm_area) + '\t' + str(self.tot_cells) + '\t' + str(self.ID) + '\t' + str(self.Pm_per_Infected) + '\t' + str(self.pmsize)

data = {}
def parse_inp(infile):
    newlines = []
    lines = infile
    del(lines[:4])
    for i in lines:
        newlines += [i.split(',')[2:]]
    return newlines

def parse_plate(infile):
    newlines = []
    lines = infile.read().split('\r')
    del(lines[0])
    for i in lines:
        newlines += [i.split('\t')[1:]]
    return newlines

percinf = parse_inp(perc_inf)
print percinf
pm_ar = parse_inp(pm_area)
totcell = parse_inp(tot_cells)
ids = parse_inp(platemap)
inuc = parse_inp(in_nuc)
pmcells = parse_inp(pmcells)
perc_LOW_ROI_B_Target_I_ObjectCount = parse_inp(perc_LOW_ROI_B_Target_I_ObjectCount)
MEAN_ROI_B_Target_I_ObjectTotalInten = parse_inp(MEAN_ROI_B_Target_I_ObjectTotalInten)

rowcount = 0
colcount = 0
order = []
for row in inuc:
    colcount = 0
    for col in row:
        pos = str(colcount+2) + chr(rowcount + ord('A'))
        if inuc[rowcount][colcount] == "":
            pass#print pos
        else:
            data[pos] = well(float(inuc[rowcount][colcount]),float(percinf[rowcount][colcount]),float(pm_ar[rowcount][colcount]),float(totcell[rowcount][colcount]),float(pmcells[rowcount][colcount]),float(perc_LOW_ROI_B_Target_I_ObjectCount[rowcount][colcount]),float(MEAN_ROI_B_Target_I_ObjectTotalInten[rowcount][colcount]),ids[rowcount][colcount])
        order.append(pos)
        colcount +=1
    rowcount += 1
print data

if len(outliers) > 0:
    for i in outliers:
       del(data[i])


outstats =open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/HW_test_stats/HW_plate_BNH.txt','w')# open(sys.argv[3],'w') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/stats_by_sample_type.txt','w')

outstats.write('Sample ID\tTotal cells average\tCV\tSEM\tPercent infected cells\tCV\tSEM\tPm Per infected cell average\tCV\tSEM\tPm Size average\tCV\tSEM\tPm cell count\tCV\tSEM\tperc_LOW_ROI_B_Target_I_ObjectCount\tCV\tSEM\tMEAN_ROI_B_Target_I_ObjectTotalInten\tCV\tSEM\tNo of replicates\n')

results = {}
for i in data:
    if data[i].ID not in results:
        results[data[i].ID] = [i]
    else:
        results[data[i].ID].append(i)

#cont = data['GAPDH']
#cont2 = data['OTP-NT']

outorder = ['1','2','3','4','5','6','7','8','9','10','WT','HK','U']



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

pmcellsav = ()
pmcellssem = ()
pmcellscv  = ()
pmcellsstd = ()

perc_LOW_ROI_B_Target_I_ObjectCountav = ()
perc_LOW_ROI_B_Target_I_ObjectCountsem = ()
perc_LOW_ROI_B_Target_I_ObjectCountcv = ()
perc_LOW_ROI_B_Target_I_ObjectCountstd = ()

MEAN_ROI_B_Target_I_ObjectTotalIntenav = ()
MEAN_ROI_B_Target_I_ObjectTotalIntensem = ()
MEAN_ROI_B_Target_I_ObjectTotalIntencv = ()
MEAN_ROI_B_Target_I_ObjectTotalIntenstd = ()
idsav = ()
# out_per_type = open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-8-14/1st_set_with_CV.txt','w')#open(sys.argv[3].replace('.txt','_per_sample_type.txt'),'w')



# inftot = []
# pmSize = []
# Totcells = []
# Infectedcells = []
# outpos = []
# for pos in results['GAPDH']:
#     inftot.append(data[pos].Pm_per_Infected)
#     pmSize.append(data[pos].pmsize)
#     Totcells.append(data[pos].tot_cells)
#     Infectedcells.append(data[pos].pct_inf)
#     outpos.append(pos)
# infcontrol = numpy.average(numpy.array(inftot))
# pmSizecontrol = numpy.average(numpy.array(pmSize))
# Totcellscontrol = numpy.average(numpy.array(Totcells))
# Infectedcellscontrol = numpy.average(numpy.array(Infectedcells))
# print Infectedcellscontrol


for i in outorder:
    print i
    inftot = []
    pmSize = []
    Totcells = []
    Infectedcells = []
    outpos = []
    perc_LOW_ROI_B_Target_I_ObjectCount = []
    MEAN_ROI_B_Target_I_ObjectTotalInten = []
    pmcells = []
    outstats.write(i + '\t')
    for pos in results[i]:
        inftot.append(data[pos].Pm_per_Infected)
        pmSize.append(data[pos].pmsize)
        Totcells.append(data[pos].tot_cells)
        Infectedcells.append(data[pos].perc_infected)
        perc_LOW_ROI_B_Target_I_ObjectCount.append(data[pos].perc_LOW_ROI_B_Target_I_ObjectCount)
        MEAN_ROI_B_Target_I_ObjectTotalInten.append(data[pos].MEAN_ROI_B_Target_I_ObjectTotalInten)
        pmcells.append(data[pos].pmcells)
        outpos.append(pos)
    inftotarr = numpy.array(inftot)

    pmSizearr = numpy.array(pmSize)

    Totcellsarr = numpy.array(Totcells)

    Infectedcellsarr = numpy.array(Infectedcells)

    perc_LOW_ROI_B_Target_I_ObjectCountarr = numpy.array(perc_LOW_ROI_B_Target_I_ObjectCount)

    MEAN_ROI_B_Target_I_ObjectTotalIntenarr = numpy.array(MEAN_ROI_B_Target_I_ObjectTotalInten)

    pmcellsarr = numpy.array(MEAN_ROI_B_Target_I_ObjectTotalInten)

    totcellsav += (numpy.average(Totcellsarr),)
    totcellsstd += (numpy.std(Totcellsarr),)
    totcellssem += (stats.sem(Totcellsarr),)
    totcellscv += (numpy.std(Totcellsarr)/totcellsav,)
    outstats.write(str(numpy.around(numpy.average(Totcellsarr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(Totcellsarr)/numpy.average(totcellsav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(Totcellsarr),3)) + '\t')

    infectcellsav += (numpy.average(Infectedcellsarr),)
    infectcellsstd += (numpy.std(Infectedcellsarr),)
    infectcellssem += (stats.sem(Infectedcellsarr),)
    infectcellscv += (numpy.std(Infectedcellsarr)/infectcellsav,)
    outstats.write(str(numpy.around(numpy.average(Infectedcellsarr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(Infectedcellsarr)/numpy.average(infectcellsav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(Infectedcellsarr),3)) + '\t')

    pmperinfav += (numpy.average(inftotarr),)
    pmperinfstd += (numpy.std(inftotarr),)
    pmperinfsem += (stats.sem(inftotarr),)
    pmperinfcv += (numpy.std(inftotarr)/pmperinfav,)
    outstats.write(str(numpy.around(numpy.average(inftotarr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(inftotarr)/numpy.average(pmperinfav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(inftotarr),3)) + '\t')

    pmsizeav += (numpy.average(pmSizearr),)
    pmsizestd += (numpy.std(pmSizearr),)
    pmsizesem += (stats.sem(pmSizearr),)
    pmsizecv += (numpy.std(pmSizearr)/pmsizeav,)
    outstats.write(str(numpy.around(numpy.average(pmSizearr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(pmSizearr)/numpy.average(pmsizeav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(pmSizearr),3)) + '\t')

    pmcellsav += (numpy.average(pmcellsarr),)
    pmcellsstd += (numpy.std(pmcellsarr),)
    pmcellssem += (stats.sem(pmcellsarr),)
    pmcellscv += (numpy.std(pmcellsarr)/pmcellsav,)
    outstats.write(str(numpy.around(numpy.average(pmcellsarr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(pmcellsarr)/numpy.average(pmcellsav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(pmcellsarr),3)) + '\t')

    perc_LOW_ROI_B_Target_I_ObjectCountav += (numpy.average(perc_LOW_ROI_B_Target_I_ObjectCountarr),)
    perc_LOW_ROI_B_Target_I_ObjectCountstd += (numpy.std(perc_LOW_ROI_B_Target_I_ObjectCountarr),)
    perc_LOW_ROI_B_Target_I_ObjectCountsem += (stats.sem(perc_LOW_ROI_B_Target_I_ObjectCountarr),)
    perc_LOW_ROI_B_Target_I_ObjectCountcv += (numpy.std(perc_LOW_ROI_B_Target_I_ObjectCountarr)/perc_LOW_ROI_B_Target_I_ObjectCountav,)
    outstats.write(str(numpy.around(numpy.average(perc_LOW_ROI_B_Target_I_ObjectCountarr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(perc_LOW_ROI_B_Target_I_ObjectCountarr)/numpy.average(perc_LOW_ROI_B_Target_I_ObjectCountav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(perc_LOW_ROI_B_Target_I_ObjectCountarr),3)) + '\t')

    MEAN_ROI_B_Target_I_ObjectTotalIntenav += (numpy.average(MEAN_ROI_B_Target_I_ObjectTotalIntenarr),)
    MEAN_ROI_B_Target_I_ObjectTotalIntenstd += (numpy.std(MEAN_ROI_B_Target_I_ObjectTotalIntenarr),)
    MEAN_ROI_B_Target_I_ObjectTotalIntensem += (stats.sem(MEAN_ROI_B_Target_I_ObjectTotalIntenarr),)
    MEAN_ROI_B_Target_I_ObjectTotalIntencv += (numpy.std(MEAN_ROI_B_Target_I_ObjectTotalIntenarr)/MEAN_ROI_B_Target_I_ObjectTotalIntenav,)
    outstats.write(str(numpy.around(numpy.average(MEAN_ROI_B_Target_I_ObjectTotalIntenarr),3)) + '\t')
    outstats.write(str(numpy.around(numpy.std(MEAN_ROI_B_Target_I_ObjectTotalIntenarr)/numpy.average(MEAN_ROI_B_Target_I_ObjectTotalIntenav),3)) + '\t')
    outstats.write(str(numpy.around(stats.sem(MEAN_ROI_B_Target_I_ObjectTotalIntenarr),3)) + '\t')

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
gs = gridspec.GridSpec(2,4,hspace = 0.5)
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

ax1 = plt.subplot(gs[0,2])
rects1 = ax1.bar(ind, pmperinfav, width, color='b', yerr=pmperinfstd, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Pm per Infected Cell', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)


ax1 = plt.subplot(gs[0,3])
rects1 = ax1.bar(ind, pmsizeav, width, color='b', yerr=pmsizestd, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Pm cell size', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

ax1 = plt.subplot(gs[1,0])
rects1 = ax1.bar(ind, pmcellsav, width, color='b', yerr=pmcellsstd, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('Pm cell count', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

ax1 = plt.subplot(gs[1,1])
rects1 = ax1.bar(ind, perc_LOW_ROI_B_Target_I_ObjectCountav, width, color='b', yerr=perc_LOW_ROI_B_Target_I_ObjectCountstd, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('ObjectCount', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

ax1 = plt.subplot(gs[1,3])
rects1 = ax1.bar(ind, MEAN_ROI_B_Target_I_ObjectTotalIntenav, width, color='b', yerr=MEAN_ROI_B_Target_I_ObjectTotalIntenstd, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('ObjectTotalInten', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

plt.savefig('/Users/mjohnpayne/Documents/PhD/THP_well_stats/HW_test_stats/HW_plate_BNH.pdf',dpi=500)

#plt.show()

outstats.close()
#out_per_type.close()
