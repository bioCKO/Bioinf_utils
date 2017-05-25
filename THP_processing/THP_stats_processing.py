import numpy
from scipy import stats
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import re

# 1. ValidObjectCount = total cells
#
# 2. Percent of cells infected is the inverse of %LOW_ROI_B_Target_I_ObjectCount (or should be, this should be the percent
# of cells with no parasites in them if it has worked correctly).
#
# 3. Pm per infected cell = MEAN_ROI_A_Target_I_ObjectCount x ValidObjectCount / percent of cells infected
#
# 4. Pm average size = MEAN_ROI_B_Target_I_ObjectTotalArea / MEAN_ROI_B_Target_I_ObjectCount
#
# 5. Mitotic index = MEAN_ROI_B_Target_I_ObjectCount / MEAN_ROI_B_Target_II_ObjectCount





indata = open(sys.argv[1],'r').read().split('\r')# open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/140121-HR-OPT-001B_ObjAbove1617forM/140121-HR-OPT-001B_ObjAbove1617forM.txt','r')
# indata2 = []

intype_lis = {}

for i in range(len(indata)):
    if 'Feature: ValidObjectCount' in indata[i]:
        intype_lis['ValidObjectCount'] = indata[i:i+19]
    if '%HIGH_ROI_B_Target_I_ObjectCount' in indata[i]:
        intype_lis['HIGH_ROI_B_Target_I'] = indata[i:i+19]
    if 'MEAN_ROI_B_Target_I_ObjectCount' in indata[i]:
        intype_lis['MEAN_ROI_B_Target_I'] = indata[i:i+19]
    if 'MEAN_ROI_B_Target_I_ObjectTotalArea' in indata[i]:
        intype_lis['MEAN_ROI_B_Target_I_area'] = indata[i:i+19]
    if 'MEAN_ROI_A_Target_II_ObjectCount' in indata[i]:
        intype_lis['MEAN_ROI_A_Target_II'] = indata[i:i+19]
    if '%LOW_ROI_B_Target_I_ObjectCount' in indata[i]:
        intype_lis['LOW_ROI_B_Target_I'] = indata[i:i+19]
    if 'MEAN_ROI_B_Target_II_ObjectCount' in indata[i]:
        intype_lis['MEAN_ROI_B_Target_II'] = indata[i:i+19]
    if 'PlateMap' in indata[i]:
        intype_lis['plateMap'] = indata[i:i+19]



class well:
    def __init__(self,tot_cells,inv_perc_inf,tot_pm_area,pm_count_b,pm_nuc_count,perc_large_pm,ID):
        self.inv_perc_inf = float(inv_perc_inf)
        self.perc_inf = 100 - self.inv_perc_inf
        self.tot_cells = float(tot_cells)
        self.tot_pm_area = float(tot_pm_area)
        self.pm_count_b = float(pm_count_b)
        self.pm_av_size = float(self.tot_pm_area)/self.pm_count_b
        self.pm_nuc_count = float(pm_nuc_count)
        self.mit_index = float(self.pm_count_b)/self.pm_nuc_count
        self.perc_large_pm = float(perc_large_pm)
        self.ID = ID
        self.pm_per_infected = float(self.pm_count_b*self.tot_cells)/self.perc_inf

data = {}
def parse_inp(infile):
    newlines = []
    del(infile[:3])
    for i in infile:
        newlines += [i.split('\t')[1:]]
    return newlines

for i in intype_lis:
    intype_lis[i] = parse_inp(intype_lis[i])



rowcount = 0
colcount = 0
order = []
for row in intype_lis['ValidObjectCount']:
    colcount = 0
    for col in row:
        pos = str(colcount+1) + chr(rowcount + ord('A'))
        if intype_lis['ValidObjectCount'][rowcount][colcount] == "":
            pass#print pos
        else:
            data[pos] = well(float(intype_lis['ValidObjectCount'][rowcount][colcount]),float(intype_lis['LOW_ROI_B_Target_I'][rowcount][colcount]),float(intype_lis['MEAN_ROI_B_Target_I_area'][rowcount][colcount]),float(intype_lis['MEAN_ROI_B_Target_I'][rowcount][colcount]),float(intype_lis['MEAN_ROI_B_Target_II'][rowcount][colcount]),float(intype_lis['HIGH_ROI_B_Target_I'][rowcount][colcount]),intype_lis['plateMap'][rowcount][colcount])
            # data[pos].tot_cells = float(intype_lis['ValidObjectCount'][rowcount][colcount])
            # data[pos].inv_perc_inf = float(intype_lis['LOW_ROI_B_Target_I'][rowcount][colcount])
            # data[pos].tot_pm_area = float(intype_lis['MEAN_ROI_B_Target_I_area'][rowcount][colcount])
            # data[pos].pm_count_b = float(intype_lis['MEAN_ROI_B_Target_I'][rowcount][colcount])
            # data[pos].pm_nuc_count = float(intype_lis['MEAN_ROI_B_Target_II'][rowcount][colcount])
            # data[pos].perc_large_pm = float(intype_lis['HIGH_ROI_B_Target_I'][rowcount][colcount])
            # data[pos].ID = intype_lis['plateMap'][rowcount][colcount]
        order.append(pos)
        colcount +=1
    rowcount += 1


outstats = open(sys.argv[1].replace('.txt','_outstats.txt'),'w') #open('/Users/mjohnpayne/Documents/PhD/THP_well_stats/27-2-14/stats_by_sample_type.txt','w')

outstats.write('Sample ID\tTotal cells Average\tCV\tSEM\tInfected Cells Average\tCV\tSEM\tPm Per infected cell Average\tCV\tSEM\tPm Size Average\tCV\tSEM\tMitotic index Average\tCV\tSEM\tPercent Large Pm cells\tCV\tSEM\tNo of replicates\n')

results = {}
for i in data:
    if data[i].ID not in results:
        results[data[i].ID] = [i]
    else:
        results[data[i].ID].append(i)

totcellsav = ()
totcellssem = ()
totcellscv = ()
infectcellsav = ()
infectcellssem = ()
infectcellscv = ()
pmperinfav = ()
pmperinfsem = ()
pmperinfcv = ()
pmsizeav = ()
pmsizesem = ()
pmsizecv = ()
mit_idxav = ()
mit_idxsem = ()
mit_idxcv = ()
perc_largeav = ()
perc_largesem = ()
perc_largecv = ()
idsav = ()
#out_per_type = open(sys.argv[1].replace('.txt','_per_sample_type.txt'),'w')

inftot = []
pmSize = []
Totcells = []
Infectedcells = []
mit_idx = []
perc_large = []
outpos = []

for pos in results['OTP-NT']:
    inftot.append(data[pos].pm_per_infected)
    pmSize.append(data[pos].pm_av_size)
    Totcells.append(data[pos].tot_cells)
    Infectedcells.append(data[pos].perc_inf)
    mit_idx.append(data[pos].mit_index)
    perc_large.append(data[pos].perc_large_pm)
    outpos.append(pos)

infcontrol = numpy.average(numpy.array(inftot))
pmSizecontrol = numpy.average(numpy.array(pmSize))
Totcellscontrol = numpy.average(numpy.array(Totcells))
Infectedcellscontrol = numpy.average(numpy.array(Infectedcells))
mit_idx_control = numpy.average(numpy.array(mit_idx))
perc_lrg_control = numpy.average(numpy.array(perc_large))




outorder = ['OTP-NT','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','GAPDH','PLK','M','']

for i in outorder:
    print i
    inftot = []
    pmSize = []
    Totcells = []
    Infectedcells = []
    mit_idx = []
    perc_large = []
    outpos = []
    outstats.write(i + '\t')
    for pos in results[i]:
        inftot.append(data[pos].pm_per_infected)
        pmSize.append(data[pos].pm_av_size)
        Totcells.append(data[pos].tot_cells)
        Infectedcells.append(data[pos].perc_inf)
        mit_idx.append(data[pos].mit_index)
        perc_large.append(data[pos].perc_large_pm)
        outpos.append(pos)
    inftotarr = numpy.array(inftot)
    #inftotarr = inftotarr/infcontrol
    pmSizearr = numpy.array(pmSize)
    #pmSizearr = pmSizearr/pmSizecontrol
    Totcellsarr = numpy.array(Totcells)
    #Totcellsarr = Totcellsarr/Totcellscontrol
    Infectedcellsarr = numpy.array(Infectedcells)
    #Infectedcellsarr = Infectedcellsarr/Infectedcellscontrol
    mit_idxarr = numpy.array(mit_idx)
    #mit_idxarr = mit_idxarr/mit_idx_control
    perc_large_arr = numpy.array(perc_large)
    #perc_large_arr = perc_large_arr/perc_lrg_control
    outstats.write(str(numpy.average(Totcellsarr)) + '\t')
    outstats.write(str(numpy.std(Totcellsarr)/numpy.average(Totcellsarr)) + '\t')
    outstats.write(str(stats.sem(Totcellsarr)) + '\t')
    totcellsav += (numpy.average(Totcellsarr),)
    totcellssem += (stats.sem(Totcellsarr),)
    totcellscv += (numpy.std(Totcellsarr)/totcellsav,)

    outstats.write(str(numpy.average(Infectedcellsarr)) + '\t')
    outstats.write(str(numpy.std(Infectedcellsarr)/numpy.average(Infectedcellsarr)) + '\t')
    outstats.write(str(stats.sem(Infectedcellsarr)) + '\t')
    infectcellsav += (numpy.average(Infectedcellsarr),)
    infectcellssem += (stats.sem(Infectedcellsarr),)
    infectcellscv += (numpy.std(Infectedcellsarr)/infectcellsav,)

    outstats.write(str(numpy.average(inftotarr)) + '\t')
    outstats.write(str(numpy.std(inftotarr)) + '\t')
    outstats.write(str(numpy.std(inftotarr)/numpy.average(inftotarr)) + '\t')
    pmperinfav += (numpy.average(inftotarr),)
    pmperinfsem += (stats.sem(inftotarr),)
    pmperinfcv += (numpy.std(inftotarr)/pmperinfav,)

    outstats.write(str(numpy.average(pmSizearr)) + '\t')
    outstats.write(str(numpy.std(pmSizearr)/numpy.average(pmSizearr)) + '\t')
    outstats.write(str(stats.sem(pmSizearr)) + '\t')
    pmsizeav += (numpy.average(pmSizearr),)
    pmsizesem += (stats.sem(pmSizearr),)
    pmsizecv += (numpy.std(pmSizearr)/pmsizeav,)

    outstats.write(str(numpy.average(mit_idxarr)) + '\t')
    outstats.write(str(numpy.std(mit_idxarr)/numpy.average(mit_idxarr)) + '\t')
    outstats.write(str(stats.sem(mit_idxarr)) + '\t')
    mit_idxav += (numpy.average(mit_idxarr),)
    mit_idxsem += (stats.sem(mit_idxarr),)
    mit_idxcv += (numpy.std(mit_idxarr)/mit_idxav,)

    outstats.write(str(numpy.average(perc_large_arr)) + '\t')
    outstats.write(str(numpy.std(perc_large_arr)/numpy.average(perc_large_arr)) + '\t')
    outstats.write(str(stats.sem(perc_large_arr)) + '\t')
    perc_largeav += (numpy.average(perc_large_arr),)
    perc_largesem += (stats.sem(perc_large_arr),)
    perc_largecv += (numpy.std(perc_large_arr)/perc_largeav,)

    outstats.write(str(len(inftot)) + '\n')
    idsav += (i,)
 #   out_per_type.write(i + '\t' + '\t'.join(outpos) + '\tAverage' + '\tSEM' + '\n')
  #  out_per_type.write('Total Cells' + '\t' + '\t'.join(map(str,Totcells)) + '\t' + str(numpy.average(Totcellsarr)) + '\t' + str(stats.sem(Totcellsarr)) + '\n')
   # out_per_type.write('Percent Infected' + '\t' + '\t'.join(map(str,Infectedcells)) + '\t' + str(numpy.average(Infectedcells)) + '\t' + str(stats.sem(Infectedcells)) + '\n')
    #out_per_type.write('Pm per Infected Cell' + '\t' + '\t'.join(map(str,inftot)) + '\t' + str(numpy.average(inftot)) + '\t' + str(stats.sem(inftot)) + '\n')
    #out_per_type.write('Pm Cell size' + '\t' + '\t'.join(map(str,pmSize)) + '\t' + str(numpy.average(pmSize)) + '\t' + str(stats.sem(pmSize)) + '\n\n')



font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 6}
matplotlib.rc('font', **font)



ind = numpy.array(range(1,len(idsav)+1))
ind += 0.3
width = 0.6
fonts = 8
xlim = len(idsav)+1.5
gs = gridspec.GridSpec(2,3,hspace = 0.3, wspace= 0.18)
ax1 = plt.subplot(gs[0,0])
rects1 = ax1.bar(ind, totcellsav, width, color='b', yerr=totcellssem, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('1-Total Cells', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)


ax1 = plt.subplot(gs[0,1])
rects1 = ax1.bar(ind, infectcellsav, width, color='b', yerr=infectcellssem, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('2-Percent Infected Cells', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.xlim = len(idsav)
ax1.set_xlim(0,xlim)
#ax1.set_ylim(0,100)

ax1 = plt.subplot(gs[0,2])
rects1 = ax1.bar(ind, pmperinfav, width, color='b', yerr=pmperinfsem, ecolor = 'k', alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('3-Pm per Infected Cell', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)


ax1 = plt.subplot(gs[1,0])
rects1 = ax1.bar(ind, pmsizeav, width, color='b', yerr=pmsizesem, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('4-Pm cell size', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

ax1 = plt.subplot(gs[1,1])
rects1 = ax1.bar(ind, mit_idxav, width, color='b', yerr=mit_idxsem, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('5-Mitotic Index', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

ax1 = plt.subplot(gs[1,2])
rects1 = ax1.bar(ind, perc_largeav, width, color='b', yerr=perc_largesem, ecolor = 'k',alpha = 0.5, linewidth = 0.5,capsize = 2)
ax1.set_title('6-Percent Large Cells', fontsize = fonts)
ax1.set_xticks(ind+width-0.3)
ax1.set_xticklabels(idsav,rotation=90)
ax1.set_xlim(0,xlim)

plt.savefig(sys.argv[-1].replace(".txt",".pdf"),dpi=500)

#plt.show()

outstats.close()
#out_per_type.close()
