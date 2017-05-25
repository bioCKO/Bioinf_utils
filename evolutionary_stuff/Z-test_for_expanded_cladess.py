__author__ = 'mjohnpayne'

import re
import math
import numpy
import glob
import time

pat = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/Distances_namura_nei/distance_complete_del_namura_nei_1-1-1-1'

indata = glob.glob(pat + '/*.meg')

#testin = ['/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/ortho_group_cds_fastas/all_with_pm/z-test_pos/z-test_out(6980).meg']

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/ortho_group_cds_fastas/all_with_pm/Z-test-stats.txt','w')

outfile.write('orthologue_cluster\tPMAA\tcomparator\tselection_statistic\tp_value\tdescription\n')

pfpml = []
pftfl = []
pmtfl = []
pftsl = []
pmtsl = []
tftsl = []

name = ['pfpm','pftf','pmtf','pfts','pmts','tfts']
#outfile.write('ID\t' + '\t'.join(name)+'\n')
ln = 0
array = []
key = []
count = 0
c = 0
for j in indata:
    array = []
    key = []
    infile = open(j,'r').readlines()
    #print j
    for i in range(len(infile)):
        if '[               ' in infile[i]:
            x = infile[i]
            for c in ['[',']','\r\n']: x = x.replace(c, "")
            #x = map(int,re.split(' *',x)[1:-1])
            #ln = x[-1]
            table = infile[i+1:i+ln+1]
            for k in table:
                k = k.replace('-214748364','1')
                k = str(k[4:])
                k = k.replace('               ',' NA ')
                for c in ['[',']','\r\n']: k = k.replace(c, "")
                k = re.split(' *',k)#[k[y:y+15].replace(' ','') for y in range(0,len(k)-1,15)]
                array.append(k[1:-1])
        elif 'No. of Taxa' in infile[i]:
            ln = int(infile[i].strip('\n\r')[-3:].replace(' ','').replace(':',''))
        elif 'Stat (above diagonal)' in infile[i]:
            key = infile[i+3:i+ln+3]
    key = [t[t.find('#')+1:].strip('\r\n') for t in key]
    xst = 0
    yen = 0
    score = {}
    pval = {}
    for i in range(xst,ln):
        for h in range(yen):
            if (h,i) == (ln-1,ln-1):
                break
            #print str(j) + ',' + str(i)
            if h-i != 0:
#                print (j,i+1)
                score[(key[h],key[i])] = array[h][i]
                pval[(key[h],key[i])] = array[i][h]
#                print (key[j],key[i+1])
        yen += 1
        xst += 1
    pairs = score.keys()
    # print len(pairs)
    # print pval
    # print score
    for i in pairs:
        # print i
        # print pval[i]
        # print score[i]
        # time.sleep(0.5)
        if 'PMAA' in i[0]:
            outfile.write(j.split('/')[-1].strip('.meg') + '\t' + i[0][4:] + '\t' + i[1][4:] + '\t' + str(score[i]) + '\t' + str(pval[i]) + '\n')
        elif 'PMAA' in i[1]:
            outfile.write(j.split('/')[-1].strip('.meg') + '\t' + i[1][4:] + '\t' + i[0][4:] + '\t' + str(score[i]) + '\t' + str(pval[i]) + '\n')
    # elif float(pval[i]) < 0.05 and 'PMAA' in i[1]:
    #     print i
    #     print pval[i]
    #     print score[i]
    #     time.sleep(0.5)

#outfile.close()