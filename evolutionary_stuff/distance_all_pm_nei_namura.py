__author__ = 'mjohnpayne'


## Modify to include all pm containing groups - perform test between Pm gene and its distance to non pm genes and distance between all non pm genes ( therefore may run multiple tests per group - one per Pm gene)

import re
import math
from scipy import stats
import glob
from time import sleep as sl

pat = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/Distances_namura_nei/distance_out_all_with_pm'

indata = glob.glob(pat + '/*.meg')

#testin = ['/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/ortho_group_cds_fastas/all_with_pm/z-test_pos/z-test_out(6980).meg']

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/Distances_namura_nei/all_with_pm_distances_ttest.txt','w')

outfile.write('orthologue_cluster\tPMAA\tt-test\tratio(+=more in pm)\n')

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

count = 0
c = 0

def check_in_pair(tup,id1,id2):
    if tup[0] == id1:
        return 'corr'
    elif tup[0] == id2:
        return 'opp'

out = []

for j in indata:
    name = j
    if count % 1000 == 0:
        print str(float(count)/len(indata)*100) + '% done'
    count += 1
    array = []
    key = []
    infile = open(j,'r').readlines()
    #print j
    for i in range(len(infile)):
        if '[            ' in infile[i]:
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
        elif 'd : Estimate' in infile[i]:
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
                pval[(key[h],key[i])] = array[i][h]

        yen += 1
        xst += 1
    pairs = pval.keys()
    pmls = []
    npmls = []
    for i in key:
        if 'PMAA' in i:
            pmls += [i]
        else:
            npmls += [i]

    npm_dists = []
    for t in npmls:
        for s in npmls:
            if (t,s) in pairs:
                npm_dists += [pval[(t,s)]]
            elif (s,t) in pairs:
                npm_dists += [pval[(s,t)]]
            else:
                continue
    for i in pmls:
        dists = []
        for j in npmls:
            if (i,j) in pairs:
                dists.append(pval[(i,j)])
            else:
                dists.append(pval[(j,i)])
        if len(dists) > 1 and len(npm_dists) > 1:
            outfile.write(name.split('/')[-1].strip('.meg') + '\t')
            cur_ttest = stats.ttest_ind(map(float,dists),map(float,npm_dists))
            outfile.write(i + '\t' + str(cur_ttest[1]) + '\t' + str(cur_ttest[0]) + '\t' + str(len(pmls)) + '\n')
            # print dists
            # print npm_dists
            # print str(cur_ttest[1]) + '\t' + str(cur_ttest[0])
            # sl(0.2)

    # print len(pairs)
    # print pval
    # print score
    # outfile.write(j.split('/')[-1].strip('.meg') + '\t' + key[1][4:] + '\t' )
    # plist = [('Pm','Tf'),('Pf','Pm'),('Pm','Ts'),('Pf','Tf'),('Tf','Ts'),('Pf','Ts')]
    # for i in pairs:
    #     for j in plist:
    #         if j[0] in i[0] and j[1] in i[1]:
    #             outfile.write(str(pval[i]) + '\t')
    # outfile.write('\n')

outfile.close()