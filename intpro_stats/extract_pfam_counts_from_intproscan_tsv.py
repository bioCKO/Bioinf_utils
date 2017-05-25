import sys
import statsmodels.sandbox.stats.multicomp as pcorrect
import scipy.stats as stats
import numpy



##adapt to extract pfams/panther/SMART/superfamilies/etc

## returns dictionary with interpo id as key and list as value in [id count,total ids count] format as well as total counts as integer
def read_to_fisher(path):
    intpro = {}
    inf = open(path,'r')
    for line in inf:
        col = line.strip('\n').split('\t')
        if col[3] == 'Pfam':
            ipr = col[4]
            pmaa = col[1]
            if ipr not in intpro:
                intpro[ipr] = [pmaa]
            else:
                if pmaa not in intpro[ipr]:
                    intpro[ipr].append(pmaa)
    inf.close()
    intproc = {}
    total = 0
    for i in intpro:
        j = len(intpro[i])
        intproc[i] = j
        total += j
    intprof = {}
    for i in intproc:
        j = [intproc[i],total - intproc[i]]
        intprof[i] = j
    return intprof, total

## iterates over input files running readtofisher and makeing dictionary of dictionaries
## also generated unique id list
## also generates dictionary of total counts

ids = []
names = {}
tots = {}

for i in range(1,len(sys.argv)-1):
    infile = sys.argv[i]
    spec = infile.split('/')[-1]
    spec = spec[:spec.index('_f')]
    intdict,tot = read_to_fisher(infile)
    names[spec] = intdict
    tots[spec] = tot
    for i in intdict:
        if i not in ids:
            ids.append(i)

outfile = open(sys.argv[-1],'w')

## adds entries to species dictionaries that have 0 counts of an interpro id where there are counts in other species

for j in names:
    for i in ids:
        if i not in names[j]:
            names[j][i] = [0,tots[j]]

pm_dict = names['Pma']

del names['Pma']

## returns sum of all species for each id and totals for each id to run in tests against pm only

def return_sums(intid,dic):
    countid = 0
    counttot = 0
    for i in dic:
        countid += dic[i][intid][0]
        counttot += dic[i][intid][1]
    return [countid,counttot]


## runs tests and produces list of pvals corresponding to ids list
## cor_pval is a list of  pvals corresponding to ids list corrected for multiple tests

pval = []
tests = {}

for i in ids:
    test = [pm_dict[i],return_sums(i,names)]
    tests[i] = [stats.fisher_exact(test)[1],stats.fisher_exact(test)[0]]
    pval.append(tests[i][0])
    
cor_pval = pcorrect.multipletests(pval,alpha=0.05, method='fdr_bh')[1]

species = []

for i in range(1,len(sys.argv)-1):
    species.append(sys.argv[i].split('/')[-1][:3])

species.remove('Pma')



outfile.write('InterPro ID\tPma\t' + '\t'.join(species) + '\tpvalue\tcorrected P value\tup/down ratio\n')

st = 0
for i in ids:
    outfile.write(i + '\t' + str(pm_dict[i][0]) + '\t')
    for j in species:
        outfile.write(str(names[j][i][0]) + '\t')
    outfile.write(str(pval[st]) + '\t' + str(cor_pval[st]) + '\t' + str(tests[i][1]) + '\n')
    st += 1

outfile.close()


##infile = sys.argv[1]
##infile2= sys.argv[2]
##infile3 = sys.argv[3]
##infile4 = sys.argv[4]
##
##outfile = open(sys.argv[-1],'w')
##
##intdict = read_to_fisher(infile)
##intdict2 = read_to_fisher(infile2)
##intdict3 = read_to_fisher(infile3)
##intdict4 = read_to_fisher(infile4)
##print intdict['IPR001969']
##print intdict2['IPR001969']
##print intdict3['IPR001969']
##print intdict4['IPR001969']
##sum_others = {}
##
##ids = []
##
##for i in intdict:
##    if i not in ids:
##        ids.append(i)
##
##for i in intdict2:
##    if i not in ids:
##        ids.append(i)
##
##for i in intdict3:
##    if i not in ids:
##        ids.append(i)
##
##for i in intdict4:
##    if i not in ids:
##        ids.append(i)
##
##pval = []
##tests = {}
##for i in ids:
##    if i not in intdict:
##        intdict[i] = [0,23393]
##    if i not in intdict2:
##        intdict2[i] = [0,34370]
##    if i not in intdict3:
##        intdict3[i] = [0,38498]       
##    if i not in intdict4:
##        intdict4[i] = [0,35409]
##for i in ids:
##    test = [intdict[i],[intdict2[i][0]+intdict3[i][0]+intdict4[i][0],intdict2[i][1]+intdict3[i][1]+intdict4[i][1]]]
##    tests[i] = [stats.fisher_exact(test)[1],stats.fisher_exact(test)[0]]
##    pval.append(tests[i][0])
##
##cor_pval = pcorrect.multipletests(pval,alpha=0.05, method='fdr_bh')[1]
##
##outfile.write('InterPro ID\tPm_No\tTs_No\tTf_no\tPf_No\tnon_Pm_No\tpvalue\tcorrected P value\tup/down ratio\n')
##st = 0
##for i in ids:
##    outfile.write(i + '\t' + str(intdict[i][0]) + '\t' + str(intdict4[i][0]) + '\t' + str(intdict3[i][0]) + '\t' + str(intdict2[i][0]) + '\t\t' + str(pval[st]) + '\t' + str(cor_pval[st]) + '\t' + str(tests[i][1]) + '\n')
##    st += 1
##    
##outfile.close()
##    
##
##### [1] = p value
####pval = stats.fisher_exact([[8,2],[1,5]])[1]
####print pval
##
####lis = [0.001,0.0023,0.234,0.00234,0.002,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02]
####
#####[1] = array of corrected p values
####print pcorrect.multipletests(lis,alpha=0.05, method='fdr_bh')[1][5]
