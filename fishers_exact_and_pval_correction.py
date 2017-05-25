import statsmodels.sandbox.stats.multicomp as pcorrect
import scipy.stats as stats

# [0] = degree i.e. 0.2 -ve corr 20 +ve corr , [1] = p value
pval = stats.fisher_exact([[1,10],[1,30]])
print pval



##lis = [0.001,0.0023,0.234,0.00234,0.002,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02]
##
###[1] = array of corrected p values
##print pcorrect.multipletests(lis,alpha=0.05, method='fdr_bh')[1][5]
