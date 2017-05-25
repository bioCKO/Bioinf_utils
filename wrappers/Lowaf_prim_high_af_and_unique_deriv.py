__author__ = 'mjohnpayne'


import sys
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Generate dictionary of SNPs - contig_pos_change:AF
# split each into low and high AF ( < or > 0.95

invp = sys.argv[1]
invcfp = open(invp,'r').readlines()
invd = sys.argv[2]
invcfd = open(invd,'r').readlines()

thresh = 0.60
#outfile = open(sys.argv[2],'w')
pdict_high = {}
pls_high = []
pvcf_high = {}
pdict_low = {}
pls_low = []
pvcf_low = {}

for line in invcfp:
    if '#' not in line:
        col = line.split('\t')
        af = float(col[-1].split(';')[1].replace('AF=',''))
        name = col[0] + "_" + col[1] + '_' + col[3] + col[4]
        if af > 0.95:
            pdict_high[name] = af
            pls_high.append(name)
            pvcf_high[name] = line
        elif af < 0.50:
            pdict_low[name] = af
            pls_low.append(name)
            pvcf_low[name] = line

ddict_high = {}
dls_high = []
dvcf_high = {}
ddict_low = {}
dls_low = []
dvcf_low = {}
for line in invcfd:
    if '#' not in line:
        col = line.split('\t')
        af = float(col[-1].split(';')[1].replace('AF=',''))
        name = col[0] + "_" + col[1] + '_' + col[3] + col[4]
        if af > 0.95:
            ddict_high[name] = af
            dls_high.append(name)
            dvcf_high[name] = line
        elif af < 0.50:
            ddict_low[name] = af
            dls_low.append(name)
            dvcf_low[name] = line


## find SNPs present in deriv high but not in primary high
deriv_high_specific = []
for i in dls_high:
    if i not in pls_high:
        deriv_high_specific.append(i)
pre_existing_snps = []

for i in deriv_high_specific:
    if i in pls_low:
        pre_existing_snps.append(i)

pre_exist = open('/'.join(invd.split('/')[:-1]) + '/' + invd.split('/')[-1].replace('.vcf','_pre_existing.vcf'),'w')

for i in pre_existing_snps:
    pre_exist.write(pvcf_low[i] + '\n')

pre_exist.close()

low_pri_snps = open('/'.join(invd.split('/')[:-1]) + '/' + invd.split('/')[-1].replace('.vcf','_low_in_pri.vcf'),'w')

for i in pdict_low:
    low_pri_snps.write(pvcf_low[i] + '\n')
low_pri_snps.close()

## Gather proportion of these SNPs that are also in primary low AF
overlap = len(pre_existing_snps)
deriv_high_spec = len(deriv_high_specific)
prim_low = len(pls_low)

fig = plt.figure
venn2(subsets = (deriv_high_spec,prim_low,overlap),set_labels = ('Unique to SD - High AF', 'Primary - Low AF'))
# plt.show()
outpdf = sys.argv[3]
plt.savefig(outpdf,dpi=800)

# values = 100*[0]
#
#
# for i in pre_existing_snps:
#     afreq = ddict_high[i]
#      afreq = round(afreq,2)
#     values[int((afreq*100)-1)] += 1