__author__ = 'mjohnpayne'


ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r').readlines()

outgff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix2.gff3','w')

for i in range(len(ingff)):
    col = ingff[i].strip('\n').split('\t')
    current = '*'
    if col[2] == 'gene':
        current = col[6]



