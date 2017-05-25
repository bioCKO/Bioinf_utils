import re
infile = open('/Users/mjohnpayne/Documents/PhD/random_python_scripts/seq_extract/Pm_from_gb_old_acc.gbk','r')
outfile = open('/Users/mjohnpayne/Documents/PhD/random_python_scripts/seq_extract/Pm_old_acc_reorder.gbk','w')

contigs = {}
c = ''
clist = []
for line in infile:
    if line[:5] == 'LOCUS':
        try:
            c = re.split(' *',line)[1]
            contigs[str(c)] = [line]
            clist.append(int(c))
        except:
            c = re.split(' *',line)[1]
            contigs[c] = [line]
            clist.append(c)           
    else:
        contigs[c].append(line)

clist = map(str,sorted(clist))

for i in clist:
    for a in contigs[i]:
        outfile.writelines(a)
outfile.close()
