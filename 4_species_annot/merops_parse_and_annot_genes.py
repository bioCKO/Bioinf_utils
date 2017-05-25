__author__ = 'mjohnpayne'

import re
import sys
from time import sleep as sl

ingenes = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

family = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/meropsdb/family.txt','r').readlines()
familysum = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/meropsdb/family_summary.txt",'r')
count = 0

outfam = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/meropsdb/family_summary_by_mp.txt",'w')

fam = {}

for i in range(len(family)):
    line = family[i]
    if "\\0family\\" in line:
        id = line[1:line.find('\\')]
        # print len(line)
        # print len(family[i+1])
        # print len(family[i+2])
        if len(family[i+1]) > 4:
            desc = family[i+1][:family[i+1].find('\\0')]
            fam[id] = [re.sub(r'[^a-zA-Z0-9\{\}\ -]','',desc),""]
        elif len(family[i+2]) > 4:
            desc = family[i+2][:family[i+2].find('\\0')]
            fam[id] = [re.sub(r'[^a-zA-Z0-9\{\}\ -]','',desc),""]
        elif len(family[i+3]) > 4:
            desc = family[i+3][:family[i+3].find('\\0')]
            fam[id] = [re.sub(r'[^a-zA-Z0-9\{\}\ -]','',desc),""]
    if "subfamily\\" in line and "external_ids" not in line:
        id = line[1:line.find('\\')]
        # print len(line)
        # print len(family[i+1])
        # print len(family[i+2])
        if len(family[i+1]) > 4:
            desc = family[i+1][:family[i+1].find('\\0')]
            fam[id] = [re.sub(r'[^a-zA-Z0-9\{\}\ -]','',desc),""]
        elif len(family[i+2]) > 4:
            desc = family[i+2][:family[i+2].find('\\0')]
            fam[id] = [re.sub(r'[^a-zA-Z0-9\{\}\ -]','',desc),""]
        elif len(family[i+3]) > 4:
            desc = family[i+3][:family[i+3].find('\\0')]
            fam[id] = [re.sub(r'[^a-zA-Z0-9\{\}\ -]','',desc),""]

for j in fam.keys():
    fam[j].append('No description')

for line in familysum:
    if "activities_and_specificities" in line:
        col = line.replace('"',"").split('\t')
        if col[0] not in fam.keys():
            fam[col[0]] = ["No example",col[2].strip('\n')]
        for j in fam.keys():
            if col[0] == j[:3]:
                fam[j][1] = col[2].strip('\n')


#
for i in fam:
    outfam.write(i + '\t' + '\t'.join(fam[i]) + '\n')
outfam.close

outfile.write("Gene ID\tMerops family name\tType protease\tFamily description\n")
pmaas = []
for line in ingenes:
    if 'query name' not in line:
        col = line.strip('\n').split('\t')
        if col[0] not in pmaas:
            # sl(0.3)
            if col[1] in fam:
                outfile.write(col[0] + '\t' + col[1] + '\t' + fam[col[1]][0] + '\t' + fam[col[1]][1] +'\n')
            elif col[1][:3] in fam:
                outfile.write(col[0] + '\t' + col[1] + '\t' + fam[col[1][:3]][0] + '\t' + fam[col[1][:3]][1] +'\n')
            pmaas.append(col[0])


#308 fams
#253 with biol functions
# extract from family,txt
## subfamily\0 to ident subfamily lines -> start of next line in description
## \0family\0 to ident family lines -> next line beginning with character is description