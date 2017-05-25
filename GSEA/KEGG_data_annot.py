__author__ = 'mjohnpayne'

import time
import itertools

inkegg = open('/Volumes/MP_HD/GSEA/kegg stuff/ko00001.txt','r')

pm_kegg = open('/Volumes/MP_HD/GSEA/kegg stuff/Pm_kegg_orthology.txt','r')

gmt_out = open('/Volumes/MP_HD/GSEA/kegg stuff/kegg_pm_assoc.gmt','w')

pm_dict = {}
for i in pm_kegg:
    col = i.strip('\n').split('\t')
    if len(col) > 1:
        if col[1] not in pm_dict:
            pm_dict[col[1]] = [col[0]]
        else:
            pm_dict[col[1]].append(col[0])

a = ''
b = ''
c = ''
d = ''
k_ids = {}
a_k_ids = {}
b_k_ids = {}
c_k_ids = {}

cinf = {}
for i in inkegg:
    if i[0] == 'A':
        a = i[i.find('>')+1:i.find('<',4)]
       # print a
    elif i[:3] == 'B  ':
        b = i[i.find('>')+1:i.find('<',4)]
       # print b
    elif i[:5] == 'C    ':
        c = ('ko' + i[5:10],i[11:].strip('\n'))
        cinf[c[0]] = c[1]
       # print c
    elif i[0] == 'D':
        d = (i[7:13],i[15:i.find('[')-1],i[i.find('[')+1:i.find(']')].replace('EC:','').split(' '))
        if d[0] in pm_dict:
            if a not in a_k_ids:
                a_k_ids[a] = [pm_dict[d[0]]]
            else:
                a_k_ids[a].append(pm_dict[d[0]])
            if b not in b_k_ids:
                b_k_ids[b] = [pm_dict[d[0]]]
            else:
                b_k_ids[b].append(pm_dict[d[0]])
            if c[0] not in c_k_ids:
                c_k_ids[c[0]] = [pm_dict[d[0]]]
            else:
                c_k_ids[c[0]].append(pm_dict[d[0]])
        # d = (i[7:13],i[15:i.find('[')-1],i[i.find('[')+1:i.find(']')].replace('EC:','').split(' '))
        # if d[0] not in k_ids:
        #     k_ids[d[0]] = [(a,b,c,d[1:])]
        #     a_k_ids[d[0]] = [a]
        #     b_k_ids[d[0]] = [b]
        #     c_k_ids[d[0]] = [c[0]]
        # else:
        #     k_ids[d[0]].append((a,b,c,d[1:]))
        #     if a not in a_k_ids[d[0]]:
        #         a_k_ids[d[0]].append(a)
        #     if b not in b_k_ids[d[0]]:
        #         b_k_ids[d[0]].append(b)
        #     if c not in c_k_ids[d[0]]:
        #         c_k_ids[d[0]].append(c[0])



def flatten(dict):
    od = {}
    for i in dict:
        out = list(set(itertools.chain.from_iterable(dict[i])))
        od[i] = out
    return od

# print k_ids['K00844']
# print a_k_ids['K00844']
# print b_k_ids['K00844']
# print c_k_ids['K00844']
# out_dict = {}
# for i in c_k_ids:
#     for j in c_k_ids[i]:
#         if j in pm_dict:
#             if i not in out_dict:
#                 out_dict[i] = []
#                 out_dict[i] = pm_dict[j]
#             else:
#                 out_dict[i] += pm_dict[j]



a_k_ids = flatten(a_k_ids)
b_k_ids = flatten(b_k_ids)
c_k_ids = flatten(c_k_ids)


for i in a_k_ids:
    # print i
    # print len(a_k_ids[i])
    gmt_out.write(i + '\ta:' + i + '\t' + '\t'.join(a_k_ids[i]) + '\n')
for i in b_k_ids:
    # print i
    # print len(b_k_ids[i])
    gmt_out.write(i + '\tb:' + i + '\t' + '\t'.join(b_k_ids[i]) + '\n')
for i in c_k_ids:
    # print i
    # print len(c_k_ids[i])
    gmt_out.write(i + '\tc:' + cinf[i] + '\t' + '\t'.join(c_k_ids[i]) + '\n')

gmt_out.close()