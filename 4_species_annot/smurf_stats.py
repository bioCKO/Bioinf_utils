__author__ = 'mjohnpayne'

import glob


inls = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/smurf/*_Secondary-Metabolite-Clusters.txt")

inbb = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/smurf/*_Backbone-genes.txt")

clusters = {}
bb = {}
for i in inbb:
    spc = i.split('/')[-1][:3]
    #print spc
    inf = open(i,"r")
    bb[spc] = {}
    for j in inf:
        col = j.strip().split()
        if len(col) > 0 and col[0] != "Backbone_gene_id":
            bb[spc][col[0]] = col[5]


for i in inls:
    name = ""
    size = 0
    spec = i.split('/')[-1][:3]
    clusters[spec] = []
    inf = open(i,"r")
    tot = 0
    no = 0
    #print spec
    for j in inf:
        if j.startswith("Cluster"):
            if size > 0:
                #print name, size-2
                tot += size-2
                no +=1
            name = j.strip('\n').replace(":","_")
            size = 0
        elif j[0] != "B":
            size +=1
            col = j.strip().split()
            if len(col) > 0 and col[2] == "0":
                clusters[spec].append(col[0])
    tot+=size-2
    no +=1
    print name, size-2
    print spec,"average",(float(tot)/no)



not_clust = {}
for t in bb:
    not_clust[t] = []
    for i in bb[t]:
        if i not in clusters[t]:
            not_clust[t].append((i,bb[t][i]))


for i in not_clust:
    print i, len(clusters[i]),len(not_clust[i]),"\n","\n".join([x[0]+"\t"+x[1] for x in not_clust[i]])

