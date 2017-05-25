__author__ = 'mjohnpayne'


infile = open("/Users/mjohnpayne/Documents/PhD/transposon insertion sites/sites_tab_file_10kbwindow_10kbstep.txt","r")


st = 0
pos = {}
no = 1
cutoff = 70
for i in infile:
    col = i.strip('\n').split('\t')
    if col[0] not in pos:
        pos[col[0]] = {}
    if int(col[3]) > cutoff and st == 0:
        st = 1
        pos[col[0]][col[0]+"_"+str(no)] = [col[1],col[2]]
    elif int(col[3]) > cutoff and st == 1:
        pos[col[0]][col[0]+"_"+str(no)][1] = col[2]
    elif int(col[3]) < cutoff and st == 1:
        st = 0
        no +=1

for i in sorted(pos.keys()):
    if len(pos[i].keys()) > 0:
        print i,pos[i]