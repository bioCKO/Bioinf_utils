import sys
import os.path
import glob

def absoluteFilePaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

class well:
    def __init__(self, cont, st, en):
        self.cont = cont
        self.st = st
        self.en = en
        self.len = int(st)-int(en)
        self.lis = []
    def __str__(self):
       return str(self.cont) + '\t' + str(self.st) + '\t' + str(self.en) + '\t' + ','.join(self.lis)

brk_dir = sys.argv[1]
outgff = open(sys.argv[2],'w')

def del_exists(deldict,DEL):
    answer = 0
    for i in deldict:
        j = deldict[i]
        if j.cont == DEL.cont:
            existr = set(range(j.st,j.en))
            newr = set(range(DEL.st,DEL.en))
            intersect = existr.intersection(newr)
            print len(intersect)
            if len(intersect)> 100:
                return i
                answer = 1
            else :
                return 'No'
                answer = 1
    if answer == 0:
       return 'No'
         
                
                

dels = {}
keylis = []
for i in absoluteFilePaths(brk_dir):
    if '_brk.txt' in i and '/.' not in i:
        strain = i.split('/')[-1][:-8]
        temp = open(i,'r')
        for line in temp:
            if line.startswith('#'):
                continue
            else:
                col = line.strip('\n').split('\t')
                if col[6] == "DEL":
                    st = col[1]
                    en = col[4]
                    uniqdel = str(col[0]) + ':' + str(st) + '-' + str(en)
                    newdel = well(col[0],int(st),int(en))
                    newdel.lis.append(strain)
                    print newdel
                    if del_exists(dels,newdel) == 'No':
                        dels[uniqdel] = newdel
#                        dels[uniqdel].lis.append(strain)
                        keylis.append(uniqdel)
                        print dels[uniqdel]
                    else:
                        existdel = del_exists(dels,newdel)
                        dels[existdel].lis.append(strain)
                        print dels[existdel]
                    print '\n\n'
        temp.close()

outgff.write('##gff-version 3\n')

for i in keylis:
   outgff.write(str(dels[i].cont) + '\tBreakdancer\tdeletion\t' + str(dels[i].st) + '\t' + str(dels[i].en) + '\t.\t*\t.\tID=' + i + ';Name=\'' + ' '.join(dels[i].lis) + '\'\n')
outgff.close()

