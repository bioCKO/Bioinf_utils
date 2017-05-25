__author__ = 'mjohnpayne'





def read_intpro_tsv(path):
    intpro = {}
    inf = open(path,'r')
    for line in inf:
        col = line.strip('\n').split('\t')
        if len(col) > 11:
            ipr = col[11]
            pmaa = col[1]
            if pmaa not in intpro:
                intpro[pmaa] = [ipr]
            else:
                if pmaa not in intpro:
                    intpro[pmaa].append(ipr)