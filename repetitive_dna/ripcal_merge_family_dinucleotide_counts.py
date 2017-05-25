__author__ = 'mjohnpayne'


import glob

#inpf = glob.glob("/Volumes/MP_HD/repDNA_data/repeatmasker/repeats_ripcal_out/pf/*")

folders = glob.glob("/Volumes/MP_HD/repDNA_data/repeatmasker/repeats_ripcal_out/*hs*")

order = ["aa","ac","ag","at","ca","cc","cg","ct","ga","gc","gg","gt","ta","tc","tg","tt"]

def returnlist_of_lists(infolder):
    inpf = glob.glob(infolder + "/*")
    dinucs = []
    used = 0
    for i in inpf:
        tmp = open(i,'r').readlines()[1].split('\t')[1:]
        if 'cannot' not in tmp[0]:
            tmp = map(float,tmp)
            used +=1
            if len(dinucs) != 0:
                dinucs = [x + y for x, y in zip(dinucs, tmp)]
            else:
                dinucs = tmp
    averages = [x/used for x in dinucs]
    return averages

outfile = open("/Volumes/MP_HD/repDNA_data/repeatmasker/repeats_ripcal_out/dinucleotide_averages_4_talaromycetes_hs.txt",'w')
outfile.write("Species" + '\t' + '\t'.join(order) + '\n')
for i in folders:
    if "txt" not in i:
        outfile.write(i[-5:-3] + '\t' + '\t'.join(map(str,returnlist_of_lists(i))) + '\n')

outfile.close()