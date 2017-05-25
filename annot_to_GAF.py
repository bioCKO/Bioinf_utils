__author__ = 'mjohnpayne'


infile = open("/Volumes/MP_HD/Blast2Go_data/B2G-Pm/final_complete_short_acc_byseq.txt",'r').readlines()
outfile = open("/Volumes/MP_HD/Blast2Go_data/B2G-Pm/final_complete_short_acc.gaf",'w')

for i in range(1,len(infile)):
    col = infile[i].strip('\n').split('\t')
    outfile.write('B2G\t' + col[0] + '\t' + col[0] + '\t\t' + col[3] + '\t' + col[0] + '\tISS\t\t' + col[2] + '\t\t\tprotein\tTalaromyces marneffei\t20140428\tmpayne\t\t\n')

outfile.close()