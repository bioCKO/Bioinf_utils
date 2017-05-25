__author__ = 'mjohnpayne'


infile = open("/Users/mjohnpayne/Documents/PhD/lab_dbs/oligodb_29_2_16",'r').read()
outfile = open("/Users/mjohnpayne/Documents/PhD/lab_dbs/oligodb_29_2_16_sc_in.txt",'w')
lines = infile.split('\r')
print len(lines)

def checker(str):
    sum = 0
    check = ['A','a','C','c','G','g','T','t']
    for i in check:
        sum += str.count(i)
    if sum == len(str):
        return 'Y'
    else:
        return 'N'

for i in lines:
    #i = i.replace('\x0b','')
    i=i.replace("\"","")
    col = i.split(',')
    if len(col) > 3:
        if len(col[3]) > 2 and checker(col[6]) == "Y":
            outfile.write(col[1] + col[0] + ',' + col[6] + ',,' + 'primer_bind,DNA,Red\n')

outfile.close()