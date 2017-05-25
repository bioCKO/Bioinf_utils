intag = raw_input('inpath: ')

infile = open(intag,'r')

outtag = intag.strip('.fasta') + '_codons.fasta'

outfile = open(outtag,'w')

count = 0
cod = []
for line in infile:
    st = 0
    en = st + 3
    if line[0] == '>':
        acc = str(line[1:].strip('\n'))
    elif line == '\n':
        continue
    else:
        line = line.strip('\n')
        while en < len(line):
            c = line[st:en]
            cod.append(c)
            st += 3
            en += 3
seq = acc + ',' + ','.join(cod)
outfile.writelines(seq)
print seq
count = count + 1

outfile.close()
