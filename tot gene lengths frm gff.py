import sys

ingff = open(sys.argv[1],'r')
totlen = 0
count = 0
for line in ingff:
    if "#" not in line:
        col = line.strip('\n').split('\t')
        if col[2] == "gene":
            length = int(col[4])-int(col[3])
            totlen += length
            count += 1
print totlen

print float(totlen)/count
        
    
