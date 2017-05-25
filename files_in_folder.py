import glob

inp = raw_input('infolder: ')

file_lis = glob.glob(inp + '/*')
file_lis.sort()

a = 0
b = 1

while b < len(file_lis):
    print '1' + file_lis[a]
    print '2' + file_lis[b]
    a+=2
    b+=2
