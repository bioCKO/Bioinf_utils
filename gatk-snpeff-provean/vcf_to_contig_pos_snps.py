__author__ = 'mjohnpayne'

inf = "/Volumes/MP_HD/CI_GENOME_SEQ/SNP_compare/012/012_samtools.vcf"
infile = open(inf,'r')
outfile = open(inf[:-4] + "short_id.txt",'w')

for line in infile:
    if line[0] == '#':
        continue
    else:
        col = line.split('\t')
        name = col[0] + '_' + col[1]
        outfile.write(name + '\n')
outfile.close()