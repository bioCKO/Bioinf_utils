
infile = open('/Volumes/MP HD/PhD/CLC_output/CI-rerun/3840_fixed_snps.fasta', "r")
outfile = open('/Volumes/MP HD/PhD/CLC_output/CI-rerun/3840_fixed_snps_old_contigs.fasta', "w")

conversions = open('/Users/mjohnpayne/Documents/PhD/wt_genome/contig name conversion pm_fix.txt','r')

DS_to_num = {}

for line in conversions:
    columns = line.strip("\n").split("\t")
    DS_to_num[columns[0]] = columns[1]

for DS_acc in DS_to_num:
#    print DS_acc
#    print DS_to_num[DS_acc]
    for line in infile:
        line.replace(DS_acc,DS_to_num[DS_acc])
#        print line
        outfile.writelines(line)


infile.close()
outfile.close()
conversions.close()
