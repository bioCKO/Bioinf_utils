__author__ = 'mjohnpayne'

from time import sleep as sl

outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/Tm_CI_orthogroups_counts.txt","w")

CI_List = ["2161","Pm1","3482","3840","3841","3871","4059","HR2","BR2","BR21","BR22","027","0271","0272","043","012","203","2031","2033","2034","702","F4"]

outfile.write("Group_name\t" + "\t".join(CI_List) + "\n")

#print len(CI_List)

for line in open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/Tm_CI_orthogroups.txt","r"):
    group = line[:line.find(":")]
    outfile.write(group)
    for i in CI_List:
        count = line.count(i+"|")
        outfile.write("\t" + str(count))
    outfile.write('\n')

outfile.close()