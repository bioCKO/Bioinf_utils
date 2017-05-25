__author__ = 'mjohnpayne'

import re
from time import sleep as sl

#inf = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_denovo.gff'
inf = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/T_flavus/vel_denovo_genome/TF_vel_pfams_para_genome_fix.gff"

#of = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_denovo_rename.gff'
of = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/T_flavus/vel_denovo_genome/TF_vel_pfams_para_genome_fix_rename.gff"

outfile = open(of,'w')
infile = open(inf,'r').read()

starts = [m.start() for m in re.finditer("# start gene ", infile)]
ends = [m.end() for m in re.finditer("# end gene ", infile)]

genes = []
for i in range(len(starts)):
    gene = infile[starts[i]:ends[i]]+'\n'
    num = gene.find('ID=g')
    num = gene[num+3:gene.find(';',num)]
    no = 6-len(num)
    id = 'TFLA_' + no*'0' + num[1:] + '0'
    gene = gene.replace(num,id)
    outfile.write(gene)

outfile.close()
