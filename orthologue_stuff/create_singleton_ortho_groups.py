__author__ = 'mjohnpayne'


specs = ["Pma","Tfl",'Pfu','Tst','Afu','Nfi','Acl','Ate','Afl','Aor','Ani','And','Pde','Pro','Pch','Pdi','Ade','Hca','Pbr','Cim','Ure','Teq','Tto']

dataorder = ['Afu','Nfi','Acl','Ate','Afl','Aor','Ani','And','Pde','Pro','Pch','Pdi','Tst','Pfu','Pma','Tfl','Ade','Hca','Pbr','Cim','Ure','Teq','Tto']

nos = [872,653,610,1130,1239,913,1312,1300,1841,2398,1959,850,1694,1489,484,2217,2446,2003,2331,2116,961,671,503]

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/singleton_counts.txt','w')

data = {}

for i in range(len(dataorder)):
    data[dataorder[i]] = nos[i]

outfile.write('ID\t' + '\t'.join(specs) + '\n')

for i in range(len(specs)):
    no = 1
    st = len(specs)*[0]
    st[i] = 1
    while no <= data[specs[i]]:
        line = 'singleton_' + specs[i] + '_' + str(no) + '\t' + '\t'.join(map(str,st)) + '\n'
        outfile.write(line)
        no += 1

outfile.close()