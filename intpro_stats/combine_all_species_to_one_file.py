__author__ = 'mjohnpayne'


nlist = ['Pma','Tfl','Pfu','Tst','Pde','Pro','Pch','Pdi','ANID','Ani','Ate','Afl','Aor','Afu','Nfi','Acl','Ade','Hca','Pbr','CIM','Ure','Teq','Tto']#names.keys()
count = 0
species = {}
for j in nlist:
    species[j] = {}
    spec_file = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/interpro_results(results spreadsheet)/all_species_expansion_tests/interpro_domain_expansions_for__' + j
    spec = open(spec_file,'r').readlines()
    for i in range(1,len(spec)):
        col = spec[i].strip().split('\t')
        species[j][col[0]] = [col[24],col[25],col[26].strip('\n')]
print species.keys()

pmfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/interpro_results(results spreadsheet)/all_species_expansion_tests/interpro_domain_expansions_for__Pma','r').readlines()

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/interpro_results(results spreadsheet)/all_species_expansion_tests/all_species_interpro_expansions.txt','w')

outfile.write('\t'.join(pmfile[0].split('\t')[:24]))

for i in nlist:
    outfile.write('\t' + i + '_pvalue\t' + i + '_corrected_pvalue\t' + i + '_up/down_ratio')

outfile.write('\n')

for i in range(1,len(pmfile)):
    col = pmfile[i].split('\t')
    outfile.write('\t'.join(col[:24]) + '\t')
    for j in nlist:
        outfile.write('\t'.join(species[j][col[0]]) + '\t')
    outfile.write('\n')

outfile.close()



