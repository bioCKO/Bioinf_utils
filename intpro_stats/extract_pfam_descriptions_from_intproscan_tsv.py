__author__ = 'mjohnpayne'


infile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/Tfl_f.fasta.tsv','r')
outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/interproscan/Tfl_intpro_assignments.txt','w')


def rn(gene):
    if 'TF' in gene:
        gene = gene[3:]
        no = 5-len(gene)
        gene = 'TFLA_' + no*'0' + gene + '0'
    elif "g" in gene:
        gene = gene[1:]
        no = 5-len(gene)
        gene = 'TFLA_' + no*'0' + gene + '0'
    elif 'Pf' in gene:
        gene = gene[3:]
        no = 5-len(gene)
        gene = 'PFUN_' + no*'0' + gene + '0'
    else:
        gene = gene[4:]
    return gene




ipro = {}
for i in infile:
    col = i.strip('\n').split('\t')
    pma = col[0][4:]
    if '\tIPR' in i:
        ipr = col[11]
        desc = col[12]
        if pma not in ipro:
            ipro[pma] = {ipr:desc}
        else:
            ipro[pma][ipr] = desc

pmas = sorted(ipro.keys())

outfile.write('pmaa\tid\tdescription\n')

for i in pmas:
    outfile.write(rn(i) + '\t' + ','.join([x + ':' + ipro[i][x] for x in ipro[i]]) + '\n')

outfile.close()