

import glob

outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_gene_coverage (generate stat for sig diff cov)/gene_coverage_doc_median_normalised/gene_coverage_norm_to_median.txt",'w')

incov = glob.glob('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/htseq_gene_cov/*.txt')

ingff = open('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_working_models_fix.gff3','r')

inmedian = open('/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/depthofcoverage_stats/CI_median_readdepths.txt','r')

genlen = {}

for i in ingff:
    if '##' not in i:
        col = i.strip('\n').split('\t')
        if col[2] == 'gene':
            len = int(col[4]) - int(col[3])
            name = col[8].split(';')[1][5:]
            genlen[name] = float(len)


median = {}
for i in inmedian:
    col = i.strip('\n').split('\t')
    median[col[0]] = float(col[1])

## normalisation is ((gene coverage / length) / median depth for strain)*100
normcov = {}
IDS = []

for file in incov:
    temp = open(file,'r').readlines()
    ID = file.split('/')[-1].replace("_GATK_processed_cov.txt","")
    IDS.append(ID)
    for i in temp[:-5]:
        col = i.split('\t')
        val = ((float(col[1])/genlen[col[0]])/median[ID])*100
        if col[0] not in normcov:
            normcov[col[0]] = [val]
        else:
            normcov[col[0]].append(val)

outfile.write("GeneID\t" + '\t'.join(IDS) +'\n')

for i in normcov:
    outfile.write(i + '\t' + '\t'.join(map(str,normcov[i])) + '\n')

outfile.close()