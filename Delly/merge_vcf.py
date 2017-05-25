__author__ = 'mjohnpayne'

import glob
import sys

invcf = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_*_TRA.vcf")
print invcf
outf = open("/Volumes/MP_HD/CI_GENOME_SEQ/CI_delly/filtered_vcf/filtered_dedup_all_TRA.vcf","w")

def remove_dupes(vcflis):
    outlis = []
    poslis = {}
    for i in vcflis:
        pres = 'N'
        col = i.split('\t')
        cont = col[0]
        pos = int(col[1])
        if cont in poslis:
            for j in poslis[cont]:
                j = float(j)
                if j-30 < pos < j+30:
                    pres = 'Y'
        if pres == 'N':
            if cont not in poslis:
                poslis[cont] = [pos]
            else:
                poslis[cont].append(pos)
            outlis.append(i)
    return outlis

outls = []
for i in invcf:
    print i
    inf = open(i,"r")
    for j in inf:
        if j[0] != "#":
            outls.append(j)
    inf.close()

outls = remove_dupes(outls)

vcf_header = """##fileformat=VCFv4.1
##fileDate=20151118
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="PE support below 3 or mapping quality below 20.">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=CONTROL,Number=1,Type=Integer,Description="Control variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods for RR,RA,AA genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts for the SV">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n"""

outf.write(vcf_header + "".join(outls))
outf.close()
