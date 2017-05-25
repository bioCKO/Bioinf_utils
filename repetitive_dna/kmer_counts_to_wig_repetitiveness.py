__author__ = 'mjohnpayne'


from Bio import SeqIO

inkmers = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_30mer_counts.txt","r")

ingenome = SeqIO.parse("/Volumes/MP_HD/repDNA_data/Pm_all_genome_remove_traling_nos.fasta","fasta")

outwig = open("/Volumes/MP_HD/repDNA_data/Pm_all_genome_rep_30mer_counts.wig","w")

"fixedStep chrom=chr3 start=400601 step=10"


kmers = {}
for i in inkmers:
    col = i.strip('\n').split('\t')
    kmers[col[0]] = col[1]

genome = {}

for i in ingenome:
    if len(i.seq) > 100000:
        outwig.write("fixedStep chrom=" + str(i.id) + " start=1 step=1\n")
        for j in range(len(str(i.seq))-30):
            win = str(i.seq)[j:j+30]
            outwig.write(kmers[win] + "\n")
        outwig.write(("1\n")*30)
outwig.close()


