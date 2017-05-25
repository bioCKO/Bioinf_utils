__author__ = 'mjohnpayne'

from Bio import SeqIO
import glob

ingenomes = glob.glob("/Volumes/MP_HD/repDNA_data/*all_genome.fasta")

outfile = open("/Volumes/MP_HD/repDNA_data/repeatmasker/repeats_ripcal_out/dinucleotide_averages_4_talaromycetes_genome.txt",'w')

order = ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]

def dinucs(genome):
    freqs = {}
    tmp = SeqIO.parse(genome,"fasta")
    tot = 0
    for i in tmp:
        print i.id
        for j in range(len(i.seq)-1):
            dn = str(i.seq[j:j+2])
            if dn not in freqs:
                freqs[dn] = 1
            else:
                freqs[dn] += 1
        tot += len(i.seq)
    new = {}
    for i in freqs:
        new[i] = float(freqs[i])/tot
    return new


outfile.write("Species\t" + "\t".join(order) + '\n')

for i in ingenomes:
    ids = i.split('/')[-1][:2]
    print ids
    tmp = dinucs(i)
    outfile.write(ids)
    for j in order:
        outfile.write('\t' + str(tmp[j]))
    outfile.write('\n')

outfile.close()