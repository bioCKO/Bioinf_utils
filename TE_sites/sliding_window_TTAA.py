import re
from Bio import SeqIO
from Bio.Seq import Seq

outfile = open("/Users/mjohnpayne/Documents/Uni/phd/transposon insertion sites/sites_tab_file_10kbwindow.txt","w")

## Make dictionary of contig numbers and sequences, contig number key

genome = {}

for record in SeqIO.parse("/Users/mjohnpayne/Documents/Uni/phd/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

for contig in genome:
    winstart = 0
    winend = 10000
    step = 1000
    window_no = 0
    if len(genome[contig]) < winend:
        sequence = genome[contig]
        count = sequence.count("TTAA")
        outfile.writelines(str(contig) + "\t" + "1" + "\t" + str(len(genome[contig])) + "\t" + str(count) + "\n")
    else:
        while winend < len(genome[contig]):
            window = genome[contig][winstart:winend]
            wincount = window.count("TTAA")
            outfile.writelines(contig + "\t" + str(winstart) + "\t" + str(winend) + "\t" + str(wincount) + "\n")
            winstart = winstart + step
            winend = winend + step
            window_no = window_no + 1

outfile.close()

