
import re
from Bio import SeqIO
from Bio.Seq import Seq

gff = open("/Users/mjohnpayne/Documents/Uni/phd/wt_genome/pmfa1_working_models.gff3","r")

outfile = open("/Users/mjohnpayne/Documents/Uni/phd/transposon insertion sites/gene_site_counts.txt","w")

#accessions = open("/Users/mjohnpayne/Documents/Uni/phd/Asp_sequences/MEGA_analysis/22_asps_acc_fix.txt","r")

#out_fasta = open("/Users/mjohnpayne/Documents/Uni/phd/Asp_sequences/500bp_5prime_alignment/500bp_5prime_seq_asps.fasta","w")

## Make dictionary of contig numbers and sequences, contig number key

genome = {}

for record in SeqIO.parse("/Users/mjohnpayne/Documents/Uni/phd/wt_genome/pmfa1_annot_scaf.fasta", "fasta"):
        genome[record.id] = record.seq

##Make Dictionaries of all genes with PMAA key and intron
##positions(as list) orientation and contig number objects resectively
accessions = []

cds = {}
gene = {}
orient = {}
contig = {}
cds_counter = 0
ID = ""
for line in gff:
        fin_gene = re.search(".Alias.",line)
        fin_cds = re.search("CDS",line)
        if fin_gene:               
                columns = line.split("\t")
                info = columns[8].split(";")
                ID = info[1].strip("Name=")
                accessions = accessions + [ID]
#                print accessions
                orient[ID] = columns[6]
                contig[ID] = columns[0]
                start = columns[3]
                end = columns[4]
                gene[ID] = [start]
                gene[ID].append(end)
#                print ID
                cds_counter = 0
                
        elif fin_cds:
                cds_counter = cds_counter + 1
                columns = line.split("\t")
                start = columns[3]
                end = columns[4]
                if cds_counter == 1:
                        cds[ID] = [start]
                        cds[ID].append(end)
                elif cds_counter > 1:
                        cds[ID].append(end)


for seq in cds:
        cds[seq] = [cds[seq][0],cds[seq][-1]]
#        print cds[seq]
        
        
gene_TTAA={}

for acc in accessions:
          start = int(gene[acc][0])
          end = int(gene[acc][1])
#          print str(start) + " " + str(end)
          sequence = str(genome[contig[acc]][start:end])
          gene_TTAA[acc] = sequence

for num in gene_TTAA:
        sequence = gene_TTAA[num]
        count = sequence.count("TTAA")
        outfile.writelines(str(num) + "\t" + str(count) + "\t" + str(len(sequence)) + "\n")
 #       print str(gene) + "\t" + str(count) + "\n"
                                  
gff.close()
outfile.close()
#accessions.close()
#out_fasta.close()


