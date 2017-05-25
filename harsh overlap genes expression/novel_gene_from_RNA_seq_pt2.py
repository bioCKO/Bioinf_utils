import random
import re
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC

######################################

# in terminal: python novel_gene_from_RNA_seq_pt2.py [output.txt from part 1] [pmfa1_annot_scaf.fasta path] [output.fasta]

######################################

indata = open(sys.argv[1],'r')
genome_fasta = open(sys.argv[2],'r')
out_fasta = open(sys.argv[3],'w')
cov = float(sys.argv[4])

##indata = open("/Users/mjohnpayne/Desktop/harsh overlap genes expression/NoPhage_1_expression_genes_25bp_window.txt",'r')
##out_fasta = open("/Users/mjohnpayne/Desktop/harsh overlap genes expression/NoPhage_1_poss_novel_>4_exp_10kb.fasta","w")
genome = {}


for record in SeqIO.parse(genome_fasta, "fasta"):
        genome[record.id] = record.seq

del genome['98']



novel = 0
counter = 0
start = 0
end = 0
gene_counter = 0
end_counter = 12
for line in indata:
    columns = re.sub('\[','',line)
    columns = re.sub("\]",'',columns)
    columns = re.sub("\,",'',columns)
    columns = columns.strip('\n').split(' ')
    contig = columns[0]
    pos = columns[1]
    coverage = float(columns[2])
    gene = int(columns[3])
    if coverage > cov and gene == 0 and novel == 0 and end_counter >= 20:
#        print "start" + " " + contig + " " + str(pos)
        novel = 1
        counter = counter + 1
        start = int(pos)
    elif gene == 1:
        end_counter = 0
        novel = 0
        counter = 0
    elif coverage > cov and gene == 0 and novel == 1:
        counter = counter + 1
    elif coverage < cov and gene == 0 and novel == 1 and counter < 20:
        novel = 0
        counter = 0
        end_counter = end_counter + 1       
    elif coverage < cov and gene == 0 and novel == 1 and counter >= 20:
#        print "start " + str(start) +" end " + contig + " " + str(pos)
        end = int(pos)
        novel_gene = genome[contig][start:end]
        gene_counter = gene_counter + 1
        out_fasta.writelines(">novel_gene_contig_" + str(contig) + "_" + str(start) + "_" + str(end) + "\n" + novel_gene + "\n")
        novel = 0
        counter = 0
        end_counter = 0
    else:
        end_counter = end_counter + 1
        counter = 0

indata.close()
out_fasta.close()
    
        
    
