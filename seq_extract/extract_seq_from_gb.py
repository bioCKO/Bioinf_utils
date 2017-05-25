
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq



###USAGE python /path/extract_seq_from_gb.py protein (or cds or gene) /path_to/genbankfile.gb /path_to/output.fasta



intype = sys.argv[1]
infile = open(sys.argv[2],'r')
gb_path = sys.argv[3]
in_gb = open(gb_path,'rU')
outfile = open(sys.argv[3],'w')


incolumns = gb_path.split('/')
last = len(incolumns)-1
col_name = incolumns[last]
#name = '_' + col_name.strip('.gb')
name = ''
gene = {}
mrna = {}
cds = {}
protein = {}


#in_gb = open('/Users/mjohnpayne/Desktop/merge_pm_embl/Pm_from_gb.gb', "rU")

sequences = SeqIO.parse(in_gb, "genbank")

pmaa = ''

for record in sequences:
    for i,feature in enumerate(record.features):
        if feature.type == 'gene':
           geneseq = feature.extract(record.seq)
           if 'locus_tag' not in feature.qualifiers:
               pmaa = feature.qualifiers['gene'][0]
           else:
               pmaa = feature.qualifiers['locus_tag'][0]
           gene[pmaa] = geneseq
           
        elif feature.type == 'mRNA':
            mrnaseq = feature.extract(record.seq)
            mrna[pmaa] = mrnaseq

        elif feature.type == 'CDS':
            cdsseq = feature.extract(record.seq)
            cds[pmaa] = cdsseq
            prot = cdsseq.translate()#feature.qualifiers['translation'][0]
            protein[pmaa] = prot

acc_list = []

for acc in infile:
    acc = acc.strip('\n')
    acc_list.append(acc)
    acc_count = acc_list.count(acc)
    if acc_count == 1:
        if intype == 'gene':
            outfile.writelines(">" + acc + name + '\n' + gene[acc] + '\n')
        elif intype == 'mrna':
            outfile.writelines(">" + acc + name + '\n' + mrna[acc] + '\n')
        elif intype == 'cds':
            if 'tRNA' not in acc:
               outfile.writelines(">" + acc + name + '\n' + cds[acc] + '\n')
        elif intype == 'protein':
            if 'tRNA' not in acc:
               outfile.writelines(">" + acc + name + '\n' + protein[acc] + '\n')
    elif acc_count > 1:
        if intype == 'gene':
            outfile.writelines(">" + acc + str(acc_count) + '_' + name + '\n' + gene[acc] + '\n')
        elif intype == 'mrna':
            outfile.writelines(">" + acc + str(acc_count) + '_' + name + '\n' + mrna[acc] + '\n')
        elif intype == 'cds':
            if 'tRNA' not in acc:
                outfile.writelines(">" + acc + str(acc_count) + '_' + name + '\n' + cds[acc] + '\n')
        elif intype == 'protein':
            if 'tRNA' not in acc:
                outfile.writelines(">" + acc + str(acc_count) + '_' + name + '\n' + protein[acc] + '\n')

outfile.close()
in_gb.close()
#infile.close()
