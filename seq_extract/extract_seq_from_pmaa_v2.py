
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
intype = raw_input('Please enter output type(gene, mrna, cds or protein):')
infile = open(raw_input('Path to accession list:').strip(' '),'r')
gb_path = raw_input('Path to genbank file:').strip(' ')
in_gb = open(gb_path,'rU')
outfile = open(raw_input('Path to output file:').strip(' '),'w')

incolumns = gb_path.split('/')
last = len(incolumns)-1
col_name = incolumns[last]
name = col_name.strip('.gb')

gene = {}
mrna = {}
cds = {}
protein = {}



#in_gb = open('/Users/mjohnpayne/Desktop/merge_pm_embl/Pm_from_gb.gb', "rU")

sequences = SeqIO.parse(in_gb, "genbank")

pmaa = ''

for record in sequences:
    contig = record.id
    for i,feature in enumerate(record.features):
        if feature.type == 'gene':
           geneseq = feature.extract(record.seq)
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

def outline(seq_object):
    outfile.writelines(">" + acc + '_' + name + '\n' + seq_object + '\n')

def ouline2(seq_object):
        outfile.writelines(">" + acc + '_' + str(acc_count) + '_' + name + '\n' + seq_object + '\n')

acc_list = []

for acc in infile:
    acc = acc.strip('\n\r')
    acc_list.append(acc)
    acc_count = acc_list.count(acc)
    if acc_count == 1:
        if intype == 'gene':
            outline(gene[acc])
        elif intype == 'mrna':
            outline(mrna[acc])
        elif intype == 'cds':
            outline(cds[acc])
        elif intype == 'protein':
            outline(protein[acc])
    elif acc_count > 1:
        if intype == 'gene':
            outline2(gene[acc])
        elif intype == 'mrna':
            outline2(mrna[acc])
        elif intype == 'cds':
            outline2(cds[acc])
        elif intype == 'protein':
            outline2(protein[acc])       

outfile.close()
in_gb.close()
infile.close()

