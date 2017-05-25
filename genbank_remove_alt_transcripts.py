__author__ = 'mjohnpayne'


from Bio import SeqIO
from Bio.Seq import Seq
from time import sleep as sl

ingb = "/Users/mjohnpayne/Documents/PhD/wt_genome/Pm_genbank_gene_containing_only.gb"

outgb = "/Users/mjohnpayne/Documents/PhD/wt_genome/Pm_genbank_gene_containing_only_no_alts.gb"

ingen = SeqIO.parse(ingb,"genbank")
newgen = []
for i in ingen:
    newi = SeqIO.SeqRecord(i.seq,id=i.id,name=i.name,description=i.description,dbxrefs=i.dbxrefs)
    for j in i.features:
        if j.type == 'mRNA':
            nt = j.qualifiers['note']
            if nt[0][-1] == 'A':
                newi.features.append(j)
        elif j.type == 'CDS':
            nt = j.qualifiers['note']
            if nt[0][-1] == 'A':
                newi.features.append(j)
        else:
            newi.features.append(j)
    newgen.append(newi)

SeqIO.write(newgen,outgb,"genbank")
