import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.CheckSum import seguid
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC

in_gb = SeqIO.parse('/Users/mjohnpayne/Desktop/merge_pm_embl/Pm_from_gb.gb', "genbank")
outfile = open('/Users/mjohnpayne/Documents/Uni/phd/wt_genome/pm_wt_dbs/pm_proteins_with_alts.fasta','w')

cds = {}
protein = {}
ids = {}
for record in in_gb:
    for i,feature in enumerate(record.features):
        if feature.type == 'CDS':
            pmaa = feature.qualifiers
            ids = pmaa['note'][0]
            accession = ids[22:34]
            prot = feature.qualifiers['translation'][0]
            protein[ids] = prot
            outfile.writelines('>' + accession + '\n' + prot + '\n')

outfile.close()
in_gb.close()



