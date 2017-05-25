__author__ = 'mjohnpayne'

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Seq import Seq
import time
import re

in_pm1 = '/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/PM_1_prot.fasta'
in_2161 = '/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/2161_prot.fasta'
# pm1dict = {}
# for i in in_pm1:
#     pm1dict[i.id] = SeqIO.SeqRecord(i.seq,i.id)
#
# for i in pm1dict:
#     print i
#     print pm1dict[i]


# inorthos = open('/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/inparanoid_PM1/table.PM_1_prot.fasta-2161_prot.fasta','r').readlines()
# outlist = open('/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/inparanoid_PM1/PM_1_PMAA_conversions.txt','w')
#
#
# for i in inorthos[1:]:
#     col = i.split('\t')
#     PMA = re.split(' *',col[2])[0]
#     PM1 = re.split(' *',col[-1])[0]
#     outlist.write(PM1 + '\t' + PMA + '\n')
#     #outseq = SeqIO.SeqRecord(in_pm1[])
#
# outlist.close()

outstats = open('/Users/mjohnpayne/Documents/PhD/wt_genome/PM1_strain_genome/inparanoid_PM1/PM_1_PMAA_conversions_stats.txt','w')

def blast_gene(seq,database):
    tempfasta = open('temp.fasta','w')
    SeqIO.write(seq,tempfasta,'fasta')
    tempfasta.close()
    run = blastp(query='temp.fasta',db=database,max_target_seqs=1,num_threads=6,outfmt=5,out='temp.xml')
    run()
    result_handle = open('temp.xml')
    result = NCBIXML.read(result_handle)
    if len(result.descriptions) > 0:
        rets = result.descriptions[0].title
        rets = [rets[rets.find('GQ'):]]
        rets.append(str(result.descriptions[0].e))
        return rets
    else:
        return ['none']
county = 0
countn = 0
for i in SeqIO.parse(in_2161,'fasta'):
    out = blast_gene(i,in_pm1)
    if out[0] != 'none':
        if float(out[1]) < float(1.0e-15):
            county += 1
            outstats.write(i.id + '\t' + '\t'.join(blast_gene(i,in_pm1)) + '\n')
    else:
        countn += 1

print 'Yes\t' + str(county)
print 'No\t' + str(countn)
