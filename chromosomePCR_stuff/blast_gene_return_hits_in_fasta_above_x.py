__author__ = 'mjohnpayne'


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
import sys
import os

indb = "/Volumes/MP_HD/pm1_2161_CI/All_strains_blast_dbs/CI_talaro_db_dedup.fasta"
in_id = "/Volumes/MP_HD/pm1_2161_CI/pm1_MADS_genes.txt"
in_eval = "1e-10"



def blast_gene(ids,eval,database,of):
    fasta_sequences = SeqIO.parse(open(database),"fasta")
    for seq in fasta_sequences:
        if seq.id == ids:
            SeqIO.write(seq,"temp.fasta", "fasta")
    run = blastp(query='temp.fasta',db=database,num_threads=6,outfmt=5,word_size=4,evalue=eval,out='temp.xml')
    run()
    result_handle = open('temp.xml')
    result = NCBIXML.read(result_handle)
    rets = []
    for i in result.descriptions:
        ttl = i.title
        e = i.e
        species = ttl.split(' ')[0]
        rets.append(species)
        rets.append(str(e))
    # for i in result.alignments:
    #     for j in i.hsps:
    #         rets.append(str(j.frame[1]))
    #         rets.append(str(j.query))
    #         rets.append(str(j.match))
    #         rets.append(str(j.sbjct_start))
    os.remove('temp.fasta')
    os.remove('temp.xml')
    genlis = []
    for i in range(0,len(rets),2):
        genlis.append(rets[i])
        print rets[i]
    fasta_sequences = SeqIO.parse(open(database),"fasta")
    seqs = []
    for seq in fasta_sequences:
        if seq.id in genlis:
            seqs.append(seq)
    SeqIO.write(seqs,of, "fasta")

for i in open(in_id,"r").readlines():
    i=i.strip('\n')
    print i
    outf = "/".join(in_id.split('/')[:-1])+"/"+i+"_blastp_hits_"+in_eval+".fasta"
    blast_gene(i,in_eval,indb,outf)
