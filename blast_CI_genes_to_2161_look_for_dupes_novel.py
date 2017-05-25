__author__ = 'mjohnpayne'

# Blast CI proteins from unnaligned contigs against Pm proteins
# if protein hits = duplicated protein
# if protein doesn't hit highly = new protein

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import subprocess

pmdb = '/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/all_pm_protein.fasta'
fungidb = '/Volumes/MP_HD/nr_fastas/All_fungi_refseq_prot_6-8-14.fasta'

##record = SeqIO.read("m_cold.fasta", format="fasta")
##result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)


## first generate list of contigs in abacas bin

def extract_contigs(cont):
    lis = []
    for line in cont: lis.append(line.replace(' \n',''))
    return lis

## second extract genes (prot seq) from those contigs from augustus output



def get_seqs(an,con):
    ## generate dictionary of contigs attached to list of gene id:seq pairs
    an = an.readlines()
    dic = {}
    cont,gene,pseq = '','',''
    for i in range(len(an)):
        line = an[i]
        ids = ''
        if '\tgene\t' in line:
            col = line.split('\t')
            cont = col[0]
            gene = col[8].strip('ID=').replace('\n','')
        st,en = 0,0
        if '= [' in line:
            st = i
            start = st
            curr = an[start]
            pseq = []
            while 'end gene' not in curr:
                pseq.append(an[start])
                start += 1
                curr = an[start]
            pseq = ''.join(pseq)
            pseq = pseq.replace('# ','').replace('protein sequence = [','').replace(']','').replace('\n','')
            pseq = Seq(pseq)
            if cont not in dic:
                dic[cont] = [SeqIO.SeqRecord(pseq,gene)]
            else:
                dic[cont] += [SeqIO.SeqRecord(pseq,gene)]
    newd = {}
    for i in con:
        if i in dic:
            newd[i] = dic[i]
    return newd


        # if 'end gene' in line:
        #     en = i
        #     print st
        #     print en
        #     pseq = ''.join(an[st:en])
        #     pseq.replace('# ','').replace('protein sequence = [','').replace(']','').strip('\n')
        #     #print pseq



## third blast these genes against ncbi NR database

def blast_gene(seq,database):
    tempfasta = open('temp.fasta','w')
    SeqIO.write(seq,tempfasta,'fasta')
    tempfasta.close()
    run = blastp(query='temp.fasta',db=database,num_descriptions=5,num_threads=6,outfmt=5,out='temp.xml')
    run()
    result_handle = open('temp.xml')
    result = NCBIXML.read(result_handle)
    rets = []
    for i in result.descriptions:
        ttl = i.title
        e = i.e
        if 'Tfl|' in ttl:
            species = 'T. flavus'
            d = ttl[ttl.find('Tfl'):]
        elif 'Pfu|' in ttl:
            species = 'P. funiculosum'
            d = ttl[ttl.find('Pfu'):]
        elif 'PMAA_' in ttl:
            species = 'T. marneffei'
            d = ttl[ttl.find('PMAA'):]
        else:
            species = ttl[ttl.find('[')+1:ttl.find(']')]
            d = ttl[ttl.find('| ')+1:ttl.find('[')-1]
        rets.append(species)
        rets.append(d)
        rets.append(str(e))
    return rets



## fourth parse blast output to determine if top hit is in Pm or another species


def main(contigs,annots,o):
    cont = open(contigs,'r')
    outfile = open(o,'w')
    conts = extract_contigs(cont)
    annot = open(annots,'r')
    seqs = get_seqs(annot,conts)
    for i in seqs:
        print i
        for j in seqs[i]:
            print j.id
            res = blast_gene(j,fungidb)
            outfile.write(i + '\t' + j.id + '\t' + '\t'.join(res) + '\n')
    outfile.close()




testseq = SeqIO.SeqRecord(Seq("MGMNINQILVESLTHLNYAFGYITPETYKIGVMPGVDASTFSDFTALKSKNSDLKTFITHLLAFMRHYGFDGVDFDWEYPGATDRQPNELNS"), id='g9865')

#print blast_gene(testseq)

testan = '/Volumes/MP_HD/CI_GENOME_SEQ/augustus_gene_finding(stats_for_fig)/velvet_assemblies_gff/012_vel_denovo.gff'
testcon = '/Volumes/MP_HD/CI_GENOME_SEQ/CI_denovo_assemblies(stats)/velvet_assemblies_scaf/abacas/012_vel_scaffolds.fasta_pmfa1_annot_scaf_concat.fasta.bin'
of = sys.argv[2][:-4] + '_unmapped_contig_blast.txt'
main(sys.argv[1],sys.argv[2],of)

