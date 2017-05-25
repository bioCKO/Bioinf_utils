__author__ = 'mjohnpayne'

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
import subprocess

## work out each piece of chromosome to ident

frags = {93:[0,1120000,1797000,1925000,2025471],94:[0,175663],95:[0,775000,1818470],96:[0,2177000,3118000,4551145,6407042],97:[0,1822000,2041000,3072962],99:[0,2873000,2932000,3229885],100:[0,1298000,2131000,4167742],102:[0,3746755],103:[0,1539000,3339384],104:[0,139459]}

# for i in frags:
#     print i
#     print frags[i]

ingenome = SeqIO.parse('/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta','fasta')


gen = {}

for i in ingenome:
    gen[i.id] = str(i.seq)
seqs = {}
for i in frags:
    seq = gen[str(i)]
    count = 1
    for j in range(len(frags[i])-1):
        seqs[str(i) + '_' + str(count)] = seq[frags[i][j]:frags[i][j+1]]
        count += 1

# for i in seqs:
#     print i
#     print len(seqs[i])

##blast primer database against each piece

def blast_gene(seq,database):
    tempfasta = open('temp.fasta','w')
    SeqIO.write(seq,tempfasta,'fasta')
    tempfasta.close()
    run = blastn(query='temp.fasta',db=database,num_descriptions=1,num_threads=6,outfmt=5,word_size=4,evalue=0.01,task="megablast",out='temp.xml')
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
    for i in result.alignments:
        for j in i.hsps:
            rets.append(str(j.frame[1]))
            rets.append(str(j.query))
            rets.append(str(j.match))
            rets.append(str(j.sbjct_start))
    return rets

def check_re(string,contig,st,en):
    sequence = gen[contig][st:en]
    num = sequence.count(string)
    return num

primer = SeqIO.SeqRecord(Seq('ACATCGATTAATAGCAACACAAAGAGCG'),'p')

infile = open("/Users/mjohnpayne/Documents/PhD/lab_dbs/oligo_db_29-4-14",'r').read()
lines = infile.split('\r')

primers = []
for i in lines:
    i = i.replace('\x0b','')
    col = i.split('\t')
    if len(col) > 3:
        if len(col[3]) > 2 and 'N' not in col[3]:
            ids = str(col[2] + col[0])
            primers += [SeqIO.SeqRecord(Seq(col[3]),ids)]


# pos = {}
# outf = open('temp_oligo_info','w')
# for i in primers:
#     res = blast_gene(i,'/Users/mjohnpayne/Documents/PhD/wt_genome/pmfa1_annot_scaf.fasta')
#     if len(res) > 0:
#         if len(i.seq) == len(res[4]):
#             if res[0] not in pos:
#                 pos[res[0]] = [(i.id,res[5])]
#             else:
#                 pos[res[0]] += [(i.id,res[5])]
#             outf.write(i.id + '\t' + str(res[0]) + '\t' + str(res[5]) + '\t' + str(res[2]) + '\n')
#             # print i.id
#             # print res[0]
#             # print res[1]
#             # print res[2]
#             # print res[5] + '\n'
# outf.close()

def re_frag_size(pattern,contig,ppos):
    repos = []
    win = len(pattern)
    st = 0
    en = st + win
    while en < len(gen[contig]):
        if gen[contig][st:en] == pattern:
            repos.append(str(st))
        st += 1
        en += 1
    down = []
    repos = [0] + repos + [len(gen[contig])]
    down = [int(i) - int(ppos) for i in repos]
    lower = max(j for j in down if j < 0)
    upper = min(k for k in down if k > 0)
    lower = repos[down.index(lower)]
    upper = repos[down.index(upper)]
    size = int(upper) - int(lower)
    return size


ininf = open('temp_oligo_info','r').readlines()

pos = {}

for j in ininf:
    i = j.strip('\n').split('\t')
    if i[1] not in pos:
        pos[i[1]] = [(i[0],i[2],i[3])]
    else:
        pos[i[1]].append((i[0],i[2],i[3]))

#re_frag_size('GCGGCCGC','93',1000000)

mot = 'GGTACC'

outfile = open('all_oligo_out.txt','w')
outfile.write('contig\tprimer1\tprimer2\tstart\tend\tsize\n')##\trefrag_size\n')
for i in sorted(pos.keys()):
    for j in pos[i]:
        for k in pos[i]:
            diff = int(j[1]) - int(k[1])
            if diff > -6000 and diff < 6000 and diff != 0:
                if diff > 200 or diff < -200:
                    if j[2] == '1' and k[2] == '-1' and int(j[1]) < int(k[1]):
                        #if check_re(mot,i,int(j[1]),int(k[1])) == 0:
                         #   resize = re_frag_size(mot,i,j[1])
                         #   if resize < 20000:
                        print i
                        print j
                        print str(k) + '\n'
                        print check_re(mot,i,int(j[1]),int(k[1]))
                        outfile.write(i + '\t' + j[0] + '\t' + k[0] + '\t' + j[1] + '\t' + k[1] + '\t' + str(int(k[1])-int(j[1])) + '\n')##+ '\t' + str(resize) + '\n')
                # elif j[2] == '-1' and k[2] == '1' and int(j[1]) > int(k[1]):
                #     print i
                #     print j
                #     print str(k) + '\n'
                #     outfile.write(i + '\t' + j[0] + '\t' + k[0] + '\t' + j[1] + '\t' + k[1] + '\t' + str(int(k[1]-j[1])) + '\n')

outfile.close()

## Pick top 2 or 3 hits

## check fragments for 8 base cutters

##output desired primers