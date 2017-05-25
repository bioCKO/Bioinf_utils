__author__ = 'mjohnpayne'
import subprocess
from Bio import SeqIO
import os
import glob
import re

orthogroup_fastas = glob.glob("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/4species_orthogroups_fastas/all_groups/*.fasta")

all_talaro = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/all_talaro_proteins.fasta"

all_talaro_cds = "/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/all_talaro_cds.fasta"

all = SeqIO.parse(all_talaro,"fasta")
all_seqs = {}

for i in all:
    all_seqs[i.id] = i

allc = SeqIO.parse(all_talaro_cds,"fasta")
all_cds_seqs = {}
counts = {}
for i in allc:
    all_cds_seqs[i.id] = i
#     if i.id not in counts:
#         counts[i.id] = 1
#     else:
#         counts[i.id] += 1
# c = 0
# for i in counts:
#     if counts[i] > 1 and "TFLA" in i:
#         print i
#         c +=1
# print c


def rn(gene):
    if 'TF' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'TFLA_' + no*'0' + gene[3:] + '0'
    elif 'Pf' in gene:
        gene = gene[4:]
        no = 5-len(gene[3:])
        gene = 'PFUN_' + no*'0' + gene[3:] + '0'
    else:
        gene = gene[4:]
    return gene






def make_blastdb(infasta):
    mkdbargs = "/usr/bin/makeblastdb -in " + infasta + " -dbtype prot"
    subprocess.Popen(mkdbargs, shell=True).wait()

def run_blast_return_top(db,fasta):
    blast_args = '/opt/local/bin/blast -d ' + db + ' -i ' + fasta + ' -m 8 -p blastp'

    blast_out = subprocess.Popen(blast_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    out = blast_out.communicate()[0]
    out = str(out)
    out = out.split('\n')
    outls = []
    for line in out[1:-1]:
        line = line.split('\t')
        outls.append(rn(line[1]))
    idstring = "".join(outls)
    if "TFLA" in idstring:
        for i in outls:
            if "TFLA" in i:
                return i
    elif "PFUN" in idstring:
        for i in outls:
            if "PFUN" in i:
                return i
    elif "TSTA" in idstring:
        for i in outls:
            if "TSTA" in i:
                return i
    else:
        return "NONE"

def process_orthogroups(ingroup):
    group = SeqIO.parse(ingroup,"fasta")
    make_blastdb(ingroup)
    for i in group:
        if "PMAA" in i.id:
            i.id = i.id[4:]
            i.description = ""
            SeqIO.write(i,"tmp.fasta","fasta")
            top = run_blast_return_top(ingroup,"tmp.fasta")
            if top == "NONE":
                top = run_blast_return_top(all_talaro,"tmp.fasta")
            os.remove("tmp.fasta")
            if top == 'NONE':
                continue
            else:
                SeqIO.write([all_cds_seqs[i.id],all_cds_seqs[top]],"/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/positive_selection_tests/genome_wide/tm_and_top_ortho/" + i.id + "_top_ortho.fasta","fasta")
    os.remove(ingroup + ".phr")
    os.remove(ingroup + ".pin")
    os.remove(ingroup + ".psq")
c = 0
for i in orthogroup_fastas:
    if c%100 == 0:
        print c*100/len(orthogroup_fastas)
    process_orthogroups(i)
    c +=1
#top = run_blast_return_top("/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pm_proteins_with_byss.fasta","/Users/mjohnpayne/Documents/PhD/ASPS/Asp_blast_fro_tree/PMAA_090410.fasta")

#print top