from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
import sys

in_tfgff = open(sys.argv[1],'r')#'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_denovo.gff','r')

in_tffas = SeqIO.parse(sys.argv[2],'fasta')#'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/PF_SSPACE.final.scaffolds.fasta','fasta')

outcds = open(sys.argv[3],'w')#'/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/P_funiculosum/Pf_vel_cds.fa','w')

#outprot = open(sys.argv[4],'w')




tf_fas = {}
for i in in_tffas:
    tf_fas[i.id] = str(i.seq)
cds = ''
name = ''
orient = '+'
scaf = ''
for line in in_tfgff:
    if 'JCVI' in line or "AUGUSTUS" in line or "Genbank" in line:
        col = line.strip('\n').split('\t')
        if col[2] == 'gene':
            inf = col[8].split(';')
            name = inf[1].strip('Name=')
            for i in tf_fas:
                if i == col[0]:
                    scaf = tf_fas[i]
            orient = col[6]
        elif col[2] == 'CDS':
            st = int(col[3])
            en = int(col[4])
            cds += scaf[st-1:en]
    elif '###' in line:
        if orient == '-':
            cds = str(Seq(cds,generic_nucleotide).reverse_complement())
        outcds.write('>' + name + '\n' + cds + '\n')
        #prot = outprot.write('>' + name + '\n' + str(Seq(cds,generic_nucleotide).translate()) + '\n')
        cds = ''
        
            
            
outcds.close()
            
