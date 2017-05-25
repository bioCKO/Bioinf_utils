__author__ = 'mjohnpayne'


import sys

import subprocess
import os

genelist = open(sys.argv[1],"r").read().split('\n')

name = sys.argv[2]

type = sys.argv[3]

wig9 = "/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/pooled/9H.Fnor.smooth.wig"

wig11 = "/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/pooled/11H.Fnor.smooth.wig"

wig15 = "/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/pooled/15H.Fnor.smooth.wig"

wig17 = "/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/pooled/17H.Fnor.smooth.wig"

inputs = [wig9,wig11,wig15,wig17]

gene_file = open("/Volumes/MP_HD/Linda_MNase_Seq/an_genome/A_nidulans_FGSC_A4_version_s10-m03-r13_danpos2_genefile.txt","r").read().split('\n')

tmp_genfile_path = "tmp_gene_files.txt"

tmp_genfile = open(tmp_genfile_path,"w")

glist = []
for i in genelist:
    glist.append(i)

for i in gene_file:
    if i.split('\t')[0] in glist:
        tmp_genfile.write(i+'\n')
tmp_genfile.close()

if type == "pair":
    profile_args = "/usr/local/bin/python /Volumes/MP_HD/Linda_MNase_Seq/danpos-2.2.2/danpos.py profile "+ ",".join(inputs) + " --genefile_paths " + tmp_genfile_path + " --plot_colors red,blue,orange,purple,skyblue,cyan,green,blue4,darkgo --periodicity 0 --heatmap 1 --wigfile_aliases 9H,11H,15H,17H --genefile_aliases all_genes --name " + name
elif type == "diff":
    profile_args = "/usr/local/bin/python /Volumes/MP_HD/Linda_MNase_Seq/danpos-2.2.2/danpos.py profile /Volumes/MP_HD/Linda_MNase_Seq/danpos_r1_out/diff/Volumes_MP_HD_Linda_MNase_Seq_bam_11H-Volumes_MP_HD_Linda_MNase_Seq_bam_15H.pois_diff.wig --genefile_paths " + tmp_genfile_path + " --plot_colors red,blue,orange,purple,skyblue,cyan,green,blue4,darkgo --periodicity 0 --heatmap 1 --wigfile_aliases 11H-15H --genefile_aliases "+name+" --name " + name+"_"+type
subprocess.Popen(profile_args, shell=True).wait()

os.remove("tmp_gene_files.txt")