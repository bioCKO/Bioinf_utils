__author__ = 'mjohnpayne'


import glob
from time import sleep as sl



inlist = glob.glob("/Users/mjohnpayne/Documents/PhD/bioinformatics_tools/4_species_genome_annot/p450_DB/*_blast_p450_out.txt")

outf = "/Users/mjohnpayne/Documents/PhD/bioinformatics_tools/4_species_genome_annot/p450_DB/"

def parse_tmhmmout(inls,outpref):
    for i in inls:
        inf = open(i,"r").readlines()
        spec = i.split('/')[-1][:-19]
        out = open(outf+spec+"p450_parsed.txt","w")
        out.write('GeneID\tcyp_top_hit\ttop_hit_species\n')
        done = []
        for j in inf:
            if j[0] == ' ' or j == "\n":
                continue
            elif "#" not in j:
                col = j.split('\t')
                id = col[0]
                if id not in done:
                    hit = col[1]
                    cyp = hit.split("_")[0]
                    cypshort = cyp.rstrip('123456789')
                    if len(hit.split("_")) > 1:
                        hitspec = "_".join(hit.split("_")[1:])
                        out.write(id+'\t'+cyp+'\t'+cypshort+'\t'+hitspec+'\n')
                    else:
                        out.write(id+'\t'+cyp+'\t'+cypshort+'\tgeneric\n')
                    done.append(id)
        out.close()


parse_tmhmmout(inlist,outf)