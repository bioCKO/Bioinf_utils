__author__ = 'mjohnpayne'

ingff = open("/Users/mjohnpayne/Documents/PhD/wt_genome/ts_wt_dbs/tsta1_working_models.gff3",'r')

outfile = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/smurf/Ts_smurf_in.txt","w")

for line in ingff:
    if line[0] != "#":
        col = line.split('\t')
        if col[2] == "gene":
            id = col[8].split(";")[1][5:]
            if col[6] == "-":
                cont = col[0]
                st = col[4]
                en = col[3]
                outfile.write(id+"\t"+cont+"\t"+st+"\t"+en+"\n")
            elif col[6] == "+":
                cont = col[0]
                st = col[3]
                en = col[4]
                outfile.write(id+"\t"+cont+"\t"+st+"\t"+en+"\n")
outfile.close()