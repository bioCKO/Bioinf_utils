import subprocess
import shlex


def run_blast_return_top(db,fasta):
    blast_args = '/opt/local/bin/blast -d ' + db + ' -i ' + fasta + ' -m 8 -p blastp'

    blast_out = subprocess.Popen(blast_args, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)

    out = blast_out.communicate()[0]
    out = str(out)
    out = out.split('\n')

    out_list = []

    for line in out:
        line = line.split('\t')
        out_list.append(line)
        #print line
    return out_list[0][1]
    # for line in out_list[:-1]:
    #    print line[1]


top = run_blast_return_top("/Users/mjohnpayne/Documents/PhD/wt_genome/pm_wt_dbs/pm_proteins_with_byss.fasta","/Users/mjohnpayne/Documents/PhD/ASPS/Asp_blast_fro_tree/PMAA_090410.fasta")

print top