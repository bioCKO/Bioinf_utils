__author__ = 'mjohnpayne'

import glob
from time import sleep as sl
import html2text


instats = glob.glob("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_indels/*.html")

outfile = open("/Volumes/MP_HD/CI_GENOME_SEQ/Clinical_isolates_sort_bam_(lofreq_in_progress)/GATK_processing/GATK_indel_calls/SNPeff_indels/indel_stats.txt","w")

outfile.write("strain\tno_ins\tno_del\tdisruptive_inframe_deletion\tdisruptive_inframe_insertion\tframeshift_variant\tinframe_deletion\tinframe_insertion\n")

for i in instats:
    spec = i.split('/')[-1].replace("_GATK_filtered_indels_pass_stats.html","")
    outfile.write(spec+"\t")
    inf = open(i,"r").read()
    print i,spec
    h = html2text.HTML2Text()
    f = h.handle(inf.decode('utf8')).split('\n')
    # for j in f[:200]:
    #     print j
    # sl(2)
    for j in f:
        if "** INS **" in j:
            no = j.replace("** INS ** |  ","").strip('   \n')
            outfile.write(no+"\t")
        elif "** DEL **" in j:
            no = j.replace("** DEL ** |  ","").strip('   \n')
            outfile.write(no+"\t")
        elif "** disruptive_inframe_deletion **" in j:
            no = j.replace("** disruptive_inframe_deletion ** |    |  ","").strip('\n')
            no = no[:no.find("|")-2]
            outfile.write(no+"\t")
        elif "** disruptive_inframe_insertion ** " in j:
            no = j.replace("** disruptive_inframe_insertion ** |    |  ","").strip('\n')
            no = no[:no.find("|")-2]
            outfile.write(no+"\t")
        elif "** frameshift_variant ** " in j:
            no = j.replace("** frameshift_variant ** |    |  ","").strip('\n')
            no = no[:no.find("|")-2]
            outfile.write(no+"\t")
        elif "** inframe_deletion **" in j:
            no = j.replace("** inframe_deletion ** |    |  ","").strip('\n')
            no = no[:no.find("|")-2]
            outfile.write(no+"\t")
        elif "** inframe_insertion ** " in j:
            no = j.replace("** inframe_insertion ** |    |  ","").strip('\n')
            no = no[:no.find("|")-2]
            outfile.write(no+"\n")

outfile.close()