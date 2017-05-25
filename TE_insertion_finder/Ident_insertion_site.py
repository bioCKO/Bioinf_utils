__author__ = 'mjohnpayne'

import glob
import subprocess
import os
import sys

from itertools import izip,islice,tee



ing = "./pmfa1_annot_scaf.fasta"

def fasta_to_dict(infasta):
    genome ={}
    gen = open(infasta,"r").read()
    gen = gen.split(">")
    for i in gen[1:]:
        j = i.split('\n')
        genome[j[0]] = "".join(j[1:])
    return genome

genome = fasta_to_dict(ing)

## function from stackoverflow to return hit with up to diffnumber misshits

def sub_findre(s,substring,diffnumber):
    sublen=len(substring)
    zip_gen=(izip(substring,islice(s,i,i+sublen)) for i in xrange(len(s)))
    for z in zip_gen:
        l,z=tee(z)
        if sum(1 for i,j in l if i==j)>=sublen-diffnumber:
            new=izip(*z)
            next(new)
            yield ''.join(next(new))

## returns contig, start and stop of hit

def run_blast_return_top(fasta):
    blast_args = './blast -d ./pmfa1_annot_scaf.fasta -i ' + fasta + ' -m 8 -p blastn -e 1'

    blast_out = subprocess.Popen(blast_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

    out = blast_out.communicate()[0]
    # print out
    # print len(out)
    if len(out) < 1:
        return "none","none","none"
    else:
        out = str(out)
        out = out.split('\n')
        hit = out[0].split('\t')
        return hit[1],str(hit[8]),str(hit[9])


inf = sys.argv[1].strip(" ") + "/*.seq"

inseqs = glob.glob(inf)

ingff = open("./pmfa1_working_models_fix.gff",'r')

## returns dictionary of contigs with subdictionaries of start positions as keys with end, pmaa and orientation as list values

def gene_pos_dict(gff):
    dict = {}
    for line in gff:
        if line[0] != "#":
            col = line.strip('\n').split('\t')
            if col[2] == 'gene':
                cont = col[0]
                st = col[3]
                en = col[4]
                pmaa = col[8].split(';')[1].replace("Name=","")
                orient = col[6]
                if cont not in dict:
                    dict[cont] = {st:[en,pmaa,orient]}
                else:
                    dict[cont][st] = [en,pmaa,orient]
    return dict

## makes tmp.fasta file of sequence read from TTAA onwards to use in blast

def make_fasta_of_site(i):
    name = i.split("/")[-1].strip('.seq')
    seq = ''
    dir = ''
    if "SS80" in i:
        dir = 5
        seq = ''
        inf = open(i,"r").read().replace('\r\n',"")
        p5_TE = "CTATCTTTCTAGGG"
        re_p5_TE = list(sub_findre(inf,p5_TE,4))[0]
        insertion = inf.find(re_p5_TE)+len(re_p5_TE)
        seq = inf[insertion:]
    elif "SS83" in i:
        dir = 3
        seq = ''
        inf = open(i,"r").read().replace('\r\n',"")
        p3_TE = "TCTTTCTAGGG"
        re_p3_TE = list(sub_findre(inf,p3_TE,4))[0]
        insertion = inf.find(re_p3_TE)+len(re_p3_TE)
        seq = inf[insertion:]
    #print seq
    outfasta = open("./tmp.fasta",'w')
    outfasta.write(">tmp\n" + seq)
    outfasta.close()
    return dir

## uses insertion position to search contig for gene start hit less than position then uses end to determine if inside gene otherwise finds orientation and position of next gene

def ident_env(cont,pos,gdict,outf):
    pos = pos+2
    site =  genome[cont][pos-2:pos+2]
    contlis = gdict[cont]
    poslis = sorted(map(int,contlis.keys()))
    for j in range(len(poslis)):
        key = str(poslis[j])
        if j == len(poslis)-1:
            key_p_1 = "terminus"
        else:
            key_p_1 = str(poslis[j+1])
        gene_inf = contlis[key]
        if int(poslis[j]) < int(pos) < int(gene_inf[0]):
            # print "contig",cont
            # print "start",poslis[j]
            # print "insertion",pos
            # print "end",gene_inf[0]
            size = int(gene_inf[0]) - poslis[j]
            lost = int(gene_inf[0]) - pos
            perc = (float(lost)/size)*100
            if contlis[key_p_1][2] == "+":
                outf.write("%s\t%s\t%s\t%s\tYes: %s\t%s\tN/A\tN/A\n" %(name,cont,pos,site,gene_inf[1],str(perc)[:5]))
                #print "%s site is in %s and disrupts %s%% of the gene" %(name,gene_inf[1],str(perc)[:5])
            elif contlis[key_p_1][2] == "-":
                perc = 100-perc
                outf.write("%s\t%s\t%s\t%s\tYes: %s\t%s\tN/A\tN/A\n" %(name,cont,pos,site,gene_inf[1],str(perc)[:5]))
                #print "%s site is in %s and disrupts %s%% of the gene" %(name,gene_inf[1],str(perc)[:5])
        elif int(gene_inf[0]) < int(pos) < int(key_p_1):
            # print "contig",cont
            # print "left gene pos",str(gene_inf[0])
            # print "insertion",pos
            # print "right gene pos",str(key_p_1)
            left_d = int(pos) - int(gene_inf[0])
            right_d = int(key_p_1) - int(pos)
            if gene_inf[2] == "-":
                if contlis[key_p_1][2] == "-":
                    outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp upstream of %s\t%sbp downstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                    #print "%s site is %sbp upstream of %s and %sbp downstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])
                elif contlis[key_p_1][2] == "+":
                    outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp upstream of %s\t%sbp upstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                    #print "%s site is %sbp upstream of %s and %sbp upstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])
            elif gene_inf[2] == "+":
                if contlis[key_p_1][2] == "-":
                    outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp downstream of %s\t%sbp downstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                    #print "%s site is %sbp downstream of %s and %sbp downstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])
                elif contlis[key_p_1][2] == "+":
                    outf.write("%s\t%s\t%s\t%s\tNo\tN/A\t%sbp downstream of %s\t%sbp upstream of %s\n" %(name,cont,pos,site,left_d,gene_inf[1],right_d,contlis[key_p_1][1]))
                    #print "%s site is %sbp downstream of %s and %sbp upstream of %s" %(name,left_d,gene_inf[1],right_d,contlis[key_p_1][1])



gene_dict = gene_pos_dict(ingff)
# genome = {}
# for j in ingenome:
#     genome[str(j.id)] = str(j.seq)

ls = []
outf = open(sys.argv[1].strip(" ") + "/insertion_site_output.tsv","w")
outf.write("ID\tcontig\tInsertion Co-ordinate\tInsertion sequence\tin gene?\tPercentage disrupted\t5 prime gene\t3 prime gene\n")

for i in inseqs:
    name = i.split("/")[-1].strip('.seq')
    #print '\n\n',name
    dir = make_fasta_of_site(i)
    contig,st,en = run_blast_return_top("./tmp.fasta")
    os.remove("./tmp.fasta")
    if contig == "none":
        outf.write(name + "\tNo Hits\n")
    else:
        start = int(st)
        end = int(en)
        coord = 0
        if start > end:
            e = start
            s = end
            site =  genome[contig][s-1:e]
            #print site,'\n'
            coord = e-4
        else:
            site =  genome[contig][start-1:end]
            #print site,'\n'
            coord = start-1
        ls.append(contig + '-' + str(coord))
        ident_env(contig,coord,gene_dict,outf)




