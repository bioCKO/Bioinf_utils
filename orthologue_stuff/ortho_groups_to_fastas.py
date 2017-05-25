import sys
import os.path
import glob
from Bio import SeqIO
import re
import time

ortho_output = open(sys.argv[1],'r')#
#ortho_output = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups.txt','r')

fastas = sys.argv[2]
#fastas = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/4_species_all_cds.fasta'
allfasta = SeqIO.parse(fastas,'fasta')
# outpath = sys.argv[3]
# st = int(sys.argv[4])
# en = int(sys.argv[5])

def absoluteFilePaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

def once_in_string(inlis,string):
    count = 0
    for i in inlis:
        if string.count(i+'|') == 1:
            count +=1
    if count == len(inlis):
        return "Yes"
    else:
        return "No"



def once_or_more_in_string(inlis,string):
    outls = []
    count = 0
    nocount = 0
    for i in inlis:
        account = string.count(i+'|')
        if account >= 1:
            count += 1
            nocount += account
            outls += [account]
    if count == len(inlis) and nocount > len(inlis):
        return ("Yes",outls)
    else:
        return "No"


def other_orthologues(inlis,string):
    outls = []
    count = 0
    nocount = 0
    for i in inlis:
        account = string.count(i+'|')
        if account >= 1:
            count += 1
            nocount += account
            outls += [account]
        else:
            outls += [account]
    if count < len(inlis) and count > 1:
        return ("Yes",outls)
    else:
        return "No"

def paralogues(inlis,string):
    outls = []
    count = 0
    nocount = 0
    for i in inlis:
        account = string.count(i+'|')
        if account >= 1:
            count += 1
            nocount += account
            outls += [account]
        else:
            outls += [account]
    if count == 1:
        return ("Yes",outls)
    else:
        return "No"

spec = {}
names = []

for i in allfasta:
    spec[i.id] = i.seq
    acc = str(i.id)
    if acc[:acc.index('|')] not in names:
        names.append(acc[:acc.index('|')])


## print total gene numbers
allfasta.close()

allfasta = open(fastas,'r').read()

gene_no = {}
for k in names:
    gene_no[k] = 0


for j in names:
    gene_no[j] = allfasta.count(j + '|')

# for i in names:
#     print i, gene_no[i]
#

### print gene numbers in groups

# ortho_output.close()
#
# ortho_output = open(sys.argv[1],'r').read()
#
# gene_no = {}
# for k in names:
#     gene_no[k] = 0
#
#
# for j in names:
#     gene_no[j] = ortho_output.count(j + '|')
#
# for i in names:
#     print gene_no[i]

####print names

##for i in absoluteFilePaths(fasta_folder):
##    if '.fasta' in i:
##        name = i[-9:-6]
##        names.append(name)
##        species[name] = SeqIO.parse(i,'fasta')
##        spec[name] = {}
##        for i in species[name]:
##            spec[name][i.id] = i.seq

### number of groups with 1 in all orthologues (ie completely conserved)



# num = 0
#
# for j in ortho_output:
#     if once_in_string(names,j) == 'Yes':
#         num +=1
# print "number 1-1-1-1 orthologue groups: " + str(num)
#
# ortho_output.close()

### number of groups with at least 1 in all orthologues (ie completely conserved with expansions)

# num2 = 0
#
# ortho_output = open(sys.argv[1],'r')
#
# speccount = [0]*22
# print speccount
#
# for j in ortho_output:
#     outs = once_or_more_in_string(names,j)
#     if outs[0] == 'Yes':
#         num2 +=1
#         for i in range(len(speccount)):
#             speccount[i] += outs[1][i]
#
# print "number N-N-N-N orthologue groups: " + str(num2)
# print "numer of genes in above species"
# for s in range(len(names)):
#     print names[s],str(speccount[s])
#
# ortho_output.close()

### Number of genes not in 1-1-1-1 or N-N-N-N but with orthologues in at least one other species


# num3 = 0
#
# ortho_output.close()
#
# ortho_output = open(sys.argv[1],'r')
#
# speccount = [0]*22
# print speccount
#
# for j in ortho_output:
#     outs = other_orthologues(names,j)
#     if outs[0] == 'Yes':
#         num3 +=1
#         for i in range(len(speccount)):
#             speccount[i] += outs[1][i]
#
# print "number other orthologue groups: " + str(num3)
# print "numer of genes in above species"
# for s in range(len(names)):
#     print names[s],str(speccount[s])
#
# ortho_output.close()

### Number of paralogue only groups

# num4 = 0
#
# ortho_output.close()
#
# ortho_output = open(sys.argv[1],'r')
#
# speccount = [0]*22
# print speccount
#
# for j in ortho_output:
#     outs = paralogues(names,j)
#     if outs[0] == 'Yes':
#         num4 +=1
#         for i in range(len(speccount)):
#             speccount[i] += outs[1][i]
#
# print "number paralogue groups: " + str(num4)
# print "numer of genes in above species"
# for s in range(len(names)):
#     print names[s],str(speccount[s])
#
# ortho_output.close()


### Extract pm genes not in any groups - Fix to add pmaas from orthos with many pmaas
# no_group_ls = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthologue_subsets/pm_in_none.txt','w')
# all_pms = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/All_species_prot_dbs/Pma.fasta','r')
# allpms2 = []
# for line in all_pms:
#     if '>' in line:
#         allpms2 += [line[5:].strip('\n')]
# pmlis = []
# for j in ortho_output:
#     if 'PMAA' in j:
#         st = [m.start() for m in re.finditer('PMAA',j)]
#         for i in st:
#             pmaa = j[i:i+11]
#             pmlis += [pmaa]
#
# no_group = list(set(allpms2) - set(pmlis))
# for i in no_group:
#     no_group_ls.write(i + '\n')
#
# no_group_ls.close()

###Number of groups with one or more from all species

# num = 0
# pmaas = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthologue_subsets/1_or_more_in_all.txt','w')
#
# for j in ortho_output:
#     lincount = 0
#     for i in names:
#         if i in j:
#             lincount += 1
#     if lincount == len(names):
#         num += 1
#         st = j.index('PMAA')
#         pmaas.write(j[st:st+11] + '\n')
# print num
#
# pmaas.close()

#### section gives % of genes in orthologue groups
        
# countin = ortho_output.read()
# print countin[3]
# goodprot = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/orthomcl/int_files/goodProteins.fasta').read()
#
# for i in names:
#   orth_count = countin.count(i+'|')
#   tot_count = goodprot.count(i+'|')
#   print i + '\t' + str(orth_count)+ '\t' + str(tot_count) + '\t' + str((float(orth_count)/tot_count)*100)


###### produces orthologous genes for each species in multifasta files where 1st gene in each file is orthologous to other first genes and so on.
# print st
# print en
#
# c=0
#
#
# for line in ortho_output:
#        if once_in_string(names,line) == "Yes" and c > st and c <= en:
#            col = line.strip('\n').split(' ')
#            col = col[1:]
#            done = []
#            for i in col:
#               for j in names:
#                  if j not in done:
#                     if j == i[:i.index('|')]:
#                        outtemp = open(str(outpath) + j + '_orthos.fasta','a')
#                        outtemp.write(">" + i + '\n' + str(spec[i]) + '\n')
#                        outtemp.close()
#                        done.append(j)
#                  else:
#                     continue
#            c += 1
#        elif once_in_string(names,line) == "Yes" and c <= st:
#            c += 1
# ortho_output.close()

# produces fasta file for each orthologue group

# c=0
#
# outpath = '/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/CI_orthogroup_fasta_files/'
# for line in ortho_output:
#       col = line.strip('\n').split(' ')
#       name = col[0][:-1]
#       col = col[1:]
#       if once_in_string(names,line) == 'Yes':
#           for i in col:
#                 outtemp = open(str(outpath) + name + '_orthos.fasta','a')
#                 if i in spec:
#                    outtemp.write(">" + i + '\n' + str(spec[i]) + '\n')
#                 else:
#                    outtemp.close()
#                    os.remove(str(outpath) + name + '_orthos.fasta')
#                    break
#                 outtemp.close()
#
# ortho_output.close()

# produces fasta file for each orthologue group with pm and removes stop and renames files

# tmpout = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/orthomcl_ortho_groups_stats_2.txt','w')
# c=0
# lis = []
# types = {}
# types2 = {}
# outpath = '/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/ortho_group_cds_fastas/all_cds_pm_orthos/'
# for line in ortho_output:
#       if 'Pma|' in line:
#           c+=1
#           col = line.strip('\n').split(' ')
#           name = line[line.index('Pma|')+4:line.index('Pma|')+15]
#           col = col[1:]
#           tfc = line.count('Tfl|')
#           pmc = line.count('Pma|')
#           tsc = line.count('Tst|')
#           pfc = line.count('Pfu|')
#           k = (pmc,tfc,pfc,tsc)
#           if k not in types:
#               types[(pmc,tfc,pfc,tsc)] = 1
#               types2[(pmc,tfc,pfc,tsc)] = [name]
#           else:
#               types[(pmc,tfc,pfc,tsc)] += 1
#               types2[(pmc,tfc,pfc,tsc)] += [name]
#           # for i in col:
#           #     if 'Pma|' in i:
#           #         # if i not in lis:
#           #         #   lis += [i]
#           #       outtemp = open(str(outpath) + name + '_orthos.fasta','a')
#           #       if i in spec:
#           #          outtemp.write(">" + i + '\n' + str(spec[i][:-3]) + '\n')
#           #       else:
#           #          outtemp.close()
#           #          os.remove(str(outpath) + name + '_orthos.fasta')
#           #          break
#           #       outtemp.close()
# # print len(lis)
# # print lis[:100]
# # print c
#
# tmpout.write('pm\ttf\tpf\tts\tcount\tids\n')
# for i in types:
#     tmpout.write('\t'.join(map(str,i)) + '\t' + str(types[i]) + '\t' + ','.join(types2[i]) + '\n')
# ortho_output.close()
# tmpout.close()

### Print names of every orthologous group that a species is in

# out = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/ortho_groups_my_species/Pfu_groups.txt','w')
#
# for line in ortho_output:
#     if 'Pfu|' in line:
#           col = line.strip('\n').split(' ')
#           out.write(col[0][:-1] + '\n')
#
# out.close()

### Make fasta file of every group 2161 is not in

#
# c=0
# lis = []
# types = {}
# types2 = {}
# outpath = '/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/CI_orthogroup_no_2161_for_annot/'
# for line in ortho_output:
#     if '2161|PMAA' not in line:
#         c+=1
#         col = line.strip('\n').split(' ')
#         groupname = col[0][:-1]
#         for i in col[1:]:
#             outtemp = open(str(outpath) + groupname + '_orthos.fasta','a')
#             if i in spec:
#                outtemp.write(">" + i + '\n' + str(spec[i]) + '\n')
#             # else:
#             #    outtemp.close()
#             #    os.remove(str(outpath) + name + '_orthos.fasta')
#             #    break
#             outtemp.close()

#### write one fasta with all genes in clusters with no 2161

c=0
lis = []
types = {}
types2 = {}
outpath = '/Volumes/MP_HD/CI_GENOME_SEQ/CI_orthomcl_data/CI_orthogroup_no_2161_for_annot/all_no_2161_genes.fasta'

outtemp = open(outpath,'w')
for line in ortho_output:
    if '2161|PMAA' not in line:
        c+=1
        col = line.strip('\n').split(' ')
        groupname = col[0][:-1]
        for i in col[1:]:
            if i in spec:
               outtemp.write(">" + i + "|" + groupname + '\n' + str(spec[i]) + '\n')
            # else:
            #    outtemp.close()
            #    os.remove(str(outpath) + name + '_orthos.fasta')
            #    break
outtemp.close()