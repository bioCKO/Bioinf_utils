__author__ = 'mjohnpayne'

from Bio import Entrez
import re
from time import sleep as sl

################### NOT WORKING ###################


# Entrez.email = "mpayne@student.unimelb.edu.au"
# handle = Entrez.efetch(db="nucleotide", id="327358302", retmode="txt",rettype="fasta")
# print handle.read()
#
#
# handle.close()


infile = open("/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/Eurotiomycetidae_proteomes/All_species_prot_dbs/Ade.fasta",'r').read()

# def change_ade_geneID(file):
#     positions = [m.start() for m in re.finditer("Ade\|", file)]
#     ids = []
#     for i in positions:
#         id = file[i+4:i+13]
#         ids.append(id)
#     Entrez.email = "mpayne@student.unimelb.edu.au"
#     handle = Entrez.efetch(db="nucleotide", id=",".join(ids), retmode="txt",rettype="gb")
#     newids = []
#     for i in handle.readlines():
#         if "/locus_tag=" in i:
#             p = i.find("/locus_tag=")
#             newid = i[p+12:p+22]
#             newids.append(newid)
#     corresp = {}
#     print len(ids),len(newids)
#     # for i in range(len(ids)):
#     #     corresp[ids[i]] = newids[i]
#     # for i in corresp:
#     #     print i,corresp[i]
#     #     sl(0.2)
#     handle.close()

def change_ade_geneID(file):
    positions = [m.start() for m in re.finditer("Ade\|", file)]
    ids = []
    for i in positions:
        id = file[i+4:i+13]
        ids.append(id)
    Entrez.email = "mpayne@student.unimelb.edu.au"
    handle = Entrez.efetch(db="nucleotide", id=",".join(ids[:10]), retmode="xml")
    records = Entrez.read(handle)
    newids = []
    for record in records:
        newid = record[u'GBSeq_other-seqids'][1].split("|")[-1][:-2]
        newids.append(newid)
    # newids = []
    # for i in handle.readlines():
    #     if "/locus_tag=" in i:
    #         p = i.find("/locus_tag=")
    #         newid = i[p+12:p+22]
    #         newids.append(newid)
    # corresp = {}
    print len(ids),len(newids)
    # for i in range(len(ids)):
    #     corresp[ids[i]] = newids[i]
    # for i in corresp:
    #     print i,corresp[i]
    #     sl(0.2)
    handle.close()

change_ade_geneID(infile)