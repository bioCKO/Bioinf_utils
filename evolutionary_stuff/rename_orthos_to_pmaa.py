__author__ = 'mjohnpayne'

import sys

infile = open(sys.argv[1],'r').read()

pos = infile.index('>Pma|')
pmaa = infile[pos+5:pos+17].strip('\n')

outfile = open('/Volumes/MP_HD/Pm_Ts_Tf_Pf_comparison/orthomcl_run/orthomcl_data/CDS_files/ortho_group_cds_fastas/cds_orthos_pm_names/' + pmaa + '_4species_orthos.fasta','w')

outfile.write(infile)

outfile.close()