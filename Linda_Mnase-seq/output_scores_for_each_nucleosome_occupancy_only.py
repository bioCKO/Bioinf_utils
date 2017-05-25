__author__ = 'mjohnpayne'

## occupancy identification and comparison
# 1 Ident differential nuc positions in all comparisons
# 2 scores for these in all other conditions even if not sig different
# 3 output nucleosome id and position relative to genes
# 4 use outputs to perform clustering


# 1 use 0.05 cutoff for occupancy changes

types = ["9H","11H","15H","17H"]

differetial_in = []

occupancy_in = []

ref_positions = "/Volumes/MP_HD/Linda_MNase_Seq/DANPOS_4_way/danpos_raw_output/reference_positions.xls"

