__author__ = 'mjohnpayne'
from Bio import SeqIO
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
import sys


max_len = 6000000
telomere_length = 40000

chr_diagram = BasicChromosome.Organism()
chr_diagram.page_size = (60*cm, 21*cm)


chr_diagram.
chr_diagram.draw('/Users/mjohnpayne/Documents/PhD/Chromosome_plots/legend_test.pdf','legend_test')