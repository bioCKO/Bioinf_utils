from Bio import SeqIO

in_gb = raw_input('input genbank file: ').strip(' ')
out_fasta = raw_input('output fasta path: ').strip(' ')
SeqIO.convert(in_gb, "genbank", out_fasta, "fasta")
