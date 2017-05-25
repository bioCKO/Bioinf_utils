'''
Author: Travis Poulsen
Date: 09 Feb. 2013

Generate some basic statistics/information about an assembly.
This script will generate a file named 'assembly_stats.txt' that will
include number of contigs, min/max/mean contig lengths, and the N50 for
the assembly. It also generates a histogram of contig lengths as a png.

Needs: numpy, scipy, matplotlib, and Biopython

Expected input: A fasta file of assembled contigs.

Options:--no_hist: generate the stats only, no histogram.
        -v, --verbose: print the stats to the screen as well as to the file.
'''
import os, sys, numpy, argparse
from Bio import SeqIO

# Main function: generates content if called as a stand alone.
def main():
        parser = argparse.ArgumentParser(description='Generate basic statistics about a collection of contigs.')
        parser.add_argument('input', help='fasta file containing the contigs for analysis.')
        parser.add_argument('-v', '--verbose', action='store_true', help = 'Print statistics to screen as well as save to an output file.')
        parser.add_argument('--no_hist', action='store_true', help = 'Generate stats but do not create histogram.')
        args = parser.parse_args()

        file_in = args.input
        with open(args.input, 'r') as seq:
                sizes = [len(record) for record in SeqIO.parse(seq, 'fasta')]
        min_contig = min(sizes)
        max_contig = max(sizes)
        avg_contig = numpy.mean(sizes)
        num_contig = len(sizes)
        with open(file_in+'_assembly_stats.txt', 'w') as handle:
                handle.write('Number of contigs:\t%i\n' % num_contig)
                handle.write('N50:\t%i\n' % int(getN50(sizes)))
                handle.write('Mean contig length:\t%.2f\n' % avg_contig)
                handle.write('Minimum contig length:\t%i\n' % min_contig)
                handle.write('Maximum contig length:\t%i\n' % max_contig)
        if args.verbose == True:
                print 'Number of contigs:\t%i' % num_contig
                print 'N50:\t%i' % int(getN50(sizes))
                print 'Mean contig length:\t%.2f' % avg_contig
                print 'Minimum contig length:\t%i' % min_contig
                print 'Maximum contig length:\t%i' % max_contig
        if args.no_hist == False:
                plot_assembly(sizes, file_in, min_contig, max_contig, avg_contig, num_contig)

# Plot the histogram of contig lengths with the specified axis labels and
# other information included in the graph.
def plot_assembly(sizes, file_in, min_contig, max_contig, avg_contig, num_contig):
        import matplotlib
        matplotlib.use('Agg')
        import pylab
        pylab.hist(sizes, bins=50)
        pylab.title('%i %s sequences\nLengths %i to %i, Average contig length: %.2f' % (num_contig, file_in, min_contig, max_contig, avg_contig))
        pylab.xlabel('Sequence length (bp)')
        pylab.ylabel('Count')
        pylab.savefig(file_in+'_contig_histogram.png')
        os.chdir('..')

# Get the N50 of the contigs. This is the sequence length at which point
# half of the bases in the entire assembly are contained in contigs of a
# smaller size.
def getN50(sizes):
        bases = []
        for read in sizes:
                for i in range(read):
                        bases.append(read)
        return numpy.median(bases)


if __name__ == '__main__':
        main()