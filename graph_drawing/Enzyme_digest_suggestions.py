# -*- coding: utf-8 -*-
#! /usr/bin/env python-32
__author__ = 'mjohnpayne'


import re

from Bio import SeqIO
from Bio.Seq import Seq
import sys
sys.path.append('/Users/mjohnpayne/PycharmProjects/Bioinf_utils/')
import simGel
import parse_xdna
from Bio import Restriction


in_RE = open('/Users/mjohnpayne/Documents/PhD/HDA_RE_list.txt', 'r').read().split()

inp = raw_input('path to xdna file: ')


min_dist = int(raw_input('minimum distance between fragments(bp): '))
min_size = int(raw_input('minimum size of fragments(bp): '))

def run(sequence):
    sequence = parse_xdna.parse_xdna(sequence)
    inseq = Seq(sequence)
    lin = False

    info = {}

    for i in Restriction.AllEnzymes:
        frags = []
        if str(i) in in_RE:
            if lin == False:
                sites = sorted(i.search(inseq, linear=False))
                frags = [sites[r] - sites[r - 1] for r in xrange(1, len(sites))]
                if len(sites) > 0:
                    if len(sites) > 1:
                        circ_frag = len(inseq) - sites[-1] + sites[0]
                        frags.append(circ_frag)
                    elif len(sites) == 1:
                        frags.append(len(str(inseq)))
                    frags = sorted(frags)
                    dist = [frags[r] - frags[r - 1] for r in xrange(1, len(frags))]
                    # print i
                    # print frags
                    # print dist
                    if len(dist) > 0:
                        info[i] = [sites, frags, dist]
                    elif len(dist) == 0:
                        info[i] = [sites, frags, [0]]
                else:
                    continue
                    # print i
                    # print 'No Sites'
                    # print len(str(inseq.seq))
            elif lin == True:
                sites = sorted(i.search(inseq, linear=True))
                frags = sorted([sites[r] - sites[r - 1] for r in xrange(1, len(sites))])
                frags = sorted([sites[0]] + frags + [len(str(inseq)) - sites[-1]])
                if len(sites) > 0:
                    dist = [frags[r] - frags[r - 1] for r in xrange(0, len(frags))]
                    if len(dist) > 0:
                        # print i
                        # print frags
                        # print dist
                        info[i] = [sites, frags, dist]
                elif len(sites) == 0:
                    frags = sorted([sites[0]] + [len(str(inseq)) - sites[0]])
                    dist = [frags[r] - frags[r - 1] for r in xrange(0, len(frags))]
                    # print i
                    # print len(str(inseq.seq))
                    # print sites
                    # print frags
                    info[i] = [sites, frags, dist]
    output = []
    count = 0
    used = []
    enzymes = ['']
    for i in sorted(info.keys()):
        if int(min(info[i][1])) >= min_size and int(min(info[i][2])) >= min_dist:
            output += [[(info[i][1][r],100) for r in range(len(info[i][1]))]]
            used.append(str(i) + ' (' + ','.join(map(str, sorted(info[i][1], reverse=True))) + ')')
            enzymes.append(str(i))
    print "Lambda Ladder  " + '  '.join(map(str,used))
    return simGel.input(output,len(output),inp.replace('.xdna','_digests.png'),enzymes)

run(inp)