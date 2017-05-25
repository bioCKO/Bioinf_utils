__author__ = 'mjohnpayne'

#### -*- coding: utf-8 -*-
####! /usr/bin/env python-32

import pygame
import math
import os
import sys
import re

from Bio import SeqIO
from Bio.Seq import Seq
import sys
from Bio import Restriction

#   version 1.2.2

# Arbitarary list of gel sizes for easy reference
# Tuples give value of height in pixels and number of lanes
GEL_SIZE = {}#{1: (360, 2),2: (360, 3),3: (360, 4),4: (360, 5),5: (360, 6),6: (360, 7),7: (360, 8),8: (360, 9)}

for i in range(1,100):
    GEL_SIZE[i] = (360,i)

# Approximate values for 1kb ladder
# Note that samples are list of tuples in the form (length , concentration)
LADDER_1KB =[
    (300, 250),
    (500, 150),
    (700, 100),
    (1000, 75),
    (1500, 50),
    (2000, 40),
    (2500, 35),
    (3000, 30),
    (4000, 25),
    (5000, 30),
    (6000, 15),
    (8000, 12),
    (10000, 15),
    (14000, 4),
    (24000, 2)]

Lambda_ladder = [(23130, 4),(9416, 15),(6557, 15),(4361, 25),(2322, 40),(2027, 42),(564, 150)]
pygame.init()
myfont = pygame.font.SysFont("arial", 11)

def gamma(x, t, k=1):
    b = 1.0 / t
    x = x * k
    c = 0.0

    for i in range(0, k):
        c += (math.exp(-b * x) * (b * x) ** i) / math.factorial(i)

    return c

def testColour(c):
    """ Convert integer to colour, which saturates if > 255 """

    if c <= 255:
        return (c, c, c)
    else:
        return (255,255,255)

def display(width, height, image):
    """ Create a Pygame image of dimensions width x height and show image """

    screen = pygame.display.set_mode((width, height))
    screen.blit(image, (60, 30))
    pygame.display.flip()

    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

class Gel:
    """ A Gel contains lanes into which DNA can be loaded and run.
        The Gel is exposured to see the location and intensity of DNA. """

    LANE_WIDTH = 30.
    LANE_HEIGHT = 4
    LANE_MARGIN = 6
    BORDER = 20

    def __init__(self, size, agarose=1):
        self.y = GEL_SIZE[size][0]
        self.x = 2 * Gel.BORDER + GEL_SIZE[size][1] * (Gel.LANE_WIDTH + Gel.LANE_MARGIN) - Gel.LANE_MARGIN

        self.agarose = agarose      # Should test whether value is reasonable (say, 0.5% - 3%)
        self.optimum_DNA_length = 2000 / agarose ** 3   # Best separate based on agarose hole size

        self.samples = []
        self.lanes = []
        self.enz = []

        for n in range(GEL_SIZE[size][1]):
            self.samples.append([])
            self.lanes.append([])

        for lane in self.lanes:
            for n in range(Gel.BORDER, self.y + 3):
                lane.append(0.0)

    def loadSample(self, lane, sample):
        """Add list containing tuple of DNA length and concentrations to lane in gel. """

        for dna in sample:
            strand = {'length': float(dna[0]), 'conc': dna[1] * 1./Gel.LANE_HEIGHT, 'position': Gel.BORDER+1}
            self.samples[lane-1].append(strand)

    def run(self, time=30.0, voltage=80.0):
        """ Move loaded DNA down the gel at a rate dependent on voltage, DNA length and agarose concentration. """

        max_dist = 0.25 * time * voltage
        for sample in self.samples:
            for dna in sample:
                g = gamma( dna['length']/20, int(self.optimum_DNA_length/20) )
                dna['position'] += max_dist * g

    def expose(self, exposure=0.1, aperture=2):
        """Returns an image of the gel with the DNA highlighted"""

        c1, c2 = exposure * 100, exposure * 200
        fill_colour = (c1, c1, c2)
        image = pygame.Surface((self.x+2, self.y+2))
        pygame.draw.rect(image, fill_colour, (1,1, self.x, self.y), 0)

        edge_colour = (140, 140, 150)
        edges = pygame.Surface((self.x+2, self.y+2))
        edges.set_colorkey((0,0,0))
        edges.set_alpha(120)
        pygame.draw.rect(image, edge_colour, (0,0, self.x+2, self.y+2), 1)

        self._findDNAConcentrations(c1)

        x = Gel.BORDER + 1
        ecount = 0
        for lane in self.lanes:

            for position in range(len(lane)):
                if lane[position] > 0:

                    brightness = lane[position] * exposure / aperture + c1
                    colour1 = testColour(brightness)
                    colour2 = testColour(brightness*0.6)
                    colour3 = testColour(brightness*0.3)

                    pygame.draw.line(image, colour3, (x-1, position+Gel.BORDER),(x+Gel.LANE_WIDTH+1, position+Gel.BORDER))
                    pygame.draw.line(image, colour2, (x, position+Gel.BORDER),(x+Gel.LANE_WIDTH, position+Gel.BORDER))
                    pygame.draw.line(image, colour1, (x+1, position+Gel.BORDER),(x+Gel.LANE_WIDTH-1, position+Gel.BORDER))
            label = myfont.render(self.enz[ecount], 1, (255,255,255))
            image.blit(label, (x+1, 2))
            pygame.draw.rect(edges, edge_colour, (x, Gel.BORDER, Gel.LANE_WIDTH, Gel.LANE_HEIGHT), 1)
            x += Gel.LANE_WIDTH + Gel.LANE_MARGIN
            ecount += 1

        image.blit(edges, (0,0))
        return image

    def _findDNAConcentrations(self, background):
        """Determines where in the concentration of DNA in every part of the gel"""

        length = len(self.lanes[0])

        for x in range(len(self.samples)):
            for dna in self.samples[x]:
                for y in range(Gel.LANE_HEIGHT-2):
                    pos = int(dna['position']) + y
                    if pos < length-4:
                        # Very crude way to create blurred line
                        self.lanes[x][pos-2] += 0.06 * dna['conc'] * dna['length']
                        self.lanes[x][pos-1] += 0.12 * dna['conc'] * dna['length']
                        self.lanes[x][pos] += 0.2 * dna['conc'] * dna['length']
                        self.lanes[x][pos+1] += 0.12 * dna['conc'] * dna['length']
                        self.lanes[x][pos+2] += 0.06 * dna['conc'] * dna['length']
#    def add_lables(self,samp_lis):


def gelrun(samples,size,path,enz):
    """ An example of how to set up and run a gel. """

    myGel = Gel(size=size+1, agarose=0.88)
    myGel.loadSample(1, Lambda_ladder)
    myGel.enz = enz
    # Samples are list of tuples in the form (length, concentration)
    for i in range(1,len(samples)+1):
        myGel.loadSample(i+1,samples[i-1])
    # sample1 = samples[0]
    # sample2 = samples[1]
    # sample3 = samples[2]

    # Set up gel, giving its size and % agarose
    # myGel = Gel(size=size, agarose=0.88)
    #
    # # Load samples
    # myGel.loadSample(1, Lambda_ladder)
    # myGel.loadSample(2, sample1)
    # myGel.loadSample(3, sample2)
    # myGel.loadSample(4, sample3)

    # Run gel for 45 minutes
    myGel.run(time=15)

    # Take a picture of your gel with a 50ms exposure
    my_image = myGel.expose(exposure=0.06)
    # Display image on pygame screen
    #display(170 + 36*size, 400, my_image)

    # Or save image as a PNG
    pygame.image.save(my_image, path)





def parse_xdna(file):
    file = open(file,'r').readlines()
    seq = file[0][112:]
    allow = ['A','T','G','C','N','a','t','c','g','n']
    end = 0
    for i in range(len(seq)):
        if seq[i] not in allow:
            end = i
            break
    seq = seq[:end]
    return seq





def run(sequence):
    sequence = parse_xdna(sequence)
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
    out_txt = open(inp.replace('.xdna','_digest_sizes.txt'),'w')
    out_txt.write("Lambda Ladder\n" + '\n'.join(map(str,used)))
    out_txt.close()
    return gelrun(output,len(output),inp.replace('.xdna','_digests.png'),enzymes)

in_RE = open('/Users/mjohnpayne/Documents/PhD/HDA_RE_list.txt', 'r').read().split()

inp = raw_input('path to xdna file: ')

min_dist = int(raw_input('minimum distance between fragments(bp): '))
min_size = int(raw_input('minimum size of fragments(bp): '))

run(inp)

while raw_input("Would you like to run another plasmid?(Y or N): ") == 'Y':
    inp = raw_input('path to xdna file: ')
    min_dist = int(raw_input('minimum distance between fragments(bp): '))
    min_size = int(raw_input('minimum size of fragments(bp): '))
    run(inp)