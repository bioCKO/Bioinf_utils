# -*- coding: utf-8 -*-
__author__ = 'mjohnpayne'

def parse_xdna(file):
    file = open(file,'r').readlines()
    seq = file[0][112:]
    allow = ['A','T','G','C']
    end = 0
    for i in range(len(seq)):
        if seq[i] not in allow:
            end = i
            break
    seq = seq[:end]
    return seq
