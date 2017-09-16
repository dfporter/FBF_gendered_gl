from __future__ import division
import HTSeq
import os
import sys
import glob

def add_a_bedgraph(ga, plus_file, minus_file):
    with open(plus_file, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
    with open(minus_file, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])
    return ga

def build_ga(file_list):
    ga = HTSeq.GenomicArray('auto', stranded=True)
    for infile in file_list:
        print "Loading {0}...".format(infile)
        ga = add_a_bedgraph(ga, infile + '_+.wig', infile + '_-.wig')
    return ga

def ave_ga(file_list):
    ga = build_ga(file_list)
    norm_ga = HTSeq.GenomicArray('auto', stranded=True)
    denom = len(file_list)
    for iv, score in ga.steps():
        norm_ga[iv] = score/denom
    return norm_ga

if __name__ == '__main__':
    file_list = ['bedgraphs_unnorm/corr/fbf1_old',
                 'bedgraphs_unnorm/corr/fbf2_old']
    norm_ga = ave_ga(file_list)
