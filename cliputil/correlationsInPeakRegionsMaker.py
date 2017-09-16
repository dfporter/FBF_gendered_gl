from __future__ import division
import sys
import glob
import os
import re
import pandas
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
import pickle
import argparse

sys.path.insert(0, '/Network/Servers/file.biochem.wisc.edu/Volumes/BioSAN/Users/dfporter/.local/lib/python2.6/site-packages')
# Need HTSeq version 0.6, hence the sys.path alteration.
sys.path.insert(
    0, '/groups/Kimble/Common/fbf_celltype/redo_fbf/analysis/src/')
import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys

import add_height
import scatterplot_correlation_by_wig


#import normalize_bedgraph
import argparse


def build_ga(file_list):
    ga = HTSeq.GenomicArray('auto', stranded=True)
    for infile in file_list:
        print "Loading {0}...".format(infile)
        ga = add_a_bedgraph(ga, infile + '_+.wig', infile + '_-.wig')
    return ga

def mk_bedgraphs(fbf1oo, fbf2oo, fbf1sp, fbf2sp, fbf1old, fbf2old):
    output_bedgraph_unnorm = 'bedgraphs_unnorm/corr/'
    output_bedgraph_norm = 'bedgraphs_norm/corr/'
    if not os.path.exists(output_bedgraph_unnorm):
        os.system('mkdir ' + output_bedgraph_unnorm)
    if not os.path.exists(output_bedgraph_norm):
        os.system('mkdir ' + output_bedgraph_norm)
    ga = {}
    ga['fbf1_oo'] = build_ga(fbf1oo)
    ga['fbf2_oo'] = build_ga(fbf2oo)
    ga['fbf1_sp'] = build_ga(fbf1sp)
    ga['fbf2_sp'] = build_ga(fbf2sp)
    ga['fbf1_old'] = build_ga(fbf1old)
    ga['fbf2_old'] = build_ga(fbf2old)
    for k in ga:
        outname = "{d}/{b}".format(
            d=output_bedgraph_unnorm,
            b=k)
        print "Creating a bedgraph {c} from {a}...".format(
            c=outname, a=k)
        outname_plus = outname + '_+.wig'
        ga[k].write_bedgraph_file(outname_plus, strand='+')
        outname_minus = outname + '_-.wig'
        ga[k].write_bedgraph_file(outname_minus, strand='-')

#    normalize_bedgraph.normalize_wig(input_bed, output_bedgraph_unnorm,
#                                     output_bedgraph_norm)

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

def get_args():
    parser = argparse.ArgumentParser(
        description='''Add sequence, has_fbe, rip_column, ect.
        info columns to a given directory -d of peaks files.
        Can also combine replicates and output the combined list.
        Will apply a normalization ratio of ~1.6 to FBF-2 reads.'''
    )
    parser.add_argument('-d', '--directory', dest='directory',
                        default='combined_filtered/',
                        help='Input directory of peaks files.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    fbf1oo = ['exp_fbf1_oo_lane2_rt1',
     'exp_fbf1_oo_lane2_rt6',
     'exp_fbf1_oo_lane2_rt9']
    fbf2oo = ['exp_fbf2_oo_lane2_rt2',
     'exp_fbf2_oo_lane2_rt11',
     'exp_fbf2_oo_lane2_rt13']
    fbf1sp = ['exp_fbf1_sp_lane1_rt1',
     'exp_fbf1_sp_lane1_rt6',
     'exp_fbf1_sp_lane1_rt9']
    fbf2sp = ['exp_fbf2_sp_lane1_rt2',
     'exp_fbf2_sp_lane1_rt13',
     'exp_fbf2_sp_lane1_rt14']
    fbf1old = ['exp_fbf1_CGGA', 'exp_fbf1_GGTT',
     'exp_fbf1_TGGC']
    fbf2old = ['exp_fbf2_CGGA', 'exp_fbf2_GGTT',
               'exp_fbf2_TGGC']        
    for x in [fbf1oo, fbf2oo, fbf1sp, fbf2sp,
              fbf1old, fbf2old,
              ]:
        for i in range(len(x)):
            x[i] = '{0}/{1}'.format('bedgraphs_norm/', x[i])
    #mk_bedgraphs(fbf1oo, fbf2oo, fbf1sp, fbf2sp, fbf1old, fbf2old)
    ga = {}
    text = {}
    ga['combined_fbf1.txt'] = add_a_bedgraph(
        HTSeq.GenomicArray('auto', stranded=True),
        'bedgraphs_unnorm/corr/fbf1_oo_+.wig',
        'bedgraphs_unnorm/corr/fbf1_oo_-.wig')
    ga['combined_fbf2.txt'] = add_a_bedgraph(
        HTSeq.GenomicArray('auto', stranded=True),
        'bedgraphs_unnorm/corr/fbf2_oo_+.wig',
        'bedgraphs_unnorm/corr/fbf2_oo_-.wig')
    print "OO FBF-1 vs FBF-2"
    text["OO FBF-1 vs FBF-2"] = scatterplot_correlation_by_wig.run(
        args, ga, a='oo_fbf1.txt', b='oo_fbf2.txt')
    ga['combined_fbf1.txt'] = add_a_bedgraph(
        HTSeq.GenomicArray('auto', stranded=True),
        'bedgraphs_unnorm/corr/fbf1_sp_+.wig',
        'bedgraphs_unnorm/corr/fbf1_sp_-.wig')
    ga['combined_fbf2.txt'] = add_a_bedgraph(
        HTSeq.GenomicArray('auto', stranded=True),
        'bedgraphs_unnorm/corr/fbf2_sp_+.wig',
        'bedgraphs_unnorm/corr/fbf2_sp_-.wig')
    print "SP FBF-1 vs FBF-2"
    text["SP FBF-1 vs FBF-2"] = scatterplot_correlation_by_wig.run(
        args, ga, a='sp_fbf1.txt', b='sp_fbf2.txt')
    ga['combined_fbf1.txt'] = add_a_bedgraph(
        HTSeq.GenomicArray('auto', stranded=True),
        'bedgraphs_unnorm/corr/fbf1_old_+.wig',
        'bedgraphs_unnorm/corr/fbf1_old_-.wig')
    ga['combined_fbf2.txt'] = add_a_bedgraph(
        HTSeq.GenomicArray('auto', stranded=True),
        'bedgraphs_unnorm/corr/fbf2_old_+.wig',
        'bedgraphs_unnorm/corr/fbf2_old_-.wig')
    print "SP FBF-1 vs FBF-2"
    text["20deg FBF-1 vs FBF-2"] = scatterplot_correlation_by_wig.run(
        args, ga, 'old_fbf1_to_fbf2_n2.txt', 'old_fbf2.txt')
    text["20deg FBF-1 vs 25deg FBF-1"] = scatterplot_correlation_by_wig.run(
        args, ga, 'old_fbf1_to_fbf2_n2.txt', 'oo_fbf1.txt')
    text["20deg FBF-2 vs 25deg FBF-2"] = scatterplot_correlation_by_wig.run(
        args, ga, 'old_fbf2.txt', 'oo_fbf2.txt')

    to_peaks = {"20deg FBF-1 vs FBF-2": (
        'combined_filtered/old_fbf1_to_fbf2_n2.txt',
        'combined_filtered/old_fbf2.txt'),
     "OO FBF-1 vs FBF-2": (
        'combined_filtered/oo_fbf1.txt',
        'combined_filtered/oo_fbf2.txt'),
     "SP FBF-1 vs FBF-2": (
         'combined_filtered/sp_fbf1.txt',
         'combined_filtered/sp_fbf2.txt'),
    "20deg FBF-1 vs 25deg FBF-1": (
        'combined_filtered/old_fbf2.txt',
        'combined_filtered/oo_fbf1.txt'),
    "20deg FBF-2 vs 25deg FBF-2": (
        'combined_filtered/old_fbf1_to_fbf2_n2.txt',
        'combined_filtered/oo_fbf2.txt'),
                }
    corrs = []
    for k in text:
        print text[k]
        spearman = re.search('spearman rho: ([^\n]*)\n', text[k])
        if spearman is not None:
            spearman = spearman.groups()[0]
        pearson = re.search('Pearson R: ([^\n]*)\n', text[k]).groups()[0]
        print spearman
        print pearson
        fbf1_targs = set(
            pandas.read_csv(to_peaks[k][0], sep='\t')['gene_name'].tolist())
        fbf2_targs = set(
            pandas.read_csv(to_peaks[k][1], sep='\t')['gene_name'].tolist())
        joint_targs = fbf1_targs & fbf2_targs
        corrs.append({'Samples':k, 'Spearman rho': spearman,
                      'Pearson R (targets of both FBF-1 and FBF-2)': pearson,
            "% of FBF-2 targets shared with FBF-1": "{0:.2f}".format(100 * len(joint_targs)/len(fbf2_targs)),
            "% of FBF-1 targets shared with FBF-2": "{0:.2f}".format(100 * len(joint_targs)/len(fbf1_targs)),
            "% of FBF-1 or FBF-2 targets shared with FBF-1 and FBF-2": \
                "{0:.2f}".format(
                    100 * len(joint_targs)/len(fbf1_targs | fbf2_targs))
        })
    corrs = pandas.DataFrame(corrs)
    corrs = corrs[['Samples', 'Spearman rho', 'Pearson R (targets of both FBF-1 and FBF-2)',
                   "% of FBF-1 targets shared with FBF-2",
                   "% of FBF-2 targets shared with FBF-1",
                   "% of FBF-1 or FBF-2 targets shared with FBF-1 and FBF-2"]]
    corrs.to_csv('tables/correlations_fbf1_and_fbf2.txt', index=False,
                 sep='\t')
    
