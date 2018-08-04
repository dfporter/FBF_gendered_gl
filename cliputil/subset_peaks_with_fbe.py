import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
#from scipy.optimize import minimize_scalar
#import pysam
import matplotlib
import sys
#from matplotlib_venn import *
import os
import HTSeq
import compare_with_ripchip


def complement(s):
    #basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)

def rc(s):
    s = s[::-1]
    s = complement(s)
    return s

def score_binding_site(peaks):
    for index, peak_row in peaks.iterrows():
        #print 'searching %s' % peaks.loc[index, 'seq']
        if re.search('tgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'has_fbe'] = 1
            peaks.loc[index, 'number_of_fbes'] = len(
                re.findall('tgt\w\w\wat', peaks.loc[index, 'seq']))
        else:
            peaks.loc[index, 'has_fbe'] = 0
            peaks.loc[index, 'number_of_fbes'] = 0
        if re.search('ctgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'minus_one_c'] = 1
        else:
            peaks.loc[index, 'minus_one_c'] = 0
        if re.search('c\wtgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'minus_two_c'] = 1
        else:
            peaks.loc[index, 'minus_two_c'] = 0
	
    if peaks.loc[index, 'minus_one_c'] == 1 or peaks.loc[index, 'minus_two_c'] == 1:
        peaks.loc[index, 'minus_one_or_two_c'] = 1
    else:
        peaks.loc[index, 'minus_one_or_two_c'] = 0

def get_sequences(combined):
    #fasta_filename = '/home/dp/Desktop/celegans_genome/wormbase_ws235/c_elegans.WS235.genomic.fa'
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        seq = sequences[chrm][start:end]
        #print "%s:%i-%i: seq %s" % (chrm, start, end, seq)
        if combined.loc[index, 'strand'] == '-':
            seq = rc(seq)
        combined.loc[index, 'seq'] = seq

def write_subset_with_fbe(peaks, top_level_dir, label):
    top_level_dir = os.path.dirname(top_level_dir)
    peaks_w_fbe = peaks[peaks['has_fbe']==1]
    if not os.path.exists('with_fbe_%s' % top_level_dir):
        os.system('mkdir with_fbe_%s' % top_level_dir)
    peaks_w_fbe.to_csv('with_fbe_%s/%s' % (top_level_dir, label), sep='\t')


if __name__ == '__main__':
    top_level_dir = sys.argv[1]
    combined = {}
    rip_targets = compare_with_ripchip.get_ripchip_targets()
    for filename in glob.glob(top_level_dir + '/combined*.txt'):
        print(filename)
        combined[filename] = pandas.read_csv(filename, sep='\t')
        get_sequences(combined[filename])
        score_binding_site(combined[filename])
        compare_with_ripchip.add_column_of_overlap(
            combined[filename], rip_targets)
        write_subset_with_fbe(
            combined[filename], top_level_dir, label=os.path.basename(filename))
