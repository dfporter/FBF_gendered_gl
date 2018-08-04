"""

Adds several columns of information to input files.
Outputs to a new directory or edits the input file.
Unclear which is better.

"""
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
import subset_peaks_with_fbe
#from peak_locations import locate_in_gene, get_distance


# Duplicated from peak_locations.py.
def locate_in_gene(gtf_sep_cols, combined):
    print("Locating peaks...")
    for index, peak_row in combined.iterrows():
        try:
            if not index % 100:
                print("Locating peak %i..." % index)
        except:
            pass
        gene = str(combined.loc[index, 'gene_name'])
        rows = gtf_sep_cols[gtf_sep_cols['gene_name']==gene]
        if len(rows) == 0:
            combined.loc[index, 'location'] = "Unknown"
            continue
        rows = rows[rows['2']=='CDS']
        if len(rows) == 0:
            combined.loc[index, 'location'] = "ncRNA"
            continue
        gene_left = min(rows['3'])
        gene_right = max(rows['4'])
        gene_strand = rows['6']
        dist = get_distance((gene_left, gene_right),
                            (combined.loc[index, 'left'],
                             combined.loc[index, 'right']))
        if combined.loc[index, 'strand'] == '-':
            dist = -1 * dist
        combined.loc[index, 'dist_to_CDS'] = dist
        if not dist:
            combined.loc[index, 'location'] = "CDS"
        if dist < 0:
            combined.loc[index, 'location'] = '''5'UTR'''
        if dist > 0:
            combined.loc[index, 'location'] = '''3'UTR'''

# Duplicated from peak_locations.py
def get_distance(geneiv, peakiv):
    if geneiv[1] < peakiv[0]:
        return peakiv[0] - geneiv[1]
    if peakiv[1] < geneiv[0]:
        return peakiv[1] - geneiv[0]
    if geneiv[0] < peakiv[0] < geneiv[1]:
        if peakiv[1] < geneiv[1]:
            return 0  # Entirely in CDS.
        else:
            in_cds = geneiv[1] - peakiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[1]
    if geneiv[0] < peakiv[1] < geneiv[1]:
        if peakiv[0] > geneiv[0]:
            return 0  # Entirely in CDS.
        else:
            in_cds = peakiv[1] - geneiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[0]
    return 0

def write_subset_of_columns(peaks, top_level_dir, label):
    top_level_dir = os.path.dirname(top_level_dir)
#    peaks_w_fbe = peaks[peaks['has_fbe']==1]
    if not os.path.exists('with_fbe_%s' % top_level_dir):
        os.system('mkdir with_fbe_%s' % top_level_dir)
    peaks_info_subset = peaks[[
        '0', '1', '2', '3', '4', '5', '6', '7', '8',
        'chrm', 'gene_name', 'height', 'left', 'right',
        'strand', 'transcript_id', 'transcript_name',
        'seq', 'has_fbe', 'number_of_fbes',
        'minus_one_c', 'minus_two_c', 'minus_one_or_two_c',
        'is_rip_chip_target', 'dist_to_CDS', 'location']]
    out_dir = 'with_info_%s' % top_level_dir
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
    peaks_info_subset.to_csv(
        'with_info_%s/fewer_columns_%s' % (top_level_dir, label), sep='\t')
    peaks.to_csv(
        'with_info_%s/%s' % (top_level_dir, label), sep='\t')

def combine_peaks(replicates):
#    print 'passed %s' % str(replicates)
    combined = {}  # Key: (chrm, left, strand) tuple.
    first_rep_name = list(replicates.keys())[0]
    as_single_df = replicates[first_rep_name]
    as_single_df = pandas.concat([replicates[rep][:] for rep in replicates], ignore_index=True)
    print("After append: size of combined df: %i." % len(as_single_df))
    for index, row in as_single_df.iterrows():
        if index % 2000 == 1:
            print("Combining %i/%i." % (index, len(as_single_df)))
        gene = as_single_df.loc[index, 'gene_name']
        iv = [as_single_df.loc[index, 'chrm'], as_single_df.loc[index, 'left'],
              as_single_df.loc[index, 'right'], as_single_df.loc[index, 'strand']]
        overlapping_peaks = overlapping(iv, as_single_df)
        if len(overlapping_peaks) > 4:
            consensus_dict, consensus_row = consensus_peak(overlapping_peaks)
            tup = (consensus_dict['chrm'],
                   consensus_dict['left'],
                   consensus_dict['strand'])
            if tup not in combined:
                combined[tup] = consensus_dict
        #if len(combined) > 10:
        #    return combined
    combined_df = pandas.DataFrame.from_dict(combined, orient='index')
    print(combined_df.head())
    return combined_df

def consensus_peak(peaks):
    rows = peaks[peaks['height']==max(peaks['height'])]
    return (rows.iloc[0].to_dict(), rows.iloc[0])

def overlapping(iv, comb_df):
#    print 'input %s' % str(comb_df.head())
    overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2]) & (comb_df['right']>=iv[1])]
    return overlap
#        print "just left"
#        print overlap.head()
#        print overlap['right']
#        overlap = overlap[overlap['right']>=iv[1]]
#        print "with right:"
#        print overlap
#    except:
#        overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2])]
#        print "Found for left with %s: %s" % (
#            str(iv), str(overlap.head()))
#        return overlap
#    overlap = comb_df[(comb_df['chrm']==iv[0]) & (comb_df['strand']==iv[3]) & (comb_df['left']<=iv[2]) & (iv[1]<=comb_df['right'])]
#    return overlap

def add_info(peaks, filename, rip_targets, gtf_sep_cols, top_level_dir):
        # Add genomic sequence in the peak.
        subset_peaks_with_fbe.get_sequences(peaks)
        # Get yes/no FBE, -1/-2 C, number of FBEs.
        subset_peaks_with_fbe.score_binding_site(peaks)
        # Overlaps with the top 1350 unique RIP-chip targets?
        compare_with_ripchip.add_column_of_overlap(
            peaks, rip_targets)
        # Locate each peak in the gene.
        locate_in_gene(gtf_sep_cols, peaks)
        subset_peaks_with_fbe.write_subset_with_fbe(
            peaks, top_level_dir, label=os.path.basename(filename))
        write_subset_of_columns(
            peaks, top_level_dir, label=os.path.basename(filename))


if __name__ == '__main__':
    top_level_dir = sys.argv[1]
    combined = {}
    #top_level_dir = 'pypeaks_fdr1_negip_local/'
    #peaks = replicates['pypeaks_fdr1_negip_local/fbf1_tccg.txt']
    rip_targets = compare_with_ripchip.get_ripchip_targets()
    gtf_sep_cols = pandas.read_csv('lib/gtf_with_names_column.txt', sep='\t')
    for filename in glob.glob(top_level_dir + '/combined*.txt'):
        continue
        print(filename)    
        combined[filename] = pandas.read_csv(filename, sep='\t')
        add_info(combined[filename], filename,
                 rip_targets, gtf_sep_cols,
                 top_level_dir)
    # Using all replicates.
    replicates = {}
    for filename in glob.glob(top_level_dir + '/fbf*'):
        replicates[filename] = pandas.read_csv(filename, sep='\t')
    combined_df = combine_peaks(replicates)
    add_info(combined_df, 'both_fbfs.txt', rip_targets, gtf_sep_cols,
             top_level_dir)
             
