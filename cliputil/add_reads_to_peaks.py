"""
Adds read information to peaks files from bedgraphs.
Edits files in place, unless -o is used.

Usage: -i <top level dir> -c config.ini

Optional: -o output directory.

For each peaks file, the config.ini must have the following:
peaks_file_basename: bedgraph_file.wig

"""

import HTSeq
import re
import sys
import os
import pandas
import numpy as np
from __init__ import *

from peaks import peaks

def up3(p):
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(p))))

def up2(p):
    return os.path.dirname(
        os.path.dirname(os.path.dirname(p)))

def guess_bedgraph(path, normed=True):
    only_replicates = True
    up = up3(path)
    if not os.path.exists(up):
        print "no %s" % up
        return False
    print "Expecting bedgraphs_unnorm folder in {0}".format(up)
    if normed:
        bgn = os.path.join(up, 'bedgraphs_norm')
    else:
        bgn = os.path.join(up, 'bedgraphs_unnorm')
    if not os.path.exists(bgn):
        print "no %s" % bgn
        return
    else:
        print "Found {0}.".format(bgn)
    bedgraphs = set()
    for g in glob.glob(bgn + '/*_+.wig'):
        print g
        guess = g.partition('_+.wig')[0]
        pos, neg = (guess + '_+.wig', guess + '_-.wig')
        if any([not os.path.exists(_x) for _x in [pos, neg]]):
            print "no %s" % g
            return
        b = os.path.basename(pos)
        if re.search('\Aall_', b): continue
        if re.search('\Acontrol_n2_[-\+].wig', b): continue
        if re.search('old_fbf', pos):
            bedgraphs.add((pos, neg))
            continue
        if only_replicates and (
            #(re.match('control_[^_]*_\+.wig', b) or \
            (
             re.match('exp_fbf\w_[^_]*_\+.wig', b))):
            print "Skipping %s" % (b)
            continue
        bedgraphs.add((pos, neg))
    return bedgraphs
        
def run(args):
    rep = {}
    peak_objs = []
    for subdir, dirs, files in os.walk(args.input):
        #print files
        file_path = os.path.join(subdir, 'null_hyp_4.txt')
        #if not re.search('old_fbf1_to_fbf2_n2', file_path): continue
        ga = {}
        if not os.path.exists(file_path):
            pass
            #print "%s does not exist" % file_path
        if os.path.exists(file_path):
            if not re.search('combined_', file_path): continue
            print file_path
            try:
                rep[file_path] = peaks(
                    file=file_path, name=os.path.basename(file_path))
                print "peaks obj\n--\n{0} --".format(rep[file_path])
            except:
                print "\n".join(['-' * 7 for x in range(7)]) 
                print "\nFailure to create peaks object from {0}.".format(
                    file_path)
                continue
            print "-------\n{0}\n---\n".format(rep[file_path])
            bedgraphs = guess_bedgraph(file_path)
            _p = peaks(file=file_path)
            if not add_reads_from_bedgraphs(_p, bedgraphs):
                print "looking for bedgraphs {0}".format(bedgraphs)
                continue
            bedgraphs_unnorm = guess_bedgraph(file_path, normed=False)
            if not add_reads_from_bedgraphs(_p, bedgraphs_unnorm, normed=False):
                continue
            _p.write_table('test.txt')
            peak_objs.append(_p)
    print "Loaded {0} peak files.".format(len(peak_objs))
    return peak_objs


def add_reads_from_bedgraphs(_p, bedgraphs, normed=True):
    ga = {}
    if (bedgraphs is None) or (
        not bedgraphs): return False
    for pos, neg in bedgraphs:
        print pos
        print neg
        k = os.path.basename(pos).partition('_+.wig')[0]
        print k
        if normed:
            ga[k] = get_a_bedgraph(pos, neg)
        else:
            ga['unnorm_' + k] = get_a_bedgraph(pos, neg)
        _p.bedgraph_filenames[k] = [pos, neg]
    for k in ga:
#        if normed:
        _p.add_reads(ga=ga[k], name=k)
#        else:
#            print 'Unnormed: ' + k
#            _p.add_reads(ga=ga[k], name='unnorm_' + k)
#            print _p.data['unnorm_' + k]
    _p.bedgraphs = list(bedgraphs)
    if normed:
        _p.set_sum(
            to_sum=[k for k in ga if re.match('exp_.*', k)],
            summed_col='exp_reads')
        _p.set_sum(
            to_sum=[k for k in ga if re.match('control_.*', k)],
            summed_col='control_reads')
        _p.set_ratio(col1='exp_reads', col2='control_reads',
                     ratio_col='ratio')
    else:
        print '\n unnormed \n'
        _p.set_sum(
            to_sum=[k for k in ga if re.match('unnorm_exp_.*', k)],
            summed_col='unnorm_exp_reads')
        _p.set_sum(
            to_sum=[k for k in ga if re.match('unnorm_control_.*', k)],
            summed_col='unnorm_control_reads')
        _p.set_ratio(col1='unnorm_exp_reads', col2='unnorm_control_reads',
                     ratio_col='unnorm_ratio')
    return True


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument(
        '-i', '--input', help='Directory of peaks files, which are searched \
recursively.')
    p.add_argument('-c', '--config', help='config.ini file',
                   default='auto.ini')
    p.add_argument('-o', '--output', help='(Optional) A directory to output \
peaks with coverages added. If not set, files are edited in place.',
                   default=None)
    args = p.parse_args()
    args.input = os.path.realpath(args.input)
    lib = config(args.config)
    run(args)