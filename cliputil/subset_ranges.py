"""
Subsets bed/bedgraphs to a genomic range.
Usage: -r "I" -b <bed folder> -o <output_dir>
# -u <bedgraph_unnorm> -n <bedgraph_norm> 
Valid ranges:
"I"
"I:99-200"
"I:99-200/+"

Defaults:
bed_folder = 'bed_collapsed/'
top = os.path.dirname(os.path.realpath(bed_folder).rstrip('/'))
bedgraph_unnorm = os.path.join(top, 'bedgraph_unnorm')
bedgraph_norm = os.path.join(top, 'bedgraph_unnorm')
output_dir = os.path.join(top, 'subset')
"""
import os
import sys
import argparse
import re
import glob
#import config


def parse_range(iv):
    print "> " + str(iv)
    iv = re.sub('\s', '', iv)
    iv = re.sub("""["']""", '', iv)
#    iv = iv.rstrip('"').lstrip('"')
#    iv = iv.rstrip("'").lstrip("'")
    iv = iv.split(':')
    def err():
        print "Invalid range: {0}".format(iv)
        sys.exit()
    if len(iv) == 0: err()
    elif len(iv) == 1:
        iv = iv[0]
        if re.match('[a-zA-Z]+', iv):
            return (iv, None, None, None)
        elif re.match('[a-zA-Z]+/[-\+]', iv):
            return (iv.split('/')[0], None, None, iv.split('/')[1])
        else: err()
    elif len(iv) == 2:
        chrom = iv[0]
        _r = iv[1].split('/')
        if len(_r) > 1:
            strand = _r[1]
        else: strand = None
        positions = _r[0].split('-')
        if len(positions) != 2: err()
        return (chrom, int(positions[0]), int(positions[1]), strand) 


def make_bed_range_check_function(iv):
    if iv[0] is None:
        print "Error: {0} bad range.".format(iv)
        sys.exit()
    if all([(x is not None) for x in iv]):
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            if s[5] != iv[-1]: return False  # Strand.
            if (iv[1] <= int(s[1]) <= iv[2]): return True
            if (iv[1] <= int(s[2]) <= iv[2]): return True
            return False
        return in_range
    elif all([(x is not None) for x in iv[:-1]]):
        # Unstranded, has positions.
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            if (iv[1] <= int(s[1]) <= iv[2]): return True
            if (iv[1] <= int(s[2]) <= iv[2]): return True
            return False
        return in_range
    elif (iv[1] is None) and (iv[2] is None) and (iv[3] is not None):
        # Stranded, no positions.
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            if s[5] != iv[-1]: return False  # Strand.
            return False
        return in_range
    elif (iv[1] is None) and (iv[2] is None) and (iv[3] is None):
        # Unstranded, no positions.
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            return False
        return in_range
    print "Error: {0} bad range.".format(iv)
    sys.exit()


def make_bedgraph_range_check_function(iv):
    if iv[0] is None:
        print "Error: {0} bad range.".format(iv)
        sys.exit()
    if all([(x is not None) for x in iv]):
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            if (iv[1] <= int(s[1]) <= iv[2]): return True
            if (iv[1] <= int(s[2]) <= iv[2]): return True
            return False
        return in_range
    elif all([(x is not None) for x in iv[:-1]]):
        # Unstranded, has positions.
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            if (iv[1] <= int(s[1]) <= iv[2]): return True
            if (iv[1] <= int(s[2]) <= iv[2]): return True
            return False
        return in_range
    elif (iv[1] is None) and (iv[2] is None) and (iv[3] is not None):
        # Stranded, no positions.
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            return False
        return in_range
    elif (iv[1] is None) and (iv[2] is None) and (iv[3] is None):
        # Unstranded, no positions.
        def in_range(s):
            if s[0] != iv[0]: return False  # Chrom
            return False
        return in_range
    print "Error: {0} bad range.".format(iv)
    sys.exit()

def subset_bed(bed_fname, iv, output_fname):
    outli = ''
    fc = make_bed_range_check_function(r)
    with open(bed_fname, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            if fc(s): outli += li
    print output_fname
    with open(output_fname, 'w') as f:
        f.write(outli)

def subset_bedgraph(bedgraph_fname, iv, output_fname):
    outli = ''
    fc = make_bedgraph_range_check_function(r)
    with open(bedgraph_fname, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            if fc(s): outli += li
    print output_fname
    with open(output_fname, 'w') as f:
        f.write(outli)

def parse_ini_file(ini_file, output_dir, out_ini_file):
    outli = ''
    for li in open(ini_file, 'r').readlines():
        s = li.rstrip('\n').split('\t')
        if re.search('\Atop:', s[0]):
            outli += 'top: {0}\n'.format(os.path.realpath(output_dir))
        else:
            outli += li
    with open(out_ini_file, 'w') as f: f.write(outli)


def _mk(output_dir):
    if not os.path.exists(output_dir): os.system('mkdir ' + output_dir)
if __name__ == '__main__':
    p = argparse.ArgumentParser(__doc__)
    p.add_argument(
        '-b', '--bed', help='Directory of bed files.')
#    p.add_argument('-c', '--config', help='config.ini file',
#                   default='auto.ini')
    p.add_argument('-r', '--range', help='Genomic range, e.g. "I:1-9/+"')
    p.add_argument('-o', '--output', help='(Optional) A directory to output \
peaks with coverages added. If not set, files are edited in place.',
                   default=None)
    args = p.parse_args()
    args.bed = os.path.realpath(args.bed)
    r = parse_range(args.range)
    bed_folder = args.bed #'bed_collapsed/'
    top = os.path.dirname(os.path.realpath(bed_folder).rstrip('/'))
    print top
    bedgraph_unnorm = os.path.join(top, 'bedgraph_unnorm')
    bedgraph_norm = os.path.join(top, 'bedgraph_unnorm')
    output_dir = os.path.join(top, 'subset')
    _mk(output_dir)
    bed_output = os.path.join(output_dir, 'bed_collapsed/')
    bedgraph_norm_output = os.path.join(output_dir, 'bedgraph_norm')
    bedgraph_unnorm_output = os.path.join(output_dir, 'bedgraph_unnorm')
    for d in [bed_output, bedgraph_norm_output, bedgraph_unnorm_output]:
        _mk(d)
    for bed in glob.glob(bed_folder + '/*.bed'):
        subset_bed(bed, r, bed_output + '/' + os.path.basename(bed))
    for bed in glob.glob(bedgraph_norm + '/*.wig')
        subset_bedgraph(
            bed, r, bedgraph_norm_output + '/' + os.path.basename(bed))
    for bed in glob.glob(bedgraph_norm + '/*.wig'):
        subset_bedgraph(
            bed, r, bedgraph_unnorm_output + '/' + os.path.basename(bed))
    ini_file = os.path.join(top, 'auto.ini')
    if not os.path.exists(ini_file):
        ini_file = glob.glob(top + '/*.ini')
        if len(ini_file) == 1: ini_file = ini_file[0]
        else: ini_file = None
    if (ini_file is not None) and (os.path.exists(ini_file)):
        parse_ini_file(
            ini_file,
            output_dir,
            os.path.join(output_dir, 'auto.ini'))
    #lib = config(args.config)

