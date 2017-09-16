import os
import ConfigParser
from peaksList import peaksList
from volcanoMaker import volcanoMaker
import gtf
import rc
from create_gtf_of_longest_txpt_per_gene import *
from sbGene import *


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([basecomplement[base] for base in list(s)])


def rc(s): return complement(s[::-1])

def mkdir(s):
    assert(type(s) == type(''))
    if os.path.isfile(s):
        s = os.path.dirname(s)
    steps = 0
    to_make = set([])
    while (steps < 5):
        if os.path.exists(s): break
        else:
            to_make.add(s)
            s = os.path.dirname(s)
    to_make = sorted(list(to_make), key=lambda x: len(x))
    for path in to_make:
        os.system('mkdir {0}'.format(path))

def get_sequences(lib=None):
    fasta_filename = '/scratch/indexes/WS235.fa'
    if lib is not None and ('fasta' in lib):
        fasta_filename = lib['fasta']
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    rc_sequences = dict(
        (p.name.split(' ')[0], rc(p.seq)) for p in HTSeq.FastaReader(fasta_filename))
    chr_lens = dict(
        [(name, len(sequences[name])) for name in sequences])
    return (sequences, rc_sequences, chr_lens)


def get_a_bedgraph(plus_file, minus_file):
    import HTSeq
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
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

def load_gtf(lib=None, add_ncrna=False, make_genes_using_peaks_file=None):
    if lib is None:
        lib = {'gtf': 'lib/gtf_with_names_column.txt',
               'gtf_one_txpt_per_gene': 'lib/gtf_one_txtp_per_gene.txt'}
    if not os.path.exists(lib['gtf_one_txpt_per_gene']):
        print "Creating a GTF-like of the longest txpt per gene..."
        create_gtf_of_longest_txpt_per_gene(
            in_filename=lib['gtf'],
            out_filename=lib['gtf_one_txpt_per_gene'])
    print "Getting sequences..."
    (sequences, rc_sequences, chr_lens) = get_sequences(lib)
    print "Getting GTF..."
    if make_genes_using_peaks_file is not None:
        g = gtf(lib['gtf'], chr_lens)
        g.sequences = sequences
        g.rc_sequences = rc_sequences
        peaks = pandas.read_csv(
            make_genes_using_peaks_file, sep='\t', index_col=False)
        peaks = g.flip_minus_strand_peak_features(peaks)
        g.flip_minus_strand_features()
        g.make_genes_using_input_peaks_file(
            peaks, sequences, rc_sequences, lib['gtf'])
    else:
        g = gtf(lib['gtf_one_txpt_per_gene'], chr_lens)
        g.sequences = sequences
        g.rc_sequences = rc_sequences
        g.flip_minus_strand_features()
        g.make_genes(sequences, rc_sequences)
    if add_ncrna:
        g.make_ncrna_sbgenes(sequences, rc_sequences)
    return g

def add_col_of_longest_txpt_containing_peaks(peaks_df, lib=None):
    if lib is None:
        lib = {'gtf': 'lib/gtf_with_names_column.txt',
               'gtf_one_txpt_per_gene': 'lib/gtf_one_txtp_per_gene.txt'}
    (sequences, rc_sequences, chr_lens) = get_sequences(lib)
    g = gtf(lib['gtf'], chr_lens)
    g.sequences = sequences
    g.rc_sequences = rc_sequences
    g.flip_minus_strand_features()
    
def printProgress (
    iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 35):
    """http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write(
        '%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")


def config(
    filepath="/groups/Kimble/Aman Prasad/redo_fbf/analysis/src/config.ini"):
    Config = ConfigParser.ConfigParser()
    Config.read(filepath)
    lib = ConfigSectionMap('library', Config)
    return lib

 
def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

