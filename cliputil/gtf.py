import HTSeq
import pandas
import os
import sys
import collections
import time

from flocs import *
from sbGene import *

class gtf(object):

    def __init__(self, filename, chr_lens, name='Unnamed', debug=False):
        self.name = name
        self.chr_lens = chr_lens
        self.minus_strand_flipped = False
        self.order = [
            '0', '1', '2', '3', '4', '5', '6', '7', '8',
            'gene_name', 'transcript_id', 'transcript_name',
            'exon_number', 'gene_id', 'biotype']
        if filename is not None:
            if os.path.exists(filename):
                self.read_csv(filename)

    def read_csv(self, in_filename):
        self.utrs = {}
        self.exon = {}
        self.cds = {}
        self.other = {}
        lt = pandas.read_csv(in_filename, sep='\t')
        ltr = lt.to_dict('records')
        to_line = collections.defaultdict(dict)
        for gene, tx, row in zip(lt['gene_name'].tolist(),
                     lt['transcript_id'].tolist(), ltr):
            if gene not in to_line:
                to_line[gene] = collections.defaultdict(list)
            to_line[gene][tx].append(row)
        for g in to_line:
            self.utrs[g] = collections.defaultdict(list)
            self.exon[g] = collections.defaultdict(list)
            self.cds[g] = collections.defaultdict(list)
            self.other[g] = collections.defaultdict(list)
            for t in to_line[g]:
                self.exon[g][t] = [x for x in to_line[g][t] if x['2'] == 'exon']
                self.utrs[g][t] = [x for x in to_line[g][t] if re.search('UTR', x['2'])]# == 'UTR']
                self.cds[g][t] = [x for x in to_line[g][t] if x['2'] == 'CDS']
                self.other[g][t] = [x for x in to_line[g][t] if \
                              x['2'] not in ['exon', 'UTR', 'CDS']]
        print "load_gtf(): created dicts of genes size: utrs: %i exons: %i cds: %i other: %i" % (
            len(self.utrs), len(self.exon), len(self.cds), len(self.other))
        print "load_gtf(): total txpt numbers are %i %i %i %i" %(
            sum([len(self.utrs[gene_name] )for gene_name in self.utrs]),
            sum([len(self.exon[gene_name]) for gene_name in self.exon]),
            sum([len(self.cds[gene_name]) for gene_name in self.cds]),
            sum([len(self.other[gene_name]) for gene_name in self.other])
            )
        self.gene_to_tx_to_line = to_line
        return (self.utrs, self.exon, self.cds, self.other)

    def as_list(self):
        return (self.utrs, self.exon, self.cds, self.other)

    def to_df(self):
        rows = []  # Of dict objects.
        for _feat in [self.utrs, self.exon, self.cds, self.other]:
            for gene in _feat:
                for tx in _feat[gene]:
                    rows.extend(_feat[gene][tx])
        self.df = pandas.DataFrame(rows)


    def to_csv(self, fname):
        self.to_df()
        self.df.to_csv(fname, sep='\t', index=False)

    def flip_minus_strand_features(self):
        chr_len = self.chr_lens
        self.minus_strand_flipped = True
        for _feat in [self.utrs, self.exon, self.cds]:
            for gene_name in _feat:
                for txpt_id in _feat[gene_name]:
                    try:
                        if _feat[gene_name][txpt_id][0]['6'] == '+':
                            continue
                        chrm = _feat[gene_name][txpt_id][0]['0']
                    except: continue
                    for _row in _feat[gene_name][txpt_id]:
                        old_right = int(_row['4'])
                        _row['4'] = chr_len[chrm] - (int(_row['3']) - 0) + 1
                        _row['3'] = chr_len[chrm] - (old_right - 0)
    # - strand features are still listed on the - strand,
    # but their genomic coordinates are relative to the 3' end
    # of the + strand now. They can be treated as + strand,
    # except that sequences should be read from rc_sequences.

    def remove_five_prime_utrs(self, minus_strand_features_flipped=False):
        if not minus_strand_features_flipped:
            print "Error: minus strand features must be flipped before",
            print " calling gtf.remove_five_prime_utrs()."
            sys.exit()
        n_five_prime, n_three_prime, no_cds = 0, 0, 0
        for gene_name in self.utrs:
            for txpt_id in self.utrs[gene_name]:
                to_remove = []
                for i, utr_row in enumerate(self.utrs[gene_name][txpt_id]):
                    if gene_name not in self.cds:
                        no_cds += 1
                        continue
                    if txpt_id not in self.cds[gene_name]:
                        no_cds += 1
                        continue
                    far_left = min([
                        int(row['3']) for row in self.cds[gene_name][txpt_id]])
                    far_right = max([
                        int(row['4']) for row in self.cds[gene_name][txpt_id]])
                    middle_cds = float(far_left + far_right)/2.
                    middle_utr = float(int(utr_row['3']) + int(utr_row['4']))/2.
                    if middle_utr > middle_cds:
                        #to_keep.append(row)
                        n_three_prime += 1
                    else:
                        to_remove.append(i)
                        n_five_prime += 1
                to_remove.sort(reverse=True)
                if len(to_remove) > 0:
                    for index in to_remove:
                        del self.utrs[gene_name][txpt_id][index]
                #self.utrs[gene_name][txpt_id] = to_keep[:]
        print """5' UTRs: {0}. 3' UTRS: {1}. UTR with no CDS: {2}""".format(
            n_five_prime, n_three_prime, no_cds)
        #sys.exit()

    def flip_minus_strand_peak_features(self, peaks):
        chr_len = self.chr_lens
        if 'chrom' in peaks.columns: chrm_label = 'chrom'
        else: chrm_label = 'chrm'
        _ivs = zip(peaks[chrm_label].tolist(), peaks.left, peaks.right,
                   peaks.strand)
        def left(_iv):
            if _iv[3] == '+': return _iv[1]
            return chr_len[_iv[0]] - _iv[2]
        def right(_iv):
            if _iv[3] == '+': return _iv[2]
            return chr_len[_iv[0]] - _iv[1] + 1
        peaks['left'] = [left(_iv) for _iv in _ivs]
        peaks['right'] = [right(_iv) for _iv in _ivs]
#        for index, row in peaks.iterrows():
#            if row['strand'] == '+': continue
#            _new_right = chr_len[row['chrm']] - row['left'] + 1
#            peaks.loc[index, 'left'] = chr_len[row['chrm']] - peaks.loc[
#                index, 'right']
#            peaks.loc[index, 'right'] = _new_right
        return peaks

    def longest_txpt(self):
        self.longest_tx = collections.defaultdict(str)
        for gene in self.exon:
            tx_len = {}
            for tx in self.exon[gene]:
                tx_len[tx] = self.len_of_a_txpt(gene, tx)
            self.longest_tx[gene] = sorted(tx_len, key=lambda x: tx_len[x])[-1]
    
    def len_of_a_txpt(self, gene, tx):
        return sum([
                    row['4'] - row['3'] for row in self.exon[gene][tx]])
                    
    def subset_to_longest_txpt(self):
        if not hasattr(self, 'longest_tx'): self.longest_txpt()
        for gene in self.exon:
            if gene in self.longest_tx:
                for _feat in [self.utrs, self.exon, self.cds, self.other]:
                    if gene not in _feat: continue
                    to_del = [tx for tx in _feat[gene] if \
                              tx != self.longest_tx[gene]]
                    for tx in to_del:
                        del _feat[gene][tx]
            else: print "No exons? %s" % gene

    def f_dict(self, _dict):
        li = ''
        for tx, _list in _dict.items():
            li += "\n-------\nTranscript:{0}:\n".format(tx)
            if _dict is None: return li
            for row in _list:
                li += '\n ~ \n'
                for k, v in row.items():
                    li += "\t{0}\t{1}\n".format(k, v)
        return li

    def inspect(self, gene):
        print "\/T\/ \n{0}:\nCDS: {1}\n UTR: {2}\n Exon: {3}\n Other: {4}".format(
            gene,
            self.f_dict(self.cds[gene] if gene in self.cds else None),
            self.f_dict(self.utrs[gene] if gene in self.utrs else None),
            self.f_dict(self.exon[gene] if gene in self.exon else None),
            self.f_dict(self.other[gene] if gene in self.other else None))
    
    def make_genes_using_input_peaks_file(
            self, peaks, sequences, rc_sequences, gtf_fname):
        self.flocs_map = {}
        self.sbgenes = {}
        if not hasattr(self, 'longest_tx'): self.longest_txpt()
        tups = zip(peaks['gene_name'].tolist(),
                   peaks['left'].tolist(), peaks['right'].tolist(),
                   peaks['chrm'].tolist(), peaks['strand'].tolist())
        peaks_by_gene = collections.defaultdict(list)
        for tup in tups:
            peaks_by_gene[tup[0]].append(tup)
        print "Making genes..."
        for gene in self.exon:
            if (gene not in self.cds) or (len(self.cds[gene]) == 0):
                print "No CDS or txpts for CDS {0}".format(gene)
                continue
            if gene in peaks_by_gene:
                n_of_peaks_in_tx = collections.defaultdict(int)
                n_of_peaks_in_coding_tx = collections.defaultdict(int)
                for peak in peaks_by_gene[gene]:
                    txpts = self.transcripts_overlapping_a_range(
                        peak, gene)
                    for tx in txpts:
                        n_of_peaks_in_tx[tx] += 1
                        if (tx in self.cds[gene]) and (len(self.cds[gene][tx])>0):
                            n_of_peaks_in_coding_tx[tx] += 1
                #print n_of_peaks_in_tx
                if len(n_of_peaks_in_tx) == 0:
                    print "no overlaping exons for {0}: error".format(gene)
                    continue
                if len(n_of_peaks_in_coding_tx) == 0:
                    pass
 #                   print "no overlapping exons for {0} in any coding tx.".format(
 #                       gene)
                else:
                    n_of_peaks_in_tx = n_of_peaks_in_coding_tx
                vals = n_of_peaks_in_tx.values(); vals.sort()
                max_peak_n = vals[-1]
                ties = [tx for tx in n_of_peaks_in_tx \
                    if n_of_peaks_in_tx[tx] == max_peak_n]
                if len(ties) == 1:
                    tx_to_use = ties[0]
                else:
                    tx_len = {}
                    for tx in ties:
                        tx_len[tx] = self.len_of_a_txpt(gene, tx)
                    tx_to_use = sorted(tx_len.keys(), key=lambda x: tx_len[x])[-1]
#                print "making gene with txpt {0}".format(tx_to_use)
                self.make_a_gene(gene, tx_to_use, sequences, rc_sequences)
            else:
                continue
                if gene not in self.longest_tx:
                    print 'Not in longest_tx.'
                    continue
                if self.longest_tx[gene] == '':
                    print 'Is empty.'
                    continue
                self.make_a_gene(gene, self.longest_tx[gene],
                                 sequences, rc_sequences)
        x = 'Y65B4BL.1'
        print "Y92H12A.4 in self.exon: {e}. In self.cds: {cds}. in self.sbgenes: {sb}.".format(
            e=(x in self.exon), cds=(x  in self.cds),
            sb=(x in self.sbgenes))

    def transcripts_overlapping_a_range(self, peak, gene):
        overlapping_tx = set()
        for tx in self.exon[gene]:
            for exon in self.exon[gene][tx]:
                #print exon
                if exon['3'] < peak[1] < exon['4']:
                    overlapping_tx.add(tx)
                    continue
                if exon['3'] < peak[2] < exon['4']:
                    overlapping_tx.add(tx)
                    continue
                if peak[1] <= exon['3'] <= peak[2]:
                    overlapping_tx.add(tx)
                    continue
                if peak[1] <= exon['4'] <= peak[2]:
                    overlapping_tx.add(tx)
                    continue
        return overlapping_tx

    def make_genes(self, sequences, rc_sequences):
        if not hasattr(self, 'longest_tx'): self.longest_txpt()
        self.flocs_map = {}
        self.sbgenes = {}
        for gene in self.cds:
            #if gene == 'nurf-1':
            #    print "#@#@"
            #    print self.cds[gene]
            #else: continue
            #if cds[0]['gene_id'] == 'WBGene00003401':
            #    print exon_dict
            #    print "^mpk-1 exon dict"
            #self.inspect(gene)
            if gene not in self.longest_tx:
                print 'not in longest_tx'
                continue
            if self.longest_tx[gene] == '':
                print 'is emp'
                continue
            self.make_a_gene(gene, self.longest_tx[gene], sequences, rc_sequences)

    def make_a_gene(self, gene, txpt, sequences, rc_sequences):
        if (txpt not in self.cds[gene]) or (
            len(self.cds[gene][txpt])==0):
            print "txpt not in cds or cds empty {0}".format(gene)
            return
        cds = sorted(
            self.cds[gene][txpt], key=lambda x: int(x['3']))
        strand, chrom = (cds[0]['6'], cds[0]['0'])
        _off = -1 if (strand == '+') else 0
        left_cds = sorted(cds, key=lambda x: x['3'])[0]['3']
        right_cds = sorted(cds, key=lambda x: x['4'])[-1]['4']
        left_utr_span = [cds[0]['3'], cds[0]['3']]
        right_utr_span = [cds[-1]['4'], cds[-1]['4']]
        given_left_utr = False
        given_right_utr = False
        if (gene in self.utrs) and (txpt in self.utrs[gene]) and \
           (len(self.utrs[gene][txpt]) > 0):
            utrs = self.utrs[gene][txpt]
            left_utr_iv = sorted(utrs, key=lambda x: x['3'])
            right_utr_iv = sorted(utrs, key=lambda x: x['4'])
            if int(left_utr_iv[0]['3']) < left_cds:
                left_utr_span = [left_utr_iv[0]['3'], left_cds]
                given_left_utr = True
            if int(right_utr_iv[-1]['4']) > right_cds:
                right_utr_span = [right_cds + 1, right_utr_iv[-1]['4']]
                given_right_utr = True
        exon_dict, exon_borders_in_seq, exon_num = ({}, {}, 1)
        # Start with the left utr region, if it exists.
        default_left = 50
        default_right = 140
        if not given_left_utr:
            left_utr_span = [left_cds - default_left, left_cds]
            left_utr_span[0] = max([left_utr_span[0], 1])
        #if not given_right_utr:
        #    right_utr_span = [right_cds, right_cds + default_right]
        #    right_utr_span[1] = min([right_utr_span[1], chr_lens[chrom]])
        left_utr_iv = [
            chrom, left_utr_span[0]+_off, left_utr_span[1]+_off, strand]
        _seq = sequences if strand == '+' else rc_sequences
        _let = 'is pos' if strand == '+' else 'is neg'
        seq = seq_from_iv([
            chrom, left_utr_span[0]+_off, left_utr_span[1]+_off, '+'], _seq)
        exon_dict[exon_num] = left_utr_iv
        exon_borders_in_seq[exon_num] = [0, len(seq)]
        for ex in cds:
            exon_num += 1
            slice_coord = [ex['0'], ex['3'] + _off, ex['4'], strand]
            init_seq_len = len(seq)
            exon_dict[exon_num] = slice_coord
            #if exon_num == 2:
            #    print slice_coord
            #    print seq_from_iv(
            #        [ex['0'], ex['3'] + _off, ex['4'], '+'], _seq)
            seq += seq_from_iv([ex['0'], ex['3'] + _off, ex['4'], '+'], _seq)
            exon_borders_in_seq[exon_num] = [init_seq_len, len(seq)]
        exon_num += 1
        # Add right utr.
        init_seq_len = len(seq)
        right_utr_iv = [
            chrom, right_utr_span[0]-1, right_utr_span[1], strand]
        seq += seq_from_iv(right_utr_iv, _seq)
        exon_dict[exon_num] = right_utr_iv
        exon_borders_in_seq[exon_num] = [init_seq_len, len(seq)]
        txpt_span = [left_utr_span[0]-1, right_utr_span[1]]
        if given_left_utr:
            txpt_span[0] = left_utr_span[0]
        if given_right_utr:
            txpt_span[1] = right_utr_span[1]
        seq = extend_to_stop(
            [left_cds, right_cds], txpt_span, seq, cds, strand,
            sequences, rc_sequences, exon_borders_in_seq, gene)
        self.flocs_map[gene] = flocs(
            txpt, chrom, txpt_span[0], txpt_span[1],
            left_cds, right_cds, strand, seq, exon_dict, gene,
            exon_borders_in_seq, given_left_utr, given_right_utr)
        self.true_strand = strand
        self.flocs_map[gene].minus_strand_flipped = self.minus_strand_flipped
        self.flocs_map[gene].given_left_utr = given_left_utr
        self.flocs_map[gene].given_right_utr = given_right_utr
#            if cds[0]['gene_id'] == 'WBGene00003401':
#                print self.flocs_map[gene].exons_dict
#                print "^mpk-1 exon dict"
        #print self.flocs_map[gene]
        self.sbgenes[gene] = sbGene(HTSeq.GenomicInterval(
            chrom, txpt_span[0], txpt_span[1], strand), name=gene)
        self.sbgenes[gene].minus_strand_flipped = self.minus_strand_flipped
        self.sbgenes[gene].true_strand = strand
        self.sbgenes[gene].gene_name = gene
        self.sbgenes[gene].chr_lens = self.chr_lens
        self.sbgenes[gene].exons_dict = exon_dict
        self.sbgenes[gene].seq = seq
    # Get start, end, left_txpt_start and right_txpt_end.
    # Also load some sequences.
    # .gtf files are 1-based. Sequence slices are 0-based.
    # The left border of slices are (of course) included.
    # The right border of some exons x, say the third A in TAA, should be found at sequences[x-1:x].
    # To say the same thing, a stop codon in the .gtf at 647-649 will be found at sequences[646:649].
    # To find a start codon on the - strand located at 2167-2169 in a .gtf, the slice is sequences[2166:2169],
    # or, in other words, nt x:x+n (not including the final nt) on either strand is [x-1:x-1+n].
    # The downstream iv should start with the nt just after the stop, and be (say) 200 nt long.
    # In the .gtf, this would be stop_right_border+1:stop_right_border+1+200.
    # As a slice, this would therefore be stop_right_border:stop_right_border+200.
    # Upstream is start_left_border-200:start_left_border in the .gtf,
    # so start_left_border-200-1:start_left_border-1 as a slice.

    def make_ncrna_sbgenes(self, sequences, rc_sequences, _pad=200):
#        print self.gene_to_tx_to_line.keys()
        no_cds_annotated = set()
        for gene in self.exon:
            if gene not in self.cds:
                no_cds_annotated.add(gene)
                continue
            has_row = False
            for tx in self.cds[gene]:
                if len(self.cds[gene][tx]) > 0:
                    has_row = True
            if not has_row:
                no_cds_annotated.add(gene)
#        no_cds_annotated = set(self.exon.keys()) - set(self.cds.keys())
        if 'linc-7' in self.exon.keys():
            print "Found linc-7"
            print self.cds['linc-7']
        print "---"
        print "r23rua;slekjf"
        nm = '''
        for gene in no_exons_cds:
            tx_len = {}
            for tx in self.gene_to_tx_to_line[gene]:
                transcript_rows = [x for x in \
                    self.gene_to_tx_to_line[gene][tx] if x['2']=='transcript']
                if len(transcript_rows) == 1: row = transcript_rows[0]
                tx_len[tx] = row['4'] - row['3']
            self.longest_tx[gene] = sorted(tx_len, key=lambda x: tx_len[x])[-1]
        '''
        for gene in no_cds_annotated:
            if gene not in self.longest_tx: continue
            if self.longest_tx[gene] == '': continue
            txpt = self.longest_tx[gene]
            rows = self.gene_to_tx_to_line[gene][txpt]
            rows = [x for x in rows if x['2']=='transcript']
            if len(rows) == 1: row = rows[0]
            else:
                print "???"
                sys.exit()
            left = max([1, int(row['3'])-_pad])
            self.sbgenes[gene] = sbGene(HTSeq.GenomicInterval(
                row['0'], left, int(row['4'])+_pad, row['6']))
            self.sbgenes[gene].minus_strand_flipped = False#self.minus_strand_flipped
            self.sbgenes[gene].true_strand = row['6']
            self.sbgenes[gene].gene_name = gene
            self.sbgenes[gene].chr_lens = self.chr_lens
            self.sbgenes[gene].exons_dict = {
                1: [row['0'], left, int(row['4'])+_pad, row['6']]}
            if gene=='linc-7': print "Added linc-7 to gtf.sbgenes"

    def add_peak_locations_to_transcripts(self, peaks):
        assert(type(peaks) == type(pandas.DataFrame()))
        peaks_d = to_dict(peaks)
        for gene_name in peaks_d:
            if gene_name != 'nurf-1': continue
            if gene_name in self.flocs_map:
                self.flocs_map[gene_name].find_peaks(peaks_d, self.chr_lens)
        return self.flocs_map


def to_dict(df):
    _li = {}
    for index, row in df.iterrows():
        _li.setdefault(row['gene_name'], [])
        _li[row['gene_name']].append(row.to_dict())
    return _li


def extend_to_stop(cds_span, txpt_span, seq, sorted_cds, strand, sequences,
                   rc_sequences,
                   exon_borders_in_seq, gene_name):
    verbose=True
    if (cds_span[1] == txpt_span[1]) and (seq[-3:] not in ['TAA', 'TAG', 'TGA']):
        #if verbose: print "Extending %s to get stop..." % gene_name
        _add = seq_from_iv(
            (sorted_cds[0]['0'], cds_span[1], cds_span[1] + 3, '+'),
            sequences if strand == '+' else rc_sequences)
        if _add in ['TAA', 'TAG', 'TGA']:
            txpt_span[1] += 3
            seq = seq + _add
            #if verbose: print "Did it work? Expected stop: %s" % seq[-3:]
            try:
                exon_borders_in_seq[len(exon_borders_in_seq)][1] = len(seq)
            except:
                loggern.warn(
                    "%s: Missing an exon? %s" % (gene_name, str(exon_borders_in_seq)))
        else:
            if verbose: print "%s: Did it work? Expected stop: %s" % (
                gene_name, _add)
    return seq


def seq_from_iv(iv, sequences):
    """Returns a slice from 0-based, slice coordinates.
    """
    if iv[2] > len(sequences[iv[0]]):
        end = len(sequences[iv[0]])
    else:
        end = iv[2]
    a_seq = sequences[iv[0]][iv[1]:end]
    if iv[3] == "-":
        a_seq = rc(a_seq)
    return a_seq


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([basecomplement[base] for base in list(s)])


def rc(s): return complement(s[::-1])

def get_sequences():
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    rc_sequences = dict(
        (p.name.split(' ')[0], rc(p.seq)) for p in HTSeq.FastaReader(fasta_filename))
    chr_lens = dict(
        [(name, len(sequences[name])) for name in sequences])
    return (sequences, rc_sequences, chr_lens)

def _time(start): return time.time() - start

if __name__ == '__main__':
    print "Debug mode for gtf.py."
    print sys.version
    start = time.time()
    seqs, rc_seqs, chr_lens = get_sequences()
    print "get_sequneces(): {0} sec".format(_time(start))
    start = time.time()
    g = gtf('/opt/lib/gtf_with_names_column.txt', debug=True)
    print "gtf init(): {0} sec".format(_time(start))
    start = time.time()
    g.flip_minus_strand_features(chr_lens)
    g.longest_txpt()
    print "gtf flip strand and longest txpt(): {0} sec".format(_time(start))
    g.make_genes(seqs, rc_seqs)
