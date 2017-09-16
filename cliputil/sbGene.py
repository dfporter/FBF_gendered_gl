import collections
import numpy as np
import HTSeq
import os
import sys
import pandas
import matplotlib.pyplot as plt
import scipy
from scipy import signal


class sbGene:

    def __init__(self, iv, name=''):
        self.iv = iv
        self.name = name
        self.n_maxima = 0
        
    def __str__(self):
        li = "sbGene instance. gene {0} ".format(self.name)
        if hasattr(self, 'exons_dict'):
            li += str(self.exons_dict)
        li += '\n'
        return li
        
    def find_subpeaks(
        self, coverage, exons_dict, rc_exons_dict_if_minus=True, debug=False):
        self.determine_coverage_accross_exons(
            coverage, exons_dict,
            rc_exons_dict_if_minus=rc_exons_dict_if_minus)
        self.n_clusters = self.find_clusters()
        #self.find_peaks_cwt()
        self.find_peaks_maxima(debug=debug)
        self.footprint = len(np.where(self.exon_coverage>1)[0])

        n = """
        self.peak_ivs_genomic_coord = [HTSeq.GenomicPosition(
            self.iv.chrom,
            self.pos_in_exon_coverage_mapped_to_genomic_pos[pos],
            self.iv.strand,
            ) for pos in self.peaks]
        self.peak_ivs_relative = [HTSeq.GenomicPosition(
            self.iv.chrom, pos, self.iv.strand,
            ) for pos in self.peaks]"""

    def find_peaks_maxima(
            self, min_coverage=2, min_frac_of_highest_peak=0.05,
            min_watershed=0.5, debug=False):
        self.n_maxima = 0
        self.peak_maxima_ga = HTSeq.GenomicArray('auto', stranded=True)
        # By whole gene.
        #np.split(self.exon_coverage,
        #         np.where(np.diff(self.exon_coverage) != 1)[0]+1)
        #hand = []
        def compress(a):
            # Allow peaks to be called at the edges by adding a zero..
            arr_as_list = list(a[:])
            positions = list(np.where(np.diff(arr_as_list) == 0)[0])
            pos_set = set(positions)
            # This set of positions lacks the index of any pos where a changes.
            cumulative_del = 0
            _to_cumulative = {}
            index_in_compressed = 0
            uncompress_index = {}
            compressed_a = []
            for i, val in enumerate(a):
                if i in pos_set:
                    pass
                else:  # There is a change in 'a' at this position.
                    uncompress_index[index_in_compressed] = i
                    index_in_compressed += 1
                    compressed_a.append(val)
            return compressed_a, uncompress_index
        compr, uncompress_index = compress(self.exon_coverage)

        indexes_of_maxima = signal.argrelextrema(np.array([0] + compr), np.greater)[0]
        indexes_of_maxima = [x-1 for x in indexes_of_maxima]
        m_indexes = list(indexes_of_maxima)
        if len(m_indexes) == 0:
            if debug:
                print "No maxima in {0}".format(self.gene_name)
            return
        # Watersheding.
        too_close = []
        hand = set()
        for i, val in enumerate(m_indexes):
            if i + 1 >= len(m_indexes):
                too_close.append(hand)
                break
            next_val = m_indexes[i+1]
            min_between = np.min([
                    compr[x] for x in range(val, next_val)])
            min_to_separate = min_watershed * np.min([
                compr[val], compr[next_val]])
            if min_between <= min_to_separate:
                too_close.append(hand)
                hand = set()
            else:
                hand.add(val)
                hand.add(next_val)
        to_del = set()
        for peak_set in too_close:
            if len(peak_set) == 0: continue
            for i in sorted(list(peak_set), key=lambda x: compr[x])[:-1]:
                to_del.add(i)
        m_indexes = [x for x in m_indexes if x not in to_del]
        # Uncompress.
        m_indexes = [uncompress_index[x] for x in m_indexes]
        # Apply height cutoffs.
        index_of_primary_peak = sorted(
                m_indexes,
                key=lambda x: self.exon_coverage[x])[-1]
        min_height_subpeak = min_frac_of_highest_peak * self.exon_coverage[
            index_of_primary_peak]
        min_cutoff = np.max([min_coverage, min_height_subpeak])
        m_indexes = [
            x for x in m_indexes if self.exon_coverage[x] >= min_cutoff]
        # Set info arrays/objects with result.
        self.indexes_in_exon_array_of_local_maxima = m_indexes
        if len(m_indexes) == 0:
            self.peak_maxima_iv = HTSeq.GenomicPosition('I', 1, '+')
            return
        self.peak_maxima_by_index_in_exon_coverage = m_indexes
        self.peak_maxima_iv = [self.to_iv(x) for x in m_indexes]
        [self.to_ga(self.peak_maxima_ga, x) for x in m_indexes]
        self.n_maxima = len(m_indexes)

    def find_peaks_cwt(self, min_coverage=2):
        self.peak_pos = [x for x in signal.find_peaks_cwt(
            self.exon_coverage, np.array([5, 250, 5]),
            noise_perc=0.05) if self.exon_coverage[x] >= min_coverage]
        self.peak_pos_iv = [self.to_iv(x) for x in self.peak_pos]
        return self.peak_pos

    def filter_by_intersection(self, list_of_ivs):
        self.filtered_peak_ivs_genomic_coord = []
        for iv in self.peak_ivs_genomic_coord:
            for alt_iv in list_of_ivs:
                if alt_iv[3] != iv.strand: continue
                if alt_iv[0] != iv.chrom: continue
                if alt_iv[1] <= iv.pos <= alt_iv[2]:
                    self.filtered_peak_ivs_genomic_coord.append(iv)
        return self.filtered_peak_ivs_genomic_coord

    def add_to_wig(self, ga):
        for iv in self.peak_ivs_genomic_coord:
            ga[iv] += 1

    def to_ga(self, ga, start, end=None):
        to_genome = self.pos_in_exon_coverage_mapped_to_genomic_pos
        if end is None:
            if start not in to_genome:
                print "{0} pos not in exon ({1} -> {2}, {3} keys())".format(
                    start, min(to_genome.keys()), max(to_genome.keys()),
                    len(to_genome))
                ga[HTSeq.GenomicPosition(self.iv.chrom, 1, self.iv.strand)] = 1
                return
            ga[HTSeq.GenomicPosition(
                self.iv.chrom,
                self.pos_in_exon_coverage_mapped_to_genomic_pos[start],
                self.iv.strand)] = self.exon_coverage[start]
            return
#        _iv = HTSeq.GenomicInterval(
#            self.iv.chrom,
#            self.pos_in_exon_coverage_mapped_to_genomic_pos[start],
#            self.pos_in_exon_coverage_mapped_to_genomic_pos[end],
#            self.iv.strand)
        for pos in range(start, end + 1):
            ga[HTSeq.GenomicPosition(
                self.iv.chrom,
                self.pos_in_exon_coverage_mapped_to_genomic_pos[start],
                self.iv.strand)] = self.exon_coverage[pos]

    def to_iv(self, start, end=None):
        to_genome = self.pos_in_exon_coverage_mapped_to_genomic_pos
        if end is None:
            if start not in to_genome:
                print "{0} pos not in exon ({1} -> {2}, {3} keys())".format(
                    start, min(to_genome.keys()), max(to_genome.keys()),
                    len(to_genome))
                return HTSeq.GenomicPosition(self.iv.chrom, 1, self.iv.strand)
            return HTSeq.GenomicPosition(
                self.iv.chrom,
                self.pos_in_exon_coverage_mapped_to_genomic_pos[start],
                self.iv.strand)
        return HTSeq.GenomicInterval(
            self.iv.chrom,
            self.pos_in_exon_coverage_mapped_to_genomic_pos[start],
            self.pos_in_exon_coverage_mapped_to_genomic_pos[end],
            self.iv.strand)

    def find_clusters(self):
        # TO DO: A significant number of objects in arr are empty ([]).
        # Not sure why. Was this figured out?
        self.clusters = []
        self.clusters_as_ga = HTSeq.GenomicArray('auto', stranded=True)
        n_clusters = 0
        nonzeroes = np.nonzero(self.exon_coverage)[0]
        clusters = self.consecutive(np.nonzero(self.exon_coverage)[0])
        self.clusters_as_indices_in_exon_coverage_array = clusters
        #print "exon_coverage= {0}\n\nnonzeroes= {1}\n\nconsecutive(nonzeroes)=\n\n{2}\
#".format(self.exon_coverage, nonzeroes, clusters)
        if (len(self.exon_coverage)==0) or (len(clusters)==0):
            n_clusters = 0
            return n_clusters
        arr = np.array([self.exon_coverage[i] for i in clusters])
        if len(arr) > 0:
            for _index_array in clusters:
                for pos in _index_array:
                    self.clusters_as_ga[self.to_iv(pos)] += 1
            nonzeros = [x for x in arr if len(x) > 0]
            #maxes = [np.max(x) for x in arr if len(x) > 0]
            self.clusters = [x for x in nonzeros if np.max(x)>1]
            n_clusters = len(self.clusters) #len([x for x in maxes if x>1])
            #n_clusters = len(filter(lambda x: np.nanmax(x)>1, arr))
        if n_clusters == 0:
            print "...Failure"     
        return n_clusters

    def consecutive(self, data):
        # type data: np.array, 1d
        # Returns np.array, 2d
        return np.split(data, np.where(np.diff(data) != 1)[0]+1)

    def rc_exons_dict_if_minus(self, exons_dict):
        if self.iv.strand == '-':
            rc_corrected_exons_dict = {}
            num_exons = len(exons_dict)
            #print 'rc -> input exons_dict len {0}'.format(num_exons)
            for exon_num in exons_dict:
                (chrm, left, right, strand) = exons_dict[exon_num]
                _left = self.chr_lens[chrm] - right
                _right = self.chr_lens[chrm] - left
                rc_corrected_exons_dict[num_exons - exon_num + 1] = (
                    chrm, _left, _right, '-')
                #print 'rc -> set {0} exon num to {1}'.format(
                #    num_exons - exon_num + 1,
                #    rc_corrected_exons_dict[num_exons - exon_num + 1])
            return rc_corrected_exons_dict
        else: return exons_dict

    def determine_coverage_accross_exons(
            self, coverage, exons_dict, rc_exons_dict_if_minus=True):
        '''
        exon    0      1       2
        genomic a      b       c
        rela    012345601234567012345
        '''
        if len(exons_dict) == 0:
            self.exon_coverage = []
            return
        if rc_exons_dict_if_minus:
            exons_dict = self.rc_exons_dict_if_minus(exons_dict)
        self.exon_coverage = np.array([])
        self.pos_in_exon_coverage_mapped_to_genomic_pos = {}
        first_exon_left = exons_dict[
            sorted(exons_dict.keys(), key=lambda x: int(x))[0]][1]
        reverse = True if (self.iv.strand == '-') else False
#        if strand == '-': reverse = True
        self.exons_dict = exons_dict
        self.set_exon_ivs(exons_dict)
        for iv in self.exon_ivs:
            #(chrm, left, right, strand) = exons_dict[exon_num]
            #iv = HTSeq.GenomicInterval(chrm, left, right, strand)
            iter_test = coverage[iv]
            this_exon_coverage = np.fromiter(coverage[iv], dtype='f')
            #print "this_exon_coverage = np.fromiter(coverage[this_exon]) = {0}, \
#len {1}".format(this_exon_coverage, len(this_exon_coverage))
            #print "{0} -> {1} iv (len = {2}). np.fromiter = {3}".format(
            #    iv.start, iv.end, iv.end - iv.start, len(this_exon_coverage))

            for pos in range(0, len(this_exon_coverage)):
                rel_pos = pos + len(self.exon_coverage)
                self.pos_in_exon_coverage_mapped_to_genomic_pos[rel_pos] \
                        = iv.start + pos
            self.exon_coverage = np.concatenate([
                self.exon_coverage, this_exon_coverage])

    def set_exon_ivs(self, exons_dict):
        self.exon_ivs = []
        reverse = True if (self.iv.strand == '-') else False
        exon_order = sorted(
            exons_dict.keys(), key=lambda x: int(x), reverse=reverse)
        for i in exon_order:
            self.exon_ivs.append(HTSeq.GenomicInterval(*exons_dict[i]))

