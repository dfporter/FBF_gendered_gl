"""
Exons: [ ] --- [ ] ---- [  ]
UTRs: [  ]
CDS: [   ]
x = coverage across mature mRNA (exons)
y_five = to_five(x) -> subsets to before start, and normalizes
"""
import sbGene

import collections
import numpy as np
import HTSeq
import os
import sys
import pandas
from __future__ import division

class normGene(sbGene):

    def __init__(self, iv):
        self.iv = iv
        self.five_len = 50
        self.cds_len = 1e3
        self.three_len = 200

    def normalize_coverage(self):
        self.norm_coverage = [
            self.pos_in_seq_to_norm(x) for x in self.exon_coverage]

    def create_norm_funcs(
        self, start_pos_in_seq, stop_pos_in_seq, txpt_end_in_seq):
        self.to_five = self.from_bounds(
            0, start_pos_in_seq, 0, 50)
        self.to_cds = self.from_bounds(
            start_pos_in_seq, stop_pos_in_seq, 50, 1050)
        self.to_three = self.from_bounds(
            stop_pos_in_seq, txpt_end_in_seq, 1050, 1250)
        def to_norm(pos):
            if pos < start_pos_in_seq: return self.to_five(pos)
            elif pos < stop_pos_in_seq: return self.to_cds(pos)
            else: return self.to_three(pos)
        self.pos_in_seq_to_norm = to_norm

    def from_bounds(self, unnorm_start, unnorm_end,
                    norm_start, norm_end):
        unnorm_len = unnorm_end - unnorm_start
        norm_len = norm_end - norm_start
        def tmp(pos):
            if (pos < unnorm_start):
                print "Out of bounds: {0} < {0}".format(pos, unnorm_start)
            if (pos < unnorm_end):
                print "Out of bounds: {0} < {0}".format(pos, unnorm_end)
            rel_pos = pos - unnorm_start
            return int(
                norm_start + int(rel_pos * norm_len/unnorm_len))
        return tmp


