
import HTSeq
import re
import sys
import os
import pandas
import numpy as np

class peaks(object):

    def __init__(
            self, data=None, name='unnamed peaks', file=None):
        
        self.name = name
        
        if data is None:
            if (file is not None):
                self.data = pandas.read_csv(file, sep='\t', index_col=False)
            else: self.data = pandas.DataFrame()
        else:
            self.data = data  # The peaks dataframe.
        
        self.file = file  # Path to peaks.txt file.
        self.bedgraph_filenames = {}

    def __str__(self):
        return "{n}:\n{f}\nRows:{d}\n".format(
            n=self.name, f=self.file,
            d=self.data_stat())

    def data_stat(self):
        if self.data is None:
            return 'None'
        return len(self.data.index)

    def absolute_read_number(self, sizes):
        for col in sizes:
            self.data['unnorm_reads_' + col] = [
                int(sizes[col] * float(x)/float(1e6)) for x in \
                self.data[col].tolist()]
        
    def add_reads(
        self, ga=HTSeq.GenomicArray('auto'), name='unnamed array'):
        if self.data is None:
            print("Asked to set reads on an empty peaks object!")
            return False
        peaks = self.data  # Shorter.
        try:
            ivs = list(zip(peaks['chrm'].tolist(), peaks['left'].tolist(),
                peaks['right'].tolist(), peaks['strand'].tolist()))
        except:
            print("Could not create ivs from %s." % self.name)
            return False
        self.data[name] = [
            get_val(ga, HTSeq.GenomicInterval(*iv)) for iv in ivs]

    def set_ratio(self, col1='fbf1_reads', col2='fbf1_n2_reads',
                  ratio_col='ratio'):
        positive_denom = [x for x in self.data[col2].tolist() if x >0]
        if len(positive_denom) == 0:
            print("All zero values for %s..." % col2)
            print("Defaulting to 1.")
            positive_denom = [1.]
        minimum_denom = min(positive_denom)
        tups = list(zip(self.data[col1].tolist(), self.data[col2].tolist()))
        self.data[ratio_col] = [
            float(x[0])/max([x[1], minimum_denom]) for x in tups]

    def get_filtered_obj(self, col='ratio', cutoff=5):
        filt_data = self.data[self.data[col]>=cutoff]
        return peaks(data=filt_data, name='filtered_' + self.name)

    def write_table(self, filename):

        if 'exp_reads' in self.data.columns:
            self.data.sort_values(by=['exp_reads'], ascending=False, inplace=True)
        
        self.data.to_csv(filename, sep='\t', index=False)

    def set_sum(self, to_sum=[], summed_col='average'):
        vals = [self.data[x].tolist() for x in to_sum]
        print(to_sum)
        print("^ to sum")
        denom = max([1, len(to_sum)])
        sums = []
        for i in range(len(self.data.index)):
            sums.append(
                float(sum([vals[x][i] for x in range(len(vals))]))/denom)
        self.data[summed_col] = sums

    def has_columns(self, cols):
        if len(set(cols) - set(self.data.columns)) > 0:
            return False
        else: return True

    def cut_columns(self, remove_those_with_patterns=[]):
        to_remove = []
        for pat in remove_those_with_patterns:
            to_remove.extend([x for x in self.data.columns if re.search(pat, x)])
        for col in set(to_remove):
            del self.data[col]


def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))
