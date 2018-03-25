from __future__ import division
import pandas
import re
import os
import sys
import glob
import numpy as np
import collections

prefix = '/groups/Kimble/Common/fbf_celltype/combined_filtered/'
prefix = './combined_filtered/'


fname_to_label = {
prefix + 'old_fbf2.txt': 'old_fbf2',
prefix + 'old_fbf1.txt': 'old_fbf1',
prefix + 'oo_both.txt': 'oo_both',
prefix + 'oo_fbf1.txt': 'oo_fbf1',
prefix + 'oo_fbf2.txt': 'oo_fbf2',
prefix + 'sp_both.txt': 'sp_both',
prefix + 'sp_fbf1.txt': 'sp_fbf1',
prefix + 'sp_fbf2.txt': 'sp_fbf2',
prefix + 'old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
}
label_to_fname = dict(zip(fname_to_label.values(), fname_to_label.keys()))


class figureMaker():

    def __init__(self, output_dir='figs/'):
        
        self.fname_to_label = fname_to_label
        self.label_to_file = label_to_fname
        
        if not os.path.exists(output_dir):
            os.system('mkdir ' + output_dir)
            
        self.df = {}
    
    def load_peaks_csv_files(self, label_to_file=None):
        """Goes through combined_filtered/ peaks files and loads into self.df as dict.
        Sets self.df (dict of dfs) and self.targs (dict of gene name sets).
        """
        
        if label_to_file is None:
            pass#self.label_to_file = self.label_to_fname
        else:
            self.label_to_file = label_to_file
        
        self.df = dict([
            (k, self.read_csv(v)) for k, v in self.label_to_file.items()])
        
        self.targs = dict([
            (k, set(self.df[k]['gene_name'].tolist())) for k in self.df])
        
        self.biotypes = {}
        for df in self.df.values():
            if 'biotype' in df.columns:
                col = 'biotype'
            elif 'Biotype' in df.columns:
                col = 'Biotype'
            else:
                continue
            
            self.biotypes.update(dict(zip(df['gene_name'], df[col])))
        
        if len(self.biotypes) > 0:
            print("Was able to read biotypes from peaks files.")
            
    def read_csv(self, fname):
        return pandas.read_csv(fname, sep='\t', index_col=False)

#    def n_gene_names(self, df):