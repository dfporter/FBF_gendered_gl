import pandas, functools, re, os, sys, glob, collections
import numpy as np

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
#prefix + 'old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
}
label_to_fname = dict(zip(fname_to_label.values(), fname_to_label.keys()))


class figureMaker():

    def __init__(self, output_dir='figs/'):
        
        self.fname_to_label = fname_to_label
        self.label_to_file = label_to_fname
        
        if not os.path.exists(output_dir):
            os.system('mkdir ' + output_dir)
            
        self.df = {}
    
    def load_peaks_from_excel_file(self, excel_fname='../tables/Table S1 Peaks.xls'):
        
        xl = pandas.ExcelFile(excel_fname)
        
        self.df, self.targs = ({}, {})
        
        for sheet_name in xl.sheet_names:
            print(sheet_name)
            df = pandas.read_excel(xl, sheetname=sheet_name)
            self.df[sheet_name] = df.copy()
            self.targs[sheet_name] = set(df['Gene name'].tolist())
            
        self.df['OO FBF (20°C, either FBF)'] = pandas.concat([
            self.df['OO FBF-1 (20°C)'], self.df['OO FBF-2 (20°C)']])
        self.targs['OO FBF (20°C, either FBF)'] = set(self.df['OO FBF (20°C, either FBF)']['Gene name'].tolist())
        
        both = 'OO FBF (20°C, both FBFs)'
        self.df[both] = self.df['OO FBF (20°C, either FBF)'].copy()
        self.df[both] = self.df[both][[
    ((name in self.targs['OO FBF-1 (20°C)']) and (name in self.targs['OO FBF-2 (20°C)'])) for name in self.df[both]['Gene name']
        ]]
        self.targs[both] = set(self.df[both]['Gene name'].tolist())
        
        self.all_targets_in_excel = functools.reduce(lambda x,y: x | y, self.targs.values())
        
    def load_peaks_csv_files(self, label_to_file=None):
        """Goes through combined_filtered/ peaks files and loads into self.df as dict.
        Sets self.df (dict of dfs) and self.targs (dict of gene name sets).
        """
        
        if label_to_file is not None:
            self.label_to_file = label_to_file
        elif hasattr(self, 'label_to_fname'):
            self.label_to_file = self.label_to_fname
        else:
            self.label_to_file = label_to_fname
        
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

