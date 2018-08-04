from __future__ import division
import pandas as pandas
import matplotlib.pyplot as plt
import seaborn as sns
import collections
import numpy as np
import glob
import pickle
import seaborn
import scipy
import scipy.stats as scs

import importlib

rl = importlib.reload

import re
from collections import defaultdict
#from scipy.spatial import distance
#from scipy.cluster import hierarchy

import figureMaker
import utils
from peaksList import peaksList
import vennMaker
from blocks import blocki, blockii, blockiii
import countsColumnsNaming
import countsFileIO

rl(countsFileIO)
rl(countsColumnsNaming)
        
class heatmapMaker(figureMaker.figureMaker,
                   utils.translator,
                   #countsColumnsNaming.countsColumnsNaming,
                   countsFileIO.countsFileIO,
                  ):

    def load_counts(self, fname='combined_counts.txt', log_scale=True,
                        only_combined_datasets=True):
        self.load_counts_file(fname=fname, log_scale=True,
                        only_combined_datasets=only_combined_datasets, style='heatmap')
        self.rm_non_target_rnas(only_combined_datasets=only_combined_datasets)
        print("""heatmapMaker.load_counts(): used {0} targets datasets to narrow  self.counts_df
        into {1} rows.""".format(len(self.targs), len(self.counts_df.index)))
        
    def rm_non_target_rnas(self, cutoff=1, only_combined_datasets=False):
        all_targets = collections.defaultdict(int)
        
        for k in self.targs:
            if only_combined_datasets:
                alt = """
                if ((not re.search('both', k) and (not re.search('old', k)))
                ):
                    continue  # Not a combined dataset
                if k[:3] == 'old':#_fbf1_to_fbf2_n2':
                    continue  # Only using old_fbf1_to_fbf2_n2
                    """
                if not (k == 'oo_both' or k == 'sp_both'):
                    continue
            if (re.search('control', k)):
                continue
            for name in self.targs[k]:
                all_targets[name] += 1
                
        print("heatmap.rm_non_target_rnas(): self.targs.keys(): ", self.targs.keys())
        print("Targets in peaks files: ", len(all_targets))

        self.counts_df = self.counts_df.loc[[
            #((x in all_targets) and (all_targets[x] >= cutoff)) for x in self.counts_df.index]]
            (x in all_targets)  for x in self.counts_df.index]]
        
        print("Kept {0} targets".format(len(self.counts_df.index)))
        here = set(self.counts_df.index)
        all_targ = self.targs['sp_both'] | self.targs['oo_both']
        
        print("all_targ ", len(all_targ))
        print("Targets here but not in self.targs: ", here - all_targ)
        print("Targets in self.targs but not here: ", all_targ - here)
        
        self.transl()  # utils.translator() loads gtf mappings.
        
        print("biotypes peaks in SP/OO (from gtf, not peaks file.):")
        self.count_biotypes(all_targ)
        print("biotypes from counts file (from gtf):")
        self.count_biotypes(here)
        
        self.cts = self.counts_df
        return self.cts
    
    def count_biotypes(self, names_list):
        biotypes = collections.defaultdict(int)
        for name in list(names_list):
            biotypes[self.name_to_biotype.get(name, '')] += 1
        print('---\nSum ', len(names_list))
        print(biotypes)
        return biotypes
    
    def add_sp_oo(self):
        from peaksList import peaksList
        self.gl_pkl = peaksList()
        self.gl_pkl.read_sp_vs_oo()
        self.counts_df['SP RPKM'] = [self.gl_pkl.name_to_sp_rpkm(x) for \
                                     x in self.counts_df.index]
        self.counts_df['OO RPKM'] = [self.gl_pkl.name_to_oo_rpkm(x) for \
                                     x in self.counts_df.index]
        #self.counts_df['SP/OO FC'] = [self.gl_pkl.name_to_fc(x) for \
        #                             x in self.counts_df.index]

    def rm_controls(self):
        for k in ['c_sp_1', 'c_n2_1', 'c_oo_1', 'ave_neg']:
            if k in self.counts_df.columns:
                del self.counts_df[k]
        combined_controls = set(['control_sp_counts.txt', 'control_oo_counts.txt',
                             'control_n2_counts.txt'])
        to_del = set([x for x in self.counts_df.columns if (\
                    re.search('control.*matching.*', x) \
                   or re.search('n2_oo_lane', x) \
                    or re.search('\Ac_', x) \
                   or re.search('control', x) # All controls removed. 
                   )])
        to_del = to_del - combined_controls
        
        for x in set(to_del):
            del self.counts_df[x]

    def heatmap_counts_file(self, fname=None):
        
        if fname is not None:
            self.load_counts_file(fname)
            
        _dict = defaultdict(dict)
        
        def _add(_d, n):
            _d[n] = 1
            
        for k in self.targs:
            for name in self.targs[k]: _dict[name][k] = 1
                
        for name in _dict:
            for k in [_ for _ in self.targs if (_ not in _dict[name])]:
                _dict[name][k] = 0

        return pandas.DataFrame(_dict.values(), index=_dict.keys())

    def pairwise(self, df=None):
        if df is None:
            df = self.counts_df
        return df.corr()


def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))
    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l
    return cluster_classes


def col_relabel():
    def math(x): return "$\mathrm{" + x + "}$"
    col_relabel_d = dict(zip(
        ['SP FBF_2', 'SP FBF_1', 'SP FBF_3', 'LT FBF1_3', 'LT FBF1_1',
         'LT FBF1_2', 'LT FBF2_3', 'LT FBF2_1', 'LT FBF2_2',
         'OO FBF_1', 'OO FBF_2', 'OO FBF_3'],
            [math(x) for x in [
    '^{25^\circ C}Sp. FBF (2)',
    '^{25^\circ C}Sp. FBF (1)',
    '^{25^\circ C}Sp. FBF (3)',
    '^{20^\circ C}Oo. FBF}$-$\mathrm{1 (3)',
    '^{20^\circ C}Oo. FBF}$-$\mathrm{1 (1)',
    '^{20^\circ C}Oo. FBF}$-$\mathrm{1 (2)',
    '^{20^\circ C}Oo. FBF}$-$\mathrm{2 (3)',
    '^{20^\circ C}Oo. FBF}$-$\mathrm{2 (1)',
    '^{20^\circ C}Oo. FBF}$-$\mathrm{2 (2)',
    '^{25^\circ C}Oo. FBF (1)',
    '^{25^\circ C}Oo. FBF (2)',
    '^{25^\circ C}Oo. FBF (3)']]))
    return col_relabel_d


def col_relabel_unicode():
    col_relabel_d = dict(zip(
        ['SP FBF_2', 'SP FBF_1', 'SP FBF_3', 'LT FBF1_3', 'LT FBF1_1',
         'LT FBF1_2', 'LT FBF2_3', 'LT FBF2_1', 'LT FBF2_2',
         'OO FBF_1', 'OO FBF_2', 'OO FBF_3'],
        [
    '25' + u'\xb0' + 'C Sp. FBF (2)',
    '25' + u'\xb0' + 'C Sp. FBF (1)',
    '25' + u'\xb0' + 'C Sp. FBF (3)',
    '20' + u'\xb0' + 'C Oo. FBF-1 (3)',
    '20' + u'\xb0' + 'C Oo. FBF-1 (1)',
    '20' + u'\xb0' + 'C Oo. FBF-1 (2)',
    '20' + u'\xb0' + 'C Oo. FBF-2 (3)',
    '20' + u'\xb0' + 'C Oo. FBF-2 (1)',
    '20' + u'\xb0' + 'C Oo. FBF-2 (2)',
    '25' + u'\xb0' + 'C Oo. FBF (1)',
    '25' + u'\xb0' + 'C Oo. FBF (2)',
    '25' + u'\xb0' + 'C Oo. FBF (3)',
    ]))
    return col_relabel_d


def recover_label():
    _dict = col_relabel()
    relab = dict(zip(_dict.values(), _dict.keys()))
    for k in relab:
        if re.search('SP FBF_', relab[k]):
            relab[k] = re.sub('SP FBF_', '25' + u'\xb0' + 'C SP FBF rep ', relab[k])
        if re.search('OO FBF_', relab[k]):
            relab[k] = re.sub('OO FBF_', '25' + u'\xb0' + 'C OO FBF rep ', relab[k])
        if re.search('LT FBF1_', relab[k]):
            relab[k] = re.sub('LT FBF1_', '20' + u'\xb0' + 'C OO FBF-1 rep ', relab[k])
        if re.search('LT FBF2_', relab[k]):
            relab[k] = re.sub('LT FBF2_', '20' + u'\xb0' + 'C OO FBF-2 rep ', relab[k])
    return relab


def make_fig(df):
    
    cmap = sns.cubehelix_palette(
        len(df.index), light=.7, dark=.2, reverse=False,
        start=1, rot=-2, as_cmap=True)
    cmap = plt.get_cmap('ocean_r')#('gray_r')#'ocean_r')
    
    plt.clf(); plt.close()
    df = df[['SP FBF_1', 'SP FBF_2', 'SP FBF_3',
             'OO FBF_1', 'OO FBF_2', 'OO FBF_3',
             'LT FBF1_1', 'LT FBF1_2', 'LT FBF1_3',
             'LT FBF2_1', 'LT FBF2_2', 'LT FBF2_3']]
    col_relabel_d = col_relabel()
    df.columns = [col_relabel_d[x] for x in df.columns]
    mat = df.as_matrix()
    #mat[abs(mat)<3e-10] = 0.001
    mat = np.clip(mat, -2, 1E6, out=mat)
#    mat = 2**mat
    #sns.set(font_scale=2)
    #plt.subplots_adjust(bottom=2, top=3)
    fig = plt.figure()
    print("Making gene by dataset heatmap (includes clustering)...")
    res = sns.clustermap(
         mat,
         yticklabels=True,
         vmin=0,
         cmap=cmap,
         #metric='correlation',  # For columns
         metric='euclidean',
         method='average',
         col_cluster=False,
         )
    print("Outputing heatmap figure...")
    plt.setp(res.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(res.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    
    pos_labels = set(
        ['gld-1','htp-1','htp-2','mpk-1','him-3','fbf-1','lip-1','syp-2',
         'fbf-2','fog-1','fem-3','syp-3','gld-3','fog-3','egl-4'])
    
    def label_if_positive(t):
        #x = df.iloc[int(t)].name
        #print(t)
        if t in pos_labels:
            return t
        else:
            return ''
        

    
    # list_of_rows_in_same_order_as_in_heatmap:
    rowlab = [df.iloc[i].name for i in res.dendrogram_row.reordered_ind]
    
    res.ax_heatmap.yaxis.set_ticklabels(rowlab)
#        [df.columns[int(t._text)] for t in \
#         res.ax_heatmap.xaxis.get_majorticklabels()])

     
    
    #print([x for x in res.ax_heatmap.yaxis.get_majorticklabels()])
    #print([x._text for x in res.ax_heatmap.yaxis.get_majorticklabels()])    
    

    new_labels = [label_if_positive(t._text) for t in res.ax_heatmap.yaxis.get_majorticklabels()]
    print(new_labels)
    res.ax_heatmap.yaxis.set_ticklabels(
        new_labels)
    
    res.ax_heatmap.xaxis.set_ticklabels(
        [df.columns[int(t._text)] for t in \
         res.ax_heatmap.xaxis.get_majorticklabels()])
    
    print('-----')

    #res.dendrogram_col.label = df.columns  # Labels diagram.
    

    # col labels in same order as heatmap
    #ylab = [df.columns[i] for i in res.dendrogram_col.reordered_ind]
    #res.dendrogram_col.label = ylab
    #df = df[ylab]
    
    cols = [t._text for t in res.ax_heatmap.xaxis.get_majorticklabels()]
    
    df = df.loc[rowlab]  
    df = df[cols]
    
    ordered = df.copy()
    plt.xticks(rotation=90)
    #print res.__dict__
    #print res.ax_heatmap.__dict__['figure'].get_size_inches()
    #print res.ax_heatmap.__dict__.keys()
    res.ax_heatmap.__dict__['figure'].set_figwidth(7)
    res.ax_heatmap.__dict__['figure'].set_figheight(10)
    #print res.ax_heatmap.__dict__['figure'].get_size_inches()
    plt.savefig('figs/Fig 3 Genes vs sample heatmap.pdf', #pad_inches=4,
                )#bbox_inches='tight')
    plt.clf(); plt.close()
    relab_dict = recover_label()
    ordered.columns = [relab_dict[x] for x in ordered.columns]
    writer = pandas.ExcelWriter('tables/File S5 Complex frequencies.xls')
    ordered.to_excel(
        writer, sheet_name='log2 reads-per-gene', index=True)
    writer.save()
    ordered.to_csv('tables/File S5 heatmap.txt', sep='\t', encoding='utf-8')
    #cl = res.dendrogram_col.linkage
    #col_clusters = scipy.cluster.hierarchy.fcluster(cl, 1)
    #res.dendrogram_col.label = col_clusters
    return ordered,


def linkages(df):
    correlations = df.corr()
    correlations_array = np.asarray(df.corr())
    row_linkage = hierarchy.linkage(
        distance.pdist(correlations_array), method='average')
    col_linkage = hierarchy.linkage(
        distance.pdist(correlations_array.T), method='average')
    ax = sns.clustermap(
        correlations, row_linkage=row_linkage, col_linkage=col_linkage,
        method="average", figsize=(13, 13))
    ax.dendrogram_col.label = df.corr().columns  # Labels diagram.
    ax.dendrogram_row.label = df.corr().index  # Labels diagram
    #print ax
    #print"^linkages return val"
    return ax


class tableMaker(object):
    
    def __init__(self, heatmapMakerObject):
        self.hm = heatmapMakerObject
        self.load_table()
    
    def load_table(self):
        self.df = self.hm.counts_df
        print("Loaded table into tableMaker object.")
        #print self.df.head(1)
        
    def define_groups(self, _dict):
        self.groups = {}
        for name, _set in _dict.items():
            self.groups[str(name)] = set(_set)
    
    def ave_reads_per_gene_in_given_cols(self, df, col_list):
        ave_reads_per_gene = []
        for dataset_name in col_list:
            ave_reads_per_gene.extend([x for x in df[dataset_name].tolist()])
        return sum(ave_reads_per_gene)/float(len(ave_reads_per_gene))

    def oo_vs_sp_reads_in_gene_as_lists(self, df):
        oo_cols = ['OO FBF_1', 'OO FBF_2', 'OO FBF_3']
        sp_cols = ['SP FBF_1', 'SP FBF_2', 'SP FBF_3']
        reads_per_gene_oo = []
        reads_per_gene_sp = []
        for dataset_name in oo_cols:
            reads_per_gene_oo.extend([
                x for x in df[dataset_name].tolist()])
        for dataset_name in sp_cols:
            reads_per_gene_sp.extend([
                x for x in df[dataset_name].tolist()])
        return reads_per_gene_oo, reads_per_gene_sp


    def oo_vs_sp_reads_in_gene_ratio_and_pvalue_vs_selfdf(self, df):
        reads_per_gene_oo, reads_per_gene_sp = \
            self.oo_vs_sp_reads_in_gene_as_lists(df)
        ratios = []
        reads_per_gene_oo_all, reads_per_gene_sp_all = \
            self.oo_vs_sp_reads_in_gene_as_lists(self.df)
        tups = zip(reads_per_gene_oo, reads_per_gene_sp,
                   reads_per_gene_oo_all, reads_per_gene_sp_all)
        df_ratio, selfdf_ratio = [], []
        df_to_selfdf_ratio_of_ratios = []
        for tup in tups:
            if (tup[0] != 0):
                df_ratio.append(tup[1]/tup[0])
            if (tup[2] != 0):
                selfdf_ratio.append(tup[3]/tup[2])
        res = scs.ttest_ind(df_ratio, selfdf_ratio)
        return np.mean(df_ratio)/np.mean(selfdf_ratio), res.pvalue

    def ttest_for_difference_in_reads_per_gene_in_given_cols(
            self, df_a, df_b, col_list):
        reads_per_gene_a = []
        reads_per_gene_b = []
        for dataset_name in col_list:
            reads_per_gene_a.extend([x for x in df_a[dataset_name].tolist()])
            reads_per_gene_b.extend([x for x in df_b[dataset_name].tolist()])
        res = scs.ttest_ind(reads_per_gene_a, reads_per_gene_b)
        return res.pvalue
    
    def fbf1_vs_fbf2_reads_in_gene_as_lists(self, df):
        fbf1_cols = ['LT FBF1_1', 'LT FBF1_2', 'LT FBF1_3']
        fbf2_cols = ['LT FBF2_1', 'LT FBF2_2', 'LT FBF2_3']
        reads_per_gene_fbf1 = []
        reads_per_gene_fbf2 = []
        for dataset_name in fbf1_cols:
            reads_per_gene_fbf1.extend([
                x for x in df[dataset_name].tolist()])
        for dataset_name in fbf2_cols:
            reads_per_gene_fbf2.extend([
                x for x in df[dataset_name].tolist()])
        return reads_per_gene_fbf1, reads_per_gene_fbf2

    def fbf1_vs_fbf2_reads_in_gene_means(self, df):
        reads_per_gene_fbf1, reads_per_gene_fbf2 = \
            self.fbf1_vs_fbf2_reads_in_gene_as_lists(df)
        return np.mean(reads_per_gene_fbf1), np.mean(reads_per_gene_fbf2)
    
    def fbf1_vs_fbf2_reads_in_gene_ratio_and_pvalue(self, df):
        reads_per_gene_fbf1, reads_per_gene_fbf2 = \
            self.fbf1_vs_fbf2_reads_in_gene_as_lists(df)
        mean_rpg_fbf1 = np.mean(reads_per_gene_fbf1)
        mean_rpg_fbf2 = np.mean(reads_per_gene_fbf2)
        if mean_rpg_fbf1 == 0:
            ratio_fbf2_vs_fbf1_rpg = 0
        else:
            ratio_fbf2_vs_fbf1_rpg = mean_rpg_fbf2/mean_rpg_fbf1
        res = scs.ttest_ind(reads_per_gene_fbf1, reads_per_gene_fbf2)
        return ratio_fbf2_vs_fbf1_rpg, res.pvalue

    def fbf1_vs_fbf2_reads_in_gene_ratio_and_pvalue_vs_selfdf(self, df):
        reads_per_gene_fbf1, reads_per_gene_fbf2 = \
            self.fbf1_vs_fbf2_reads_in_gene_as_lists(df)
        ratios = []
        mean_rpg_fbf1 = np.mean(reads_per_gene_fbf1)
        mean_rpg_fbf2 = np.mean(reads_per_gene_fbf2)
        reads_per_gene_fbf1_all, reads_per_gene_fbf2_all = \
            self.fbf1_vs_fbf2_reads_in_gene_as_lists(self.df)
        tups = zip(reads_per_gene_fbf1, reads_per_gene_fbf2,
                   reads_per_gene_fbf1_all, reads_per_gene_fbf2_all)
        df_ratio, selfdf_ratio = [], []
        df_to_selfdf_ratio_of_ratios = []
        for tup in tups:
            if (tup[0] != 0):
                df_ratio.append(tup[1]/tup[0])
            if (tup[2] != 0):
                selfdf_ratio.append(tup[3]/tup[2])
        res = scs.ttest_ind(df_ratio, selfdf_ratio)
        return np.mean(df_ratio)/np.mean(selfdf_ratio), res.pvalue

    def make_table(self, output_dir='tables/'):
        self.define_programs()
        self.df['Program'] = [self.program.get(x, '') for x in self.df['Wormbase ID'].tolist()]
        known_items = set(self.df.index)
        # Subset to the names we have in the heatmapMaker object.
        for name, _set in self.groups.items():
            self.groups[name] = self.groups[name] & known_items
        self.groups['All {0} RNAs'.format(len(known_items))] = set(known_items)
        def as_int(x):
            if x:
                return 1
            else:
                return 0
        # Fill in the info for the table.
        # Programs: "Oogenic only", "Spermatogenic only",
        # and 'Oogenic and Spermatogenic'
        self.define_targets()
        table = collections.defaultdict(dict)
        for block_name, _set in self.groups.items():
            table[block_name] = {}
            table[block_name]['# RNAs'] = len(_set)
            #print dict(self.df['Program'].value_counts())
            self.df['in_block'] = [(x in _set) for x in self.df.index]
            this_block = self.df[self.df['in_block']].copy()
            fbf1_mean, fbf2_mean = self.fbf1_vs_fbf2_reads_in_gene_means(this_block)
            table[block_name]['Mean reads/gene 20C FBF-1'] = fbf1_mean
            table[block_name]['Mean reads/gene 20C FBF-2'] = fbf2_mean
            ratio, p = self.fbf1_vs_fbf2_reads_in_gene_ratio_and_pvalue(this_block)
            table[block_name]['Reads/gene ratio FBF-2/FBF-1'] = ratio
            table[block_name]['Reads/gene FBF-1 vs FBF-2 t-test p value'] = p
            ratio, p = self.fbf1_vs_fbf2_reads_in_gene_ratio_and_pvalue_vs_selfdf(
                this_block)
            table[block_name][
'FBF-2/FBF-1 RPG ratio vs all RNAs in heatmap'] = ratio
            table[block_name][
'FBF-2/FBF-1 RPG ratio vs all RNAs in heatmap p value'] = p            


            table[block_name]['Mean reads/gene in oogenic CLIP'] = \
                self.ave_reads_per_gene_in_given_cols(
                this_block, ['OO FBF_1', 'OO FBF_2', 'OO FBF_3'])
            table[block_name]['Mean reads/gene in spermatogenic CLIP'] = \
                self.ave_reads_per_gene_in_given_cols(
                this_block, ['SP FBF_1', 'SP FBF_2', 'SP FBF_3'])

            ratio, p = self.oo_vs_sp_reads_in_gene_ratio_and_pvalue_vs_selfdf(
                this_block)
            table[block_name]['OO/SP RPG ratio vs all RNAs in heatmap'] = ratio
            table[block_name]['OO/SP RPG ratio vs all RNAs in heatmap p value'] = p
            
            table[block_name][
'P value (t-test) for reads/gene in oogenic CLIP vs all RNAs in heatmap'] = \
                self.ttest_for_difference_in_reads_per_gene_in_given_cols(
                self.df, this_block, ['OO FBF_1', 'OO FBF_2', 'OO FBF_3'])
            table[block_name][
'P value (t-test) for reads/gene in spermatogenic CLIP vs all RNAs in heatmap'] = \
                self.ttest_for_difference_in_reads_per_gene_in_given_cols(
                self.df, this_block, ['SP FBF_1', 'SP FBF_2', 'SP FBF_3'])
            table[block_name][
'P value (t-test) for reads/gene in all 25C CLIP vs all RNAs in heatmap'] = \
                self.ttest_for_difference_in_reads_per_gene_in_given_cols(
                self.df, this_block, ['SP FBF_1', 'SP FBF_2', 'SP FBF_3',
                                      'OO FBF_1', 'OO FBF_2', 'OO FBF_3'])
            program_counts = dict(this_block['Program'].value_counts())
            perc_program = collections.defaultdict(float)
            for k, v in program_counts.items():
                perc_program[k] = float(int(
                    10000 * float(v)/float(sum(program_counts.values()))))/100.0
#            n_oo = [as_int(x in self.program['Oogenic only']) \
#                for x in _set]
#            n_sp = [as_int(x in self.program['Spermatogenic only']) \
#                for x in _set]
#            n_ne = [as_int(x in self.program['Oogenic and Spermatogenic']) \
#                for x in _set]
            table[block_name]['% only in oogenic program'] = \
                perc_program['Oogenic only']
            table[block_name]['% only in spermatogenic program'] = \
                perc_program['Spermatogenic only']
            table[block_name]['% in oogenic and spermatogenic programs'] = \
                perc_program['Oogenic and Spermatogenic']
            with_sp_peaks_in_block = set(self.peaks['sp_both'].df.gene_name) & _set
            with_oo_peaks_in_block = set(self.peaks['oo_both'].df.gene_name) & _set
            self.peaks['sp_both'].df['in_block'] = [(x in _set) for x in 
                self.peaks['sp_both'].df.gene_name]
            self.peaks['oo_both'].df['in_block'] = [(x in _set) for x in 
                self.peaks['oo_both'].df.gene_name]
            self.sp_peaks_in_block = \
                self.peaks['sp_both'].df[self.peaks['sp_both'].df['in_block']].copy()
            self.oo_peaks_in_block = \
                self.peaks['oo_both'].df[self.peaks['oo_both'].df['in_block']].copy()
            table[block_name]['% with peaks in 25C sp. clip'] = float(int(
                10000 * len(with_sp_peaks_in_block)/len(_set)
                )/100.)
            table[block_name]['% with peaks in 25C oo. clip'] = float(int(
                10000 * len(with_oo_peaks_in_block)/len(_set)
                )/100.)
            both_sp_and_oo_peaks = set(self.peaks['sp_both'].df.gene_name) & \
                set(self.peaks['oo_both'].df.gene_name)
            table[block_name]['% with peaks in 25C oo. and 25 sp. clip'] = float(int(
                10000 * len(both_sp_and_oo_peaks & _set)/len(_set)
                )/100.)
            table[block_name]['Average peak height (sp)'] = \
                self.sp_peaks_in_block['exp_reads'].mean()
            table[block_name]['Average peak height (oo)'] = \
                self.oo_peaks_in_block['exp_reads'].mean()
#        table['All sp. targets']['Average peak height (sp)'] = \
#            self.peaks['sp_both'].df['exp_reads'].mean()
#        table['All oo. targets']['Average peak height (oo)'] = \
#            self.peaks['oo_both'].df['exp_reads'].mean()
#        table['All {0} RNAs'.format(len(self.df.index))][
#                'Mean reads/gene in oogenic CLIP'] = \
#                self.ave_reads_per_gene_in_given_cols(
#                this_block, ['OO FBF_1', 'OO FBF_2', 'OO FBF_3'])
#        table[][
#                'Mean reads/gene in spermatogenic CLIP'] = \
#                self.ave_reads_per_gene_in_given_cols(
#                this_block, ['SP FBF_1', 'SP FBF_2', 'SP FBF_3'])
#for k, _dict in table.items():
#            print "***"
#            print k
            #for k, v in _dict.items():
            #    print k, v
        self.tabledf = pandas.DataFrame(table)
        print(self.tabledf)
        self.tabledf.to_csv('tables/Table_1_block_stats.txt', sep='\t', index=False)
        writer = pandas.ExcelWriter('tables/Table 1 block stats.xls')
        self.tabledf.to_excel(writer)
        writer.save()

    def define_targets(self):
        self.peaks = {}
        from vennMaker import fname_to_label, label_to_fname
        for dataset in [
            'oo_both', 'sp_both', 'old_fbf1_to_fbf2_n2', 'old_fbf2']:
            fname = label_to_fname[dataset]
            self.peaks[dataset] = peaksList()
            self.peaks[dataset].read_csv(fname)

    def define_programs(self):
        self.original_df = self.df.copy()
        self.df['gene_name'] = self.df.index
        self.pkl = peaksList(dataframe=self.df, gene_name_col='gene_name')
#        self.volcano.gl_rnaseq()
#        self.volcano.read_clip_deseq_csv('tables/lt_fbf1_and_2_vs_ht_fbf_deseq.txt')
        self.pkl.read_sp_vs_oo_as_programs()
        self.df = self.pkl.df.copy()
        self.program = self.pkl.program
        
#fname_to_label = {
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/old_fbf2.txt': 'old_fbf2',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/old_fbf1.txt': 'old_fbf1',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/oo_both.txt': 'oo_both',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/oo_fbf1.txt': 'oo_fbf1',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/oo_fbf2.txt': 'oo_fbf2',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/sp_both.txt': 'sp_both',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/sp_fbf1.txt': 'sp_fbf1',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/sp_fbf2.txt': 'sp_fbf2',
#'/groups/Kimble/Common/fbf_celltype/combined_filtered/old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
#}
#label_to_fname = dict(zip(fname_to_label.values(), fname_to_label.keys()))

if __name__ == '__main__':
    print("Creating heatmapMaker() object()...")
    v = heatmapMaker()
    print("Created heatmapMaker() object(). Loading peaks files...")
    v.load(label_to_fname)  # Global from vennMaker used.
    plt.clf()
    plt.rc('font', size=1)
    fig, ax = plt.subplots()
    df = v.load_counts_file()
    tableM = tableMaker(v)
    tableM.define_groups({'I': blocki, 'II': blockii, 'III': blockiii})
    tableM.make_table()
    #sys.exit()
    # Create the heatmap figure.
    with sns.plotting_context("paper", font_scale=0.8):
        plt.clf()
        plt.close()
        g = make_fig(df)
        plt.clf()
        plt.close()
    # Create the database vs database correlation block figure.
    with sns.plotting_context("paper", font_scale=1):
        #df = v.load_counts_file()
        df = df[['SP FBF_1', 'SP FBF_2', 'SP FBF_3',
             'OO FBF_1', 'OO FBF_2', 'OO FBF_3',
             'LT FBF1_1', 'LT FBF1_2', 'LT FBF1_3',
             'LT FBF2_1', 'LT FBF2_2', 'LT FBF2_3']]
        col_relabel_d = col_relabel()
        df.columns = [col_relabel_d[x] for x in df.columns]
        pwdf = df.corr('spearman')
        mask = np.zeros_like(pwdf)
        mask[np.triu_indices_from(mask)] = True
        np.place(mask, mask<1, 2)
        np.place(mask, mask==1, 0)
        np.place(mask, mask==2, 1)
        cmap = sns.cubehelix_palette(
            len(df.index), light=.7, dark=.2, reverse=False,
            start=1, rot=-2, as_cmap=True)
        cmap = plt.get_cmap('Greys')
        with sns.axes_style("white"):
            ax = sns.heatmap(pwdf, cmap=cmap, annot=True,
                            #col_cluster=False,
                            #row_cluster=False,
                            mask=mask)
#        ax.set_xticklabels(rotation=0)
#        plt.yticks(rotation=45)
#        plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
#        locs, labels = plt.yticks()
#        plt.setp(labels, rotation=45)
        print('--')
        plt.savefig('figs/Fig 2 Spearman corr.pdf')
#    mat = df.as_matrix()
    nah_dog = '''
    # clustering
    import sklearn
    from sklearn import cluster
    km = cluster.KMeans(n_clusters=8)
    km.fit(mat)
    # Get cluster assignment labels
    labels = km.labels_
    # Format results as a DataFrame
    results = pandas.DataFrame([df.index,labels]).T
    print results
    print results.columns
    df['cluster'] = results[results.columns[1]].tolist()
    df.sort(columns=['cluster'], inplace=True)
    del df['cluster']
    mat = df.as_matrix()
    '''
    # output

#    ann_df = df.copy()
#    g = seaborn.clustermap(df)#, method='average')
    # _____ Annotate df with clusters, compute linkages ______
    # rlinks is a n by 4 array constituting a linkage matrix.
    # Each len=4 row represents the combination of the two branches in
    # cols 0 and 1, at distance in col 2, combining the total number
    # of primary observations in col 3. Branches with index<=len(df.columns)
    # are the index of a primary observation.
#    rlinks = scipy.cluster.hierarchy.linkage(
#        distance.pdist(df.corr()), method='average')
    # den is a dict with info to depict a linkage matrix as a dendrogram.
    # color_list holds labels based off of distance thresholding.
#    rden = scipy.cluster.hierarchy.dendrogram(
#        rlinks, labels=df.columns)
#    rden['color_list'].append(rden['color_list'][-1])  # Not sure what's going on here.
    # Dict of column label to cluster color.
#    rto_color = dict(zip(rden['ivl'], rden['color_list']))
#    ann_df.loc['col cluster'] = [rto_color[x] for x in df.columns]
    # Now rows.
#    clinks = scipy.cluster.hierarchy.linkage(
#        distance.pdist(df.corr().T), method='average')
#    cden = sp.cluster.hierarchy.dendrogram(
#        clinks, labels=df.index)
    # Dict of column label to cluster color.
#    cto_color = dict(zip(cden['ivl'], cden['color_list']))
#    ann_df['row cluster'] = ['k'] + [to_color[x] for x in df.index]
    # Alt.
#    col_lut = dict(zip(map(str, df.columns), col_pal))
    #cols = df.columns.get_level_values("network")
#    cols_colors = pd.Series(df.columns).map(col_lut)
#    cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    #pal = sns.palplot(sns.diverging_palette(150, 275, s=80, l=55, n=9))
    # _____ Make figure ______

