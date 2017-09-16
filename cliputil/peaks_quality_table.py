from __future__ import division
import scipy.stats as scs
from peaksList import *
from collections import defaultdict
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import glob
#from volcanoMaker import *
import os
import sys
import re
from blocks import *


def get_fasta(df, top_n=200):
    fasta = ''
    top = df.copy()
    top.sort(columns='exp_reads', ascending=False, inplace=True)
    top.drop_duplicates(subset='exp_reads', inplace=True)
    top = top.head(top_n)
    for t in zip(top.gene_name, top.left, top.seq):
        fasta += '>{0}_{1}\n{2}\n'.format(t[0], t[1], t[2])
    return fasta


def is_target(x):
    if (x == ''): return False
    if type(x) == type(''):
        return True
    return False


def make_bar(pkl_list):
    #import seaborn as sns
    #sns.set_context('paper')
    #sns.set_style('ticks')
    cell_text_l = []
    scores_df = []
    common_targets = set()
    def frac_pum2(df):
        df['overlap'] = [
            is_target(x) for x in df['PUM2 target HGNC symbol'].tolist()]
        overlapdf = df[df['overlap']].copy()
        overlap = set(overlapdf.gene_name)
        return 100 * len(overlap)/len(df.index)
    positive_control_genes = set(
"""gld-1,htp-1,htp-2,mpk-1,him-3,fbf-1,lip-1,syp-2,fbf-2,fog-1,fem-3,syp-3,gld-3,fog-3,egl-4""".split(','))
    def frac_pos(genes):
        return len(set(genes) & positive_control_genes)/len(positive_control_genes)
    for pkl in pkl_list:
        _score = get_score(pkl.df)
        print pkl.df
        top = get_top_n(pkl.df, n=500)
        perc_fbe_top_n = get_score(top)
        has_pos = set(pkl.df.gene_name) & positive_control_genes
        has_pos_top = set(top.gene_name) & positive_control_genes
        scores_df.append(
            {'Name': pkl.name,
            '% FBE': _score,
            '% Verified targets': frac_pos(pkl.df.gene_name),
            '% PUM2 targets': frac_pum2(pkl.df)
            })
        scores_df.append({
            'Name': u"{0} (top 500)".format(pkl.name),
            '% FBE': perc_fbe_top_n,
            '% Verified targets': frac_pos(top.gene_name),
            '% PUM2 targets': frac_pum2(top)
            })
    oo = [x for x in pkl_list if x.name==to_name['oo_both']]
    oo = oo[0].df
    sp = [x for x in pkl_list if x.name==to_name['sp_both']]
    sp = sp[0].df
    joint = set(oo['gene_name'].tolist()) & set(sp['gene_name'].tolist())
    oo['oo_and_sp'] = [(x in joint) for x in oo['gene_name'].tolist()]
    jointdf = oo[oo['oo_and_sp']]
    scores_df.append(
        {'Name': "OO/SP common targets",
         "% FBE": get_score(jointdf),
         "% Verified targets": frac_pos(jointdf.gene_name),
         "% PUM2 targets": frac_pum2(jointdf)})
    top = get_top_n(jointdf, n=500)
    scores_df.append({
            'Name': u"OO/SP common targets (top 500)",
            '% FBE': get_score(top),
            '% Verified targets': frac_pos(top.gene_name),
            '% PUM2 targets': frac_pum2(top)
            })
    
    w_blocks = scores_df[:]#dict(zip(scores_df.keys(), scores_df.values()))
    sp['Block I'] = [(x in blocki) for x in sp.gene_name]
    sp['Block II'] = [(x in blockii) for x in sp.gene_name]
    sp['Block III'] = [(x in blockiii) for x in sp.gene_name]
    blocki_df = sp[sp['Block I']].copy()
    blockii_df = sp[sp['Block II']].copy()
    blockiii_df = sp[sp['Block III']].copy()
    w_blocks.append({
        'Name': 'Block I',
        '% PUM2 targets': frac_pum2(blocki_df)})
    w_blocks.append({
        'Name': 'Block II',
        '% PUM2 targets': frac_pum2(blockii_df)})
    w_blocks.append({
        'Name': 'Block III',
        '% PUM2 targets': frac_pum2(blockiii_df)})
    w_blocks = pandas.DataFrame(w_blocks)
    w_blocks.index = w_blocks.Name
    w_blocks = w_blocks.reindex([

u'oo FBF-1 (20' + u'\xb0' + u'C)',
u'oo FBF-1 (20' + u'\xb0' + u'C) (top 500)',
u'oo FBF-2 (20' + u'\xb0' + u'C)',
u'oo FBF-2 (20' + u'\xb0' + u'C) (top 500)',
u'sp FBF (25' + u'\xb0' + u'C)',
u'sp FBF (25' + u'\xb0' + u'C) (top 500)',
u'oo FBF (25' + u'\xb0' + u'C)',
u'oo FBF (25' + u'\xb0' + u'C) (top 500)',
u'OO/SP common targets',
u'OO/SP common targets (top 500)',
'Block I', 'Block II', 'Block III'
])
    scores_df = pandas.DataFrame(scores_df)
    scores_df.index = scores_df.Name
    scores_df = scores_df.reindex([
u'oo FBF-1 (20' + u'\xb0' + u'C)',
u'oo FBF-1 (20' + u'\xb0' + u'C) (top 500)',
u'oo FBF-2 (20' + u'\xb0' + u'C)',
u'oo FBF-2 (20' + u'\xb0' + u'C) (top 500)',
u'sp FBF (25' + u'\xb0' + u'C)',
u'sp FBF (25' + u'\xb0' + u'C) (top 500)',
u'oo FBF (25' + u'\xb0' + u'C)',
u'oo FBF (25' + u'\xb0' + u'C) (top 500)',
u'OO/SP common targets',
u'OO/SP common targets (top 500)'

#u'oo FBF-1 (25' + u'\xb0' + u'C)',
#u'oo FBF-2 (25' + u'\xb0' + u'C)',
])
    
    plot_col(
        scores_df, '% FBE', fname='figs/peak_quality_perc_fbe_grouped.pdf')
    plot_col(
        scores_df, '% Verified targets',
        fname='figs/peak_quality_perc_pos_grouped.pdf')
    plot_col(
        scores_df, '% PUM2 targets',
        fname='figs/peak_quality_perc_pum2_grouped.pdf')
    plot_col(
        w_blocks, '% PUM2 targets',
        fname='figs/peak_quality_perc_pum2_w_blocks_grouped.pdf')
    #row_labels = []
    for pkl in pkl_list:
        #row_labels.append(pkl.name)
        tup = build_row(pkl, return_text=False)
        cell_text_l.append(dict(zip(tup[0], tup[1])))
        cell_text_l[-1]['Name'] = pkl.name
    df = pandas.DataFrame(cell_text_l)
    df.index = df.Name
    all_col = df.copy()
    df = df.reindex([
u'oo FBF-1 (20' + u'\xb0' + u'C)',
u'oo FBF-2 (20' + u'\xb0' + u'C)',
u'sp FBF (25' + u'\xb0' + u'C)',
u'oo FBF (25' + u'\xb0' + u'C)',
#u'oo FBF-1 (25' + u'\xb0' + u'C)',
#u'oo FBF-2 (25' + u'\xb0' + u'C)',
])
    sp = all_col.reindex([
u'sp FBF-1 (25' + u'\xb0' + u'C)',
u'sp FBF-2 (25' + u'\xb0' + u'C)',
])
    plot_col(
        df, '# RNAs', fname='figs/peak_quality_num_rnas.pdf')
    plot_col(
        df, '% FBE', fname='figs/peak_quality_perc_fbe.pdf')
    plot_col(
        df, '% Verified targets', fname='figs/peak_quality_perc_pos.pdf')
    plot_col(
        df, '% FBE (top 500 targets)',
        fname='figs/peak_quality_perc_fbe_top.pdf')
    plot_col(
        df, '% Verified targets',
        fname='figs/peak_quality_perc_pos_control.pdf')
    
    plot_col(
        sp, '# RNAs', fname='figs/peak_quality_sp_num_rnas.pdf')
    plot_col(
        sp, '% FBE', fname='figs/peak_quality_sp_perc_fbe.pdf')
    plot_col(
        sp, '% FBE (top 500 targets)',
        fname='figs/peak_quality_sp_perc_fbe_top.pdf')
    plot_col(
        sp, '% Verified targets',
        fname='figs/peak_quality_sp_perc_pos_control.pdf')


def plot_col(df, col, ylabel=None,
             fname='figs/peak_quality_bar.pdf'):
    plt.clf(); plt.close()
    width = 0.7
    if ylabel is None: ylabel = col
    fig, ax = plt.subplots(figsize=(2,2))
    if len(df.index) == 4:
        xpos = [1, 2, 4, 5, ]#7, 8]
    elif len(df.index) == 8:
        xpos = [1,2, 3.6,4.6, 6.2,7.2, 8.8,9.8]
    elif len(df.index) == 10:
        xpos = [1,2, 3.6,4.6, 6.2,7.2, 8.8,9.8, 11.2,12.2]
    elif len(df.index) == 13:
        xpos = [1,2, 3.6,4.6, 6.2,7.2, 8.8,9.8, 11.2,12.2, 13.6,14.6,15.6]
    elif len(df.index) == 2:
        xpos = [1.3,1.7]
        width = 0.3
    else: xpos = np.arange(len(df.index))
    plt.bar(xpos, df[col].tolist(),
            width=width,
            align='center',
            facecolor='grey',
            edgecolor='black')
    print '---'
#    labels = [from_name[x] for x in df.index]
    #if len(df.index) == 6:
    ax.set_xticks(xpos)
    ax.set_xticklabels(tuple(df.index))
    ax.set_ylabel(col)
    ax = despine(ax)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    #ax.despine()
    plt.savefig(fname, bbox_inches='tight')
    plt.clf(); plt.close()

def despine(ax1):
    #plt.tick_params(which='both', direction='out')
    #plt.tick_params(which='both', right=False)
    ax1.get_yaxis().set_tick_params(which='both', direction='out')
    ax1.get_yaxis().set_tick_params(which='both', right=False)
    ax1.get_xaxis().set_tick_params(which='both', direction='out')
    ax1.get_xaxis().set_tick_params(which='both', top=False)
#    ax1.set_ylabel(ylab)
#    ax1.set_xlabel(xlab)
    ax1.spines['top'].set_visible(False)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    return ax1

def make_table(pkl_list):
    cell_text_d = {}
    #row_labels = []
    for pkl in pkl_list:
        #row_labels.append(pkl.name)
        cell_text_d[pkl.name] = build_row(pkl)[1]
    cell_text = [cell_text_d[x] for x in datasets if x in cell_text_d]
    row_labels = [x for x in datasets if x in cell_text_d]
    col_labels = ['# RNAs', '% Verified targets',
                  '% FBE', '% FBE (top 500 targets)']
    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    _table = ax.table(
        cellText=cell_text,
        rowLabels=row_labels,
        colLabels=col_labels,
        )
    plt.savefig('figs/peak_quality_table.pdf',
            bbox_inches='tight')
    plt.clf(); plt.close()
    print cell_text


def get_score(_df):
    motif = re.compile('TGT\w\w\wAT', re.IGNORECASE)
    scores = []
    for seq in _df.seq:
        if motif.search(seq):
            scores.append(1)
        else:
            scores.append(0)
    perc_fbe = 100 * sum(scores)/len(scores)
    return perc_fbe


def get_top_n(_df, n=500):
    top = _df.copy()
    top.sort(columns='exp_reads', ascending=False, inplace=True)
    top.drop_duplicates(subset='exp_reads', inplace=True)
    top = top.head(n).copy()
    return top

def build_row(pkl, return_text=True):    
    perc_fbe = get_score(pkl.df)
    top = get_top_n(pkl.df, n=500)
    perc_fbe_top_n = get_score(top)
    header = ['# RNAs', '% Verified targets',
              '% FBE', '% FBE (top 500 targets)']
    positive_control_genes = set(
        """gld-1,htp-1,htp-2,mpk-1,him-3,fbf-1,lip-1,syp-2,fbf-2,fog-1,fem-3,syp-3,gld-3,fog-3,egl-4""".split(','))
    has_pos = set(pkl.df.gene_name) & positive_control_genes
    if return_text:
        cell_text = [
            str(len(set(pkl.df.gene_name))),
            '%.3f' % float(len(has_pos)/len(positive_control_genes)),
            '%.1f' % perc_fbe,
            '%.1f' % perc_fbe_top_n]
    else:
        cell_text = [
            len(set(pkl.df.gene_name)),
            len(has_pos)/len(positive_control_genes),
            perc_fbe,
            perc_fbe_top_n]
    return header, cell_text
        
to_name = {
    'old_fbf1': 'Discard',
    'old_fbf1_to_fbf2_n2': u'oo FBF-1 (20' + u'\xb0' + u'C)',
    'old_fbf2': u'oo FBF-2 (20' + u'\xb0' + u'C)',
    'oo_both': u'oo FBF (25' + u'\xb0' + u'C)',
    'oo_fbf1': u'oo FBF-1 (25' + u'\xb0' + u'C)',
    'oo_fbf2': u'oo FBF-2 (25' + u'\xb0' + u'C)',
    'sp_both': u'sp FBF (25' + u'\xb0' + u'C)',
    'sp_fbf1': u'sp FBF-1 (25' + u'\xb0' + u'C)',
    'sp_fbf2': u'sp FBF-2 (25' + u'\xb0' + u'C)',
}
from_name = dict(zip(to_name.values(), to_name.keys()))
datasets = [
u'oo FBF-1 (20' + u'\xb0' + u'C)',
u'oo FBF-2 (20' + u'\xb0' + u'C)',
u'sp FBF (25' + u'\xb0' + u'C)',
u'oo FBF (25' + u'\xb0' + u'C)',
u'sp FBF-1 (25' + u'\xb0' + u'C)',
u'sp FBF-2 (25' + u'\xb0' + u'C)',
u'oo FBF-1 (25' + u'\xb0' + u'C)',
u'oo FBF-2 (25' + u'\xb0' + u'C)',
]

from __init__ import mkdir
pkl_list = []
for f in glob.glob('combined_filtered/ortho_info_added/*txt'):
    bname = os.path.basename(f).partition('.txt')[0]
    if bname in  to_name: bname = to_name[bname]
    #print ">{0}".format(bname)
    p = peaksList(name=bname)
    p.read_csv(f)
    pkl_list.append(p)
    mkdir('data/fasta/')
    top_n = 500
    f = open('data/fasta/top{0}_{1}'.format(
        top_n,
        os.path.basename(f).partition('.txt')[0] + '.fa'), 'w')
    #f.encode('utf8')
    f.write(get_fasta(p.df, top_n=top_n))
    f.close()
make_table(pkl_list)
make_bar(pkl_list)
