from __future__ import division
import scipy.stats as scs
from peaksList import *
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

def gc(pkl):
    pkl.df = pkl.df.head(2000)
    # pkl = cliputil.peaksList
    def per_nt(a):
        a = a.upper()
        d = defaultdict(int)
        for n in list(a): d[n] += 1
        for n in d: d[n] = d[n]/float(len(a))
        return d
    dists = []
    tups = []
    for seq, height in zip(
        pkl.df['seq'].tolist(), pkl.df['exp_reads'].tolist()):
        freqs = per_nt(seq)
        tups.append((np.mean([freqs['G'], freqs['C']]), np.log10(height)))
        dists.append(freqs)
    per_gc = [tup[0] for tup in tups]
    df = pd.DataFrame(dists)
    df = df[['A', 'T', 'G', 'C']]
    tups = pandas.DataFrame(tups, columns=['Percent GC', 'Height'])
    print tups
    #print df
    print "% GC mean {0}".format(np.mean(per_gc))
    return per_gc, tups, df
sk = '''
ga['ht_both'] = cliputil.get_a_bedgraph(
    'individual_clip/oo_both/bedgraph_norm/control_oo_+.wig',
    'individual_clip/oo_both/bedgraph_norm/control_oo_-.wig')
#ga['lt_fbf1'] = cliputil.get_a_bedgraph(
#    'individual_clip/old_fbf1_to_fbf2_n2/bedgraph_norm/control_n2_+.wig',
#    'individual_clip/old_fbf1_to_fbf2_n2/bedgraph_norm/control_n2_-.wig')
ga['ht_fbf1'] = cliputil.get_a_bedgraph(
    'individual_clip/oo_fbf1/bedgraph_norm/control_oo_+.wig',
    'individual_clip/oo_fbf1/bedgraph_norm/control_oo_-.wig')
    '''
per_gcs = {}
tups = {}
p_oo = peaksList(name='ht_both')
p_oo.read_csv('combined_filtered/oo_both.txt')
per_gcs['oo_both'], tups['oo_both'], _ = gc(p_oo)
#p.read_csv('combined_filtered/oo_fbf1.txt')
#per_gcs['oo_fbf1'], _ = gc(p)
#p.read_csv('combined_filtered/oo_fbf2.txt')
#per_gcs['oo_fbf2'], _ = gc(p)
p_fbf1 = peaksList(name='FBF1')
p_fbf1.read_csv('combined_filtered/old_fbf1_to_fbf2_n2.txt')
per_gcs['lt_fbf1'], tups['lt_fbf1'], _ = gc(p_fbf1)
p_fbf2 = peaksList(name='FBF2')
p_fbf2.read_csv('combined_filtered/old_fbf2.txt')
per_gcs['lt_fbf2'], tups['lt_fbf2'], _ = gc(p_fbf2)

print p_fbf1.df.columns
low_temp = set(p_fbf1.df['gene_name'].tolist()) | set(p_fbf2.df['gene_name'].tolist())
new_targs = set(p_oo.df['gene_name'].tolist()) - low_temp
p_oo.df['keep'] = [(x in new_targs) for x in p_oo.df.gene_name]
p_oo.df = p_oo.df[p_oo.df['keep']].copy()
per_gcs['only_high'], tups['only_high'], _ = gc(p_oo)

for k1, v1 in per_gcs.items():
    for k2, v2 in per_gcs.items():
        if k1 >= k2: continue
        print k1, ' vs ', k2,
        print " Means {0:.3g} vs {1:.3}".format(np.mean(v1), np.mean(v2)),
        print "n={0}, n={1}".format(len(v1), len(v2))
        print scs.mannwhitneyu(v1, v2)
with sns.plotting_context('paper'):
    sns.set_style('ticks')
    mpl.rc("figure", figsize=(12, 4))
    fig, ax = plt.subplots(ncols=3)
    ax[0] = sns.kdeplot(
        tups['oo_both']['Percent GC'], tups['oo_both']['Height'],
        cmap='Reds', shade=False, shade_lowest=False, alpha=1,
        ax=ax[0])#, aspect=2, size=2)
    ax[0].set_title(
        '25C Oo. FBF-1/FBF-2 (n={0})'.format(len(per_gcs['oo_both'])))
    ax[0].set_ylabel('Peak height (log10 reads)')
    ax[1] = sns.kdeplot(
        tups['lt_fbf1']['Percent GC'], tups['lt_fbf1']['Height'],
        cmap='Blues', shade=False, shade_lowest=False, alpha=1,
        ax=ax[1])
    ax[1].set_title(
        '20C Oo. FBF-1 (n={0})'.format(len(per_gcs['lt_fbf1'])))
    ax[1].set_ylabel('Peak height (log10 reads)')

    ax[2] = sns.kdeplot(
        tups['lt_fbf2']['Percent GC'], tups['lt_fbf2']['Height'],
        cmap='Greens', shade=False, shade_lowest=False, alpha=1,
        ax=ax[2])
    ax[2].set_title(
        '20C Oo. FBF-2 (n={0})'.format(len(per_gcs['lt_fbf2'])))
    ax[2].set_ylabel('Peak height (log10 reads)')
    print "Spearman corr (gc %, height):"
    print "\t 25deg OO {0}".format(
        scs.spearmanr(tups['oo_both']['Percent GC'], tups['oo_both']['Height']))
    print "\t 20deg FBF-1 {0}".format(
        scs.spearmanr(tups['lt_fbf1']['Percent GC'], tups['lt_fbf1']['Height']))
    print "\t 20deg FBF-2 {0}".format(
        scs.spearmanr(tups['lt_fbf2']['Percent GC'], tups['lt_fbf2']['Height']))
    print "Pearson corr (gc %, height):"
    print "\t 25deg OO {0}".format(
        scs.pearsonr(tups['oo_both']['Percent GC'], tups['oo_both']['Height']))
    print "\t 20deg FBF-1 {0}".format(
        scs.pearsonr(tups['lt_fbf1']['Percent GC'], tups['lt_fbf1']['Height']))
    print "\t 20deg FBF-2 {0}".format(
        scs.pearsonr(tups['lt_fbf2']['Percent GC'], tups['lt_fbf2']['Height']))

#        
#    ax = sns.distplot(
#        [x[0] for x in tups['lt_fbf1']],
##        [x[1] for x in tups['lt_fbf1']],
 #       label='20C Oo. FBF-1 (n={0})'.format(len(per_gcs['lt_fbf1'])))
#    ax.set_xlabel('% GC in peak')
#    ax.set_ylabel('Peak height (max reads)')
    sns.despine()
    plt.legend()
    plt.savefig('figs/Fig S3 gc_richness_distplot_height_v_gc_fbf1_v_ht_fbf.pdf',
                bbox_inches='tight')
    plt.clf()
    plt.close()

with sns.plotting_context('paper'):
    sns.set_style('ticks')
    ax = sns.distplot(
        per_gcs['oo_both'],
        hist=False,
        label='25C Oo. FBF-1/FBF-2 (n={0})'.format(len(per_gcs['oo_both'])))
    ax = sns.distplot(
        per_gcs['lt_fbf1'],
        hist=False,

        label='20C Oo. FBF-1 (n={0})'.format(len(per_gcs['lt_fbf1'])))
    ax = sns.distplot(
        per_gcs['lt_fbf2'],
        hist=False,

        label='20C Oo. FBF-2 (n={0})'.format(len(per_gcs['lt_fbf2'])))
    ax.set_xlabel('% GC in peak')
    ax.set_ylabel('# Peaks (density)')
    sns.despine()
    plt.legend()
    plt.savefig('figs/Fig S2 gc_richness_distplot_both_ht_fbf.pdf')
    plt.clf()
    plt.close()


    
    sys.exit()

    
    sns.set_style('ticks')
    ax = sns.distplot(
        per_gcs['oo_fbf1'],
        label='25C Oo. FBF-1 (n={0})'.format(len(per_gcs['oo_fbf1'])))
    ax = sns.distplot(
        per_gcs['lt_fbf1'],
        label='20C Oo. FBF-1 (n={0})'.format(len(per_gcs['lt_fbf1'])))
    ax.set_xlabel('% GC in peak')
    ax.set_ylabel('# Peaks (density)')
    sns.despine()
    plt.legend()
    plt.savefig('figs/gc_richness_distplot_fbf1.pdf')
    plt.clf()
    plt.close()
    sns.set_style('ticks')
    ax = sns.distplot(
        per_gcs['oo_fbf2'],
        label='25C Oo. FBF-2 (n={0})'.format(len(per_gcs['oo_fbf2'])))
    ax = sns.distplot(
        per_gcs['lt_fbf2'],
        label='20C Oo. FBF-2 (n={0})'.format(len(per_gcs['lt_fbf2'])))
    ax.set_xlabel('% GC in peak')
    ax.set_ylabel('# Peaks (density)')
    sns.despine()
    plt.legend()
    plt.savefig('figs/gc_richness_distplot_fbf2.pdf')
    plt.clf()
    plt.close()

#p.add_reads(ga['ht_both'], name='ht_both')
#p.add_reads(ga['ht_fbf1'], name='ht_fbf1')
