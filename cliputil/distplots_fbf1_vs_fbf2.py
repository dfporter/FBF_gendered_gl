import pandas
import scipy as sp
import numpy as np
import scipy.stats as scs
import seaborn as sns
import matplotlib.pyplot as plt

def deseq_volcano():
    plt.clf(); plt.close()
    clip_deseq = pandas.read_csv('tables/6_reps_lt_fbf1_vs_fbf2.txt', sep='\t')
    x_vals = [x for x in clip_deseq['log2FoldChange'].tolist()]
    y_vals = [(1-x) for x in clip_deseq['pvalue'].tolist()]
    #fig, ax = plt.subplots()
    #plt.rc({'markeredgewidth':0.0})
    #matplotlib.rcParams['lines.markeredgewidth'] = 0.0
    plt.rcParams['lines.markeredgewidth'] = 0.0
    sns.set_style("ticks")
    s = plt.scatter(x_vals, y_vals, color='k',#color=cs,
                    alpha=0.2,
                marker='o',
                s=20)
    plt.axhline(y=(0.95), c='k', linestyle='--', alpha=0.5, lw=1)
    sns.despine()
    s.set_lw(0)
    plt.savefig('figs/lt_fbf1_vs_fbf2_deseq_volcano.pdf')
    plt.clf(); plt.close()

deseq_volcano()
    
df1 = pandas.read_csv(
    'combined_filtered/old_fbf1_to_fbf2_n2.txt', sep='\t', index_col=False)
df2 = pandas.read_csv(
    'combined_filtered/old_fbf2.txt', sep='\t', index_col=False)
#df1 = df1[df1['exp_reads']>50]
#df2 = df2[df2['exp_reads']>50]
df1_unique = set(df1['gene_name'].tolist()) - set(df2['gene_name'].tolist())
df2_unique = set(df2['gene_name'].tolist()) - set(df1['gene_name'].tolist())
both = set(df1['gene_name'].tolist()) & set(df2['gene_name'].tolist())

def which_cat(x):
    if x in df1_unique: return 'FBF-1'
    if x in df2_unique: return 'FBF-2'
    if x in both: return 'Both'
    return ''

df1['cat'] = [which_cat(x) for x in df1['gene_name'].tolist()]
df2['cat'] = [which_cat(x) for x in df2['gene_name'].tolist()]
df1_unique_df = df1[df1['cat']=='FBF-1']
df1_unique_df.sort(columns='exp_reads', ascending=False, inplace=True)
df2_unique_df = df2[df2['cat']=='FBF-1']
both_df = df1[df1['cat']=='Both']

df1x = df1_unique_df['exp_reads'].tolist()
df2x = df2_unique_df['exp_reads'].tolist()
both_dfx = both_df['exp_reads'].tolist()
with sns.plotting_context("paper", font_scale=0.8):
    #sns.set_style("white")
    sns.set_style("ticks")
#    print df1x
#    print both_dfx
    #f, ax = plt.subplots(2)
    ax = sns.distplot(np.log10(df1x), label='FBF-1 unique targets')
    #sns.distplot(df2x)
    ax = sns.distplot(np.log10(both_dfx), label='FBF-1 and FBF-2 targets')
    ax.set_xlabel('Peak height')
    ax.set_ylabel('# Peaks (density)')
    sns.despine()#offset=1, ax=ax, trim=True)
    plt.legend()
    plt.savefig('figs/dist_plot_heights_lt_fbf.pdf')
    print scs.mannwhitneyu(df1x, both_dfx)
    print df1_unique_df.loc[154]
    print len(df1_unique_df.index)
    print len(df2_unique_df.index)
    print len(both_df.index)
    plt.clf()
    plt.close()


