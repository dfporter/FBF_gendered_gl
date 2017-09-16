import cliputil
import pandas
import os
import sys
import collections
import glob
import re

df = pandas.read_csv('tables/sp_block_from_table_s4.txt', header=None)
print df
print df.columns
print df[0]

from peaksList import *
sp = peaksList()
sp.read_csv('combined_filtered/sp_both.txt')
oo = peaksList()
oo.read_csv('combined_filtered/oo_both.txt')

sp_block = set(df[0].tolist())
sp_targets = set(sp.df['gene_name'].tolist())
oo_targets = set(oo.df['gene_name'].tolist())

sp_height = sp.dict_by_first_instance(
    zip(sp.df.gene_name, sp.df['exp_reads'].tolist())
    )
sp_block_only_sp = (sp_block & sp_targets) - oo_targets
sp_y = [
    (x, sp_height[x]) for x in (sp_block & sp_targets)]
sp_block_only_sp_high = [x for x in sp_y if (
    (x[1] > 50) and (x[0] in sp_block_only_sp))]
print """SP block size {0}.\nSP targets {1}.\nOO targets {2}
SP block and SP targets {3}
SP block and SP targets - OO targets {4}
SP block and SP targets - OO targets, above 50 reads/mill {5}""".format(
    len(sp_block), len(sp_targets), len(oo_targets),
    len(sp_block & sp_targets), len(sp_block_only_sp),
    len(sp_block_only_sp_high))

sp_block_df = pandas.DataFrame(
    sp_block_only_sp_high, columns=['Gene name', 'Peak height (Sp. FBF)'])
sp_block_df.sort(
    columns='Peak height (Sp. FBF)', ascending=False, inplace=True)
sp_block_df.to_csv('tables/sp_block_genes.txt', sep='\t', index=False)
print sp_block_only_sp_high

import seaborn as sns
import matplotlib.pyplot as plt
with sns.plotting_context("paper", font_scale=0.8):
    #sns.set_style("white")
    sns.set_style("ticks")
    ax = sns.distplot([x[1] for x in sp_y], label='SP block, SP FBF targets')
    ax = sns.distplot([x[1] for x in sp_y if (x[0] in sp_block_only_sp)], label='SP block, SP FBF unique targets (not OO)')
    sns.despine()
    ax.set_xlabel('FBF peak height in Sp. GLs')
    ax.set_ylabel('# Peaks (density)')
    plt.legend()
    plt.savefig('figs/peak_height_for_sp_block_genes.pdf')
    plt.clf()
    plt.close()

