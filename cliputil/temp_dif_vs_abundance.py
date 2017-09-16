from __future__ import division
import scipy.stats as scs
from peaksList import *
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from volcanoMaker import *
#$^\circ$ Text

def compare_temp_changes_vs_abundance(v):
    ortiz_df = pandas.read_csv(
        '/opt/lib/ortiz/DESeq_genes_in_gonad.txt', sep='\t')
    print v.df.columns
    #tups = zip(#v.df['gene_name'].tolist(),
    v.df = v.df[v.df['baseMean']>=50]
#    v.df = v.df[v.df['padj']<0.01]
    y = v.df['log2FoldChange'].tolist()
    x = [v.name_to_oo_rpkm(name) for name in v.df.gene_name]
    z = v.df['padj'].tolist()
    #tups = np.array(zip(x, y, z))
    df = pandas.DataFrame(zip(x, y, z), columns=['x', 'y', 'padj'])
    df['x'] = [np.log10(_x) for _x in df.x]
    df = df.replace([np.inf, -np.inf], np.nan)
    df.dropna(inplace=True)
    #df = df.head(100)
    print df
    #x = v.df['Expression in  fog-2 gonads (RPKM)'].tolist()
    with sns.plotting_context('paper'):
        sns.set_style('ticks')
        mpl.rc("figure", figsize=(4, 4))
        fig, ax = plt.subplots(nrows=1)
        xy = df[['x', 'y']]
        ax = sns.kdeplot(
            xy, shade=True, shade_lowest=False, alpha=1,
            ax=ax)#, aspect=2, size=2)
        ax.set_title(
            'Temperature enrichment as a function of RNA abundance ' + \
            '(n={0}).'.format(len(x)))
        ax.set_ylabel('25' + u'\xb0' + 'C/20' + u'\xb0' +'FBF reads/gene')
        ax.set_xlabel('RNA abundance (RPKM)')
        df = df[df['padj']<0.01]
        ax = sns.kdeplot(
            df[['x', 'y']], cmap='Reds',
            shade=True, shade_lowest=False, alpha=0.5,
            ax=ax)#, aspect=2, size=2)
        plt.savefig('figs/25deg enrichment vs RNA abundance kdeplot.pdf',
                    bbox_inches='tight')
        plt.clf()
        plt.close()

def run():
    v = volcanoMaker()
    v.gl_rnaseq()
    v.read_clip_deseq_csv('tables/lt_fbf1_and_2_vs_ht_fbf_deseq.txt')
    v.read_programs()
    df = v.df[v.df['has_ortiz']].copy()
    compare_temp_changes_vs_abundance(v)
    sys.exit()
    df = df[df['padj']<0.01]
    df_ht_down = df[df['log2FoldChange']>1]
    df_ht_up = df[df['log2FoldChange']<-1]
    df_ht_down.to_csv('tables/deseq_down_at_ht.txt', sep='\t')
    df_ht_up.to_csv('tables/deseq_up_at_ht.txt', sep='\t')
    v.table_of_stats()
    v.volcano_plot(
        reverse_x=True,
        xlabel='FBF binding: high temperature/low temperature (log2)',
        output_name='Fig Sx volcano_of_temp_effect.pdf')
    v = volcanoMaker()
    v.gl_rnaseq()
    v.read_clip_deseq_csv('tables/6_reps_sp_vs_oo.txt')
    v.read_programs()
    df = v.df[v.df['has_ortiz']].copy()
    df = df[df['padj']<0.01]
    df_sp_down = df[df['log2FoldChange']<-1]
    df_sp_up = df[df['log2FoldChange']>1]
    df_sp_down.to_csv('tables/deseq_down_at_sp.txt', sep='\t')
    df_sp_up.to_csv('tables/deseq_up_at_sp.txt', sep='\t')
    v.table_of_stats()
    v.volcano_plot(xlabel='FBF spermatogenic/oogenic (log2)',
                   output_name='Fig 3A volcano.pdf')


if __name__ == '__main__':
    run()
