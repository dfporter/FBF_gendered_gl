
import scipy.stats as scs
from .peaksList import *
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from .volcanoMaker import *
table = """gld-1
larp-1
ppw-2
wago-4
prg-1
T07A9.14
ifet-1
cgh-1
cpb-3
ncl-1
car-1
lin-41
spat-2
sip-1
egg-6
trcs-1
cpg-1
ima-3
kin-19
gck-1
plk-3
htp-1
smk-1""".split('\n')


iib = """gld-1
lin-41
cpg-1
plk-3
cgh-1
wago-4
ifet-1
cpb-3
ppw-2
ima-3
car-1
trcs-1
egg-6
spat-2
ncl-1
T07A9.14
htp-1
larp-1
gck-1
kin-19
sip-1
prg-1
smk-1""".split('\n')
from .block_iii import blockiii
#iib = blockiii
print(iib)
def sort_drop_dups(_df):
    _df.sort(columns=['exp_reads'], ascending=False, inplace=True)
    _df.drop_duplicates(subset='gene_name', inplace=True)
    return _df

def compare_temp_changes_vs_abundance(v):
    ortiz_df = pandas.read_csv(
        '/opt/lib/ortiz/DESeq_genes_in_gonad.txt', sep='\t')
    print(v.df.columns)
#    v.df['keep'] = [(x in iib) for x in v.df['gene_name'].tolist()]
    #tups = zip(#v.df['gene_name'].tolist(),
#    v.df = v.df[v.df['baseMean']>=50]
    print("len iib {0} len found {1}. missing: {2}".format(
        len(iib), len(v.df.gene_name),
        set(iib) - set(v.df.gene_name)))
    def get_xy_df(xyz, name_set):
        print(str(xyz)[:100])
        sublist = [
            (name, x, y, z) for (name, x, y, z) in xyz if name in name_set]
        df = pandas.DataFrame(sublist, columns=['gene_name', 'x', 'y', 'padj'])
        print(df[df['gene_name']=='larp-1'])
        print('---')
        df['x'] = [np.log10(_x) for _x in df.x]
        df = df.replace([np.inf, -np.inf], np.nan)
        df.dropna(inplace=True)
        _x = np.array(df.x)
        abundances = _x[np.isfinite(_x)]
        return df, abundances
    names = v.df['gene_name'].tolist()
    print(v.df.head())
#    print names
    x = [v.name_to_oo_rpkm(name) for name in v.df.gene_name]
    y = v.df['log2FoldChange'].tolist()
    z = v.df['padj'].tolist()
    xyz = list(zip(names, x, y, z))
    name_set = set(v.df['gene_name'].tolist())
    df,  abundances_all = get_xy_df(xyz, name_set)
    ii, abundances_ii = get_xy_df(xyz, set(iib))
    print("all genes abundances median {0} block abundances median {1}".format(
        np.median([10**x for x in abundances_all]),
        np.median([10**x for x in abundances_ii])))
    print("all genes sp/oo median {0} block sp/oo median {1}".format(
        np.median([2**a for a in df['y'].tolist()]),
        np.median([2**a for a in ii['y'].tolist()])))

    oo_df = pandas.read_csv('combined_filtered/oo_both.txt', sep='\t')
    oo_df = sort_drop_dups(oo_df)
    oo_df['keep'] = [(n in set(iib)) for n in oo_df.gene_name]
    ii_peaks = oo_df[oo_df['keep']].copy()
    ii_peaks.to_csv('tables/ii_peaks.txt', sep='\t')
    print("ii not in targs: {0}".format(set(iib) - set(
        oo_df.gene_name)))
    sp_df = pandas.read_csv('combined_filtered/sp_both.txt', sep='\t')
    sp_df = sort_drop_dups(sp_df)
    #print oo_df[oo_df['gene_name']=='larp-1']
    def top_n(_df):
        _df = sort_drop_dups(_df.copy())
        return _df.head(100).copy()
    targs_set = set(oo_df.gene_name)# & set(oo_df.gene_name)
    targs_xy_df, abundances_targs = get_xy_df(xyz, targs_set)
    top_oo = top_n(oo_df)
    #top500_sp = top500(sp_df)
    top_oo_df, abundances_top_oo = get_xy_df(xyz, set(top_oo.gene_name))
    #print oo_df[oo_df['gene_name']=='larp-1']


    oo_d = oo_df.to_dict('records')
    by_name = {}
    for row in oo_d:
        by_name[row['gene_name']] = row
    for gene in set(ii.gene_name):
        if gene not in set(by_name):
            #print gene
            #print "/"
            continue
        else:
            pass
            #print row
    #x = v.df['Expression in  fog-2 gonads (RPKM)'].tolist()
    with sns.plotting_context('paper'):
        sns.set_style('ticks')
        mpl.rc("figure", figsize=(4, 4))
        fig, ax = plt.subplots(nrows=1)
        xy_all = df[['x', 'y']]
        plt.scatter(df['x'].tolist(), df['y'].tolist(), color='k', alpha=0.05)
        plt.scatter(ii['x'].tolist(), ii['y'].tolist(), color='r',
                    edgecolor=None)
        ax.set_ylim([-2, 2])
        #ax = sns.kdeplot(
        #    xy, shade=True, shade_lowest=False, alpha=1,
        #    ax=ax)#, aspect=2, size=2)
        #ax.set_title(
        #    'Temperature enrichment as a function of RNA abundance ' + \
        #    '(n={0}).'.format(len(x)))
        #ax.set_ylabel('25' + u'\xb0' + 'C/20' + u'\xb0' +'FBF reads/gene')
        #ax.set_xlabel('RNA abundance (RPKM)')
        #df = df[df['padj']<0.05]
        #ax = sns.kdeplot(
        #    df[['x', 'y']], cmap='Reds',
        #    shade=True, shade_lowest=False, alpha=0.5,
        #    ax=ax)#, aspect=2, size=2)
        plt.savefig('figs/gld1_block 25deg enrichment vs RNA abundance scatter.pdf',
                    bbox_inches='tight')
        plt.clf()
        plt.close()
        fig, ax = plt.subplots(nrows=1)
        to_pad = len(abundances_all) - len(abundances_ii)
#        np.stack()
        arr = np.vstack([abundances_all,
             np.concatenate((abundances_ii,[np.nan] * to_pad))]).T
        adf = pandas.DataFrame(
            arr,
            columns=['all', 'ii'])
        adf = adf.apply(np.log10)
        adf = adf.replace([np.inf, -np.inf], np.nan)
        #adf.dropna(inplace=True)
        tups = []
        import random
        for r in random.sample(abundances_all, 500):
            tups.append(['Oo. transcriptome', r])
        for r in abundances_targs:
            tups.append(['Oo. FBF targets', r])
        for r in abundances_ii:
            tups.append(['Block II', r])
        for r in abundances_top_oo:
            tups.append(['Top 100 Oo. FBF targets', r]) 
        adf = pandas.DataFrame(
            tups, columns=[' ', 'RNA abundance in Oo. GL (log10 RPKM)'])
        adf.dropna(inplace=True)
        #adf['val'] = [np.log10(x) for x in adf.val]
        #adf = adf.replace([np.inf, -np.inf], np.nan)
        #print adf
        sns.set_style('ticks')
        ax = sns.violinplot(
            x=' ',
            y='RNA abundance in Oo. GL (log10 RPKM)',
            data=adf,
            palette="Set3",
            #color='k',
            ax=ax)#, aspect=2, size=2)
        sns.despine()
        #ax = sns.kdeplot(
        #    abundances_ii, ax=ax, color='r')
        #ax.set_xticklabels(rotation=30)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig('figs/gld1_block RNA abundance hist.pdf',
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

    
    df = df[df['padj']<0.05]
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
    df = df[df['padj']<0.05]
    df_sp_down = df[df['log2FoldChange']<-1]
    df_sp_up = df[df['log2FoldChange']>1]
    df_sp_down.to_csv('tables/deseq_down_at_sp.txt', sep='\t')
    df_sp_up.to_csv('tables/deseq_up_at_sp.txt', sep='\t')
    v.table_of_stats()
    v.volcano_plot(xlabel='FBF spermatogenic/oogenic (log2)',
                   output_name='Fig 3A volcano.pdf')


if __name__ == '__main__':
    run()
