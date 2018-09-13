#from bokeh.models import HoverTool
#from bokeh.plotting import ColumnDataSource, figure, show, output_file
#import bokeh

import os, sys, glob, pandas, collections, random, matplotlib
import numpy as np
import matplotlib.pyplot as plt
import peaksList

import scipy.stats as scs

from blocks import blocki, blockii, blockiii
#from cliputil import *


def id_to_gl_deseq_val(wbid, gl_deseq):
    if wbid not in gl_deseq:
        return 0
    val = gl_deseq[wbid]['log2 fold change (of normalized reads)']
    try:
        val=float(val)
        return val
    except:
        return 0


class volcanoMaker(peaksList.peaksList):

    def __init__(self, dataframe=pandas.DataFrame(), name="Unnamed",
                 gene_name_col='gene_name'):
        self.name = name
        self.df = dataframe
        self.gene_name_col = gene_name_col
        self.translator()
        df = pandas.read_csv(
            '/opt/lib/name_wbid_map', sep='\t', index_col=False)
        self.name_to_wbid = dict(list(zip(
            df.gene_name, df.gene_id)))
        self.wbid_to_name = dict(list(zip(
            df.gene_id, df.gene_name)))
    
    def read_counts_files(self, dirname='counts_6_reps/'):
        
        self.counts = {}
        self.combined_counts = collections.defaultdict(int)

        for fname in glob.glob(dirname + '/*'):
            name = os.path.basename(fname)
            df = pandas.read_csv(fname, sep='\t', header=None, index_col=False)
            
            self.counts[name] = dict(list(zip(df[0].tolist(), df[1].tolist())))
            
            for gene, val in list(self.counts[name].items()):
                self.combined_counts[gene] += val
                
        for gene in self.combined_counts:
            self.combined_counts[gene] = self.combined_counts[gene]/len(self.counts)

    def filter_read_counts(self, lower_cutoff_average=10):
        
        if not(hasattr(self, 'combined_counts')):
            self.read_counts_files()  # Filter by raw read count.
            
        in_len = len(self.clipdf.index)
        
        above = [(self.combined_counts.get(x, 0)>=lower_cutoff_average) 
            for x in self.clipdf['gene_name'].tolist()]
        
        def one_if_true(_):
            if _: return 1
            return 0
        
        print(("{0}/{1} genes in clipdf above {2} average raw counts.".format(
            sum([one_if_true(n) for n in above]), len(self.clipdf.index), 
            lower_cutoff_average)))
        
        self.clipdf = self.clipdf[above].copy()
        print(("Input clipdf len: {0} output {1}".format(
            in_len, len(self.clipdf.index))))
    
    def translator(self):
        if os.path.exists('/opt/lib/name_wbid_map'): return True
        with open('/opt/lib/name_wbid_map', 'w') as f:
            f.write('gene_name\tgene_id\n')
            f.write("\n".join(
                ["{}\t{}".format(x, y) for x, y in list(self.name_to_wbid.items())]
                ))

    def gl_rnaseq(self, fname='/opt/lib/ortiz/DESeq_genes_in_gonad.txt'):
        """ Define the df and some useful dicts of GL RNA-seq."""
        # This function sets self.ortiz as the DataFrame of fname,
        # and defines variety of dicts for (gene name -> rpkm), 
        # (wb id -> rpkm), and (gene name -> log2 fold change). 
        self.read_sp_vs_oo()
        
        self.gl_deseq_df = self.ortiz  # Make an AKA for the same thing.
        
        self.gl_deseq = {}
        for row in self.ortiz.to_dict('records'):
            self.gl_deseq[row['Gene name']] = row
            
        # Should really get these next two dicts from a gtf.
        self.wbid240_to_name = dict(list(zip(self.ortiz['WormBase ID (WS240)'],
                                self.ortiz['Gene name'])))
        self.name_to_wbid240 = dict(list(zip(list(self.wbid_to_name.values()),
                                     list(self.wbid_to_name.keys()))))

    def read_programs(self):
        self.read_sp_vs_oo_as_programs()
        # The above sets self.program as (wb id -> program) dict.

    def read_clip_deseq_csv(self, fname, header_ok=False):
        
        if not header_ok:
            self.fix_header(fname)
            
        if not hasattr(self, 'gl_deseq'):
            self.gl_rnaseq()
            
        self.clipdf = pandas.read_csv(fname, sep='\t', index_col=False)
        
        self.clipdf = self.clipdf[self.clipdf['gene_name'] != '_no_feature']
        self.clipdf = self.clipdf[self.clipdf['gene_name'] != '_ambiguous']
        
        self.clipdf['gene_id'] = [self.name_to_wbid[x] for x in self.clipdf.gene_name]
        self.clipdf['Wormbase ID'] = self.clipdf['gene_id']
        self.clipdf['has_ortiz'] = [id_to_gl_deseq_val(x, self.gl_deseq)!=0 for \
                               x in self.clipdf['gene_name']]
        
    def read_clip_deseq_excel(self, fname):

        if not hasattr(self, 'gl_deseq'):
            self.gl_rnaseq()
            
        self.clipdf = pandas.read_excel(fname)
        self.clipdf['gene_name'] = self.clipdf.index
        
        self.clipdf = self.clipdf[self.clipdf['gene_name'] != '_no_feature']
        self.clipdf = self.clipdf[self.clipdf['gene_name'] != '_ambiguous']
        
        self.clipdf['gene_id'] = [self.name_to_wbid[x] for x in self.clipdf.gene_name]
        self.clipdf['Wormbase ID'] = self.clipdf['gene_id']
        self.clipdf['has_ortiz'] = [id_to_gl_deseq_val(x, self.gl_deseq)!=0 for \
                               x in self.clipdf['gene_name']]

    def fix_header(self, fname):
        l1r = ''
        l2r = ''
        outli = ''
        with open(fname) as f:
            l1r = next(f)
            l2r = next(f)
            for li in f:
                outli += li
        if len(l1r.split('\t')) == (len(l2r.split('\t')) -1):
            l1r = 'gene_name\t' + l1r
            with open(fname, 'w') as f:
                f.write(l1r)
                f.write(l2r)
                f.write(outli)
            print(("Fixed header to %s" % l1r))

    def table_of_stats(self, outf='tables/volcano_stats.txt'):
        
        positive_control_genes = set(
            "gld-1,htp-1,htp-2,mpk-1,him-3,fbf-1,lip-1,syp-2,fbf-2,fog-1,fem-3,syp-3,gld-3,fog-3,egl-4".split(','))
        
        print('-------table_of_stats() self.clipdf.index = ', len(self.clipdf.index))
        clip_deseq = self.clipdf[self.clipdf['has_ortiz']].copy()
        print('-------table_of_stats() has_ortiz self.clipdf.index = ', len(clip_deseq.index))
        
        sig_dif = clip_deseq[clip_deseq['padj']<0.01].copy()
        sig_dif['sig_up'] = [(x>=1) for x in sig_dif.log2FoldChange]
        sig_dif['sig_down'] = [(x<=-1) for x in sig_dif.log2FoldChange]
        sig_dif['sig_2_fold'] = [
            ((x>=1) or (x<=-1)) for x in sig_dif.log2FoldChange]
        
        n_gl = float(len(clip_deseq.index))
        
        def as_p(x, y):
            if type(x) == type(pandas.DataFrame()):
                x = len(x.index)
            if y == 0: return 0
            if int(float(100*x)/float(y)) == 0:
                return "{} (0.{}%)".format(x, int(float(1000*x)/float(y)))
            return "{} ({}%)".format(x, int(float(100*x)/float(y)))
        
        sp_enriched_sig = sig_dif[sig_dif['sig_up']].copy()
        oo_enriched_sig = sig_dif[sig_dif['sig_down']].copy()
        
        st = {  # Stats table.
'Total genes': len(self.clipdf.index),
'With data in RNA-seq of germlines (Ortiz, et. al)': len(clip_deseq.index),
'Oogenic': as_p(clip_deseq[clip_deseq['Program']=='Oogenic only'], n_gl),
'Spermatogenic': as_p(clip_deseq[clip_deseq['Program']=='Spermatogenic only'], n_gl),
'Oogenic and Spermatogenic': as_p(
    clip_deseq[clip_deseq['Program']=='Oogenic and Spermatogenic'], n_gl),
'Significantly different between SP and OO CLIP (P<0.01)': \
as_p(sig_dif, n_gl),
'Out of {0} genes, Significantly different between SP and OO CLIP (P<0.01) and 2-fold different'.format(n_gl):\
as_p(sig_dif['sig_2_fold'].sum(), n_gl),
'Out of {0} genes, significantly different between SP and OO CLIP (P<0.01) and 2-fold SP enriched'.format(n_gl):\
as_p(sig_dif['sig_up'].sum(), n_gl),
'Significantly different between SP and OO CLIP (P<0.01) and 2-fold OO enriched':\
as_p(sig_dif['sig_down'].sum(), n_gl),
'Out of {0} genes, Significantly 2-fold enriched in SP and in only spermatogenic program'.format(len(sp_enriched_sig.index)):\
as_p(sp_enriched_sig[sp_enriched_sig['Program']=='Spermatogenic only'], len(sp_enriched_sig.index)),
'Out of {0} genes, Significantly 2-fold enriched in SP and in only oogenic program'.format(len(sp_enriched_sig.index)):\
as_p(sp_enriched_sig[sp_enriched_sig['Program']=='Oogenic only'], len(sp_enriched_sig.index)),
'Out of {0} genes, Significantly 2-fold enriched in OO and in only spermatogenic program'.format(len(oo_enriched_sig.index)):\
as_p(oo_enriched_sig[oo_enriched_sig['Program']=='Spermatogenic only'], len(oo_enriched_sig.index)),
'Out of {0} genes, Significantly 2-fold enriched in OO and in only oogenic program'.format(len(oo_enriched_sig.index)):\
as_p(oo_enriched_sig[oo_enriched_sig['Program']=='Oogenic only'], len(oo_enriched_sig.index)),

}
        sig_records = sig_dif.to_dict('records')
        gl_records = clip_deseq.to_dict('records')
        sig_ds = {}
        for r in sig_records: sig_ds[r['gene_name']] = r
        gl_ds = {}
        for r in gl_records: gl_ds[r['gene_name']] = r
        positives = collections.defaultdict(dict)
        clip_deseq_genes = set(self.clipdf['gene_name'].tolist())
        ortiz_and_clip_genes = set(clip_deseq['gene_name'].tolist())
        for k in positive_control_genes:
            # Defaults
            positives[k]["In CLIP DESeq2"] = False
            positives[k]['In GL RNA-seq'] = False
            positives[k]['Sig dif in RNA-seq'] = False
            positives[k]['padj'] = 0.
            positives[k]['log2FoldChange'] = 0. 
            positives[k]['Program'] = ""
            if k in clip_deseq_genes:
                positives[k]['In CLIP DESeq2'] = True
                if k in ortiz_and_clip_genes:
                    positives[k]['In GL RNA-seq'] = True
                    positives[k]['padj'] = float(gl_ds[k]['padj'])
                    positives[k]['log2FoldChange'] = float(gl_ds[k]['log2FoldChange'])
                    positives[k]['Program'] = gl_ds[k]['Program']
                    if k in sig_ds:
                        positives[k]['Sig dif in RNA-seq'] = True
        positives[k]['Dataset'] = ''
        nest = ['In CLIP DESeq2', 'In GL RNA-seq']
        for v in list(positives.values()):
            for dataset in nest:
                if v[dataset]:
                    v['Dataset'] = dataset
        st['Positive controls'] = "\n\tName\tP<0.01?\tlog2FC\tP\tDataset\tProgram\n"
        for k in positives:
            st['Positive controls'] += "\t{name}\t{sig}\t{lg:.2}\t{pv:.1e}\t{d}\t{pg}\n".format(
                name=k, d=positives[k]['Dataset'], sig=positives[k]['Sig dif in RNA-seq'],
                lg=positives[k]['log2FoldChange'], pv=positives[k]['padj'],
                pg=positives[k]['Program'])
        print('--Begin stats table--')
        for k in st:
            print(("{}\t{}".format(k, st[k])))
        print('--End of stats table--')

    def id_to_color(self, x):
        if x not in self.program: return '#000000' #'#75968f'
        prog = self.program[x]
        if prog not in set(['Oogenic only', 'Spermatogenic only',
                           'Oogenic and Spermatogenic']):
            return '#000000'# 'gray'#'#75968f'
        else:
            return {'Oogenic only': "#F0326E",
                      'Spermatogenic only': "#0087B9", # "#6600ff",
                    'Spermatogenic onl': "#0087B9",
                      'Oogenic and Spermatogenic': '#000000'}[prog]

    @staticmethod
    def color_if_sig_dif(x, y, x_cutoff=1, y_cutoff=0.01):
        if abs(x) >= x_cutoff and y<=y_cutoff:
            return 'r'
        return 'k'

    def volcano_plot(
            self, xlabel='FBF spermatogenic/oogenic (log2)',
            output_name='figs/Fig 2D volcano.pdf',
            reverse_x=False, ylim=(0, 20), xlim=(-9, 9),
            temp_cf_settings=False):
        
        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce",
              "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
        
        clip_deseq = self.clipdf[self.clipdf['has_ortiz']]
        print('-------volcano_plot() self.clipdf.index = ', len(self.clipdf.index))
        
        cs = [self.id_to_color(x) for x in clip_deseq['gene_id'].tolist()]
        if temp_cf_settings:
            cs = [self.color_if_sig_dif(x, y) for x, y in zip(
                clip_deseq['log2FoldChange'].tolist(), clip_deseq['padj'].tolist())]

        x_vals = clip_deseq['log2FoldChange'].tolist()
        if reverse_x:
            x_vals = [-x for x in x_vals]
            
        y_vals = [-np.log10(x) for x in clip_deseq['padj'].tolist()]
        
        r_val, p_val = scs.pearsonr([2**w for w in x_vals], [2**w for w in y_vals])
        print("Pearson R for FoldChange vs FoldChange: {0} (R^2 {1} P value {2})".format(r_val, r_val**2, p_val))
        r_val, p_val = scs.pearsonr(x_vals, y_vals)
        print("Pearson R for log2FoldChange vs log2FoldChange: {0} (R^2 {1} P value {2})".format(r_val, r_val**2, p_val))
        
        pltclose()
        plt.rcParams['lines.markeredgewidth'] = 0.0
        fig = plt.figure()
        
        s = plt.scatter(
            x_vals, y_vals, color=cs, alpha=0.4, marker='o', s=5)
        
        s.set_lw(0)
        to_label = ['sygl-1', 'lst-1', 'fog-1', 'fog-3', 'fem-3', 'mpk-1',
                    ]#'linc-7', 'linc-29', 'linc-4']
        
        labels = clip_deseq['gene_name'].tolist()

        for label, x, y in zip(labels, x_vals, y_vals):
            if not(temp_cf_settings) and (label in to_label):

                plt.annotate(
                    label, xy=(x, y), xytext=(10, 5),
                    textcoords= 'offset points',
                    fontsize=6,
                    arrowprops=dict(arrowstyle='->'))
                
        plt.ylim(*ylim)
        plt.xlim(*xlim)
        
        plt.axhline(y=-np.log10(0.01), c='k', linestyle='--', alpha=0.5, lw=1)
        
        plt.xlabel(xlabel)
        plt.ylabel('-log10(p value) for differential binding')
        
        fig.set_figwidth(3)
        fig.set_figheight(3)
        plt.savefig(output_name)
        pltclose()

    def volcano_plot_blocks(
        self, xlabel='FBF spermatogenic/oogenic (log2)',
        output_name='figs/Fig 3A volcano blocks.pdf',
        reverse_x=False,
        ylim=(0, 50), xlim=(-7, 7)):
        
        
        print(('(((())))\n' * 7))

        all_blocks = blocki | blockii | blockiii
        
        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce",
              "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
                   
        df = self.clipdf[self.clipdf['has_ortiz']].copy()

        block_i_df = df[[(name in blocki) for name in df.gene_name]].copy()
        block_ii_df = df[[(name in blockii) for name in df.gene_name]].copy()
        block_iii_df = df[[(name in blockiii) for name in df.gene_name]].copy()
        
        not_in_block_df = df[[(name not in all_blocks) for name in df.gene_name]].copy()
        pltclose()

        plt.rcParams['lines.markeredgewidth'] = 0.0

        fig = plt.figure()#figsize=(3, 3))
        ax = fig.add_subplot(111)
        
        if reverse_x:
            x_vals = [-x for x in x_vals]

        plt.rc({'markeredgewidth':0.0})
        matplotlib.rcParams['lines.markeredgewidth'] = 0.0
        plt.rcParams['lines.markeredgewidth'] = 0.0
        
        points = list(zip(not_in_block_df['log2FoldChange'].tolist(),
                         [-np.log10(x) for x in not_in_block_df['padj'].tolist()]))
        
        points = np.array(points, dtype=float)
        points = points[~np.isnan(points).any(axis=1)]
        
        #background = random.sample(background, int(len(background)/10))
        
        outside = []
        mid = []
        
        sb = ax.scatter(
            [x[0] for x in points],
            [x[1] for x in points],
            c='grey',#[to_grey(x[0]) for x in background],
            marker='o', alpha=0.2, s=5)
        
        sb.set_lw(0)
        
        def add_plot(_df, color, size, alpha):
            return ax.scatter(
                _df['log2FoldChange'].tolist(),
                [-np.log10(x) for x in _df['padj'].tolist()],
                color=color, alpha=alpha,
                    marker='o', s=size)
        
        
        s1 = add_plot(block_i_df, '#0087B9', 20, 0.5)
        s2 = add_plot(block_ii_df, 'green', 20, 0.25)
        s3 = add_plot(block_iii_df, 'purple', 20, 0.5)
        
        s1.set_lw(0); s2.set_lw(0); s3.set_lw(0);
        
        to_label = ['sygl-1', 'lst-1', 'fog-1',  'fem-3', 'mpk-1',
                    ]#'linc-7', 'linc-29', 'linc-4']
        #def _label(x): return x if x in to_label else ''
        labels = df['gene_name'].tolist()
        #for label, x, y in zip(labels, x_vals, y_vals):
        #    if label in to_label:
        #        ax.annotate(
        #            label, xy=(x, y), xytext=(10, 5),
        #            textcoords= 'offset points',
        #            fontsize=8,
        #            arrowprops=dict(arrowstyle='->'))

        plt.ylim(*ylim)
        plt.xlim(*xlim)

        plt.xlabel(xlabel, fontsize=10)
        plt.ylabel('-log10(p value) for differential binding', fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.tick_params(axis='both', which='minor', labelsize=6)

        fig.set_figwidth(3)
        fig.set_figheight(3)

        print((plt.rcParams["figure.figsize"]))
        print('---')
        plt.savefig(output_name)
        pltclose()

        
def pltclose():
    plt.clf()
    plt.close()
#cf = pandas.read_csv(, sep='\t')

def write_excel_of_deseq(df_in, header='SPvOO'):
    
    df = df_in.copy()
    print('-------write_excel_of_deseq() clipdf.index = ', len(df.index))
    df = df[df['has_ortiz']].copy()
    print('-------write_excel_of_deseq() has_ortiz clipdf.index = ', len(df.index))
    df = df[['gene_name','Wormbase ID', 'log2FoldChange', 
                'pvalue', 'padj', 'Program']]
    

    
    def by_pval(_df, pval):
        df_s = _df[_df['BH-corrected p value']<pval].copy()
        if 'log2 fold enrichment SP/OO (FBF iCLIP)' in df_s.columns:
            log2col = 'log2 fold enrichment SP/OO (FBF iCLIP)'
        if 'log2 fold enrichment 25/20 degrees (FBF iCLIP)' in df_s.columns:
            log2col = 'log2 fold enrichment 25/20 degrees (FBF iCLIP)'
        df_sp = df_s[df_s[log2col]>1]
        df_oo = df_s[df_s[log2col]<-1]
        return df_sp, df_oo
    
    if header == 'SPvOO':
        df.columns = [
            'Gene name', 'Wormbase ID',
            'log2 fold enrichment SP/OO (FBF iCLIP)',
            'Unadjusted p value', 'BH-corrected p value',
            'Program (Noble et al. 2016)']
        
        writer = pandas.ExcelWriter(
            'tables/Table S4 Gender differences DESeq2.xlsx')
        
        df_sp, df_oo = by_pval(df, 0.01)
        
        df_sp.to_excel(
            writer, sheet_name='SP enriched P<0.01', index=False)
        df_oo.to_excel(
            writer, sheet_name='OO enriched P<0.01', index=False)
        
        df_sp, df_oo = by_pval(df, 0.05)
        
        df_sp.to_excel(
            writer, sheet_name='SP enriched P<0.05', index=False)
        df_oo.to_excel(
            writer, sheet_name='OO enriched P<0.05', index=False)
        df.to_excel(
            writer, sheet_name='All DESeq2 results', index=False)
        
        writer.save()
        
    elif header == 'HTvLT':
        df.columns = [
            'Gene name', 'Wormbase ID',
            'log2 fold enrichment 25/20 degrees (FBF iCLIP)',
            'Unadjusted p value', 'BH-corrected p value',
            'Program (Noble et al. 2016)']
        writer = pandas.ExcelWriter(
            'tables/Table Sx Temperature differences DESeq2.xlsx')
        df_ht, df_lt = by_pval(df, 0.01)
        df_ht.to_excel(
            writer, sheet_name='25 enriched P<0.01', index=False)
        df_lt.to_excel(
            writer, sheet_name='20 enriched P<0.01', index=False)
        df_ht, df_lt = by_pval(df, 0.05)
        df_ht.to_excel(
            writer, sheet_name='25 enriched P<0.05', index=False)
        df_lt.to_excel(
            writer, sheet_name='20 enriched P<0.05', index=False)
        df.to_excel(
            writer, sheet_name='All DESeq2 results', index=False)
        writer.save()
    else:
        print(("Failed to parse this header argument: {0}".format(header)))


def run():
    # Volcano plot for LT vs HT (20C vs 25C, oogenic).
    v = volcanoMaker()
    v.gl_rnaseq()
    v.read_clip_deseq_csv('tables/lt_fbf1_and_2_vs_ht_fbf_deseq.txt')
    v.filter_read_counts(lower_cutoff_average=20)
    v.read_sp_vs_oo_as_programs()
    write_excel_of_deseq(v.clipdf.copy(), header='HTvLT')
    v.table_of_stats()
    v.volcano_plot(
        reverse_x=True, temp_cf_settings=True, xlim=(-3, 3),
        xlabel='FBF binding: high temperature/low temperature (log2)',
        output_name='figs/Fig Sx volcano_of_temp_effect.pdf')

    # Volcano plot for OO vs SP (25C).
    print("\n\n\n***SP vs OO***\n\n\n")
    v = volcanoMaker()
    v.gl_rnaseq()
    v.read_clip_deseq_csv('tables/6_reps_sp_vs_oo.txt')
    v.filter_read_counts(lower_cutoff_average=20)
    v.read_sp_vs_oo_as_programs()
    write_excel_of_deseq(v.clipdf.copy(), header='SPvOO')
    v.table_of_stats()
    v.volcano_plot(
        xlabel='FBF spermatogenic/oogenic (log2)',
        output_name='figs/Fig 2D volcano.pdf')
    v.volcano_plot_blocks(
        xlabel='FBF spermatogenic/oogenic (log2)',
        output_name='figs/Fig 3B volcano blocks.pdf')


if __name__ == '__main__':
    in_file = '/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/tables/File S6 Blocks.xlsx'

    #import programs_in_blocks 
    run()
