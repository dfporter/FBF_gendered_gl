
import argparse
import sys
import os
import pandas
import collections
import numpy as np
import matplotlib.pyplot as plt
#import cliputil
import scipy.stats as scs
import scipy
import matplotlib

import volcanoMaker
from volcanoMaker import  pltclose


def id_to_gl_deseq_val(wbid, gl_deseq):
    if wbid not in gl_deseq:
        return 0
    val = gl_deseq[wbid]['log2 fold change (of normalized reads)']
    try:
        val=float(val)
        return val
    except:
        return 0

class scatterplotMaker(volcanoMaker.volcanoMaker):
    """Make scatterplots.
    Inherits from volcanoMaker, which inherits from peaksList.
    """

    def scatterplot(
            self, clip_fname='tables/6_reps_sp_vs_oo.txt',
            outfname='figs/Fig 2C Scatterplot.pdf'):

        print("Making a scatterplot for {}.".format(clip_fname))

        # Set self.gl_deseq as dict of (name -> row of deseq info), which calls
        # peaksList.read_sp_vs_oo() to load RNA-seq programs.
        self.gl_rnaseq()
        
        # Set self.clipdf (as clip DESeq2 data).
        # volcanoMaker function. Calls volcanoMaker.gl_rnaseq() if it
        # wasn't called before. Sets a has_ortiz column in self.clipdf based on
        # the output of self.gl_rnaseq().
        self.read_clip_deseq_csv(clip_fname)
        
        # Set self.program (dict by wbid) and adds Program column to self.df,
        # and, if self.clipdf exists, to self.clipdf.
        # This is a peaksList function.
        self.read_sp_vs_oo_as_programs()
        
        # Now the figure.        
        # Only keep genes with Ortiz germline RNA-seq data.
        clip_deseq = self.clipdf[self.clipdf['has_ortiz']]

        print("Making a scatterplot using {} RNAs.".format(clip_deseq.shape[0]))

        #clip_deseq = clip_deseq[clip_deseq['baseMean']>=50]
        gl_deseq_vals = [id_to_gl_deseq_val(x, self.gl_deseq) \
                         for x in clip_deseq['gene_name'].tolist()]

        clip_deseq['gl_deseq'] = gl_deseq_vals
        cs = [self.id_to_color(x) for x in clip_deseq.gene_id]

        _y = clip_deseq['log2FoldChange'].tolist()

        pltclose()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.scatter(
            gl_deseq_vals, _y, linewidth=0,
            color=cs, alpha=0.2, s=5)
        #plt.plot([-20, 20], [-21, 19], 'r--')
        #plt.plot([-20, 20], [-19, 21], 'r--')
        plt.plot([-20, 20], [-20, 20], 'k--', alpha=0.2)
        plt.plot([-20, 20], [0, 0], 'k:', alpha=0.5)
        
        plt.ylim([-10, 15])
        plt.xlim([-10, 15])
        plt.axes().set_aspect(1./plt.axes().get_data_ratio())
        #plt.axhline(y=0, linestyle='--', c='k', alpha=0.5)
        
        plt.ylabel('FBF spermatogenic/oogenic (log2)', fontsize=10)
        plt.xlabel('RNA-seq spermatogenic/oogenic (log2)', fontsize=10)
        
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.tick_params(axis='both', which='minor', labelsize=6)
        
        fig.set_figheight(3)
        fig.set_figwidth(3)
        
        plt.savefig(outfname, format='pdf')
        pltclose()
        
        
        pear_r, _ = scs.pearsonr(gl_deseq_vals, _y)
        spear_r, _ = scs.spearmanr(gl_deseq_vals, _y)
        
        print("Correlation between CLIP and RNA abundance log2 fold change vs log2 fold change:")
        print("\tPearson R {0}\n\tSpearman rho {1}".format(
            pear_r, spear_r))

        pear_r, _ = scs.pearsonr([2**x for x in gl_deseq_vals], [2**x for x in _y])
        spear_r, _ = scs.spearmanr([2**x for x in gl_deseq_vals], [2**x for x in _y])
        print("Correlation between CLIP and RNA abundance nonlog fold change vs nonlog fold change:")
        print("\tPearson R {}\n\tSpearman rho {}".format(pear_r, spear_r))
        
        
        sig = clip_deseq[clip_deseq['padj']<0.05].copy()
        sp = sig[sig['log2FoldChange']>1].copy()
        oo = sig[sig['log2FoldChange']<-1].copy()
        
        print(sig['Program'].value_counts())
        
        def pr_vc(vc, label='value counts'):
            vc = dict(vc['Program'].value_counts())
            total = sum(vc.values())
            print(label)
            for k in vc: print('{0}: {1} ({2:.2f}%)'.format(
                k, vc[k], 100* vc[k]/total))
                
        pr_vc(sig, label='genes < 0.05 p value')
        pr_vc(sp, label='sp')
        pr_vc(oo, label='oo')
        
        return

    
        # Abundance in Oo. transcriptomes vs correlation.
        # Subset to well-bound genes.
        clip_deseq['RPKM in Oo. GL'] = [np.log10(self.name_to_oo_rpkm(x)) for \
                                        x in clip_deseq.gene_id]
        clip_deseq['RPKM in Sp. GL'] = [np.log10(self.name_to_sp_rpkm(x)) for \
                                        x in clip_deseq.gene_id]
        clip_deseq['log10_baseMean'] = [np.log10(x) for x in clip_deseq.baseMean]
        oo_deseq = clip_deseq[clip_deseq['Program']!='Spermatogenic only'].copy()
        sp_deseq = clip_deseq[clip_deseq['Program']!='Oogenic only'].copy()
        both_deseq = clip_deseq[clip_deseq['Program']=='Oogenic and Spermatogenic'].copy()
        def dist_to_diagonal(df):
            tups = list(zip(df['log2FoldChange'].tolist(), df['gl_deseq'].tolist()))
            return [(a - b) for a, b in tups]
        f, ax = plt.subplots(3, 2, figsize=(8, 12))

        # Oo.
        ax[0][0].hexbin(
            oo_deseq['RPKM in Oo. GL'].tolist(),
             #oo_deseq['RPKM in Oo. GL'].tolist(),
             dist_to_diagonal(oo_deseq), #linewidth=0,
            #color="#933b41",
            #alpha=0.2
            )
        ax[0][0].axhline(y=0, linestyle='--', c='k', alpha=0.5)
        ax[0][1].hexbin(
            oo_deseq['log10_baseMean'].tolist(),
             #oo_deseq['RPKM in Oo. GL'].tolist(),
             dist_to_diagonal(oo_deseq), #linewidth=0,
            #color="#933b41", alpha=0.2
            )
        ax[0][1].axhline(y=0, linestyle='--', c='k', alpha=0.5)

        # Sp.
        ax[1][0].hexbin(
            sp_deseq['RPKM in Sp. GL'].tolist(),
            dist_to_diagonal(sp_deseq),
#            linewidth=0,
#            color="#6600ff",
#            alpha=0.2
            )
        ax[1][0].axhline(y=0, linestyle='--', c='k', alpha=0.5)
        ax[1][1].hexbin(
            sp_deseq['log10_baseMean'].tolist(),
             #oo_deseq['RPKM in Oo. GL'].tolist(),
             dist_to_diagonal(sp_deseq),# linewidth=0,
            #color="#933b41", alpha=0.2)
            )
        ax[1][1].axhline(y=0, linestyle='--', c='k', alpha=1)

        # Both.
        t = list(zip(both_deseq['RPKM in Sp. GL'].tolist(),
                both_deseq['RPKM in Oo. GL'].tolist()))
        ax[2][0].hexbin(
            [np.mean([a, b]) for a,b in t],
            dist_to_diagonal(both_deseq),
            #linewidth=0,
            #color='#75968f',
            #alpha=0.2
            )
        ax[2][0].axhline(y=0, linestyle='--', c='k', alpha=1)
        ax[2][1].hexbin(
            both_deseq['log10_baseMean'].tolist(),
             #oo_deseq['RPKM in Oo. GL'].tolist(),
             dist_to_diagonal(both_deseq), #linewidth=0,
            #color="#933b41",
            #alpha=0.2
            )
        ax[2][1].axhline(y=0, linestyle='--', c='k', alpha=1)
        for i, _ in enumerate(ax):
            #ax[i][0].set_xlim([-4, 5])
            #ax[i][0].set_ylim([-15, 15])
            ax[i][0].set_ylabel(
                'SP/OO enrichment, log2(FBF/RNA-seq)', fontsize=8)
            ax[i][1].set_ylabel(
                'SP/OO enrichment, log2(FBF/RNA-seq)', fontsize=8)
            ax[i][0].set_xlabel(
'RNA abundance, log10(RPKM, ' + ['Oo. GL)', 'Sp. GL)', 'Oo./Sp. GL average)'][i],
fontsize=8)
            ax[i][1].set_xlabel(
'FBF binding, log10 AU ' + ['(Oo. program)', 'Sp. program)', 'Oo. and Sp. program)'][i],
 fontsize=8)
            
        plt.tight_layout()
        #Y47A7.2'
        plt.savefig('figs/Fig S2 Distance_to_45degree_by_abundance.pdf')
        pltclose()


if __name__ == '__main__':
    s = scatterplotMaker()
    s.scatterplot(clip_fname='tables/6_reps_sp_vs_oo.txt')
