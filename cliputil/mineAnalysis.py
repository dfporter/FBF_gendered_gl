
import re, os, pandas, matplotlib, collections
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import statsmodels.stats.multitest
from blocks import blocki, blockii, blockiii, blockiv

class mine(object):
    
    def __init__(self, fname='blockii.txt', col_to_use='RNAi Phenotype Observed',
                 protein_coding=None, gene_subset=None):

        self.df = pandas.read_csv(fname, sep='\t')
        self.df.drop_duplicates(subset=col_to_use, inplace=True)

        if protein_coding is not None:
            self.df = self.df[[(x in protein_coding) for x in self.df['Gene'].tolist()]]

        if gene_subset is not None:
            self.df = self.df[[(x in gene_subset) for x in self.df['Public Name'].tolist()]]

        rnai = [set(str(x).split(',')) for x in \
            self.df[col_to_use].tolist()]# if not(set(['N.A.']) & set(str(x).split(',')))]

        pheno_freq = collections.defaultdict(int)
        n_genes = len(rnai)
        pheno_freq['all apoptosis'] = 0

        for pheno_set in rnai:
            for pheno in pheno_set:
                pheno_freq[pheno] += 1

        self.n_with_term = dict(list(pheno_freq.items()))
        self.n_without_term = collections.defaultdict(int)

        for key in pheno_freq:
            self.n_without_term[key] = n_genes - pheno_freq[key]
            pheno_freq[key] = pheno_freq[key]/n_genes

        self.pheno_freq = pheno_freq

        print('Created mine object from {} RNAs.'.format(len(self.df.index)))
    
class plotTerms(object):

    def __init__(self, protein_coding, all_targs, msp, moo, m1, m2, m3, m4, fname='figs/terms.pdf',
                 include_re_cat=False, include_na=False):
        self.mrna, self.all_targs = protein_coding, all_targs
        self.m1, self.m2, self.m3, self.m4 = m1, m2, m3, m4
        freqs = {
            'Genome': protein_coding.pheno_freq,
            'SP or OO targets': all_targs.pheno_freq,
            'SP': msp.pheno_freq,
            'OO': moo.pheno_freq,
            'I': m1.pheno_freq, 
            'II': m2.pheno_freq, 
            'III': m3.pheno_freq,
            'IV': m4.pheno_freq}
        name_order = ['Genome', 'SP or OO targets', 'SP', 'OO', 'I', 'II', 'III', 'IV']
        top_terms = set()
        levels = [-1, -1, -1, -1, -1]

        cutoff = 6
        if not include_na:
            for x in list(freqs.values()):
                if 'N.A.' in x:
                    del x['N.A.']
        if not include_re_cat:
            for x in list(freqs.values()):
                if 'all apoptosis' in x:
                    del x['all apoptosis']
                if 'germ' in x:
                    del x['germ']
#        for pheno_freq in [freq1, freq2, freq3]:
#            for k in ['sterile', 'embryonic lethal', 'maternal sterile', 'N.A.']:
#                del pheno_freq[k]
        pvals = collections.defaultdict(dict)

        for name, mine_object in zip(
            name_order,
            [protein_coding, all_targs, msp, moo, m1, m2, m3]):
            
            mine_object.name = name
            for term in mine_object.pheno_freq:
                
                #if mine_object.pheno_freq[term] < 0.01: continue
                pval = self.fisher(term, protein_coding, mine_object)
                pval = -np.log10(pval)
                pvals[name][term] = pval

        while len(top_terms) < cutoff:
            for level, mine_obj in zip(
                    [0, 0, 1], [all_targs, all_targs, m3]):
                _sorted = sorted(
                    mine_obj.pheno_freq, key=lambda x: pvals[mine_obj.name][x])
                top_terms.add(_sorted[levels[level]])
                levels[level] -= 1

        #sys.exit()
        to_df = []
        to_fold_df = []
        pvals_to_corr = []
        pvals_to_corr_fold_df = []
        for term in sorted(top_terms, key=lambda x: pvals[all_targs.name][x]):

            for name, block, mine_object in zip(
                name_order, 
                [freqs[k] for k in name_order],
                [protein_coding, all_targs, msp, moo, m1, m2, m3, m4]):

                print(term, name)
                pval = self.fisher(term, protein_coding, mine_object)
                log10_pval = -np.log10(pval)
                print(log10_pval)

                to_df.append(
                    {'term': term, 'block': name, 'freq': block[term], 'log10 pval': log10_pval,
                             'pval': pval})

                pvals_to_corr.append(pval)

                if freqs['Genome'][term]>0:

                    to_fold_df.append(
                        {'term': term, 'block': name,
                        'freq': block[term]/freqs['Genome'][term],
                        'log10 pval': log10_pval})

                    pvals_to_corr_fold_df.append(pval)


        rej_array, cor_pvals, a_sidak, a_ = statsmodels.stats.multitest.multipletests(
            pvals=pvals_to_corr, method='fdr_bh')
        for n in range(len(to_df)):
            to_df[n]['log10 pval'] = -np.log10(cor_pvals[n])

        rej_array, cor_pvals_fold_df, a_sidak, a_ = statsmodels.stats.multitest.multipletests(
            pvals=pvals_to_corr_fold_df, method='fdr_bh')

        for n in range(len(to_fold_df)):
                to_fold_df[n]['log10 pval'] = cor_pvals_fold_df[n]# -np.log10(cor_pvals_fold_df[n])

        print('+=++=\nHypotheses tested {0}'.format(len(to_df)))
        df = pandas.DataFrame(to_df)
        fold_df = pandas.DataFrame(to_fold_df)

        sns.set_context('paper')
        sns.set(style='whitegrid')
        matplotlib.rcParams.update({'font.size': 6})

        fig = plt.figure()

        ax = sns.barplot(y="block", x="log10 pval", hue="term", data=df,
                  palette=['#82239B', '#455CFA', '#BCF8FF', '#C3D8A3', '#E3BE74', '#5E5E5E'])

        ax.tick_params(labelsize=6)
        sns.despine()

        plt.rc('legend', fontsize=6)
        leg = plt.legend()

        fig.set_figwidth(2)
        fig.set_figheight(2)

        ax.set_xlabel('-log10(p value)')

        if not os.path.exists(os.path.dirname(fname)):
            os.system('mkdir ' + os.path.dirname(fname))

        plt.savefig(fname)
        plt.clf(); plt.close()
        
        fig = plt.figure()

        ax = sns.barplot(y="block", x="freq", hue="term", data=fold_df,
                palette=['#82239B', '#455CFA', '#BCF8FF', '#C3D8A3', '#E3BE74', '#5E5E5E'])
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(6)
        plt.rc('legend', fontsize=6)
        leg = plt.legend()
#        print leg.__dict__
#        print ax.__dict__['figure'].get_size_inches()
        #ax.set_xlim([0, 10])
        print(ax.__dict__['figure'].get_size_inches())
        #ax.__dict__['figure'].c(3)
        #ax.__dict__['figure'].set_figheight(3)
        fig.set_figwidth(3)
        fig.set_figheight(3)
        print(ax.__dict__['figure'].get_size_inches())
        plt.savefig(re.sub('\.pdf', '_fold.pdf', fname))
        plt.clf(); plt.close()
        
    def fisher(self, term, m1, m2):
        # m2 is assumed to be a subset of m1
        m1_len = len(set(m1.df['Gene'].tolist()))
        m2_len = len(set(m2.df['Gene'].tolist()))

        if term not in m2.n_with_term:
            return 1.

#        if term in ['germ', 'antibody staining reduced', 'transgene induced cosuppression variant']:
#            return 1.

        # 2x2 table: Has term or not vs in m2 or not
        table =[ # First row:In m2 
            [m2.n_with_term[term],  # Has term
             m2.n_without_term[term]],  # Does not have term
             # Second row: not in m2 
            [(m1.n_with_term[term]) - (m2.n_with_term[term]), 
             m1.n_without_term[term] - m2.n_without_term[term]]
        ]
        for x in table:
            for y in x:
                if y<0:
                    print(m1.name)
                    print(m2.name)
                    print(term, table)
        #print table
        odds, pval = sp.stats.fisher_exact(table, alternative='greater')
        return pval


if __name__ == '__main__':

    phenotypes_file = os.path.dirname(__file__) + '/phenotypes_all_protein_coding_genes_ws235.txt'
    if not os.path.exists(phenotypes_file):
        raise IOError("Looked for a file of phenotypes for each protein coding gene at {}, but file does not exist.".format(
            phenotypes_file))

    protein_coding = set(
        pandas.read_csv(phenotypes_file, 
                sep='\t', index_col=False)['Gene'].tolist())

    for col in ['RNAi Phenotype Observed', 'Allele Phenotype Observed',
                'Expr_pattern Tissue', 'Genomic Study Tissue',
                'Expr_pattern LifeStage',
                'Genomic Study LifeStage'
                    ]:

        if col != 'RNAi Phenotype Observed':
            continue

        mi = mine(
            fname=phenotypes_file, gene_subset=blocki,
            col_to_use=col, protein_coding=protein_coding)
        mii = mine(
            fname=phenotypes_file, gene_subset=blockii,
            col_to_use=col, protein_coding=protein_coding)
        miii = mine(
            fname=phenotypes_file, gene_subset=blockiii,
            col_to_use=col, protein_coding=protein_coding)
        miv = mine(
            fname=phenotypes_file, gene_subset=blockiv,
            col_to_use=col, protein_coding=protein_coding)

        all_genome = mine(fname=phenotypes_file, col_to_use=col)

        sp_genes =  set(
                pandas.read_csv('combined_filtered/sp_both.txt', 
                sep='\t', index_col=False)['gene_name'].tolist())   
        
        oo_genes =  set(
                pandas.read_csv('combined_filtered/oo_both.txt', 
                sep='\t', index_col=False)['gene_name'].tolist())

        only_sp = sp_genes - oo_genes
        only_oo = oo_genes - sp_genes
        both = sp_genes & oo_genes
        
        all_spoo = mine(fname=phenotypes_file,
            gene_subset=both,
            col_to_use=col, protein_coding=protein_coding)
        
        msp = mine(
            fname=phenotypes_file,
            gene_subset=sp_genes, col_to_use=col, 
            protein_coding=protein_coding)
        
        moo = mine(
            fname=phenotypes_file,
            gene_subset=oo_genes, col_to_use=col, 
            protein_coding=protein_coding)
        
        if col == 'RNAi Phenotype Observed':
            include_re_cat = True
        else:
            include_re_cat = True

        pt = plotTerms(all_genome, all_spoo, msp, moo, mi, mii, miii, miv,
                       fname='figs/{0}.pdf'.format(col),
                include_re_cat=include_re_cat)
