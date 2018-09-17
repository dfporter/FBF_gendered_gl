import pandas, re, os, collections, sys, glob, matplotlib, matplotlib_venn, importlib
import numpy as np

import matplotlib.pyplot as plt
import scipy.stats as scs
import matplotlib.lines as mlines
from matplotlib_venn import venn2, venn3, venn3_circles

import utils
import peaksList
import figureMaker
from peaksList import *

importlib.reload(figureMaker)

class vennPainter():
    """All functions return a venn object and do nothing else.
    """
    def old_fbf1_vs_fbf2(self, ax):
        a = self.targs['old_fbf1']
        b = self.targs['old_fbf2']
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{1}$',
                '$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{2}$')
        v_plt = venn2([a, b], ax=ax, set_labels=_lab)
        v_plt.get_patch_by_id('10').set_color("green")
        v_plt.get_patch_by_id('10')._linewidth = 0
        v_plt.get_patch_by_id('01').set_color('#933b41')
        v_plt.get_patch_by_id('01')._linewidth = 0
        return v_plt

    def new_oo_vs_old_fbf(self, ax):
        
        if 'OO FBF (25°C)' in self.targs:
            a = self.targs['OO FBF-1 (20°C)']
            b = self.targs['OO FBF (25°C)']
            c = self.targs['OO FBF-2 (20°C)']
        else:
            a = self.targs['old_fbf1']
            b = self.targs['oo_both']
            c = self.targs['old_fbf2']
            
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{1}$',
                '$\mathrm{^{25^\circ C}Oo. FBF}$',
                '$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{2}$')
        v_plt = venn3(subsets=(a, b, c), ax=ax, set_labels=_lab)
        v_plt.get_patch_by_id('100').set_color("g")#("#933b41")
        v_plt.get_patch_by_id('100')._linewidth = 0
        v_plt.get_patch_by_id('010').set_color("r")
        v_plt.get_patch_by_id('010')._linewidth = 0
        v_plt.get_patch_by_id('001').set_color('#933b41')
        v_plt.get_patch_by_id('001')._linewidth = 0
        cir = venn3_circles(subsets=(a, b, c), ax=ax, linestyle='solid')
        cir[0].set_lw(0)
        cir[1].set_lw(0)
        cir[2].set_lw(0)
        return v_plt

    def new_oo_vs_old_fbf1(self, ax):
        a = self.targs['old_fbf1']
        b = self.targs['oo_both']
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{1}$',
                '$\mathrm{^{25^\circ C}Oo. FBF}$')
        v_plt = venn2(subsets=(a, b), ax=ax, set_labels=_lab)
        return self.set_venn2_red(v_plt)

    def set_venn2_red(self, v_plt):
        v_plt.get_patch_by_id('10').set_color("#F4848F")#("#933b41")
        v_plt.get_patch_by_id('10')._linewidth = 1
        v_plt.get_patch_by_id('10')._edgecolor = (0.9137254901960784, 0.11372549019607843, 0.18823529411764706, 0.40000000000000002)
        v_plt.get_patch_by_id('01').set_color("#EA1A43")
        v_plt.get_patch_by_id('01')._linewidth = 1
        v_plt.get_patch_by_id('01')._edgecolor = (0.9137254901960784, 0.11372549019607843, 0.18823529411764706, 0.40000000000000002)
        v_plt.get_patch_by_id('11').set_color("#D1202B")
        v_plt.get_patch_by_id('11')._linewidth = 1
        v_plt.get_patch_by_id('11')._edgecolor = (0.9137254901960784, 0.11372549019607843, 0.18823529411764706, 0.40000000000000002)
        return v_plt
    
    def new_oo_vs_old_fbf2(self, ax):
        a = self.targs['old_fbf2']
        b = self.targs['oo_both']
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{2}$',
                '$\mathrm{^{25^\circ C}Oo. FBF}$')
        v_plt = venn2(subsets=(a, b), ax=ax, set_labels=_lab)
        return self.set_venn2_red(v_plt)

    def new_oo_vs_old_both_fbfs(self, ax):
        if 'OO FBF (25°C)' in self.targs:
            #a_1 = self.targs['OO FBF-1 (20°C)']
            #a_2 = self.targs['OO FBF-2 (20°C)']
            #a = a_1 & a_2
            a = self.targs['OO FBF (20°C, both FBFs)']
            b = self.targs['OO FBF (25°C)']            
        else:
            a_1 = self.targs['old_fbf1']
            a_2 = self.targs['old_fbf2']
            a = a_1 & a_2
            b = self.targs['oo_both']
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$',
                '$\mathrm{^{25^\circ C}Oo. FBF}$')
        v_plt = venn2(subsets=(a, b), ax=ax, set_labels=_lab)
        return self.set_venn2_red(v_plt)
    
    def new_sp_vs_old_fbf(self, ax):
        a = self.targs['old_fbf1']
        b = self.targs['sp_both']
        c = self.targs['old_fbf2']
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{1}$',
                '$\mathrm{^{25^\circ C}Sp. FBF}$',
                '$\mathrm{^{20^\circ C}Oo. FBF}$-$\mathrm{2}$')
        v_plt = venn3(subsets=(a, b, c), ax=ax, set_labels=_lab)
        v_plt.get_patch_by_id('100').set_color("green")#("#933b41")
        v_plt.get_patch_by_id('100')._linewidth = 0
        v_plt.get_patch_by_id('010').set_color("blue") ##75968f
        v_plt.get_patch_by_id('010')._linewidth = 0
        v_plt.get_patch_by_id('001').set_color('#933b41')
        v_plt.get_patch_by_id('001')._linewidth = 0
        return v_plt
        
    def new_oo_fbf1_vs_fbf2(self, ax):
        a = self.targs['oo_fbf1']
        b = self.targs['oo_fbf2']
        _lab = ('$\mathrm{^{25^\circ C}Oo. FBF}$-$\mathrm{1}$',
                '$\mathrm{^{25^\circ C}Oo. FBF}$-$\mathrm{2}$')
        v_plt = venn2([a, b], ax=ax, set_labels=_lab)
        v_plt.get_patch_by_id('10').set_color("green")
        v_plt.get_patch_by_id('01').set_color('#933b41') 
        v_plt.get_patch_by_id('10')._linewidth = 0
        v_plt.get_patch_by_id('01')._linewidth = 0
        return v_plt

    def new_sp_fbf1_vs_fbf2(self, ax):
        a = self.targs['sp_fbf1']
        b = self.targs['sp_fbf2']
        _lab = ('$\mathrm{^{25^\circ C}Sp. FBF}$-$\mathrm{1}$',
                '$\mathrm{^{25^\circ C}Sp. FBF}$-$\mathrm{2}$')
        v_plt = venn2([a, b], ax=ax, set_labels=_lab)
        v_plt.get_patch_by_id('10').set_color("green")
        v_plt.get_patch_by_id('01').set_color('#933b41') 
        v_plt.get_patch_by_id('10')._linewidth = 0
        v_plt.get_patch_by_id('01')._linewidth = 0
        return v_plt

    def venn_combined_programs_vs_sp_oo(self, ax, programs):
        assert(type(programs)==type(dict()))
       # all_targ = self.targs['sp_both'] | self.targs['oo_both']
        _lab = ['SP/OO programs', 'Oo FBF targets', 'Sp FBF targets']
        v = venn3([
            set(programs.keys()), self.targs['oo_both'], self.targs['sp_both'], 
            ], ax=ax, set_labels=_lab)
        return v
    
class vennMaker(figureMaker.figureMaker, vennPainter, utils.translator):
    
    def top_n_ranks(self, n=500):
        
        for list_name, df in self.df.items():
            print("input {0} {1}".format( list_name, len(df.index)))

            df.drop_duplicates(subset='Gene name', keep='first', inplace=True)
            self.df[list_name] = df.head(n).copy()
            print(">", len(self.df[list_name].index))
            self.targs[list_name] = set(self.df[list_name]['Gene name'].tolist())
            
    def compare_biotypes_assignments(self):
        all_targs = self.targs['oo_both'] | self.targs['sp_both']
        self.transl()  # utils.translator. Sets self.name_to_biotype dict.
        
        agree, disagree = (set(), set())
        
        for name in list(all_targs):
            gtf_biotype = self.name_to_biotype.get(name, '')
            peaks_file_biotype = self.biotypes.get(name, '')
            if gtf_biotype == peaks_file_biotype:
                agree.add(name)
            else:
                disagree.add(name)
        
        print("The peaks file and gtf agreed on the biotypes of {0} RNAs, and disagreed on {1}".format(
            len(agree), len(disagree)))
        
        print("The disagreements were: ")
        for name in disagree:
            print("{0}: peaks: {1} gtf {2}".format(name, self.biotypes.get(name, ''), self.name_to_biotype.get(name, '')))
            
    def oo_vs_sp(self, ax, make_figure=True, look_at_biotypes=True, **kwargs):
        
        if ('a' not in kwargs) and ('b' not in kwargs):
            a, b = (self.targs['oo_both'], self.targs['sp_both'])
        else:
            a, b = (kwargs['a'], kwargs['b'])
            
        print("---\n{0} RNAs were targets in SP or OO germlines.".format(len(a | b)))

        num_total = len(a | b)

        print("Among those, {} ({:.3}%) were OO only, {} ({:.3}%) were SP only, and {} ({:.3}%) were shared.\n".format(
            len(a - b), 100* len(a - b)/num_total,
            len(b - a), 100* len(b - a)/num_total,
            len(a & b), 100* len(a & b)/num_total))
        
        if look_at_biotypes:
            mrna_a = set([name for name in a if self.biotypes.get(name, '') == 'protein_coding'])
            mrna_b = set([name for name in b if self.biotypes.get(name, '') == 'protein_coding'])
            print("For mRNA:")
            self.oo_vs_sp(None, make_figure=False, look_at_biotypes=False, a=mrna_a, b=mrna_b)
            print('')
            print("For ncRNA:")
            self.oo_vs_sp(None, make_figure=False, look_at_biotypes=False, a=a-mrna_a, b=b-mrna_b)
            print('')
            
        if make_figure:
            
            _lab = ('$\mathrm{^{25^\circ C}Oo. FBF}$',
                    '$\mathrm{^{25^\circ C}Sp. FBF}$')

            v_plt = venn2([a, b], ax=ax, set_labels=_lab)
            v_plt.get_patch_by_id('10').set_color("r")
            v_plt.get_patch_by_id('10')._linewidth = 0
            v_plt.get_patch_by_id('01').set_color('blue') ##75968f
            v_plt.get_patch_by_id('01')._linewidth = 0

            return v_plt
        
    def programs_vs_combined_sp_oo(self, ax, programs):
        assert(type(programs)==type(dict()))
        all_targ = self.targs['sp_both'] | self.targs['oo_both']
        oo_prog = set([x for x in programs if re.search('[Oo]ogenic', programs[x])])
        sp_prog = set(
            [x for x in programs if re.search('[Ss]permatogenic', programs[x])])

        _lab = ['SP/OO FBF targets', 'Sp program', 'Oo program']
        v = venn3([all_targ, sp_prog, oo_prog, ], ax=ax, set_labels=_lab)
        gl_rnas = set(programs.keys())
        #v = venn2([all_targ, oo_prog], ax=ax, set_labels=_lab)

        print("""There are {0} targets in SP or OO germlines. There are
 {1} RNAs in a SP or OO program. {2} ({3}%) of those are targets in the SP or
 OO FBF datasets.
        """.format(
            len(all_targ), len(programs), len(all_targ & gl_rnas),
            int(1000*len(all_targ & gl_rnas)/len(gl_rnas))/10.
            ))

        print("""only targ {0}. only oo prog {1} only sp prog {2}. 
 targ+oo {3} targ+sp {4} not targ but sp+oo {5}. targ+sp+oo {6}""".format(
        len(all_targ - oo_prog - sp_prog), len(oo_prog - all_targ -sp_prog),
        len(sp_prog - all_targ - oo_prog), len(all_targ & oo_prog - sp_prog),
        len(all_targ & sp_prog - oo_prog), len(sp_prog & oo_prog - all_targ),
        len(all_targ & oo_prog & sp_prog)))
        return v


    def heatmap_binary(self):
        _dict = collections.defaultdict(dict)
        def _add(_d, n): _d[n] = 1
        for k in self.targs:
            for name in self.targs[k]: _dict[name][k] = 1
        for name in _dict:
            for k in [_ for _ in self.targs if (_ not in _dict[name])]:
                _dict[name][k] = 0
            #_dict[name]['gene_name'] = name
        df = pandas.DataFrame(_dict.values(), index=_dict.keys())
        return df


if __name__ == '__main__':
    
    v = vennMaker()
    v.load_peaks_csv_files()
    
    # Use the peaksList object just to get sp/oo programs.
    pk = peaksList(dataframe=pandas.read_csv(
        'combined_filtered/oo_both.txt', sep='\t', index_col=False))
    pk.read_sp_vs_oo_as_programs()
    programs = pk.programs_as_public_names()
    
    fig, ax = plt.subplots()
    v.programs_vs_combined_sp_oo(ax, programs)
    
#    v.combined_programs_vs_sp_oo(ax[1], programs)
#    v.oo_vs_sp(ax[2])
    plt.tight_layout()
    plt.savefig('figs/Programs vs targets venn.pdf')
    plt.clf()
    
    
    plt.rc('font', size=10)
    fig, ax = plt.subplots(2, 2)
    v.new_oo_vs_old_fbf1(ax[0][0])
    v.new_oo_vs_old_fbf2(ax[0][1])
    v.new_oo_vs_old_both_fbfs(ax[1][0])
    v.oo_vs_sp(ax[1][1])
    plt.tight_layout()
    plt.savefig('figs/Fig 1 venns.pdf')
    plt.clf()
    plt.close()
    
    plt.rc('font', size=10)
    fig, ax = plt.subplots(2, 3)
    v.old_fbf1_vs_fbf2(ax[0][0])
    v.new_oo_fbf1_vs_fbf2(ax[0][1])
    v.new_sp_fbf1_vs_fbf2(ax[0][2])
    v.new_oo_vs_old_fbf(ax[1][0])
    v.new_sp_vs_old_fbf(ax[1][2])
    v.oo_vs_sp(ax[1][1])
    plt.tight_layout()
    plt.savefig('figs/Fig S1 venns.pdf')
    plt.clf()
    plt.close()
    #df = v.heatmap_binary()
    sys.exit()

    # Convert DataFrame to matrix
    mat = df.as_matrix()
    # Using sklearn
    import sklearn
    from sklearn import cluster
    km = cluster.KMeans(n_clusters=8)
    km.fit(mat)
    # Get cluster assignment labels
    labels = km.labels_
    # Format results as a DataFrame
    results = pandas.DataFrame([df.index,labels]).T
    df['cluster'] = results[results.columns[1]].tolist()
    df.sort(columns=['cluster'], inplace=True)
    del df['cluster'] 
    #results = pandas.DataFrame(data=labels, columns=['cluster'], index=collapsed.index)

    plt.pcolor(df)
    #plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
    plt.savefig('figs/heatmap.pdf')
     

    import numpy as np 
    from pandas import DataFrame
    import matplotlib.pyplot as plt

    a = '''
    Index= ['aaa', 'bbb', 'ccc', 'ddd', 'eee']
    Cols = ['A', 'B', 'C', 'D']
    df = DataFrame(abs(np.random.randn(5, 4)), index=Index, columns=Cols)

    plt.pcolor(df)
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
    plt.show()'''
