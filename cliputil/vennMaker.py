from __future__ import division
import pandas
import re
import os
import sys
import glob
import numpy as np
import collections
import peaksList
    
from peaksList import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scs
import matplotlib.lines as mlines
import matplotlib
from matplotlib_venn import venn2, venn3, venn3_circles
import matplotlib_venn
print(matplotlib_venn.__file__)

prefix = '/groups/Kimble/Common/fbf_celltype/combined_filtered/'
prefix = './combined_filtered/'


fname_to_label = {
prefix + 'old_fbf2.txt': 'old_fbf2',
prefix + 'old_fbf1.txt': 'old_fbf1',
prefix + 'oo_both.txt': 'oo_both',
prefix + 'oo_fbf1.txt': 'oo_fbf1',
prefix + 'oo_fbf2.txt': 'oo_fbf2',
prefix + 'sp_both.txt': 'sp_both',
prefix + 'sp_fbf1.txt': 'sp_fbf1',
prefix + 'sp_fbf2.txt': 'sp_fbf2',
prefix + 'old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
}
label_to_fname = dict(zip(fname_to_label.values(), fname_to_label.keys()))

class vennMaker(object):

    def __init__(self, output_dir='figs/'):
        if not os.path.exists(output_dir): os.system('mkdir ' + output_dir)
        self.df = {}
        
    def load(self, label_to_file):
        self.df = dict([
            (k, self.read_csv(v)) for k, v in label_to_file.items()])
        self.label_to_file = label_to_file
        self.targs = dict([
            (k, set(self.df[k]['gene_name'].tolist())) for k in self.df])

    def read_csv(self, fname):
        return pandas.read_csv(fname, sep='\t', index_col=False)

    def old_fbf1_vs_fbf2(self, ax):
        a = self.targs['old_fbf1_to_fbf2_n2']
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
        a = self.targs['old_fbf1_to_fbf2_n2']
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
        a = self.targs['old_fbf1_to_fbf2_n2']
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
        a_1 = self.targs['old_fbf1_to_fbf2_n2']
        a_2 = self.targs['old_fbf2']
        a = a_1 & a_2
        b = self.targs['oo_both']
        _lab = ('$\mathrm{^{20^\circ C}Oo. FBF}$',
                '$\mathrm{^{25^\circ C}Oo. FBF}$')
        v_plt = venn2(subsets=(a, b), ax=ax, set_labels=_lab)
        return self.set_venn2_red(v_plt)

    def oo_vs_sp(self, ax):
        a = self.targs['oo_both']
        b = self.targs['sp_both']
        _lab = ('$\mathrm{^{25^\circ C}Oo. FBF}$',
                '$\mathrm{^{25^\circ C}Sp. FBF}$')
        v_plt = venn2([a, b], ax=ax, set_labels=_lab)
        v_plt.get_patch_by_id('10').set_color("r")
        v_plt.get_patch_by_id('10')._linewidth = 0
        v_plt.get_patch_by_id('01').set_color('blue') ##75968f
        v_plt.get_patch_by_id('01')._linewidth = 0
        return v_plt
        
    def new_sp_vs_old_fbf(self, ax):
        a = self.targs['old_fbf1_to_fbf2_n2']
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

    def combined_programs_vs_sp_oo(self, ax, programs):
        assert(type(programs)==type(dict()))
       # all_targ = self.targs['sp_both'] | self.targs['oo_both']
        _lab = ['SP/OO programs', 'Oo FBF targets', 'Sp FBF targets']
        v = venn3([
            set(programs.keys()), self.targs['oo_both'], self.targs['sp_both'], 
            ], ax=ax, set_labels=_lab)
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
    v.load(label_to_fname)
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
    plt.savefig('Fig 1 venns.pdf')
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
    plt.savefig('Fig S1 venns.pdf')
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
    plt.savefig('heatmap.pdf')
     

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
