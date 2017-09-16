from __future__ import division
from peaksList import peaksList 
import annotate_peaks
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas
import sys
import os
import scipy.stats as scs

combined_bed_size = {
    'exp_fbf1': 9177044,
    'exp_fbf2': 4787043,
    'exp_fbf_oo': 3545033,
    'exp_fbf_sp': 4119309,
}

def threeD(m):
    # Argument is peaksList
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(
        m.df['unnorm_old_fbf1'].tolist(),
        m.df['unnorm_old_fbf2'].tolist(),
        m.df['unnorm_oo_both'].tolist(),
               c='k', marker='.')
    ax.set_xlabel('20 deg FBF-1')
    ax.set_ylabel('20 deg FBF-2')
    ax.set_zlabel('25 deg FBF')
    ax = fig.add_subplot(111, projection='3d')
    plt.savefig('figs/3d_scatter_old_1_and_2_vs_new_oo.pdf')
    plt.clf(); plt.close()

def get_table_of_merged(filename_a, filename_b, filename_c=None,
                        filename_d=None):
    a = peaksList()
    a.read_csv(filename_a)
    b = peaksList()
    b.read_csv(filename_b)

    m = a.merge_ranges_with_other_peaksList(b, use_name='gene_name')
    df = pandas.DataFrame(m, columns=['chrm', 'left', 'right', 'strand',
                                      'gene_name'])
    print df
    m = peaksList(dataframe=df)
    if filename_c is not None:
        c = peaksList()
        c.read_csv(filename_c)
        _m = c.merge_ranges_with_other_peaksList(m, use_name='gene_name')
        df = pandas.DataFrame(_m, columns=['chrm', 'left', 'right', 'strand',
                                      'gene_name'])
        m = peaksList(dataframe=df)
    if filename_d is not None:
        d = peaksList()
        d.read_csv(filename_d)
        _m = d.merge_ranges_with_other_peaksList(m, use_name='gene_name')
        df = pandas.DataFrame(_m, columns=['chrm', 'left', 'right', 'strand',
                                      'gene_name'])
        m = peaksList(dataframe=df)

    ga = {}
    ga['old_fbf2'] = annotate_peaks.load_bedgraph(
        'bedgraphs_norm/exp_fbf1.wig')
    ga['old_fbf1'] = annotate_peaks.load_bedgraph(
        'bedgraphs_norm/exp_fbf2.wig')
    m.add_reads(ga=ga['old_fbf2'], name='norm_old_fbf2')
    m.add_reads(ga=ga['old_fbf1'], name='norm_old_fbf1')
#    m.df['norm_old_fbf1'] = [
#        1e6 * x/combined_bed_size['exp_fbf1'] for x in m.df.unnorm_old_fbf1]
#    m.df['norm_old_fbf2'] = [
#        1e6 * x/combined_bed_size['exp_fbf2'] for x in m.df.unnorm_old_fbf2]
    if filename_c is not None:
        ga['oo_both'] = annotate_peaks.load_bedgraph(
            'bedgraphs_norm/exp_fbf_oo.wig',#all_exp.wig'
        )
        m.add_reads(ga=ga['oo_both'], name='norm_oo_both')
#        m.df['norm_oo_both'] = [
#            1e6 * x/combined_bed_size['exp_fbf_oo'] for x in m.df.unnorm_oo_both]
    if filename_d is not None:
        ga['sp_both'] = annotate_peaks.load_bedgraph(
            'bedgraphs_norm/exp_fbf_sp.wig'
            )
        m.add_reads(ga=ga['sp_both'], name='norm_sp_both')
#        m.df['norm_sp_both'] = [
#            1e6 * x/combined_bed_size['exp_fbf_sp'] for x in m.df.unnorm_sp_both]
    #tups = zip(m.df.unnorm_old_fbf1, m.df.unnorm_old_fbf2)
    #m.df['unnorm_old_fbf'] = [sum(t)/2. for t in tups]
    tups = zip(m.df.norm_old_fbf1, m.df.norm_old_fbf2)
    m.df['norm_old_fbf'] = [sum(t)/2. for t in tups]
    def log_or_nan(n):
        if n <= 0: return np.nan
        else: return np.log(n)
    for name in [#'unnorm_old_fbf1', 'unnorm_old_fbf2', 'unnorm_oo_both',
                 'norm_old_fbf1', 'norm_old_fbf2', 'norm_oo_both',
                 #'unnorm_old_fbf',
                 'norm_old_fbf']:
        m.df[name] = [log_or_nan(x) for x in m.df[name].tolist()]
    if filename_d is not None:
        for name in ['norm_sp_both', ]:
            m.df[name] = [log_or_nan(x) for x in m.df[name].tolist()]
    return m

def sp_and_oo_mode():
    if raw_input('Remake merged oo.+sp. regions [y/n]?').upper()[0] == 'Y':
        m = get_table_of_merged(
            'combined_filtered/old_fbf1_to_fbf2_n2.txt',
            'combined_filtered/old_fbf2.txt',
            filename_c='combined_filtered/oo_both.txt',
            filename_d='combined_filtered/sp_both.txt'
            )
        if not os.path.exists('data/'): os.system('mkdir data')
        m.to_csv('data/sp_and_oo_merged_peak_regions_with_reads_in_peak.txt')
    m = peaksList()
    m.read_csv('data/sp_and_oo_merged_peak_regions_with_reads_in_peak.txt')
    df = m.df[[#'unnorm_oo_both', 'unnorm_old_fbf1', 'unnorm_old_fbf2',
           #'unnorm_sp_both',
               'norm_oo_both', 'norm_old_fbf1', 'norm_old_fbf2', 'norm_sp_both',
               #'unnorm_old_fbf',
               'norm_old_fbf'
               ]].copy()
    return df

def oo_mode():
    if raw_input('Remake merged oo. regions [y/n]?').upper()[0] == 'Y':
        m = get_table_of_merged(
            'combined_filtered/old_fbf1_to_fbf2_n2.txt',
            'combined_filtered/old_fbf2.txt',
            'combined_filtered/oo_both.txt'
            )
        if not os.path.exists('data/'): os.system('mkdir data')
        m.to_csv('data/oo_merged_peak_regions_with_reads_in_peak.txt')
    m = peaksList()
    m.read_csv('data/oo_merged_peak_regions_with_reads_in_peak.txt')
    df = m.df[[#'unnorm_oo_both', 'unnorm_old_fbf1', 'unnorm_old_fbf2',
               'norm_oo_both', 'norm_old_fbf1', 'norm_old_fbf2',
               #'unnorm_old_fbf',
               'norm_old_fbf']].copy()
    return df

#df = sp_and_oo_mode()
df = oo_mode()
df['hue'] = 'k'

print "Added reads."

df.dropna(axis=0, how='any', inplace=True)
print df
xcol, ycol = 'norm_old_fbf', 'norm_oo_both'
xcol, ycol = 'norm_old_fbf1', 'norm_old_fbf2'
#xcol, ycol = 'norm_old_fbf1', 'norm_old_fbf2'
#xcol, ycol = 'norm_sp_both', 'norm_oo_both'
sk = '''
g = sns.lmplot(
    x=xcol, y=ycol,
    data=df,
    hue='hue',
    palette=dict(k="k"),
    scatter_kws=dict(alpha='0.1'),
    fit_reg=False,
    )'''
plt.clf(); plt.close()
x = df[xcol].tolist()
y = df[ycol].tolist()
xmin, ymin = 0, 0
xmax, ymax = 8, 8
extent = [xmin, xmax, ymin, ymax]
gs = matplotlib.gridspec.GridSpec(1,1)
ax = plt.subplot(gs[0, 0], aspect='equal')
#ax.set_xlim(xmin, xmax)
#ax.set_ylim(ymin, ymax)
ax.hexbin(x, y, cmap='Spectral', mincnt=1,
          extent=extent)#gridsize=30)# c='k', alpha=0.1, s=4)
ax.axis(extent)
slope, intercept, r_value, p_value, std_err = scs.linregress(x,y)
for x in [slope, intercept, r_value, p_value, std_err]:
    x_name = [k for k,v in locals().items() if v is x]
    print x_name, x
print slope, intercept, r_value, p_value, std_err
def t(x): return slope * x + intercept

#ax.plot([xmin, xmax], [ymin, ymax], 'k-', alpha=0.5)
#ax.plot([xmin, xmax], [t(ymin), t(ymax)], 'r--', alpha=.5)
#ax.plot([xmin, xmax], [ymin-np.log10(2), ymax-np.log10(2)], 'r:')
#ax.plot([xmin, xmax], [ymin+np.log10(2), ymax+np.log10(2)], 'r:')
sk = '''
g.set(xlim=(-2, 10))
g.set(ylim=(-2, 10))
def t(x): return slope * x + intercept
g.ax.plot([-2, 10], [-2, 10], 'k-')
g.ax.plot([-2, 10], [t(-2), t(10)], 'r--')
g.ax.plot([-2, 10], [-2-np.log10(2), 10-np.log10(2)], 'r:')
g.ax.plot([-2, 10], [-2+np.log10(2), 10+np.log10(2)], 'r:')
'''

#print plt.axes().get_data_ratio()
#print plt.axes().__dict__
#plt.axes().set_aspect(1./plt.axes().get_data_ratio())

plt.savefig('figs/2d_scatter_old_1_and_2_vs_new_oo.pdf')
plt.clf(); plt.close()


import seaborn as sns
sns.set_context('paper', font_scale=0.8)
sns.set(style="ticks")
ax = sns.regplot(
    x='norm_old_fbf', y='norm_oo_both', data=df)

#ax[0][1].sns.regplot(
#    x='unnorm_old_fbf2', y='oo_both', data=m.df)
#ax[1][1].sns.regplot(
#    x='unnorm_old_fbf2', y='oo_both', data=m.df)
#ax[0][0].set_xlabel('20 deg FBF-1')
#ax[0][0].set_ylabel('20 deg FBF-2')
#ax[1][0].set_ylabel('25 deg FBF')
plt.savefig('figs/2d_regression_old_1_and_2_vs_new_oo.pdf')
plt.clf(); plt.close()
