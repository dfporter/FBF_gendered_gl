import glob, pickle, seaborn, scipy, re, collections, pandas, importlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as scs


rl = importlib.reload

from collections import defaultdict
#from scipy.spatial import distance
#from scipy.cluster import hierarchy
import figureMaker
import utils
from peaksList import peaksList
import vennMaker
from blocks import blocki, blockii, blockiii, blockiv
import countsColumnsNaming
import countsFileIO
import figureMaker

rl(figureMaker)
rl(countsFileIO)
rl(countsColumnsNaming)


class spearmanCorrelationsMaker(
	countsFileIO.countsFileIO,
	figureMaker.figureMaker):

	def __init__(self):
		pass

	@staticmethod
	def top_n_genes(df, n=4000):
		df['sum'] = df.copy().sum(axis=1)
		df.sort_values(by='sum', inplace=True, ascending=False)
		del df['sum']
		print(df.head())

		return df.head(n).copy()

	def subset_counts_to_targets(self):
		all_targets = set()

		for target_set_name, targ_set in self.targs.items():
			all_targets |= targ_set

		self.counts_df = self.counts_df[[
			(x in all_targets) for x in self.counts_df.index]]

	def make_spearman_correlations_figure(self, columns=[], figname='default.pdf', only_20C=False):

		df = self.counts_df

		print('counts columns: {}'.format(df.columns))

		with sns.plotting_context("paper", font_scale=0.4):
			#df = v.load_counts_file()
			df = df[columns].copy()
			#col_relabel_d = col_relabel()
			#df.columns = [col_relabel_d[x] for x in df.columns]

			#df = self.top_n_genes(df)

			print("Making Spearman correlations from {} RNAs.".format(df.shape[0]))

			# Create correlation matrix.
			pwdf = df.corr('spearman')

			# Make it a diagonal, not a square.
			mask = np.zeros_like(pwdf)
			mask[np.triu_indices_from(mask)] = True
			np.place(mask, mask<1, 2)
			np.place(mask, mask==1, 0)
			np.place(mask, mask==2, 1)

			# Define color map.
			cmap = sns.cubehelix_palette(
				len(df.index), light=.7, dark=.2, reverse=False,
				start=1, rot=-2, as_cmap=True)
			cmap = plt.get_cmap('Greys')

			# Make the spearman rho correlations figure.
			with sns.axes_style("white"):
				ax = sns.heatmap(pwdf, cmap=cmap, annot=True,
				#col_cluster=False,
				#row_cluster=False,
				mask=mask)

			plt.savefig(figname)
			plt.clf()
			plt.close()


		self.fisher_test(pwdf, only_20C=only_20C)


	@staticmethod
	def fisher_test(df, only_20C=False):
		try:
			sp_vs_sp = [
				df.loc['SP FBF_1', 'SP FBF_2'],
				df.loc['SP FBF_1', 'SP FBF_3'],
				df.loc['SP FBF_2', 'SP FBF_3'],
				]

			oo_vs_oo = [
				df.loc['OO FBF_1', 'OO FBF_2'],
				df.loc['OO FBF_1', 'OO FBF_3'],
				df.loc['OO FBF_2', 'OO FBF_3'],
			]
			sp_vs_oo = [
				df.loc['SP FBF_1', 'OO FBF_1'],
				df.loc['SP FBF_1', 'OO FBF_2'],
				df.loc['SP FBF_1', 'OO FBF_3'],
				df.loc['SP FBF_2', 'OO FBF_1'],
				df.loc['SP FBF_2', 'OO FBF_2'],
				df.loc['SP FBF_2', 'OO FBF_3'],
				df.loc['SP FBF_3', 'OO FBF_1'],
				df.loc['SP FBF_3', 'OO FBF_2'],
				df.loc['SP FBF_3', 'OO FBF_3'],

			]

			lt_fbf1_vs_fbf2 = [
				df.loc['LT FBF1_1', 'LT FBF2_1'],
				df.loc['LT FBF1_1', 'LT FBF2_2'],
				df.loc['LT FBF1_1', 'LT FBF2_3'],
				df.loc['LT FBF1_2', 'LT FBF2_1'],
				df.loc['LT FBF1_2', 'LT FBF2_2'],
				df.loc['LT FBF1_2', 'LT FBF2_3'],
				df.loc['LT FBF1_3', 'LT FBF2_1'],
				df.loc['LT FBF1_3', 'LT FBF2_2'],
				df.loc['LT FBF1_3', 'LT FBF2_3'],
			]

			lt_vs_ht = [
				df.loc['LT FBF1_1', 'OO FBF_1'],
				df.loc['LT FBF1_1', 'OO FBF_2'],
				df.loc['LT FBF1_1', 'OO FBF_3'],
				df.loc['LT FBF1_2', 'OO FBF_1'],
				df.loc['LT FBF1_2', 'OO FBF_2'],
				df.loc['LT FBF1_2', 'OO FBF_3'],
				df.loc['LT FBF1_3', 'OO FBF_1'],
				df.loc['LT FBF1_3', 'OO FBF_2'],
				df.loc['LT FBF1_3', 'OO FBF_3'],
				df.loc['LT FBF2_1', 'OO FBF_1'],
				df.loc['LT FBF2_1', 'OO FBF_2'],
				df.loc['LT FBF2_1', 'OO FBF_3'],
				df.loc['LT FBF2_2', 'OO FBF_1'],
				df.loc['LT FBF2_2', 'OO FBF_2'],
				df.loc['LT FBF2_2', 'OO FBF_3'],
				df.loc['LT FBF2_3', 'OO FBF_1'],
				df.loc['LT FBF2_3', 'OO FBF_2'],
				df.loc['LT FBF2_3', 'OO FBF_3'],


			]

			same_gender = oo_vs_oo + sp_vs_sp

			print("25C SP vs SP correlations: {}\n\tAverage: {}".format(sp_vs_sp, np.mean(sp_vs_sp)))
			print("25C OO vs OO correlations: {}\n\tAverage: {}".format(oo_vs_oo, np.mean(oo_vs_oo)))
			print("25C SP vs OO correlations: {}\n\tAverage: {}".format(sp_vs_oo, np.mean(sp_vs_oo)))
			print("25C OO vs 20C OO correlations: {}\n\tAverage: {}".format(lt_vs_ht, np.mean(lt_vs_ht)))
			
			res = scs.ttest_ind(same_gender, sp_vs_oo, )
			print("SP vs SP and OO vs OO rho values cf. SP vs OO rho values t-test: {}".format(res))

			res = scs.ttest_ind(oo_vs_oo, sp_vs_oo)
			print("OO vs OO rho values cf. SP vs OO rho values t-test: {}".format(res))

			res = scs.ttest_ind(sp_vs_sp, sp_vs_oo)
			print("SP vs SP rho values cf. SP vs OO rho values t-test: {}".format(res))

			res = scs.ttest_ind(sp_vs_oo, lt_fbf1_vs_fbf2)
			print("OO vs SP and 20C FBF-1 vs FBF-2 rho values t-test: {}".format(res))

			res = scs.ttest_ind(sp_vs_oo, lt_vs_ht)
			print("OO vs SP and 20C vs 25C rho values t-test: {}".format(res))

		except:
			pass

		if only_20C:
			print(df.columns, df.index)
			lt_vs_ht = [
				df.loc['LT FBF1_1', 'LT FBF2_1'],
				df.loc['LT FBF1_1', 'LT FBF2_2'],
				df.loc['LT FBF1_1', 'LT FBF2_3'],
				df.loc['LT FBF1_2', 'LT FBF2_1'],
				df.loc['LT FBF1_2', 'LT FBF2_2'],
				df.loc['LT FBF1_2', 'LT FBF2_3'],
				df.loc['LT FBF1_3', 'LT FBF2_1'],
				df.loc['LT FBF1_3', 'LT FBF2_2'],
				df.loc['LT FBF1_3', 'LT FBF2_3'],
			]

			print("20C FBF-1 vs 20C FBF-2 correlations: {}\n\tAverage: {}".format(lt_vs_ht, np.mean(lt_vs_ht)))



if __name__ == '__main__':
	m = spearmanCorrelationsMaker()
	m.load_counts_file(log_scale=False, style='including_25C_single_replicates')
	print(m.counts_df.head())
	m.load_peaks_csv_files()
	m.subset_counts_to_targets()

	m.make_spearman_correlations_figure(columns=[
		'SP FBF1_1', 'SP FBF1_2', 'SP FBF1_3',
		'SP FBF2_1', 'SP FBF2_2', 'SP FBF2_3'],
		figname='figs/Spearman corr SP FBF1 vs FBF2.pdf')

	m.make_spearman_correlations_figure(columns=[
		'OO FBF1_1', 'OO FBF1_2', 'OO FBF1_3',
		'OO FBF2_1', 'OO FBF2_2', 'OO FBF2_3',],
		figname='figs/Spearman corr OO FBF1 vs FBF2.pdf')

	m.make_spearman_correlations_figure(columns=[
		'SP FBF_1', 'SP FBF_2', 'SP FBF_3',
		'OO FBF_1', 'OO FBF_2', 'OO FBF_3',
		'LT FBF1_1', 'LT FBF1_2', 'LT FBF1_3',
		'LT FBF2_1', 'LT FBF2_2', 'LT FBF2_3'],
		figname='figs/Figure S2 Spearman corr.pdf')

	m.make_spearman_correlations_figure(columns=[
		'LT FBF1_1', 'LT FBF1_2', 'LT FBF1_3',
		'LT FBF2_1', 'LT FBF2_2', 'LT FBF2_3'],
		figname='figs/20C Spearman corr.pdf', only_20C=True)
	