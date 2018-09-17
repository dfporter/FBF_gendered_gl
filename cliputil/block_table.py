

import peaksList
from collections import defaultdict
import pandas as pd


import numpy as np

from volcanoMaker import *



from blocks import blocki, blockii, blockiii, blockiv


def sort_drop_dups(_df):
    _df.sort_values(by=['exp_reads'], ascending=False, inplace=True)
    _df.drop_duplicates(subset='gene_name', inplace=True)
    return _df


def compare_temp_changes_vs_abundance():

    oo_pk = peaksList.peaksList()
    oo_pk.read_csv('combined_filtered/oo_both.txt')
    oo_pk.add_wbid()
    oo_df = oo_pk.df.copy()

    sp_pk = peaksList.peaksList()
    sp_pk.read_csv('combined_filtered/sp_both.txt')
    sp_pk.add_wbid()
    sp_df = sp_pk.df.copy()

#    oo_df = pandas.read_csv('combined_filtered/oo_both.txt', sep='\t')
#    sp_df = pandas.read_csv('combined_filtered/sp_both.txt', sep='\t')

    oo_df = sort_drop_dups(oo_df)
    oo_targs = set(oo_df['gene_name'].tolist())

    sp_df = sort_drop_dups(sp_df)
    sp_targs = set(sp_df['gene_name'].tolist())

    i = sp_df

    def in_block(name, a_set):
        if name in a_set:
            return True
        return False

    i_peaks = sp_df[[in_block(_name, set(blocki)) for _name in sp_df.gene_name]].copy()

    robust_i_peaks = sp_df[[((n in set(blocki)) and (n in sp_targs) and \
                      (n not in oo_targs) \
                  and (exp >= 25)) for n, exp in \
                  zip(sp_df.gene_name, sp_df['exp_reads'].tolist())]].copy()


    only_sp_df = sp_df[[((n in set(blocki)) and (n in sp_targs) and \
                      (n not in oo_targs)) for n, exp in \
                  zip(sp_df.gene_name, sp_df['exp_reads'].tolist())]].copy()

    print("""Of the {} RNAs in Block I, {} were only targets of FBF in \
spermatogenic germlines, and {} of those were robust peaks.""".format(
    len(blocki), len(only_sp_df.index), len(robust_i_peaks.index)))

    common_df = oo_df[[((n in oo_targs) and (n in sp_targs)) for n in \
                     oo_df['gene_name'].tolist()]].copy()

    to_sp_exp = dict(list(zip(sp_df.gene_name, sp_df.exp_reads)))
    to_oo_exp = dict(list(zip(oo_df.gene_name, oo_df.exp_reads)))

    common_df['keep'] = [(n in set(blockii) and (n in sp_targs) \
                  and (n in oo_targs) and (to_oo_exp[n] >= 50) and \
                    (to_sp_exp[n] >= 50)) for n in common_df.gene_name]

#    ii_peaks = common_df[common_df['keep']].copy()
    ii_peaks = oo_df[[in_block(_name, set(blockii)) for _name in oo_df.gene_name]].copy()
    sp_only_ii = set(blockii) - set(ii_peaks.gene_name)

    sp_ii_peaks = sp_df[[in_block(_name, sp_only_ii) for _name in sp_df.gene_name]].copy()
    ii_peaks = pd.concat([ii_peaks, sp_ii_peaks])

#    common_df['keep'] = [(n in set(blockii) and (n in sp_targs) \
#                  and (n in oo_targs) and (to_oo_exp[n] >= 50)) \
#                         for n in common_df.gene_name]
    common_ii_df = common_df[common_df['keep']].copy()

    print("""Of the {0} RNAs in Block II, {1} were targets of FBF in \
spermatogenic and oogenic germlines, and had robust peaks.""".format(
        len(blockii), len(common_ii_df.index)))

    iii_peaks = oo_df[[in_block(_name, set(blockiii)) for _name in oo_df.gene_name]].copy()
    sp_only_iii = set(blockiii) - set(iii_peaks.gene_name)
    sp_iii_peaks = sp_df[[in_block(_name, sp_only_iii) for _name in sp_df.gene_name]].copy()
    iii_peaks = pd.concat([iii_peaks, sp_iii_peaks])
    
    common_df['keep'] = [(n in set(blockiii) and (n in sp_targs) \
                  and (n in oo_targs) and (to_oo_exp[n] >= 50) and \
                    (to_sp_exp[n] >= 50)) for n in common_df.gene_name]
    robust_iii_peaks = common_df[common_df['keep']].copy()

    common_df['keep'] = [(n in set(blockiii) and (n in sp_targs) \
                  and (n in oo_targs)) \
                         for n in common_df.gene_name]
    common_iii_df = common_df[common_df['keep']].copy()

    print("""Of the {0} RNAs in Block III, {1} were  targets of FBF in \
spermatogenic and oogenic germlines, and {2} were robust peaks.""".format(
    len(blockiii), len(common_iii_df.index), len(robust_iii_peaks.index)))

    dont_write = """
    writer =  pandas.ExcelWriter(
        'tables/Table S6 Blocks.xls')

    i_b = pandas.DataFrame(
        list(zip(list(range(len(blocki))), blocki)), columns=['index', 'Gene name'] )
    i_b.to_excel(
        writer, sheet_name='Block I', columns=['Gene name'], index=False)
    ii_b = pandas.DataFrame(list(zip(list(range(len(blockii))), blockii)), columns=['index', 'Gene name'] )
    ii_b.to_excel(
        writer, sheet_name='Block II', columns=['Gene name'], index=False)
    iii_b = pandas.DataFrame(list(zip(list(range(len(blockiii))), blockiii)), columns=['index', 'Gene name'] )
    iii_b.to_excel(
        writer, sheet_name='Block III', columns=['Gene name'], index=False)
    i_peaks.to_csv('tables/i_peaks.txt', sep='\t')
    ii_peaks.to_csv('tables/ii_peaks.txt', sep='\t')
    iii_peaks.to_csv('tables/iii_peaks.txt', sep='\t')


    i_peaks.to_excel(
        writer, sheet_name='Block I subnetwork', #columns=['gene_name'],
        index=False)
    ii_peaks.to_excel(
        writer, sheet_name='Block II subnetwork', #columns=['gene_name'],
        index=False)
    iii_peaks.to_excel(
        writer, sheet_name='Block III subnetwork', #columns=['gene_name'],
        index=False)

    writer.save()        """



def run():
    #v = volcanoMaker()
    compare_temp_changes_vs_abundance()


if __name__ == '__main__':
    run()
