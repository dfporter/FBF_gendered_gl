import pandas
import csv
import re


def get_ripchip_targets():
    ripchip_filename = 'lib/ripchip/sd01.txt'
    ripchip = pandas.read_csv(
        ripchip_filename, sep='\t')
    top_1350 = list()
    for index, row in ripchip.iterrows():
        if row['Gene public name'] not in top_1350:
            top_1350.append(row['Gene public name'])
        if len(list(set(top_1350))) >= 1350:
            break
    #top_1350_rows = ripchip[0:1349]  # There is a row zero.
    #rip_targets = list(top_1350['Gene public name'])
    print("Top targets: {t}.\nLast targets: {b}.".format(
        t=str(top_1350[0:2]), b=str(top_1350[-2:])))
    return top_1350

def add_column_of_overlap(df, rip_targets):
    rip_targets = set(rip_targets)
    df['is_rip_chip_target'] = 0
    for index, row in df.iterrows():
        if df.loc[index, 'gene_name'] in rip_targets:
            df.loc[index, 'is_rip_chip_target'] = 1
    
