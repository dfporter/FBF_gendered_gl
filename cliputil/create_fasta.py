import xlwt
import glob
import sys
import pandas
import collections
import re
import os
import numpy as np


def _mk(t):
    if os.path.isfile(t):
        t = os.path.dirname(t)
    if not os.path.exists(t):
        os.system('mkdir ' + t)
_mk('fasta/')
_mk('fasta/top500')

dfs = {}
seqs = {}
topseqs = {}
for fname in glob.glob('combined_filtered/*'):
    f = re.sub('\.txt', '', os.path.basename(fname))
    dfs[f] = pandas.read_csv(fname, sep='\t')
    dfs[f] = dfs[f][dfs[f]['biotype']=='protein_coding']
    dfs[f]['height'] = dfs[f]['exp_reads']
    dfs[f].sort(columns=['height'], inplace=True, ascending=False)
    tups = zip(dfs[f]['gene_name'], dfs[f].index, dfs[f]['seq'].tolist())
    seqs[f] = [">{0}__{1}\n{2}\n".format(a, b, c) for a, b, c in tups]
    topseqs[f] = seqs[f][:min(len(seqs[f]), 500)]
    seqs[f] = "".join(seqs[f])
    topseqs[f] = "".join(topseqs[f])
    # Write.
    print topseqs[f]
    with open('fasta/{0}.fa'.format(f), 'w') as fh:
        fh.write(seqs[f])
    with open('fasta/top500/{0}.fa'.format(f), 'w') as fh:
        fh.write(topseqs[f])

