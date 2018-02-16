
import sys
import os
import pandas
import collections
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
import random
import volcanoMaker

from blocks import blocki, blockii, blockiii

#blocki = set(blocks.blocki)
#blockii = set(blocks.blockii)
#blockiii = set(blocks.blockiii)
all_blocks = blocki | blockii | blockiii

print("I {i} II {ii} III {iii}".format(
    i=len(blocki), ii=len(blockii), iii=len(blockiii)))
p = volcanoMaker.volcanoMaker()
p.read_clip_deseq_csv('tables/6_reps_sp_vs_oo.txt')
p.read_programs()
df = p.df[p.df['has_ortiz']].copy()


def as_p(x, y):
    if type(x) == type(pandas.DataFrame()):
        x = len(x.index)
    if y == 0: return 0
    if int(float(100*x)/float(y)) == 0:
        return "{} (0.{}%)".format(x, int(float(1000*x)/float(y)))
    return "{} ({}%)".format(x, int(float(100*x)/float(y)))


block_i_df = df[[(name in blocki) for name in df.gene_name]].copy()
block_ii_df = df[[(name in blockii) for name in df.gene_name]].copy()
block_iii_df = df[[(name in blockiii) for name in df.gene_name]].copy()


def row(_df, name):
    n = len(_df.index)
    d = {'Total genes': len(_df.index),
            'Oogenic only': as_p(_df[_df['Program']=='Oogenic only'], n),
'Spermatogenic only': as_p(_df[_df['Program']=='Spermatogenic only'], n),
'Oogenic and Spermatogenic': as_p(
    _df[_df['Program']=='Oogenic and Spermatogenic'], n),}
    print(("""{name} (n={n})
\t{oo}\t{oo_n}
\t{sp}\t{sp_n}
\t{oo_sp}\t{oo_sp_n}""".format(
    name=name, n=n, oo='Oogenic only', oo_n=d['Oogenic only'],
    sp='Spermatogenic only', sp_n=d['Spermatogenic only'],
    oo_sp='Oogenic and Spermatogenic',
    oo_sp_n=d['Oogenic and Spermatogenic'])))

for name, _df in [('Block I', block_i_df), ('Block II', block_ii_df),
             ('Block III', block_iii_df)]:
    row(_df, name)
