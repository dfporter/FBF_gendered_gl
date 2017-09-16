from __future__ import division
import scipy.stats as scs
from peaksList import *
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from volcanoMaker import *
from block_iii import blockiii

print blockiii


def run():
    v = volcanoMaker()
    v.gl_rnaseq()
    v.read_clip_deseq_csv('tables/6_reps_sp_vs_oo.txt')
    v.read_programs()
    df = v.df[v.df['has_ortiz']].copy()
    compare_temp_changes_vs_abundance(v)
    sys.exit()



if __name__ == '__main__':
    run()
