import os
import sys
import re
import collections
import pandas
from read_yeast_worm_ortho_file import *


def read_yeast_pufs_as_ensembl_ids(fname):
    targs = collections.defaultdict(set)
    with open(fname, 'r') as f:
        _ = next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            targs['puf3'].add(s[0])
            if len(s) > 1:
                targs['puf4'].add(s[1])
            if len(s) > 2:
                targs['puf5'].add(s[2])
    return targs

def get_puf3_targets():
    """Return a set. Tagging."""
    return read_yeast_pufs_as_ensembl_ids(
        'YeastPUF-TaggingTargets.txt')['puf3']

def get_puf5_targets():
    """Return a set. CLIP."""
    return(
        set([x.rstrip('\n') for x in open('puf5_targets.txt')]) - set([''])
        )

def get_set_of_all_yeast_tagging_targets():
    # Yeast targets
    yeast_targs = read_yeast_pufs_as_ensembl_ids('YeastPUF-TaggingTargets.txt')
    all_yeast_targs = set()
    for a_set in [yeast_targs[protein] for protein in yeast_targs]:
        all_yeast_targs |= a_set  
    return all_yeast_targs
