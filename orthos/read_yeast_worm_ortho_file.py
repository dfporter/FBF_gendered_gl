"""
The file ortho_map has a crazy format.
This script parses it.
"""
import os
import sys
import re
import collections
import pandas


def get_translation_table(
        fname, delim='------------------------------------------------------------------------',
        use_star=True):
    transl_yeast_to_worm = {}
    transl_worm_to_yeast = {}
    with open(fname, 'r') as f:
        for pair in olines(f, delim):
            process_pair(pair, transl_yeast_to_worm,
                         transl_worm_to_yeast, use_star=use_star)
    return transl_worm_to_yeast, transl_yeast_to_worm


def olines(f, newline):
    buf = ''
    while True:
        while newline in buf:
            pos = buf.index(newline)
            yield buf[:pos]
            buf = buf[pos + len(newline):]
        chunk = f.read()
        if not chunk:
            yield buf
            break
        buf += chunk


def process_pair(pair, transl_yeast_to_worm, transl_worm_to_yeast,
                 use_star=True):
    if use_star:
        if '*' not in pair:
            print 'no *'
            return
        star_pos = pair.index('*')
    elif 'Score diff' in pair:
        score_pos = pair.index('Score diff')
        for pos in range(score_pos, len(pair)):
            if pair[pos] == '\n':
                star_pos = pos + 1
                break
    else:
        print "Couldn't find gene."
        return
    boot_pos = pair.index('Bootstrap')
    chunk = pair[star_pos:boot_pos]
    left = set()
    right = set()
    for li in chunk.split('\n'):
        li = re.sub(' ', '', li)
        li = re.sub('\*', '', li)
        s = li.split('\t')
        left.add(s[0])
        if len(s) > 3:
            right.add(s[3])
    for gene in right:
        transl_yeast_to_worm.setdefault(gene, set())
        transl_yeast_to_worm[gene] |= left
        transl_yeast_to_worm[gene] -= set([''])
    for gene in left:
        transl_worm_to_yeast.setdefault(gene, set())
        transl_worm_to_yeast[gene] |= right
        transl_worm_to_yeast[gene] -= set([''])
