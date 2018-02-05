import pandas
from language import *
from annotater import *
import write_excel
import os
import re
import sys
import numpy as np
from fbfTargetSet import fbfTargetSet
from human_gene_ensmbl_map import human_gene_ensmbl_map


def read_ortholist_file(fname, max_orthologs=2):
    df = pandas.read_csv(fname, sep=',')
    cols = ['Ensembl Compara', 'InParanoid', 'Homologene', 'OrthoMCL']
    from_locus_id = dict((re.sub(' ', '', row['Locus ID']), row) \
                         for row in df.to_dict('records'))
    human_ensembl_to_locus_id = {}
    cmap = complexTranslation('Locus ID', 'Ensmbl', {},
                              max_orthologs=max_orthologs)
    for i, k in enumerate(from_locus_id):
        if not(i % 2e2): print i
        aset = set()
        for vals in [from_locus_id[k][x] for x in cols]:
            if isinstance(vals, type('a')):
                ensg_ids = [re.sub(' ', '', x) for x in re.split('[ ,]', vals)]
            else:
                continue
            aset |= set(ensg_ids)
            for ensg in ensg_ids:
                human_ensembl_to_locus_id[ensg] = k
            aset = aset - set([]) - set(['NULL', '', '-'])
            if len(aset) > 0:
                cmap.add_mapping(frozenset([k]), frozenset(aset))
        from_locus_id[k] = aset - set([]) - set(['NULL', '', '-'])
    print cmap.info()
    return from_locus_id, human_ensembl_to_locus_id, cmap


def pum2_targets_as_list_of_ensmbl_sets():
#    ortho_locus_to_ensmbl, human_ensembl_to_locus_id, ortho_cmap = orthos
    human_gene_symbol_to_human_ensmbl, human_ensmbl_to_human_gene_symbol = \
                                       human_gene_ensmbl_map()
    pum2_targets = set([x.rstrip('\n') for x in open('pum2_targets.txt', 'r')])
#    pum2_targets = set([x.rstrip('\n') for x in open('lib/twenty_lines_pum2_targets.txt', 'r')])
    pum2_targets -= set(['NULL', '', '-'])
    symb_to_ens_tr = translation(
        'Gene symbol', 'Ensmbl', human_gene_symbol_to_human_ensmbl)
    pum2_in_ensmbl = targetSet(pum2_targets, 'Gene symbol').translate(
        symb_to_ens_tr, 'Ensmbl')
#    print "PUM2 target lines: {0}. Gene symbols: {1}. Ensmbl
    return pum2_in_ensmbl


def pum2_targets_as_ensembl(orthos):
    ortho_locus_to_ensmbl, human_ensembl_to_locus_id, ortho_cmap = orthos
    # Human ensmbl <-> worm locus
    print "Total number of loci in ortholist: {v}. \
    With translation to ESNGxx: {z}.".format(
        v=len(ortho_locus_to_ensmbl),
        z=len([x for x in ortho_locus_to_ensmbl if len(ortho_locus_to_ensmbl[x])]))
    locus_to_ensmbl_tr = translation(
        'Loci', 'Ensmbl', ortho_locus_to_ensmbl)
    # PUM2 targets.
    pum2_in_ensmbl = pum2_targets_as_list_of_ensmbl_sets()
    print "PUM2 complexTargetSet:"
    pum2_ensembl_ts = complexTargetSet(pum2_in_ensmbl, 'Ensmbl', name='PUM2')
    pum2_ensembl_ts.stats()
    pum2_ensembl_ts.remove_nontranslatable_cgenes(
        ortho_cmap.get_reverse_translator())
    return pum2_ensembl_ts


if __name__ == '__main__':
    if not os.path.exists('output/'): os.system('mkdir output/')
    orthos = read_ortholist_file('lib/OrthoList_20130531.csv', max_orthologs=2)
    ortho_locus_to_ensmbl, human_ensembl_to_locus_id, ortho_cmap = orthos
    pum2_ensembl_ts = pum2_targets_as_ensembl(orthos)
    print "Genes: {0}".format(pum2_ensembl_ts.num())
#    print pum2_ensembl_ts.cgenes
    #pum2_ann = annotater(pandas.DataFrame(pum2), name='PUM2 targets')
    # FBF targets.
    fbf = fbfTargetSet(#input_filename='lib/ii_peaks.txt')
        input_filename='lib/ii_peaks.txt')
#        input_filename='lib/combined_filtered/old_fbf1_to_fbf2_n2.txt')
#        input_filename='lib/combined_filtered/old_fbf2.txt')
    fbf.fbf_cts.remove_nontranslatable_cgenes(ortho_cmap)
    fbf.fbf_cts.overlap_with_complexTargetSet(pum2_ensembl_ts, ortho_cmap)
    pum2_in_ensmbl = pum2_targets_as_list_of_ensmbl_sets()
    fbf_targs = fbf.annotate_fbf_targs_df(pum2_in_ensmbl, orthos)
    ann = annotater(
        pandas.read_csv('lib/OrthoList_20130531.csv', sep=','),
        name='Ortholist')
    ann.create_mapping('Locus ID')
    fbf_targs = ann.annotate_a_dataframe(fbf_targs, key='Locus ID')
    write_excel.write_excel(fbf_targs, tabs=[3, 5])
    fbf_targs.to_csv('output/fbf_targs_with_overlap.txt', index=False, sep='\t')



