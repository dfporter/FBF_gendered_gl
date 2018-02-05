import pandas
from language import *
from annotater import *
import write_excel
import os
import re
import sys
import numpy as np
import glob
from fbfTargetSet import fbfTargetSet
from human_gene_ensmbl_map import human_gene_ensmbl_map
from pum2_vs_fbf import *
import fbf_from_peaks

def trunc_locus(df):
    _fbf_targs = df.transcript_id
    fbf_targs = []
    for t in _fbf_targs:
        s = t.split('.')
        if len(s) > 1:
            s[1] = re.sub('[a-z]$', '', s[1])
            fbf_targs.append('.'.join(s[0:2]))
        else:
            fbf_targs.append(t)
    df['Locus ID'] = fbf_targs
    return df
            
if __name__ == '__main__':
    if not os.path.exists('output/'): os.system('mkdir output/')
    anndir_tr = 'combined_filtered/translatable/'
    if not os.path.exists(anndir_tr): os.system('mkdir ' + anndir_tr)
    anndir = 'combined_filtered/ortho_info_added/'
    if not os.path.exists(anndir): os.system('mkdir ' + anndir)
    orthos = read_ortholist_file(
        'lib/OrthoList_20130531.csv', max_orthologs=5)
    ortho_locus_to_ensmbl, human_ensembl_to_locus_id, ortho_cmap = orthos
    pum2_ensembl_ts = pum2_targets_as_ensembl(orthos)
    print("Genes: {v}".format(v=pum2_ensembl_ts.num()))
#    print pum2_ensembl_ts.cgenes
    #pum2_ann = annotater(pandas.DataFrame(pum2), name='PUM2 targets')
    # FBF targets.
#    fbf = fbfTargetSet(#input_filename='lib/ii_peaks.txt')
#        input_filename='lib/ii_peaks.txt')
#        input_filename='lib/combined_filtered/old_fbf1_to_fbf2_n2.txt')
#        input_filename='lib/combined_filtered/old_fbf2.txt')
#    fbf.fbf_cts.remove_nontranslatable_cgenes(ortho_cmap)
#    fbf.fbf_cts.overlap_with_complexTargetSet(pum2_ensembl_ts, ortho_cmap)
    pum2_in_ensmbl = pum2_targets_as_list_of_ensmbl_sets()
#    fbf_targs = fbf.annotate_fbf_targs_df(pum2_in_ensmbl, orthos)
    ann = annotater(
        pandas.read_csv('lib/OrthoList_20130531.csv', sep=','),
        name='Ortholist')
    ann.create_mapping('Locus ID')
#    fbf_targs = ann.annotate_a_dataframe(fbf_targs, key='Locus ID')
    for fname in glob.glob('lib/combined_filtered/*txt'):
        pk = fbfTargetSet(
            input_filename=fname)
        #fbf_targs = pk.annotate_fbf_targs_df(pum2_in_ensmbl, orthos)
        df = pandas.read_csv(fname, sep='\t', index_col=False)
        df = trunc_locus(df)
        fbf_targs = ann.annotate_a_dataframe(df, key='Locus ID')
        print(fbf_targs)
        f = fbfTargetSet(fname)
        f.fbf_cts.remove_nontranslatable_cgenes(ortho_cmap)
        f.fbf_cts.overlap_with_complexTargetSet(pum2_ensembl_ts, ortho_cmap)
        fbf_targs_tr = f.annotate_fbf_targs_df(pum2_in_ensmbl, orthos)
        fbf_targs_tr.to_csv(
            anndir_tr + os.path.basename(fname), sep='\t', index=False)
        ann_tr = annotater(
            pandas.read_csv(
                'combined_filtered/translatable/' + os.path.basename(fname),
                sep='\t', index_col=False), name='tr')
        ann_tr.create_mapping('Locus ID')
        fbf_targs = ann_tr.annotate_a_dataframe(fbf_targs, key='Locus ID')
        fbf_targs.to_csv(
            anndir + os.path.basename(fname), sep='\t', index=False)
        
