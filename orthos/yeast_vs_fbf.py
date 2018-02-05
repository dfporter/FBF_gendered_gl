import os
import sys
import re
import collections
import pandas
from read_yeast_worm_ortho_file import *
from get_yeast_target_lists import *
from language import *
from annotater import *
import fbf_from_peaks
import write_excel

# Yeast targets
all_yeast_targs = get_set_of_all_yeast_tagging_targets()
all_yeast_targs = set([frozenset([x]) for x in all_yeast_targs])
yeast_tagging_ts = complexTargetSet(
    all_yeast_targs, 'Yeast Ensembl', 'Yeast tagging targets')
print yeast_tagging_ts.info()
# FBF.
fbf_ts = fbf_from_peaks.fbf_from_peaks()  # TargetSet.
fbf_ts = complexTargetSet(
    set([frozenset([x]) for x in fbf_ts.genes]), 'Worm transcript', name='FBF')
print fbf_ts.info()
# Translation.
worm_to_yeast_df = pandas.read_csv('lib/worm_yeast_map.txt', sep=',')
wty = worm_to_yeast_df.to_dict('records')
worm_yeast = complexTranslation('Worm transcript', 'Yeast Ensembl', {})
for row in wty:
    if pandas.isnull(
        row['Saccharomyces cerevisiae Ensembl Gene ID']):
        continue
    s = row['Ensembl Transcript ID'].split('.')
    s[1] = re.sub('[a-z]$', '', s[1])
    locus = '.'.join(s[0:2])
    worm_yeast.add_mapping(
        set([locus]),
        set([row['Saccharomyces cerevisiae Ensembl Gene ID']]))
worm_yeast.refresh()
print worm_yeast.info()

fbf_ts.remove_nontranslatable_cgenes(worm_yeast)
yeast_tagging_ts.remove_nontranslatable_cgenes(
    worm_yeast.get_reverse_translator())
fbf_ts.overlap_with_complexTargetSet(yeast_tagging_ts, worm_yeast)
fbf_yeast_df = fbf_ts.to_dataframe()
fbf_yeast_df.to_csv('fbf_yeas_df', sep='\t', index=False)
fbf_yeast_df['Locus ID'] = fbf_yeast_df['a_names']
yeast_ens_to_yeast_name = dict(zip(
  worm_to_yeast_df['Saccharomyces cerevisiae Ensembl Gene ID'].tolist(),
  worm_to_yeast_df['Saccharomyces cerevisiae associated gene name'].tolist()
  ))

fbf_pum2 = pandas.read_csv('fbf_targs_with_overlap.txt', sep='\t')
fbf_yeast_df['Yeast Ensembl ID'] = fbf_yeast_df['translations']
import numpy as np
ann = annotater(fbf_pum2, keys=['Ensembl IDs'], name='FBF PUM2')
df = ann.annotate_a_dataframe(fbf_yeast_df, key='Locus ID')
df['Worm locus ID'] = df['Locus ID']
def tr_if_can(sgdid):
    if type(sgdid) != type(''):
        print "???"
        print sgdid
        return ''
    y = set([yeast_ens_to_yeast_name[x] for x in sgdid.split(', ')])
    y -= set([np.nan])
    #print y
    if len(y) == 0: return ''
    return ", ".join(y)

df['Yeast gene name'] = [
    tr_if_can(sgdid) for sgdid in df['Yeast Ensembl ID'].tolist()]
df['yeast shared'] = df['shared']
df.drop(['shared'], axis=1, inplace=True)

if not os.path.exists('output/'): os.system('mkdir output/')
df.to_csv('output/fbf_yeast_with_pum2_if_present', sep='\t', index=False)
write_excel.write_excel(df, tabs=[2, 3, 4, 5])

#ann = annotater(fbf_yeast_df, keys=['Ensembl IDs'], name='FBF PUM2')
#df = ann.annotate_a_dataframe(fbf_pum2, key='Locus ID')
#df.to_csv('fbf_pum2_with_yeast_if_present', sep='\t', index=False)


#print df
sys.exit()
import scipy.stats as sps
print """There are {i} Puf3/5p targets. {a} of these have a C. elegans homolog.
There are {x} FBF targets. The expected number of shared targets is {sh}.
There are {b} shared FBF and Puf3/5p targets.
There are {pm2} PUM2 targets.
There are {pm2fbf} targets shared between FBF and PUM2.
There expected number of shared targets between FBF, PUM2 and Puf3/5p is {ethree}.
There are {obsthree} shared targets between FBF, PUM2 and Puf3/5p.
Targets in all three are:\n{shared_list}.""".format(
    i=len(puf35_targets), a=len(puf35_targets_in_elegans),
    b=len(puf35_fbf_targets), x=len(fbf_targs),
    sh=len(puf35_targets_in_elegans) * float(len(fbf_targs))/20000.,
    pm2=len(pum2_targets),
    pm2fbf=len(conserved_pum2_fbf_targs_worm_locus_id),
    ethree=float(len(conserved_pum2_fbf_targs_worm_locus_id)) * float(len(puf35_targets_in_elegans))/20000.,
    obsthree=len(all_three),
    shared_list='\n'.join([
        '\t'.join([locusid, common_names[locusid],
                   str(mart[locusid]['Yeast Ensembl Gene ID']) if locusid in mart else 'NA',
                   str(mart[locusid]['Human Ensembl Gene ID']) if locusid in mart else 'NA',
                ]) for locusid in all_three])
)
