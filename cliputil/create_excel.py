import xlwt
import glob
import sys
import pandas
import collections
import re
import os
import numpy as np

def get_ripchip_targets(get_sam_rank=False):
    ripchip_filename = '/opt/lib/ripchip/sd01.txt'
    ripchip = pandas.read_csv(
        ripchip_filename, sep='\t')
    top_1350 = []
    by_sam_rank = {}
    for index, row in ripchip.iterrows():
        if row['Gene public name'] not in top_1350:
            top_1350.append(row['Gene public name'])
            by_sam_rank[row['Gene public name']] = row['SAM rank']
        if len(list(set(top_1350))) >= 1350:
            break
    #top_1350_rows = ripchip[0:1349]  # There is a row zero.
    #rip_targets = list(top_1350['Gene public name'])
    print "Top targets: {t}.\nLast targets: {b}.".format(
        t=str(top_1350[0:2]), b=str(top_1350[-2:]))
    if get_sam_rank:
        return set(top_1350), by_sam_rank
    return set(top_1350)

'''
chrm	gene_name   left    	control_sp_lane3_rt15	control_sp_lane3_rt16	exp_fbf2_sp_lane1_rt14	exp_fbf2_sp_lane1_rt13	exp_fbf2_sp_lane1_rt2	control_sp_lane3_rt3
exp_reads	control_reads	ratio	unnorm_control_sp_lane3_rt15
unnorm_control_sp_lane3_rt16	unnorm_exp_fbf2_sp_lane1_rt14	unnorm_exp_fbf2_sp_lane1_rt13	unnorm_control_sp_lane3_rt3	unnorm_exp_fbf2_sp_lane1_rt2
unnorm_exp_reads    unnorm_control_reads	unnorm_ratio	biotype	seq	has_fbe	number_of_fbes_fbe	minus_one_c	minus_two_c	minus_one_or_two_c
neg_ip_local_norm	right	strand	
transcript_name	

transcript_id
iv height max_bin
'''
def name_cols(df):
    df['Rank'] = range(1, len(df.index)+1)
    df['Chrm'] = df['chrm']
    df['Nucleotide location left end of peak'] = df['left']
    df['Nucleotide location right end of peak'] = df['right']
    df['Gene name'] = df['gene_name']
    df['Peak height (reads/million)'] = df['height']
    df['Strand'] = df['strand']
    df['Transcript'] = df['transcript_name']
    df['Biotype'] = df['biotype']
    #df['Location'] = df['location']
    df['Has UGUNNNAU (FBE)?'] = df['has_fbe']
    df['Has -1C FBE?'] = df['minus_one_c']
    df['Has -2C FBE?'] = df['minus_two_c']
    df['Has -1 or -2 C FBE?'] = df['minus_one_or_two_c']
    #df['FBF-1 max coverage in peak (reads/million)'] = df['fbf1_reads']
    #df['FBF-2 max coverage in peak (reads/million)'] = df['fbf2_reads']
    #df['N2 control for FBF-1 max coverage in peak (reads/million)'] = df['fbf1_n2_reads']
    #df['N2 control for FBF-2 max coverage in peak (reads/million)'] = df['fbf2_n2_reads']
    #df['RNA-seq modencode (Acc. 4594) max coverage in peak (reads/million)'] = df['rna_seq_modencode']
    #df['RNA-seq Ortiz et al., oogenic germlines max coverage in peak (reads/million)'] = df['rna_seq_oo']
    #df['RNA-seq Ortiz et al., oogenic germlines (RPKM)'] = df['RNA abundance']
    df['Sequence under peak'] = df['seq']
    for x in [y for y in df.columns if re.search('_lane', y)]:
        del df[x]
#    df.columns = [re.sub('_lane.*', '_rep', x) for x in df.columns]
    df.columns = [re.sub('ratio', 'enrichment_over_control', x) for x in df.columns]
    for x in ['iv', 'max_bin']:
        if x in df.columns: del df[x]
    return df


def _mk(t):
    if os.path.isfile(t):
        t = os.path.dirname(t)
    if not os.path.exists(t):
        os.system('mkdir ' + t)


def normalize_to_rna_seq(df):
    ortiz_min = min([x for x in df.rna_seq_oo if x > 0])
    rpkm_min = min([x for x in df['RNA abundance'].tolist() if x > 0])
    modencode_min = min([x for x in df.rna_seq_modencode if x > 0])
    ortiz_norm = [float(t[0])/float(max([t[1], ortiz_min])) for t in zip(df.height, df.rna_seq_oo)]
    ortiz_rpkm_norm = [float(t[0])/float(max([t[1], rpkm_min])) for t in zip(df.height, df['RNA abundance'].tolist())]
    modencode_norm = [float(t[0])/float(max([t[1], modencode_min])) for t in zip(df.height, df.rna_seq_modencode)]
    df['Peak height/RNA abundance (RPKM for oogenic germlines, Ortiz et al.)'] = ortiz_rpkm_norm
    df['Peak height/RNA abundance (Max RNA coverage for oogenic germlines, Ortiz et al.)'] = ortiz_norm
    df['Peak height/RNA abundance (Max RNA coverage for whole worms, modencode (Acc. 4594))'] = modencode_norm

dfs = {}
for fname in glob.glob('combined_filtered/*'):
    print '-->', fname
    f = re.sub('\.txt', '', os.path.basename(fname))
    try:
        dfs[f] = pandas.read_csv(fname, sep='\t')
    except:
        print 'Not reading in {0} due to being unable to read'.format(fname)
        continue
    dfs[f]['height'] = dfs[f]['exp_reads']
    dfs[f].sort(columns=['height'], inplace=True, ascending=False)

reads = []
for fname in glob.glob('all_bed_collapsed/*bed'):
    f = re.sub('\.bed', '', os.path.basename(fname))
    if re.match('$n2_oo.*', f): continue
    reads.append({
'Sample': f,
'Number of unique reads mapping uniquely': len(open(fname, 'r').readlines())})
reads = sorted(reads, key=lambda x: x['Sample'], reverse=False)
reads = pandas.DataFrame(reads)

print reads
targset = get_ripchip_targets()

def isin(z, aset):
    if z in aset: return 1
    else: return 0

for x in dfs.values():
    x['is_ripchip'] = [isin(y, targset) for y in x.gene_name]
    x['Is a target by RIP-chip (Kershner et al.)?'] = x['is_ripchip']
    name_cols(x)
#    normalize_to_rna_seq(x)

columns = [
    'Rank', 'Chrm', 'Nucleotide location left end of peak',
    'Nucleotide location right end of peak', 'Gene name',
    'Peak height (reads/million)',
    'Strand', 'Transcript', 'Biotype', 'Location',
    'Has UGUNNNAU (FBE)?', 'Has -1C FBE?',
    'Has -1 or -2 C FBE?', 'Has -2C FBE?',
    'FBF-1 max coverage in peak (reads/million)',
    'FBF-2 max coverage in peak (reads/million)',
    'N2 control for FBF-1 max coverage in peak (reads/million)',
    'N2 control for FBF-2 max coverage in peak (reads/million)',
    'RNA-seq modencode (Acc. 4594) max coverage in peak (reads/million)',
    'RNA-seq Ortiz et al., oogenic germlines max coverage in peak (reads/million)',
    'RNA-seq Ortiz et al., oogenic germlines (RPKM)',
    'Peak height/RNA abundance (RPKM for oogenic germlines, Ortiz et al.)',
    'Peak height/RNA abundance (Max RNA coverage for oogenic germlines, Ortiz et al.)',
    'Peak height/RNA abundance (Max RNA coverage for whole worms, modencode (Acc. 4594))',
    'Is a target by RIP-chip (Kershner et al.)?',
    'Sequence under peak']
raw_columns = [
    'chrm', 'left', 'right', 'gene_name',
    'strand', 'ratio', 'biotype', 'location',
    'height', 'seq', 'control', 'fbf1_n2',
    'fbf2_reads', 'fbf1_n2_reads', 'exp_fbf1_GGTT_reads', 'iv',
    'control_reads', 'minus_one_or_two_c', 'fbf1_rep_2', 'fbf1_rep_1',
    'exp_fbf1_TGGC_reads', 'fbf1_reads', 'fbf2_rep_3', 'exp_fbf1_CGGA_reads',
    'fbf2_rep_1', 'exp', 'n2_from_fbf2_reads', 'minus_two_c',  'transcript_id',
    'fbf2_n2_reads', 'exp_fbf1_CGGA', 'number_of_fbes_fbe', 'fbf1_rep_3',
    'has_fbe', 'exp_fbf1_TGGC','minus_one_c', 'exp_fbf1_GGTT', 'Rank',
    'fbf2_rep_2',
    'n2_from_fbf1_reads', 'transcript_name', 'fbf2_n2', 'RNA abundance',
    'rna_seq_modencode', 'rna_seq_oo'
    ]

_mk('tables/')
def in_a_not_in_b(a, b, min_height=50):
    robust_a = a[a['height']>min_height].copy()
    #print a['height'].value_counts()
    print "{0}/{1} peaks above {2} height...".format(
        len(robust_a.index), len(a.index), min_height)
    robust_a['not_in_b'] = [
        (x not in set(b.gene_name)) for x in robust_a.gene_name]
    a_only = robust_a[robust_a['not_in_b']].copy()
    del a_only['not_in_b']
    print "\t->{0}/{1} unique ({2}%)".format(
        len(a_only.index), len(robust_a.index),
        100*len(a_only.index)/len(robust_a.index))
    return a_only.copy()
dfs['Robust, unique FBF-1'] = in_a_not_in_b(
    dfs['old_fbf1_to_fbf2_n2'], dfs['old_fbf2'])
dfs['Robust, unique FBF-2'] = in_a_not_in_b(
    dfs['old_fbf2'], dfs['old_fbf1_to_fbf2_n2'])

########
# Table S1: Peaks.
########
to_sheet = {
    'old_fbf1': 'Discard',
    'old_fbf1_to_fbf2_n2': u'oo FBF-1 (20' + u'\xb0' + u')',
    'old_fbf2': u'oo FBF-2 (20' + u'\xb0' + u'C)',
    'oo_fbf1': u'oo FBF-1 (25' + u'\xb0' + u'C)',
    'oo_fbf2': u'oo FBF-2 (25' + u'\xb0' + u'C)',
    'sp_fbf1': u'sp FBF-1 (25' + u'\xb0' + u'C)',
    'sp_fbf2': u'sp FBF-2 (25' + u'\xb0' + u'C)',
    'oo_both': u'oo FBF (25' + u'\xb0' + u'C)',
    'sp_both': u'sp FBF (25' + u'\xb0' + u'C)',
    'Robust, unique FBF-1': u'R.U. FBF-1 (20' + u'\xb0' + u'C)',
    'Robust, unique FBF-2': u'R.U. FBF-2 (20' + u'\xb0' + u'C)',
}
writer = pandas.ExcelWriter('tables/Table S1 Peaks.xls')
# Tab 1.
for f in dfs:
    if f == 'old_fbf1': continue
    new_cols = [x for x in dfs[f].columns if x in columns]
    if 'Is a target by RIP-chip (Kershner et al.)?' == new_cols[0]:
        del new_cols[0]
    new_cols.append('Is a target by RIP-chip (Kershner et al.)?')
    old_cols = [x for x in dfs[f].columns if x not in new_cols]
    order_columns = new_cols + old_cols
    dfs[f] = dfs[f][order_columns]
    # degree symbol = u'\xb0'
    dfs[f].to_excel(
        writer, sheet_name=to_sheet[f], index=False)
writer.save()
del dfs['Robust, unique FBF-1']
del dfs['Robust, unique FBF-2']


########
# Table S2: Stats.
########
import collections
def count_coding(x):
    if x == 'protein_coding': return 1
    else: return 0
def count_fbe(x):
    try:
        if x == 1: return 1
        else: return 0
    except:
        print "error: {0}".format(x)
        return 0
res = []#collections.defaultdict(dict)#pandas.Dataframe()
for f in dfs:
    top = dfs[f][dfs[f]['Rank']<201].copy()
    res.append({
        'Name': to_sheet[f], 'Total # peaks': len(dfs[f].index),
        'Total # target RNAs': len(set(dfs[f]['gene_name'].tolist())),
        '% peaks with FBE': int(np.mean([count_fbe(x) for x in dfs[f]['has_fbe'].tolist()]) * 100.),
        '% protein coding': int(np.mean([count_coding(x) for x in dfs[f]['biotype'].tolist()]) * 100.),
        '% peaks with FBE (top 200)': int(np.mean(top['has_fbe'].tolist()) * 100.),
        '% protein coding (top 200)': int(np.mean([count_coding(x) for x in dfs[f]['biotype'].tolist()]) * 100.),
        })
res = pandas.DataFrame(res)
first = ['Name', 'Total # peaks', 'Total # target RNAs']
order = first + [x for x in res.columns if x not in first]
res = res[order]
index_order = [
'exp_fbf1_CGGA',
'exp_fbf1_GGTT',
'exp_fbf1_TGGC',
'exp_fbf2_CGGA',
'exp_fbf2_GGTT',
'exp_fbf2_TGGC',
'control_n2_matching_fbf2_CCGG',
'control_n2_matching_fbf2_GGCA',
'control_n2_matching_fbf2_TTGT',
'exp_fbf1_oo_lane2_rt1',
'exp_fbf1_oo_lane2_rt6',
'exp_fbf1_oo_lane2_rt9',
'exp_fbf1_sp_lane1_rt1',
'exp_fbf1_sp_lane1_rt6',
'exp_fbf1_sp_lane1_rt9',
'control_sp_lane3_rt15',
'control_sp_lane3_rt16',
'control_sp_lane3_rt3',
'exp_fbf2_oo_lane2_rt11',
'exp_fbf2_oo_lane2_rt13',
'exp_fbf2_oo_lane2_rt2',
'exp_fbf2_sp_lane1_rt13',
'exp_fbf2_sp_lane1_rt14',
'exp_fbf2_sp_lane1_rt2',
'control_oo_lane1_rt15',
'control_oo_lane1_rt16',
'control_oo_lane1_rt3']
index_to_exp = {
'exp_fbf1_CGGA': u'oo FBF-1 (20' + u'\xb0' + u') barcode CGGA',
'exp_fbf1_GGTT': u'oo FBF-1 (20' + u'\xb0' + u') barcode GGTT',
'exp_fbf1_TGGC': u'oo FBF-1 (20' + u'\xb0' + u') barcode TGGC',
'exp_fbf2_CGGA':u'oo FBF-2 (20' + u'\xb0' + u') barcode CGGA',
'exp_fbf2_GGTT':u'oo FBF-2 (20' + u'\xb0' + u') barcode GGTT',
'exp_fbf2_TGGC':u'oo FBF-2 (20' + u'\xb0' + u') barcode TGGC',
'control_n2_matching_fbf2_CCGG':u'N2 control for oo FBF (20' + u'\xb0' + u') barcode CCGG',
'control_n2_matching_fbf2_GGCA':u'N2 control for oo FBF (20' + u'\xb0' + u') barcode GGCA',
'control_n2_matching_fbf2_TTGT':u'N2 control for oo FBF (20' + u'\xb0' + u') barcode TTGT',
'exp_fbf1_oo_lane2_rt1': u'oo FBF-1 (25' + u'\xb0' + u') rep 1',
'exp_fbf1_oo_lane2_rt6': u'oo FBF-1 (25' + u'\xb0' + u') rep 2',
'exp_fbf1_oo_lane2_rt9': u'oo FBF-1 (25' + u'\xb0' + u') rep 3',
'exp_fbf1_sp_lane1_rt1': u'sp FBF-1 (25' + u'\xb0' + u') rep 1',
'exp_fbf1_sp_lane1_rt6': u'sp FBF-1 (25' + u'\xb0' + u') rep 2',
'exp_fbf1_sp_lane1_rt9': u'sp FBF-1 (25' + u'\xb0' + u') rep 3',
'control_sp_lane3_rt15': u'control for sp FBF (25' + u'\xb0' + u') rep 1',
'control_sp_lane3_rt16': u'control for sp FBF (25' + u'\xb0' + u') rep 2',
'control_sp_lane3_rt3':  u'control for sp FBF (25' + u'\xb0' + u') rep 3',
'exp_fbf2_oo_lane2_rt11': u'oo FBF-2 (25' + u'\xb0' + u') rep 1',
'exp_fbf2_oo_lane2_rt13': u'oo FBF-2 (25' + u'\xb0' + u') rep 2',
'exp_fbf2_oo_lane2_rt2': u'oo FBF-2 (25' + u'\xb0' + u') rep 3',
'exp_fbf2_sp_lane1_rt13':u'sp FBF-2 (25' + u'\xb0' + u') rep 1',
'exp_fbf2_sp_lane1_rt14': u'sp FBF-2 (25' + u'\xb0' + u') rep 2',
'exp_fbf2_sp_lane1_rt2':u'sp FBF-2 (25' + u'\xb0' + u') rep 3',
'control_oo_lane1_rt15':u'control for oo FBF (25' + u'\xb0' + u') rep 1',
'control_oo_lane1_rt16': u'control for oo FBF (25' + u'\xb0' + u') rep 2',
'control_oo_lane1_rt3': u'control for oo FBF (25' + u'\xb0' + u') rep 3',
}
reads['keep'] = [(x in index_order) for x in reads['Sample'].tolist()]
reads = reads[reads['keep']]
del reads['keep']
reads['Sample'] = [index_to_exp[x] for x in reads['Sample']]
writer = pandas.ExcelWriter('tables/Table S2 Peak stats.xls')
res.to_excel(
        writer, sheet_name="Peak stats", index=False)
reads = reads[['Sample', 'Number of unique reads mapping uniquely']]
reads.to_excel(writer, sheet_name='Unique read counts', index=False)
# Correlations in peak regions:
corrs = pandas.read_csv(
    'tables/correlations_fbf1_and_fbf2.txt', sep='\t',
    index_col=False)
corrs.to_excel(writer, sheet_name='Corr FBF-1 v -2', index=False)
writer.save()

