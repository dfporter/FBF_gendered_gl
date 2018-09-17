import glob, HTSeq, sys, pandas, collections, re, os, xlwt
import numpy as np

import excelHelper, __init__


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


def in_a_not_in_b(a, b, min_height=50):
	
    robust_a = a[a['height']>min_height].copy()
    #print a['height'].value_counts()
    print("{0}/{1} peaks above {2} height...".format(
                                                     len(robust_a.index), len(a.index), min_height))
    robust_a['not_in_b'] = [
                     (x not in set(b.gene_name)) for x in robust_a.gene_name]
    a_only = robust_a[robust_a['not_in_b']].copy()
    del a_only['not_in_b']
    print("\t->{0}/{1} unique ({2}%)".format(
                                      len(a_only.index), len(robust_a.index),
                                      100*len(a_only.index)/len(robust_a.index)))
    return a_only.copy()
#dfs['Robust, unique FBF-1'] = in_a_not_in_b(
#dfs['old_fbf1_to_fbf2_n2'], dfs['old_fbf2'])
#dfs['Robust, unique FBF-2'] = in_a_not_in_b(
#   dfs['old_fbf2'], dfs['old_fbf1_to_fbf2_n2'])

def isin(obj, aset):
    if obj in aset:
        return 1
    return 0

def correlations_in_peak_regions(dfs):
    ga = {}
    for g in glob.glob('bedgraph_norm/*_+.wig'):

        print(g)

        guess = g.partition('_+.wig')[0]
        pos, neg = (guess + '_+.wig', guess + '_-.wig')
        
        ga[os.path.basename(guess)] = __init__.get_a_bedgraph(pos, neg)

        print("Read {}".format(os.path.basename(guess)))

    def max_in_region(_ga, _iv):
        max_val = 0
        for stepiv, val in _ga[_iv].steps():
            if val > max_val:
                max_val = val
        return max_val

    def vals_in_peak_regions(peak_sheet_names, ga_names):
        ivs = collections.defaultdict(list)
        for f in peak_sheet_names:

            dfs[f].sort_values(by=['Chrm', 'Nucleotide location left end of peak'], inplace=True, ascending=False)

            for (chrom, left, right, strand) in zip(
                dfs[f]['Chrm'], dfs[f]['Nucleotide location left end of peak'],
                dfs[f]['Nucleotide location right end of peak'], dfs[f]['Strand']):

                ivs[(chrom, strand)].append((left, right))

        overlapping_2 = collections.defaultdict(list)
        replacement = collections.defaultdict(list)
        for chrom_strand in ivs:
            for (left1, right1) in ivs[chrom_strand]:
                for (left2, right2) in ivs[chrom_strand]:
                    if (left1 <= left2 <= right1) or (
                        left2 <= left1 <= right2) or (
                        ):
                        overlapping_2[chrom_strand].append((left1, right1))
                        overlapping_2[chrom_strand].append((left2, right2))
                        replacement[chrom_strand].append(
                            (min([left1, left2]), max([right1, right2]))
                            )

        for chrom_strand in ivs:
            ivs[chrom_strand] = set(ivs[chrom_strand])

            ivs[chrom_strand] -= set(overlapping_2[chrom_strand])

            ivs[chrom_strand] = ivs[chrom_strand] | set(replacement[chrom_strand])

        print("Intervals used in comparison: {}".format(sum([len(x) for x in ivs.values()])))

        vals = []
        for chrom_strand in ivs:
            for (left, right) in ivs[chrom_strand]:
                _iv = HTSeq.GenomicInterval(chrom_strand[0], left, right, chrom_strand[1])
                vals.append([max_in_region(ga[x], _iv) for x in ga_names])

        print(len(vals))
        return vals

    fbfs_25C =[
        'exp_fbf1_sp_1', 'exp_fbf1_sp_2', 'exp_fbf1_sp_3',
        'exp_fbf2_sp_1', 'exp_fbf2_sp_2', 'exp_fbf2_sp_3',
        'exp_fbf1_oo_1', 'exp_fbf1_oo_2', 'exp_fbf1_oo_3',
        'exp_fbf2_oo_1', 'exp_fbf2_oo_2', 'exp_fbf2_oo_3',
        ]

    vals = vals_in_peak_regions(['sp_fbf1', 'sp_fbf2', 'oo_fbf1', 'oo_fbf2'], fbfs_25C)

    with open('tables/individual_replicate_binding_normalized_25C.txt', 'w') as f:
        f.write('\t'.join(fbfs_25C) + '\n')
        for row in vals:
            f.write('\t'.join([str(x) for x in row]) + '\n')


    fbfs_20C = [ 
        'exp_fbf1_CGGA', 'exp_fbf1_GGTT', 'exp_fbf1_TGGC', 
        'exp_fbf2_CGGA', 'exp_fbf2_GGTT', 'exp_fbf2_TGGC',
        ]

    vals = vals_in_peak_regions(['old_fbf1', 'old_fbf2'], fbfs_20C)

    with open('tables/individual_replicate_binding_normalized_20C.txt', 'w') as f:
        f.write('\t'.join(fbfs_20C) + '\n')

        for row in vals:
            f.write('\t'.join([str(x) for x in row]) + '\n')

    print("...Finished determining correlations in peak height in peak regions.")

if __name__ == '__main__':

    helper = excelHelper.excelHelper()

    dfs = helper.read_combined_filtered_directory_into_dict_of_dataframes()

    # Add RIP-chip data and gene descriptions.
    rip_targset = helper.get_ripchip_targets()

    for x in dfs.values():
        x['is_ripchip'] = [isin(y, rip_targset) for y in x.gene_name]
        x['Is a target by RIP-chip (Kershner et al.)?'] = x['is_ripchip']
        helper.name_cols(x)
        helper.add_gene_description(x)
#    normalize_to_rna_seq(x)

    _mk('tables/')

    ########
    # Table S1: Peaks.
    ########
    to_sheet = {
        'old_fbf1': u'OO FBF-1 (20' + u'\xb0' + u'C)',
        'old_fbf1_to_fbf2_n2': 'Discard',
        'old_fbf2': u'OO FBF-2 (20' + u'\xb0' + u'C)',
        'oo_fbf1': u'OO FBF-1 (25' + u'\xb0' + u'C)',
        'oo_fbf2': u'OO FBF-2 (25' + u'\xb0' + u'C)',
        'sp_fbf1': u'SP FBF-1 (25' + u'\xb0' + u'C)',
        'sp_fbf2': u'SP FBF-2 (25' + u'\xb0' + u'C)',
        'oo_both': u'OO FBF (25' + u'\xb0' + u'C)',
        'oo_fbf': u'OO FBF (25' + u'\xb0' + u'C)',
        'sp_both': u'SP FBF (25' + u'\xb0' + u'C)',
        'sp_fbf': u'SP FBF (25' + u'\xb0' + u'C)',
        'Robust, unique FBF-1': u'R.U. FBF-1 (20' + u'\xb0' + u'C)',
        'Robust, unique FBF-2': u'R.U. FBF-2 (20' + u'\xb0' + u'C)',
    }

    #if os.path.exists(
    #    'tables/Table S1 Peaks.xls'):
    #    os.system('rm tables/Table\ S1\ Peaks.xls')

    writer = pandas.ExcelWriter('tables/Table S1 Peaks.xls')

    for f in dfs:
        
        new_cols = [x for x in dfs[f].columns if x in helper.columns]
        
        if 'Is a target by RIP-chip (Kershner et al.)?' == new_cols[0]:  # Move to the end.
            del new_cols[0]
        new_cols.append('Is a target by RIP-chip (Kershner et al.)?')

        old_cols = [x for x in dfs[f].columns if x not in new_cols]

        dfs[f] = dfs[f][new_cols + old_cols]  # Rearrange.

        out_sheet = dfs[f].copy()

        if 'transcript_id' in out_sheet.columns:
            del out_sheet['transcript_id']

        out_sheet.to_excel(
            writer, sheet_name=to_sheet[f], index=False)

    writer.save()


    ########
    # Table S2: Stats.
    ########

    print("Generating metrics for peaks...")

    ########
    # S2: FBF-1 vs FBF-2 correlations.
    #correlations_in_peak_regions(dfs)

    import spearmanCorrelationsMaker
    sp = spearmanCorrelationsMaker.spearmanCorrelationsMaker()
    sp.load_counts_file(log_scale=False, style='including_25C_single_replicates')
    sp.load_peaks_csv_files()
    sp.subset_counts_to_targets()

    def correlations(columns):
        df = sp.counts_df.copy()
        df = df[columns]

        print("Making Spearman correlations from {} RNAs.".format(df.shape[0]))

        # Create correlation matrix.
        cor = df.corr('pearson')

        return cor

    cor = correlations([
        'SP FBF1_1', 'SP FBF1_2', 'SP FBF1_3',
        'SP FBF2_1', 'SP FBF2_2', 'SP FBF2_3',
        'OO FBF1_1', 'OO FBF1_2', 'OO FBF1_3',
        'OO FBF2_1', 'OO FBF2_2', 'OO FBF2_3',])

    sp_fbf1_vs_fbf2 = [
            cor.loc['SP FBF1_1', 'SP FBF2_1'],
            cor.loc['SP FBF1_1', 'SP FBF2_2'],
            cor.loc['SP FBF1_1', 'SP FBF2_3'],
            cor.loc['SP FBF1_2', 'SP FBF2_1'],
            cor.loc['SP FBF1_2', 'SP FBF2_2'],
            cor.loc['SP FBF1_2', 'SP FBF2_3'],
            cor.loc['SP FBF1_3', 'SP FBF2_1'],
            cor.loc['SP FBF1_3', 'SP FBF2_2'],
            cor.loc['SP FBF1_3', 'SP FBF2_3'],
        ]

    oo_fbf1_vs_fbf2 = [
            cor.loc['OO FBF1_1', 'OO FBF2_1'],
            cor.loc['OO FBF1_1', 'OO FBF2_2'],
            cor.loc['OO FBF1_1', 'OO FBF2_3'],
            cor.loc['OO FBF1_2', 'OO FBF2_1'],
            cor.loc['OO FBF1_2', 'OO FBF2_2'],
            cor.loc['OO FBF1_2', 'OO FBF2_3'],
            cor.loc['OO FBF1_3', 'OO FBF2_1'],
            cor.loc['OO FBF1_3', 'OO FBF2_2'],
            cor.loc['OO FBF1_3', 'OO FBF2_3'],
        ]

    ########
    # S2: Dataset sizes.

    reads = []
    for fname in glob.glob('bed_collapsed/*bed'):

        f = re.sub('\.bed', '', os.path.basename(fname))

        if re.match('$n2_oo.*', f):
            continue

        reads.append({
                     'Sample': f,
                     'Number of unique reads mapping uniquely': len(open(fname, 'r').readlines())})

    print("Finished counting dataset sizes (used the directory {})".format('bed_collapsed/'))

    reads = sorted(reads, key=lambda x: x['Sample'], reverse=False)
    reads = pandas.DataFrame(reads)


    if len(reads.index) == 0:
        print("Found no bed_collapsed/*bed files, and will not output read stats as a result.")
        sys.exit()

    print(reads)

    ########
    # S2: FBE, total RNA number, biotype

    def count_coding(x):

        if x == 'protein_coding':
            return 1
        else:
            return 0

    def count_fbe(x):

        try:
            if x == 1: return 1
            else: return 0
        except:
            print("error: {0}".format(x))
            return 0


    res = []#collections.defaultdict(dict)#pandas.Dataframe()

    for f in dfs:
        top = dfs[f][dfs[f]['Rank']<501].copy()
        res.append({
            'Name': to_sheet[f], 'Total # peaks': len(dfs[f].index),
            'Total # target RNAs': len(set(dfs[f]['Gene name'].tolist())),
            '% peaks with FBE': int(np.mean([count_fbe(x) for x in dfs[f]['Has UGUNNNAU (FBE)?'].tolist()]) * 100.),
            '% protein coding': int(np.mean([count_coding(x) for x in dfs[f]['Biotype'].tolist()]) * 100.),
            '% peaks with FBE (top 500)': int(np.mean(top['Has UGUNNNAU (FBE)?'].tolist()) * 100.),
            '% protein coding (top 500)': int(np.mean([count_coding(x) for x in dfs[f]['Biotype'].tolist()]) * 100.),
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

    #reads = reads[[(x in index_order) for x in reads['Sample'].tolist()]]

    #reads['Sample'] = [index_to_exp[x] for x in reads['Sample']]
    
    print("Outputing:\n{}".format(reads))
    writer = pandas.ExcelWriter('tables/Table S2 Peak stats.xls')

    res.to_excel(
            writer, sheet_name="Peak stats", index=False)
    
    reads = reads[['Sample', 'Number of unique reads mapping uniquely']]
    reads.to_excel(writer, sheet_name='Unique read counts', index=False)

    cor.to_excel(writer, sheet_name='Corr FBF-1 vs FBF-2')
    # Correlations in peak regions:
    #corrs = pandas.read_csv(
    #    'tables/correlations_fbf1_and_fbf2.txt', sep='\t',
    #    index_col=False)
    #corrs.to_excel(writer, sheet_name='Corr FBF-1 v -2', index=False)
    writer.save()

