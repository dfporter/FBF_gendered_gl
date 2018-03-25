import pandas, re
import numpy as np

import countsColumnsNaming
import importlib

rl = importlib.reload
rl(countsColumnsNaming)
to_len={
'exp_fbf_oo_rt9_and_1_counts.txt': 2108384,
'exp_fbf_oo_rt11_and_6_counts.txt': 1084505,
'exp_fbf_oo_rt2_and_13_counts.txt': 1104292,

    
'fbf2_oo_lane2_rt2_counts.txt': 593974,
'fbf1_oo_lane2_rt6_counts.txt': 813150,
'fbf2_oo_lane2_rt13_counts.txt': 510318,
'fbf1_oo_lane2_rt9_counts.txt': 2087441,
'exp_fbf2_oo_lane2_rt11_counts.txt': 271355,
'exp_fbf1_oo_lane2_rt6_counts.txt': 813150,
'exp_fbf1_oo_lane2_rt1_counts.txt': 20943,
'fbf2_oo_lane2_rt11_counts.txt': 271355,
'fbf1_oo_lane2_rt1_counts.txt': 20943,


'control_oo_counts.txt': 669244,
'control_sp_lane3_rt15_counts.txt': 72787,
'control_oo_lane1_rt15_counts.txt': 337373,
'control_oo_lane1_rt16_counts.txt': 165648,
'control_oo_lane1_rt3_counts.txt': 166223,


'exp_fbf1_CGGA_counts.txt': 1341859,
'exp_fbf2_TGGC_counts.txt': 1932659,
'exp_fbf2_GGTT_counts.txt': 348788,
'exp_fbf2_CGGA_counts.txt': 2505596,
'exp_fbf1_GGTT_counts.txt': 3558177,
'exp_fbf1_TGGC_counts.txt': 4277008,

'fbf2_sp_lane1_rt13_counts.txt': 1298522,
'exp_fbf2_sp_lane1_rt14_counts.txt': 999598,
'exp_fbf2_sp_lane1_rt2_counts.txt': 762755,
'fbf1_sp_lane1_rt6_counts.txt': 597193,
'fbf1_sp_lane1_rt1_counts.txt': 1128000,
'exp_fbf1_sp_lane1_rt1_counts.txt': 1128000,
'fbf1_sp_lane1_rt9_counts.txt': 1386301,
'exp_fbf1_sp_lane1_rt6_counts.txt': 597193,

'exp_fbf_sp_rt1_and_14_counts.txt': 2127598,
'exp_fbf_sp_rt6_and_9_counts.txt': 1983494,
'exp_fbf_sp_rt2_and_13_counts.txt': 2061277,
    
'exp_fbf1_oo_lane2_rt9_counts.txt': 2087441,


'control_sp_lane3_rt3_counts.txt': 139152,
'control_sp_counts.txt': 390956,
'control_sp_lane3_rt16_counts.txt': 179017,

'control_n2_counts.txt': 268768,
'control_n2_matching_fbf2_CCGG_counts.txt': 13314,
'control_n2_matching_fbf2_TTGT_counts.txt': 40194,
'control_n2_matching_fbf2_GGCA_counts.txt': 215260,

'exp_fbf2_sp_lane1_rt13_counts.txt': 1298522,
'exp_fbf1_sp_lane1_rt9_counts.txt': 1386301,
'fbf2_sp_lane1_rt2_counts.txt': 762755,
'fbf2_sp_lane1_rt14_counts.txt': 999598,



'n2_oo_lane1_rt16_counts.txt': 165648,
'n2_oo_lane1_rt15_counts.txt': 337373,
'n2_oo_lane1_rt3_counts.txt': 166223,

'exp_fbf2_oo_lane2_rt13_counts.txt': 510318,
'exp_fbf2_oo_lane2_rt2_counts.txt': 593974,
'SP RPKM': 1e6,
'OO RPKM': 1e6}


class countsFileIO(countsColumnsNaming.countsColumnsNaming):
    
    def load_counts_file(self, fname='combined_counts.txt', log_scale=True,
                        only_combined_datasets=True, style='heatmap'):
        
        self.to_len = to_len
        
        self.counts_df = pandas.read_csv(
            fname, sep='\t', index_col=False)#.head(1000)
        
        # We only keep the LT FBF replicates, the combined SP/OO FBF
        # replicates, and the uncombined LT and HT controls for normalization.
        
        if style == 'heatmap':

            self.counts_df = self.counts_df[[x for x in self.counts_df.columns if (
                not(self.name_is_one_of_11_ht_reps(x)) or x=='gene' or re.search('control', x))]]
            
            is_rrna = [x for x in self.counts_df.gene if re.match('rrn-\d+.*', x)]
            print("RNA genes: {0}\n{1}\n".format(len(is_rrna), is_rrna))

            self.counts_df.set_index('gene', inplace=True)
            self.counts_df.drop(is_rrna, inplace=True)
            self.counts_df = self.to_reads_per_mil(self.counts_df)
            self.shorten_names_and_subset_columns()
            self.counts_df = self.subtract_controls(self.counts_df)

            if log_scale:
                self.counts_df = self.scale_columns(self.counts_df)

            self.rm_controls()

        return self.counts_df
    
    def to_reads_per_mil(self, df):
        for col in [x for x in df.columns if (x != 'gene')]:
            df[col] = [1e6 * x/self.to_len[col] for x in df[col].tolist()]
        return df

    def scale_columns(self, df):
#        def fix(x):
#            return re.sub('all_bed_collapsed/', '',
#                        re.sub('\.bed', '_counts.txt', x))
        def _log(x):
            if (x>=0.25): return np.log2(x)
            else: return -2
        
        def ceiling(x):
            if x < 1E3:
                return x
            else:
                return 1E3
        
        for col in df.columns:
            if col == 'gene': continue
            df[col] = [_log(y) for y in df[col].tolist()]
            #df[col].replace([np.inf, -np.inf, float(np.log2(0))], 0)
        return df

    def subtract_controls(self, df):
        
        def subtract(exp_cols, control_col, df):
            for k in exp_cols:
                df[k] = [a - b for a, b in \
                         zip(df[k].tolist(), df[control_col].tolist())]
            return df
        
        tups = zip(
            df['c_sp_1'].tolist(), df['c_oo_1'].tolist(),
            df['c_n2_1'].tolist())
        
        df['ave_neg'] = [np.sum(t)/3. for t in tups]
        
        df = subtract(('SP FBF_1', 'SP FBF_2', 'SP FBF_3'), 'ave_neg', df)
        df = subtract(('OO FBF_1', 'OO FBF_2', 'OO FBF_3'), 'ave_neg', df)
        df = subtract(('LT FBF1_1', 'LT FBF1_2', 'LT FBF1_3'), 'ave_neg', df)
        df = subtract(('LT FBF2_1', 'LT FBF2_2', 'LT FBF2_3'), 'ave_neg', df)
        
        return df

    def normalize_columns(self, df):
        cols_to_norm = [x for x in df.columns if x != 'gene']
        df[cols_to_norm] = df[cols_to_norm].apply(
            lambda x: (x - x.mean()) / (x.max() - x.min()))
        return df

