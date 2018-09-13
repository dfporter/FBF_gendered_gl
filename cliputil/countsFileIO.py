import pandas, re
import numpy as np

import countsColumnsNaming
import importlib

rl = importlib.reload
rl(countsColumnsNaming)

nope = """
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
'OO RPKM': 1e6,

"control_oo_counts.txt": 393784,
"control_oo_1_counts.txt": 97070,
"control_oo_2_counts.txt": 200654,
"control_oo_3_counts.txt": 96060,
"control_sp_counts.txt": 359370,
"control_sp_1_counts.txt": 128443,
"control_sp_2_counts.txt": 65914,
"control_sp_3_counts.txt": 165013,
"exp_fbf1_oo_1_counts.txt": 13486,
"exp_fbf1_oo_2_counts.txt": 681860,
"exp_fbf1_oo_3_counts.txt": 1698217,
"exp_fbf1_sp_1_counts.txt": 852744,
"exp_fbf1_sp_2_counts.txt": 475910,
"exp_fbf1_sp_3_counts.txt": 1064502,
"exp_fbf2_oo_1_counts.txt": 426458,
"exp_fbf2_oo_2_counts.txt": 186319,
"exp_fbf2_oo_3_counts.txt": 354981,
"exp_fbf2_sp_1_counts.txt": 569792,
"exp_fbf2_sp_2_counts.txt": 1001336,
"exp_fbf2_sp_3_counts.txt": 740220,
"exp_fbf_oo_1_counts.txt": 439944,
"exp_fbf_oo_2_counts.txt": 1036841,
"exp_fbf_oo_3_counts.txt": 1884536,
"exp_fbf_sp_1_counts.txt": 1477246,
"exp_fbf_sp_2_counts.txt": 1592964,
"exp_fbf_sp_3_counts.txt": 1634294,
}
"""

class countsFileIO(countsColumnsNaming.countsColumnsNaming):
    
    def load_counts_file(self, fname='combined_counts.txt', log_scale=True,
                        only_combined_datasets=True, style='heatmap',
                        total_read_numbers='total_read_numbers.txt'):
        
        self.to_len = {}

        with open(total_read_numbers) as f:
            for li in f:
                (_fname, read_count) = li.rstrip('\n').split('\t')
                self.to_len[_fname.split('/')[-1]] = int(read_count)
        
        #fname = 'combined_counts.txt'
        
        print("Loading {} as combined counts file.".format(fname))
        self.counts_df = pandas.read_csv(
            fname, sep='\t', header=0, index_col=False)#.head(1000)

        if 'gene' not in self.counts_df.columns:
            self.counts_df['gene'] = self.counts_df.index
            print("Set gene column in counts_file {}".format(self.counts_df['gene'].head(2)))

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

            self.counts_df = self.counts_df.loc[[
                (x not in ['_no_feature', '_ambiguous']) for x in self.counts_df.index]]
        
        else:
            
            self.counts_df = self.counts_df[[x for x in self.counts_df.columns if (
                not(self.name_is_one_of_11_ht_reps(x)) or x=='gene' or re.search('control', x))]]

            is_rrna = [x for x in self.counts_df.gene if re.match('rrn-\d+.*', x)]
            print("RNA genes: {0}\n{1}\n".format(len(is_rrna), is_rrna))

            self.counts_df.set_index('gene', inplace=True)
            self.counts_df.drop(is_rrna, inplace=True)
            self.shorten_names_and_subset_columns()

            self.rm_controls()

            self.counts_df = self.counts_df.loc[[
                (x not in ['_no_feature', '_ambiguous']) for x in self.counts_df.index]]
            

        return self.counts_df
    
    def rm_controls(self):
        for k in ['c_sp_1', 'c_n2_1', 'c_oo_1', 'ave_neg']:
            if k in self.counts_df.columns:
                del self.counts_df[k]
        combined_controls = set(['control_sp_counts.txt', 'control_oo_counts.txt',
                             'control_n2_counts.txt'])
        to_del = set([x for x in self.counts_df.columns if (\
                    re.search('control.*matching.*', x) \
                   or re.search('n2_oo_lane', x) \
                    or re.search('\Ac_', x) \
                   or re.search('control', x) # All controls removed. 
                   )])
        to_del = to_del - combined_controls
        
        for x in set(to_del):
            del self.counts_df[x]

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
            df['control_sp'].tolist(), df['control_oo'].tolist(),
            #df['c_n2_1'].tolist()
            )
        
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

