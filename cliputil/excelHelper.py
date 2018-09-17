import xlwt
import glob
import sys
import pandas
import collections
import re
import os
import numpy as np

class excelHelper():
    
    def __init__(self):
        self.columns = [
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
    
    @classmethod
    def get_ripchip_targets(cls, get_sam_rank=False):
        ripchip_filename = '/opt/lib/ripchip/sd01.txt'
        ripchip = pandas.read_csv(ripchip_filename, sep='\t')
        
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
        print("Top targets: {t}.\nLast targets: {b}.".format(
            t=str(top_1350[0:2]), b=str(top_1350[-2:])))
        
        if get_sam_rank:
            return set(top_1350), by_sam_rank
        
        return set(top_1350)
    
    @classmethod
    def name_cols(cls, df):
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
        df['P value vs negative IP, using a local normal distribution'] = df['neg_ip_local_norm']
        df['Enrichment over negative IP (reads/million)'] = df['ratio']
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

        df.columns = [re.sub('ratio', 'enrichment_over_control', x) for x in df.columns]
        for x in ['iv', 'max_bin', 'chrm', 'left', 'right', 'height', 'strand', 'transcript_name',
                 'biotype', 'has_fbe', 'minus_one_c', 'minus_two_c', 'minus_one_or_two_c',
                 'seq', 'gene_name', 'neg_ip_local_norm', 'ratio', 'enrichment_over_control']:
            if x in df.columns:
                del df[x]
                
        col_to_del_final = [  # Duplicate col names.
            #                 'transcript_id', 
            'transcript_name', # Duplicate col names.
                            '\Aexp_.*',
                            '\Acontrol_.*', 'unnorm_.*',
                  'minus_one_or_two_c', 'is_ripchip',  # Duplicate col names.
                            'number_of_fbes_fbe',
            
                  ]
        for pat in col_to_del_final:
            for x in [y for y in df.columns if re.search(pat, y)]:
                del df[x]

        return df
    

    
    @classmethod
    def read_combined_filtered_directory_into_dict_of_dataframes(cls):
        
        dfs = {}
        
        for fname in glob.glob('combined_filtered/*txt'):
            print('-->', fname)
            
            if re.search('block', fname):
                continue

            f = re.sub('\.txt', '', os.path.basename(fname))
            
            try:
                dfs[f] = pandas.read_csv(fname, sep='\t')
            except:
                print('Not reading in {0} due to being unable to read'.format(fname))
                continue
                
            dfs[f]['height'] = dfs[f]['exp_reads']
            dfs[f].sort_values(by=['height'], axis=0, inplace=True, ascending=False)
           
        return dfs
    
    @staticmethod
    def simplify_locus_id(locus_id):

        txpt = locus_id.split('.')
        if len(txpt) > 1:
            locus_id = '.'.join(txpt[:2])

        return re.sub('[a-z]$', '', locus_id)
    
    def read_mart_export(self):
        self.trxpt_to_desc = {}
        with open('/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/tables/mart_export.txt') as f:
            for li in f:
                s = li.split('\t')
                self.trxpt_to_desc[self.simplify_locus_id(s[0])] = s[1].rstrip('\n').split(' [')[0]

    def add_gene_description(self, df, use='transcript_id'):
        if not hasattr(self, 'trxpt_to_desc'):
            self.read_mart_export()
        
        if 'transcript_id' in df.columns:
            df['Description'] = [self.trxpt_to_desc.get(name, name) for name in df['transcript_id']]
            
        elif ('Gene name' in df.columns) or ('gene_name' in df.columns):
            gtf = pandas.read_csv('/opt/lib/gtf_with_names_column.txt', sep='\t')
            name_to_txpt = dict(zip(gtf['gene_name'].tolist(), gtf['transcript_id'].tolist()))
            
            if 'Gene name' in df.columns:
                df['transcript_id'] = [self.simplify_locus_id(name_to_txpt.get(x, '')) for x in df['Gene name'].tolist()]
            else:
                df['transcript_id'] = [self.simplify_locus_id(name_to_txpt.get(x, '')) for x in df['gene_name'].tolist()]
                
            df['Description'] = [self.trxpt_to_desc.get(name, name) for name in df['transcript_id']]  
            
        return df
    
    