import pandas, importlib, os

import simple
import simple.targets

import fbfTargetHelper
importlib.reload(simple.targets)

class fbfTargets(simple.targets.targets):

    def read_fbf_targets(self):
        top_dir = '/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/'

        overlap_table = []

        all_fbf_targets = set()
        
        file_sheet_name_tups = [
            (top_dir + '/tables/Table S1 Peaks.xls', 'OO FBF (25°C)'),
            (top_dir + '/tables/Table S1 Peaks.xls', 'SP FBF (25°C)'),
            (top_dir + '/tables/Table S1 Peaks.xls', 'OO FBF-1 (20°C)'),  
            (top_dir + '/tables/Table S1 Peaks.xls', 'OO FBF-2 (20°C)'),
            (top_dir + '/tables/Table S6 Blocks.xls', 'Block I'),
            (top_dir + '/tables/Table S6 Blocks.xls', 'Block II'),
            (top_dir + '/tables/Table S6 Blocks.xls', 'Block III'),
            ]
        
        self.read_excel(file_sheet_name_tups)
        
        ##################
        # Worm Gene name <-> Locus ID translator.
        worm_translator = fbfTargetHelper.fbfTargetHelper().get_gene_name_locus_id_mapper(
            '/opt/lib/worm_txpt_id_gene_name_map.txt')

        # Set the transcript_id column with the ENSG.
        # (That's all this does.)
        for fname, sheetname in file_sheet_name_tups:

            df = self.target_lists_dfs[sheetname]
            
            if ('transcript_id' not in df.columns) and ('Transcript' in df.columns):
                df['transcript_id'] = df['Transcript']

            elif ('transcript_id' not in df.columns) and ('Gene name' in df.columns):
                print('Using Gene name...')
                df['transcript_id'] = [
                    worm_translator.translate(
                        input_id=name, return_str=True,
                        multiple_homologs_in_native_language_possible=False
                    ) for name in df['Gene name']
                ]

        self.all_targets = set()
        for sheetname in self.target_lists_dfs:
            self.all_targets |= set(self.target_lists_dfs[sheetname]['transcript_id'])
    def find_overlaps(self, otherTargets):
        
        overlap_table = []
        
        for fname, sheetname in file_sheet_name_tups:

            df = self.target_lists_dfs[sheetname]
            #all_fbf_targets |= set(df.transcript_id)

            print("Loaded file, finding overlaps.")
            overlap_table.append(find_and_write_overlap_with_pum2(
                pum2_ensembl_ts, orthos, sheetname, writer, fbf_df=df.copy()))

            # Top 500.
            if re.search('S1', fname):
                df.sort_index(by='Rank', inplace=True)
                
            overlap_table.append(find_and_write_overlap_with_pum2(
                pum2_ensembl_ts, orthos, sheetname + ' top 500', writer, fbf_df=top_500_genes(df)))
