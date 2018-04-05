import functools, pandas

class overlapDictInfo():
    
    def __init__(self):
        pass
    
    @staticmethod
    def union_of_column_of_sets(df, col):
        try:
            union = functools.reduce(lambda x,y: set.union(x, y), df[col])
        except:
            if type(df[col].tolist()[0]) == type(''):
                _col = [set([x]) for x in df[col].tolist()]
            else:
                _col = [set(x) for x in df[col].tolist()]
            union = functools.reduce(lambda x,y: set.union(x, y), _col)
        return union
    
    @staticmethod
    def to_type(x, y):
        if x and not y:
            return 'FBF only'
        elif (not x) and y:
            return 'PUM only'
        elif x and y:
            return 'Shared'
        elif (not x) and not y:
            return 'Neither'
        
    def add_info(self, name, df):
        df['has FBF targets'] = [len(x)>0 for x in df['FBF targets (Worm public name)']]
        df['has PUM targets'] = [len(x)>0 for x in df['PUM targets (human gene symbols)']]
        df['Shared'] = [self.to_type(x, y) for x,y in zip(df['has FBF targets'], df['has PUM targets'])]
        
        results = dict(df['Shared'].value_counts())
        results['FBF translatable targets (homologs collapsed)'] = sum([1 for x in df['has FBF targets'] if x])
        results['PUM translatable targets (homologs collapsed)'] = sum([1 for x in df['has PUM targets'] if x])
        
        flat_fbf_targets_public_name = self.union_of_column_of_sets(
            df, 'FBF targets (Worm public name)')
        flat_pum_targets_gene_symbols = self.union_of_column_of_sets(
            df, 'PUM targets (human gene symbols)')
        
        def _d(x, y):
            if x:
                return y
            return set()

 #       df['fbf_if_shared'] = [_d(len(p), f) for f,p in zip(
 #           df['FBF targets (Worm public name)'], df['PUM targets (human gene symbols)'])]
        
 #       flat_fbf_targets_public_name_if_shared = self.union_of_column_of_sets(
 #           df, 'fbf_if_shared')

#        df['pum_if_shared'] = [_d(len(f), p) for f,p in zip(
#            df['FBF targets (Worm public name)'], df['PUM targets (human gene symbols)'])]
        
#        flat_pum_targets_gene_symbols_if_shared = self.union_of_column_of_sets(
#            df, 'pum_if_shared')
        
 #       results['# translatable FBF targets (Worm public name)'] = len(flat_fbf_targets_public_name)
 #       results['# translatable PUM targets (human gene symbols)'] = len(flat_pum_targets_gene_symbols)
 #       results['# translatable FBF targets shared (Worm public name)'
 #             ] = len(flat_fbf_targets_public_name_if_shared)
        
        results['% translatable FBF targets shared (homologs collapsed)'] = int(100 * (
            results['Shared']/results['FBF translatable targets (homologs collapsed)']))
        results['% translatable PUM targets shared (homologs collapsed)'] = int(100 * (
            results['Shared']/results['PUM translatable targets (homologs collapsed)']))
        results['% translatable FBF targets shared (not collapsing homologs)'] = int(100 * (
            len(flat_fbf_targets_public_name_if_shared)/len(flat_fbf_targets_public_name)))
        results['% translatable PUM targets shared (not collapsing homologs)'] = int(100 * (
            len(flat_pum_targets_gene_symbols_if_shared)/len(flat_pum_targets_gene_symbols)))
        
        return df, results
        