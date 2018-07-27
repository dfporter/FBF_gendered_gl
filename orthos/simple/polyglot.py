import pandas, functools, collections, re, importlib, random

import orthoListHelper
import human_gene_ensmbl_map

importlib.reload(human_gene_ensmbl_map)
importlib.reload(orthoListHelper)
    
class polyglot():
    
    def __init__(self, df=None):
        self.df = df
    
    def add_a_column_of_targets(self, targ_set, language, col_name='targs'):
        new_col = []
        for _id in self.df[language].tolist():
            new_col.append(_id & targ_set)
        self.df[col_name] = new_col

    def add_a_column_of_random_targets(self, language, n_to_pick=100, col_name='random'):
        
        genes = self.df[language].tolist()
        n_genes = len(genes)
        
        picked_genes = set()
        while len(picked_genes) <= n_to_pick:
            sample = random.sample(range(0, n_genes), 1)[0]
            picked_genes |= genes[sample]
        
        new_col = []
        for _id in genes:
            new_col.append(_id & picked_genes)
            
        self.df[col_name] = new_col
        
    @staticmethod
    def make_polyglot_file(config):

        print("Making polyglot file.")
        
        # This is a map of worm loci <-> human ENSG ids of orthologs.
        orthos = orthoListHelper.orthoListHelper().read_ortholist_file(
            config.s['ortholist'], max_orthologs=10)

        # This is a map of human gene symbol <-> human ENSG.
        human_ensmbl_gene = human_gene_ensmbl_map.human_gene_ensmbl_map()
        _polyglot = []
        ensg_added = set()

        for n, (worm_names, ensg_set) in enumerate(orthos.transl.items(), start=1):
            
            #if n < 20:
            #    break

            gene_symbol = human_ensmbl_gene.translate_list(list(ensg_set), reverse=True)
            if '' in gene_symbol:
                gene_symbol.remove('')
            
            if len(gene_symbol):
                gene_symbol = functools.reduce(lambda x, y: set.union(x, y), gene_symbol)
            else:
                gene_symbol = set()

            _polyglot.append({
                'Worm public name': worm_names,
                'Human ENSG': ensg_set,
                'Human gene symbol': gene_symbol})

            if not(n % 1000):
                print('On ensg set {0}/{1}.'.format(n, len(list(orthos.transl.values()))))
                print(_polyglot[-1])

            ensg_added |= ensg_set 
            
        for n, (human_ensmbl, gene_symbol) in enumerate(human_ensmbl_gene.transl.items(), start=1):
            break
            #worm_names = orthos.translate_list(list(human_ensmbl), reverse=True)
            #worm_names -= set([''])
            
            #if len(worm_names):
            #    pass
                #worm_names = functools.reduce(lambda x, y: set.union(x, y), worm_name)
            #else:
            #    worm_names = set()
            if len(human_ensmbl) != 1:
                print("ERROR")
                print(human_enmbl)
                raise ValueError('asdfasdf')
                
            if human_ensmbl not in ensg_added:
                _polyglot.append({
                    'Worm public name': set(),
                    'Human ENSG': human_ensmbl,
                    'Human gene symbol': gene_symbol})
                
            if not(n % 1000):
                print('(rev) On ensg set {0}/{1}.'.format(n, len(list(human_ensmbl_gene.transl.values()))))
                print(_polyglot[-1])
                
        df = pandas.DataFrame(_polyglot)

        for col in df.columns:
            df[col] = [x - set(['']) for x in df[col]]

        print(df.head())

        return df
    
    @staticmethod
    def to_nums(tup):
        if len(tup[0]):
            if len(tup[1]):
                return [1, 1]
            else:
                return [1, 0]
        if len(tup[1]):
            return [0, 1]
        else:
            return [0, 0]