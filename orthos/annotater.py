import collections
import scipy.stats
import pandas

class annotater(object):

    def __init__(self, df, keys=None, name="Unnamed"):
        assert(type(df) == type(pandas.DataFrame()))
        
        self.df = df
        self.name = name
        self.df_list = df.to_dict('records')
        self.by = collections.defaultdict(dict)
        
        if keys is not None:
            for k in keys:
        
                if (type(k) == type([])) or type(k) == type(set()):
                    k = frozenset(k)
                
                self.create_mapping(k)

    def create_mapping(self, key):
        
        if key not in self.df.columns:
            print("%s not found in columns for %s." % (
                key, self.name))
            print("\tColumns are %s" % str(self.df.columns))
            return
        
        self.by[key] = dict([(self.df_list[i][key], self.df_list[i]) \
                             for i in range(len(self.df_list))])

    def annotate(self, name, key='gene_name'):
        self.can_map(key)
        return self.by[key]

    def add_info(self, info_dict, key='gene_name'):
        self.can_map(key)
        df_list = df.to_dict('records')
        for k in self.by[key]:
            if k not in info_dict:
                info_dict[k] = {}
        for i, row in enumerate(df_list):
            k = df_list[i][key]
            df_list[i].update(info_dict[k])
        self.df = pandas.DataFrame(df_list)

    def can_map(self, key):
        if key not in self.by:
            print("%s: Creating mapping by %s if possible..." % (self.name, key))
            self.create_mapping(key)

    def annotate_a_dataframe(self, df, key='gene_name'):
        print("Annotating by %s..." % key)
        
        self.can_map(key)
        
        df_list = df.to_dict('records')
     
        for i, row in enumerate(df_list):

            row.update(dict([(k, '') for k in self.df.columns if k not in row]))
            row.update(self.by[key].get(row[key], {}))
                
        return pandas.DataFrame(df_list)
