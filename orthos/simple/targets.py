import pandas, os, re, warnings

class targets():
    
    def __init__(self):
        pass
    
    @staticmethod
    def instructed(text, config):
        if text in config and config[text]:
            return True
        return False
    
    def read_input(self, config, verbose=False):
        
        self.fname, self.language, self.species = (
            config['filename'], config['language'], config['species'])
        
        if verbose:
            print('---\ntargets(): Loading with configurations:')
            print(config)
            print("---")
            
        if self.instructed('read_csv', config):
            
            if self.instructed('read_every_line_as_gene', config):
                self.targets = [x.rstrip('\n') for x in open(self.fname, 'r')]
                
        if self.instructed('load_sheets_of_excel', config):
            
            self.sheetnames = [x.strip(' ') for x in config['sheetnames'].split(',')]
            self.targets = set()
            _dfs = {}
            for sheetname in self.sheetnames:
                df = pandas.read_excel(
                    self.fname, index_col=False, sheetname=sheetname)
                
                if self.instructed('subset', config):
                    print("Before running rule ", len(df.index))
                    
                    rule = {}  # This is a workaround for the weird namespace/exec handling.
                    exec("rule['f'] = "  + config['subset_function'])
                    df = df[rule['f'](df)]
    
                    print("After rule ", len(df.index))
                    
                self.targets |= set(df[config['column_name']].tolist())
        
        if self.instructed('capitalize', config) or self.instructed('lowercase', config):
            
            if self.instructed('capitalize', config):
                func = lambda x: x.upper()
            else:
                func = lambda x: x.lower()
            _new = set()
            
            for target in list(self.targets):
                try:
                    _new.add(func(target))
                except:
                    warnings.warn(
                        "Non-string data type {0} for target {1} in {2}. Skipping.".format(
                            type(target), target, self.fname))
                    
            self.targets = _new
        
        self.targets = set(self.targets) - set(['NULL', '', '-'])
        
    def read_excel(self, file_sheet_name_tups):

        self.target_lists_dfs = {}
        
        for fname, sheetname in file_sheet_name_tups:
            print("\n---\nLoading {0} - {1}".format(fname, sheetname))

            # Load file.
            if os.path.basename(fname).split('.')[-1] in ['xls', 'xlsx']:
                self.target_lists_dfs[sheetname] = pandas.read_excel(
                    fname, index_col=False, sheetname=sheetname)
    
    def make_cgenes(self, translator):
        
        print("Making cgenes attribute as a list of frozensets.")
        
        if translator is None:
            self.cgenes = [frozenset([target]) for target in self.targets]
            return
        
        tups = []
        for target in self.targets:
            translates_to = translator.translate(target)
            tups.append([
                translates_to,
                translator.translate(translates_to, reverse=True)
            ])
                    
        tups = translator.collapse_list_of_paired_sets(tups)
        self.cgenes = [tup[0] for tup in tups]
        self.cgenes_translations = [tup[1] for tup in tups]
        
        print("Finished making cgenes.")