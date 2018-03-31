

class ManyToManyTranslator(object):
    """
    self.transl (dict): {frozenset a} -> {frozenset b}
    self.reverse_transl (dict): {frozenset b} -> {frozenset a}
    """
    
    def __init__(self, input_name, output_name,
                 name='Unnamed ManyToManyTranslater object'):
        
        self.input_name = input_name
        self.output_name = output_name
        self.name = name
        #self.max_orthologs = max_orthologs
        self.transl = {}
        self.reverse_transl = {}

    def reverse(self):
        (self.input_name, self.output_name) = (self.output_name, self.input_name)
        self.transl = self.reverse_transl
        self.reverse_transl = self.transl
    
    def define_mappings_from_one_to_one_list(
            self, list_of_paired_ids):
        for (seta, setb) in list_of_paired_ids:
            self.transl[frozenset([seta])] = frozenset([setb])
            self.reverse_transl[frozenset([setb])] = frozenset([seta])

    def define_mappings_from_list_of_paired_sets(
            self, list_of_paired_sets,
            max_orthologs=5):
        
        for (seta, setb) in list_of_paired_sets:
            if (len(seta)>max_orthologs) or (len(setb)>max_orthologs):
                continue
            self.transl[frozenset(seta)] = frozenset(setb)
            self.reverse_transl[frozenset(setb)] = frozenset(seta)
    
    def combine_list_of_sets_based_on_homology(self, _list, reverse_transl=False):
        
        if reverse_transl:
            groupings = list(self.reverse_transl.keys())
        else:
            groupings = list(self.transl.keys())
            
        print('combine_list_of_sets_based_on_homology():')
        
        expanded = self.expand_list_of_sets_to_include_all_in_group(
            _list, groupings)
        
        print('Expanded input list of {0} items (=output length {1}).'.format(
            len(_list), len(expanded)))
        
        collapsed = list(set(expanded))#self.collapse_list_of_sets(expanded)
        print("Collapsed a list of sets to length {0}.".format(len(collapsed)))
        return collapsed
    
    @staticmethod
    def expand_list_of_sets_to_include_all_in_group(_t, groupings):
        expanded = [set() for x in _t]
        
        for n, _a in enumerate(_t):
            
            for _b in groupings:
                
                if _a & _b:
                    expanded[n] |= _b
        
        for n in range(len(expanded)):
            expanded[n] = frozenset(expanded[n])
            
        return expanded
    
    @staticmethod
    def collapse_list_of_sets(_t):
        print("Collapsing a list of sets of length {0}".format(len(_t)))
        collapsed = []
        
        def remove_a_row(_list):
            for n1, _a in enumerate(_t):
                for n2, _b in enumerate(_t[n1+1:]):
                    if _a & _b:
                        return (True, n1, n2)
            return (False, -1, -1)
        
        has_overlaps = True
        
        while has_overlaps:
            has_overlaps, n1, n2 = remove_a_row(_t)
            if has_overlaps:
                _t[n1] |= _t[n2]
                del _t[n2]
                
        print("Collapsed a list of sets to length {0}.".format(len(_t)))
        return _t
    
    @staticmethod
    def collapse_list_of_paired_sets(_t):
        collapsed = []
        
        def update_row(a, b, row):
            if a & row[0]:
                collapsed[n][0] |= a
                collapsed[n][1] |= b
                return True
            elif b & row[1]:
                collapsed[n][0] |= a
                collapsed[n][1] |= b
                return True
            return False
                
        for _a, _b in _t:
            
            found = False
            
            for n, _row in enumerate(collapsed):
                if (not found) and update_row(_a, _b, _row):
                    found = True
                if found:
                    break
            
            if not found:
                collapsed.append([_a, _b])
            
        return collapsed
    
    def translate_list(self, _list, **kwargs):
        if ('reverse' not in kwargs) or (not kwargs['reverse']):
            print("Translating list of length {0} from {1} to {2}".format(
                len(_list), self.input_name, self.output_name))
        else:
            print("Translating list of length {0} from {1} to {2}".format(
                len(_list), self.output_name, self.input_name))
        
        output = []
        for n, item in enumerate(_list, start=1):
            
            output.append(self.translate(item, **kwargs))
            if not (n % 200):
                print("Translated item {0} ({1}) to {2}".format(n, item, output[-1]))
            
        return output
    
    def translate(
        self, input_id, verbose=False, reverse=False,
        multiple_homologs_in_native_language_possible=True,
        return_str=False):
        
        # input_id is frozenset
        if type(input_id) == type(''):
            input_id = frozenset([input_id])
            
        if verbose:
            print('looking for ', input_id)
            
        translates_to = set()
        if not reverse:
            
            if multiple_homologs_in_native_language_possible:
                for k in self.transl:
                    if len(k & input_id) > 0:
                        translates_to |= self.transl[k]
            else:
                translates_to = self.transl[input_id]
            
        else:
            for k in self.reverse_transl:
                if len(k & input_id) > 0:
                    translates_to |= self.reverse_transl[k]
                    
        if return_str:
            _list = list(translates_to)
            
            if len(_list) == 1:
                return str(_list[0])
            
            elif len(_list) == 0:
                
                if len(list(input_id)) == 1:
                    return str(list(input_id)[0])
                else:
                    return str(input_id)
                
            else:
                return str(translates_to)
            
        return translates_to