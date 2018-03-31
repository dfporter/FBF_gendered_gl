import collections
import scipy.stats
import pandas



class polyglot(object):

    def __init__(self):
        self.transl = {}

    def add_translation(self, trans):
        from_to = (trans.input_name, trans.output_name)
        if from_to not in self.transl:
            self.transl[from_to] = trans
        else:
            self.transl[from_to] = self.merge_translations(
                trans, self.transl[from_to])

    def translate(self, input_id, input_name, output_name):
        from_to = (input_name, output_name)
        if from_to in self.transl:
            return self.transl[from_to][input_id]
        else:
            return ''

    def merge_translations(self, a, b):
        c = collections.defaultdict(set)
        for input_name in a.map:
            if input_name in b.map:
                c[input_name] = a.map[input_name] | b.map[input_name]
            else:
                c[input_name] = a.map[input_name]
        for input_name in b.map:
            if input_name not in c:
                c[input_name] = b.map[input_name]
        return translation(a.input_name, a.output_name, c)


class translation(object):

    def __init__(self, input_name, output_name, mapping):
        self.input_name = input_name
        self.output_name = output_name
        self.transl = mapping  # Dict: {str} -> set(str, str, ...)
        for k in self.transl:
            if type(self.transl[k]) == type(''):
                self.transl[k] = set([self.transl[k]])
            if type(self.transl[k]) == type([]):
                self.transl[k] = set(self.transl[k])
        self.map = self.transl

    def translate(self, input_id, verbose=False):
        if verbose:
            print(("Looking for {a} in a dict with keys like {b}...".format(
            a=input_id, b=list(self.transl.keys())[:3])))
        
        if type(input_id) == type(set()):
            for _id in list(input_id):
                if _id in self.map:
                    return self.transl[_id]
            return set([''])
    
        elif input_id in self.map:
            if verbose: print('found')
            return self.transl[input_id]
        else:
            if verbose: print('did not find')
            return set([''])


def is_set_or_list(k):
    assert((type(k) == type(set([]))) or (type(k) == type([])) \
           or (type(k) == type(frozenset([]))))


    
class complexTranslation(object):
    """
    self.transl (dict): {frozenset a} -> {frozenset b}
    self.reverse_transl (dict): {frozenset b} -> {frozenset a}
    """
    
    def __init__(self, input_name, output_name, mapping,
                 name='Unnamed complexTranslation object',
                 max_orthologs=1e4):
        
        self.input_name = input_name
        self.output_name = output_name
        self.name = name
        self.max_orthologs = max_orthologs
        self.too_many_orthologs_a = set()
        self.too_many_orthologs_b = set()
        
        # mapping should be a dict of {frozenset} -> {frozenset}
        assert(type(mapping) == type({}))
        
        for k in mapping:
            is_set_or_list(k)
            is_set_or_list(mapping[k])
            
        self.transl = collections.defaultdict(frozenset)
        self.transl.update(dict([
            (frozenset(k), frozenset(mapping[k])) for k in mapping
        ]))
        
        for k in self.transl:
            is_set_or_list(self.transl[k])
            self.transl[k] = frozenset(
                self.transl[k] - set(['NULL', '', '-']))
            
        self.remove_untranslatable()
        
        self.map = self.transl
        self.reverse_transl = dict([
            (self.transl[k], k) for k in self.transl])
        
        self._info = {'gene sets a': len(self.transl),
                      'gene sets b': len(self.reverse_transl),
                      'name': self.name}
        
        if len(self.transl) > 1:
            self.qc()

    def get_reverse_translator(self):
        return complexTranslation(
            self.output_name, self.input_name, self.reverse_transl)

    def remove_untranslatable(self):
        for x in set([x for x in self.transl if (len(self.transl[x]) == 0)]):
            del self.transl[x]

    def check_if_too_many_orthologs(self, frozen_set_lang_a, seta_to_add, setb_to_add, del_from_a):
        not_too_many = True
        # Are there two many homologs in language A to include this mapping?
        # Specifically, are there too many locus IDs mapping to an overlapping
        # set of ENSGs?
        if len(frozen_set_lang_a | seta_to_add) > self.max_orthologs:
            self.too_many_orthologs_a |= seta_to_add
            self.too_many_orthologs_b |= setb_to_add
            if len(frozen_set_lang_a):
                del_from_a[frozen_set_lang_a] = ''
            not_too_many = False

        # Are there two many homologs in language B to include this mapping?
        # Specifically, are there too many ENSGs mapped to the same locus ID?
        if len(frozenset(setb_to_add | self.transl[frozen_set_lang_a])) > self.max_orthologs:
            self.too_many_orthologs_b |= seta_to_add
            self.too_many_orthologs_a |= setb_to_add
            if len(frozen_set_lang_a):
                del_from_a[frozen_set_lang_a] = ''
            not_too_many = False
        
        return not_too_many
    
    def define_mappings_from_list_of_paired_sets(
        self, list_of_paired_sets, max_orthologs=5):
        for (seta, setb) in list_of_paired_sets:
            if (len(seta)>max_orthologs) or (len(setb)>max_orthologs):
                continue
            self.transl[frozenset(seta)] = frozenset(setb)
            self.reverse_transl[frozenset(setb)] = frozenset(seta)
            
    def add_mapping(
            self, set_translate_from, set_translate_to, 
            multiple_homologs_in_native_language_possible=True, verbose=False):
        
        is_set_or_list(set_translate_from)
        is_set_or_list(set_translate_to)
        seta, setb = (frozenset(set_translate_from), frozenset(set_translate_to))

        too_many = False  
        if (len(seta & self.too_many_orthologs_a) > 0):
            self.too_many_orthologs_a |= seta
            too_many = True
        if (len(setb & self.too_many_orthologs_b) > 0):
            self.too_many_orthologs_b |= setb
            too_many = True
            
        if too_many:
            if verbose:
                print("{0} has too many orthologs, ignoring...".format(seta))
            return#'''
        
        found = False
        del_from_a = {}
        del_from_b = {}
        
        if not multiple_homologs_in_native_language_possible:
            if self.check_if_too_many_orthologs(frozenset(), seta, setb, del_from_a):
                self.transl[seta] |= setb
                    
        else:
            for frozen_set_lang_a in self.transl:

                if len(frozen_set_lang_a & seta) > 0:
                    print("The same language A ID {0} has been observed before: {1}: {2}".format(
                        seta, frozen_set_lang_a, self.transl[frozen_set_lang_a]))

                    if len(frozen_set_lang_a) > 1:
                        print("???: {0} {1}".format(frozen_set_lang_a, seta))
                        input()

                    if self.check_if_too_many_orthologs(frozen_set_lang_a, seta, setb, del_from_a):
                        self.transl[frozen_set_lang_a] = frozenset(setb | self.transl[frozen_set_lang_a])

                    found = True

            if not found:
                if self.check_if_too_many_orthologs(frozenset(), seta, setb, del_from_a):
                    self.transl[seta] = setb

        for frozen_set_lang_a in del_from_a:
            del self.transl[frozen_set_lang_a]

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
    def refresh(self):
        self.remove_untranslatable()
        
        self._info = {
            'Name': self.name,
            'gene sets a': len(self.transl),
            'gene sets b': len(self.reverse_transl)}

    def info(self):
        self.refresh()
        self.qc()
        return ' '.join(['language.complexTranslation info:', str(self._info)])

    def qc(self):
        known = set()
        in_two = set()
        
        for k in self.transl:
            if len(known & self.transl[k]) > 0:
                #print "Have an ID in two gene sets..."
                in_two |= known & self.transl[k]

            known |= self.transl[k]
            
        if len(in_two) > 0:
            print("Have {i} {b} language IDs in two gene sets...".format(
                i=len(in_two), b=self.output_name))
            
        else:
            print("complexTranslation.qc(): All IDs are in one gene set.")
    
    def translate(
        self, input_id, verbose=False, reverse=False,
        multiple_homologs_in_native_language_possible=True,
        return_str=False):
        
        # input_id is frozenset
        if type(input_id) == type(''):
            input_id = frozenset([input_id])
            
        if verbose:
            print('looking for ', input_id)
            
        translates_to = set([])
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




        
        
