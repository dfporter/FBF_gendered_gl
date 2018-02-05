import collections
import scipy.stats
import pandas, random
import targetSet
import functools

class complexTargetSet(targetSet.targetSet):
    """Holds a set of frozensets of target genes.
    Arguments:
        cgenes: a set of sets, or other iterable of iterables. The iterables are ids.
        language: the language of the ids, known to this object as its native language.
    
    This object will calculate overlaps with another complexTargetSet object. If
    the two complexTargetSets are in different languages, a translator must be provided.
    """
    
    def __init__(self, cgenes, language, name='Unnamed'):
        
        self.cgenes = cgenes  # Set of sets.
        self.name, self.language = (name, language)
        self.overlaps_info = {}
        
        if type(self.cgenes) == type([]):
            self.cgenes = set([frozenset(x) for x in self.cgenes])
            
        if (type(self.cgenes) == type(set([]))) or \
           (type(self.cgenes) == type(frozenset([]))):
            self.cgenes = set([frozenset(x) for x in cgenes])
            
        self.cgenes -= set(['NULL', '', '-'])
        
        self.refresh()
        
        print("\n...Created a complexTargetSet (%s) in language %s:" % (name, language))

    def stats(self):
        
        hist = collections.defaultdict(int)
        for x in self.cgenes:
            hist[len(x)] += 1
            
        df = pandas.DataFrame([(size, freq) for size, freq in list(hist.items())])
        
        print("{0}".format(self._info))
        #print('--\n\n')

    def num(self):
        return len(self.cgenes)

    def refresh(self):
        self._info = {
            'Number of gene sets': len(self.cgenes),
            'Native language': self.language,
            'Name': self.name}        

    def set_of_all_targets(self, language=None, translator=None):
        if (language == self.language) or (language == None):
            return set(functools.reduce(lambda a, b: a | b, self.cgenes))
        else:
            translated = [translator.translate(x) for x in self.cgenes if len(translator.translate(x))>0]
            return set(functools.reduce(lambda a, b: a| b, translated))
                       
    def translatable_cgenes(self, translator):
        targets_with_translations = [
            x for x in self.cgenes if (len(translator.translate(x))>0)]
        
        return set(targets_with_translations)

    def remove_nontranslatable_cgenes(self, translator):
        """Sets self.cgenes as targets with translations."""
        
        print("Removing non-translatable genes: input {0}".format(len(self.cgenes)))
        self.cgenes = set([
            x for x in self.cgenes if (len(translator.translate(x))>0)])
        print("...{0} were translatable.".format(len(self.cgenes)))
        
        self.refresh()
        
    def combine_targets_based_on_homology(self, translator):
        
        new_cgenes = [set([]) for x in range(len(self.cgenes))]
      
        #print('{0} len cgens ----->'.format(len(self.cgenes)))
        #print(new_cgenes)

        for index1, genes1 in enumerate(list(self.cgenes)):
            transl1 = translator.translate(genes1, multiple_homologs_in_native_language_possible=False)
            #print('-----> index1 {0} genes1 {1}'.format(index1, genes1))
            #print(new_cgenes[:100])
            
            #print('Before processing index1: new_cgenes[index1 + 1 = {0}] = {1}'.format(index1 + 1, new_cgenes[index1 + 1]))
            
            #print('before |=, new_cgenes[index1 = {0}] = {1}'.format(index1, new_cgenes[index1]))
            #print('type of new_cgenes=', type(new_cgenes))
            new_cgenes[index1] = new_cgenes[index1] | set(genes1)
            #print(new_cgenes[:100])
            #print('after |=, new_cgenes[index1 = {0}] = {1}'.format(index1, new_cgenes[index1]))
            
            #print('After |= processing index1: new_cgenes[index1 + 1 = {0}] = {1}'.format(index1 + 1, new_cgenes[index1 + 1]))
            #print('After |= processing index1: new_cgenes[index1 + 2 = {0}] = {1}'.format(index1 + 2, new_cgenes[index1 + 2]))
            if len(transl1) == 0:
                continue
                
            
            for index2, genes2 in enumerate(list(self.cgenes)):
                transl2 = translator.translate(genes2, multiple_homologs_in_native_language_possible=False)
                
                if transl1 & transl2:
                    #print('transl1 {0} transl2 {1}'.format(transl1, transl2))
                    new_cgenes[index1] |= genes2
                    new_cgenes[index2] |= genes2 | genes1
                    #print('{0}, {1} indexes together, len {2}'.format(index1, index2, len(new_cgenes[index1])))
                    #print('new_cgenes[index1 = {0}] = {1}'.format(index1, new_cgenes[index1]))
                    #print('new_cgenes[index2 = {0}] = {1}'.format(index2, new_cgenes[index2]))
                    
            #print('After processing index1: new_cgenes[index1 + 1 = {0}] = {1}'.format(index1 + 1, new_cgenes[index1 + 1]))
            
        print(str(new_cgenes)[:1000])
        print("Combining targets in language {0} based on overlapping translations {1}->{2}".format(
            self.language, translator.input_name, translator.output_name))
        
        print("Input is {0}  IDs (unflattened), {1} (flattened)".format(
            len(self.cgenes), len(set.union(*[set(x) for x in self.cgenes]))))
        
        self.cgenes = set([frozenset(x) for x in new_cgenes])
        
        # Reduce to one example gene.
        self.cgenes = set([frozenset(random.sample(x, 1)) for x in self.cgenes])
        
        hm = list(self.cgenes)[0]
        print('hm len = {0}'.format(len(hm)))
        print("After collapse: {0}  IDs (unflattened), {1} (flattened)".format(
            len(self.cgenes), len(set.union(*[set(x) for x in self.cgenes]))))
        
    def overlap_with_complexTargetSet(self, cts, translator, verbose=False):
        """Translate set of frozensets <-> set of frozensets.
        translator is complexTranslation object.
        """
        translator.refresh()
        
        if verbose:
            print("This complexTargetSet has {0} IDs of the form: {1}".format(
                len(self.cgenes), random.sample(self.cgenes, 1)[0]))

            print("It's info is:")
            print(self.info())

        a_targets_with_translations = [x for x in self.cgenes if (len(translator.translate(x))>0)]
        b_targets_with_translations = [x for x in cts.cgenes if (len(translator.translate(x, reverse=True))>0)]
        
        all_a_targets_with_translations_in_language_a = set.union(*[set(x) for x in a_targets_with_translations])
        all_b_target_ids_in_lang_b = set.union(*[set(x) for x in b_targets_with_translations])            
        
        def display(iterable):
            return (len(iterable), random.sample(iterable, 1))
        
        print("""
a_language_ids_with_translation (superset of a_targets): {0}, ex: "{1}"
b_language_ids_with_translation (superset of b_targets): {2}, ex: "{3}"
a_targets_with_translations to lang b (unflattened): {4}, ex: "{5}"
b_targets_with_translations to lang a (unflattened): {6}, ex: "{7}"
a_targets_with_translations to lang b (flattened): {8}, ex: "{9}"
all_a_target_ids in lang b after many<->many translation (flattened): {10}, ex: "{11}"
all_b_target_ids in lang b (flattened): {12}, ex: "{13}"
""".format(
    *display(list(translator.transl.keys())),
    *display(list(translator.reverse_transl.keys())),
    *display(a_targets_with_translations),
    *display(b_targets_with_translations),
    *display(all_a_targets_with_translations_in_language_a),  
    *display(set.union(*[translator.translate(x) for x in self.cgenes])),
    *display(all_b_target_ids_in_lang_b),
    ))
        
        self.annotated_cgenes = collections.defaultdict(dict)
        table = collections.defaultdict(int)
        
        # For each set of a_language_ids_with_translations (superset of targets).
        for frznset_lang_a, frznset_lang_b in list(translator.transl.items()):  
            
            targs_in_set = frznset_lang_a & all_a_targets_with_translations_in_language_a
    
            self.annotated_cgenes[frznset_lang_a]['translations'] = frznset_lang_b
            self.annotated_cgenes[frznset_lang_a]['a_names'] = targs_in_set
            
            if len(targs_in_set) > 0:
                
                if len(frznset_lang_b & all_b_target_ids_in_lang_b) > 0:
                    table['shared_targets'] += 1
                    self.annotated_cgenes[frznset_lang_a]['a_target'] = 1
                    self.annotated_cgenes[frznset_lang_a]['shared'] = 1
                    
                else:
                    table['a_target_only'] += 1
                    self.annotated_cgenes[frznset_lang_a]['a_target'] = 1
                    self.annotated_cgenes[frznset_lang_a]['shared'] = 0
                    
            else:
                
                if len(frznset_lang_b & all_b_target_ids_in_lang_b) > 0:
                    table['b_target_only'] += 1
                    self.annotated_cgenes[frznset_lang_a]['a_target'] = 0
                    self.annotated_cgenes[frznset_lang_a]['shared'] = 0
                    
                else:
                    table['not_a_or_b_target'] += 1
                    self.annotated_cgenes[frznset_lang_a]['a_target'] = 0
                    self.annotated_cgenes[frznset_lang_a]['shared'] = 0
        
        print(table)
        self.fisher_table(table)
        return table
        nope = '''
        for frznset in b_language_ids_with_translations:
            #a_set = translator.translate(frznset, reverse=True)
            if len(frznset & all_a_target_ids) > 0:
                if len(frznset & all_b_target_ids) > 0:
                    n['shared_targets'] += 1
                else:
                    n['a_target_only'] += 1
            else:
                if len(frznset & all_b_target_ids) > 0:
                    n['b_target_only'] += 1
                else:
                    n['not_a_or_b_target'] += 1
        print n
        self.fisher_table(n)'''

    def to_dataframe(self):
        #ld = [self.annotated_cgenes[x] for x in self.annotated_cgenes]
        df = pandas.DataFrame(list(self.annotated_cgenes.values()))
        df['translations'] = [", ".join(x) for x in df.translations]
        df['a_names'] = [", ".join(x) for x in df.a_names]
        return df
        
    def fisher_table(self, n):
        table = [
            [n['not_a_or_b_target'], n['b_target_only']],
            [n['a_target_only'], n['shared_targets']]]
        print("Fisher 2x2 table=", table)
        odds, pval = scipy.stats.fisher_exact(table)
        print("Fisher 2x2 pval=", pval)
        
    def overlap_with_targetSet(self, ts_b, translator):
        tup_a = (self.language, ts_b.language)
        if (self.language != translator.input_name) or (
            ts_b.language != translator.output_name):
            print("No translator for {a} > {b}.".format(
                a=self.language, b=ts_b.language))
            return
        has_translation = collections.defaultdict(set)
        no_translation = collections.defaultdict(set)
        first = 1
        for aset in self.cgenes:
            found_translation = False
            if first:
                print(aset)
                print(type(aset))
            if type(aset) == type(''):
                aset = frozenset([aset])
            if first:
                print(aset)
                print(type(aset))
            for x in aset:
                y = translator.translate(x, verbose=False)
                if y != set(['']):
                    has_translation[aset] |= y
                    found_translation = True
                else:
                    pass
            if not found_translation:
                no_translation[aset] = set()
            first = 0
        hits = ts_b.overlap_with_list_of_name_sets(
            list(has_translation.values()))
        print("""Has translation {has_tra}. No translation {no_trans}.
{i}/{all_genes} genes in this set (names {z}) \
have a name in {y}. {a}/{i} of this gene set is shared.""".format(
    has_tra=len(has_translation), no_trans=len(no_translation),
    i=len(has_translation), all_genes=len(self.cgenes), z=self.language,
    y=ts_b.language, a=len(hits)))
        self.overlaps_info.update({
'Have translation': set(has_translation.keys()),
'Have translation, are shared': set(hits),
'Have translation, are not shared': set(has_translation.keys()) - set(hits),
'No translation': self.cgenes - set(has_translation.keys()) })
        self._info.update({
'Have translation': len(set(has_translation.keys())),
'Have translation, are shared': len(set(hits)),
'Have translation, are not shared': len(set(has_translation.keys()) - set(hits)),
'No translation': len(self.cgenes - set(has_translation.keys())) })

    def info(self):
        self.refresh()
        return ' '.join(["complexTargetSet info: ", str(self._info)])

    def overlap_with_list_of_name_sets(self, list_of_name_sets):
        list_of_name_sets = set([frozenset(x) for x in list_of_name_sets])
        overlappers = collections.defaultdict(list)
        non_overlappers = collections.defaultdict(list)
        n = 0
        for aset in self.cgenes:
            aset = frozenset(aset)
            found = False
            for bset in list_of_name_sets:
                if len(aset & bset) > 0:
                    found = True
                    overlappers[aset].append(bset)
            if not found:
                non_overlappers[aset] = []
            n += 1
        print("As a set of frozensets, list to overlap with, length={a}".format(
            a=len(list_of_name_sets)))
        num_w_overlap = len(list(overlappers.keys()))
        num_wo_overlap = len(list(non_overlappers.keys()))
        print("{a}/{total} w overlap. {b} without.".format(
            a=num_w_overlap, total=len(self.cgenes), b=num_wo_overlap,
            ))
        print(n)
        return overlappers
    