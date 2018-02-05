import collections
import scipy.stats
import pandas


class targetSet(object):

    def __init__(self, genes, language):
        assert(type(genes) == type(set([])))
        self.genes = genes  # Set of strings, one per gene.
        self.language = language  # String
    
    def translate(self, translator, output_language):
        
        tup = (self.language, output_language)
        
        if self.language != translator.input_name:
            print(("Wrong input language for translator ({a} vs {b}).".format(
                a=self.language, b=translator.input_name)))
            
        if output_language != translator.output_name:
            print(("Wrong output language for translator ({a} vs {b}).".format(
                a=self.language, b=translator.output_name)))
            
        return [
            translator.translate(x) for x in self.genes]

    def overlap_with(self, targset_b, translator):
        tup_a = (self.language, targset_b.language)
        if tup_a not in translator.transl:
            return "No translator for {a} > {b}.".format(
                a=self.language, b=targset_b.language)
        else:
            a_to_b = self.translate(translator, targset_b.language)
            have_transl = [x for x in a_to_b if x != '']
            shared_hits = targset_b.overlap_with_list_of_name_sets(have_transl)
            print(("""{i}/{all_genes} genes in this set (names {z}) \
have a name in {y}. {a}/{i} of this gene set is shared.""".format(
    i=len(have_transl), all_genes=len(self.genes), z=self.language,
    y=targset_b.language, a=len(shared_hits))))

    def overlap_with_list_of_name_sets(self, list_of_name_sets):
        overlapers = set()
        non_overlappers = set()
        for aset in list_of_name_sets:
            if len(aset & self.genes) > 0:
                overlappers.add(aset)
            else:
                non_overlappers.add(aset)
        return overlappers