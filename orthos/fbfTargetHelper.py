import re
import language
import targetSet
import complexTargetSet
import human_gene_ensmbl_map

class fbfTargetHelper():
    
    @classmethod
    def get_gene_name_locus_id_mapper(cls, fname='/opt/lib/worm_txpt_id_gene_name_map.txt'):
        
        worm_translator = language.complexTranslation('Gene name', 'Locus ID', {}, max_orthologs=10000, name='worm IDs')
        
        with open(fname) as f:
            """Format:
            Gene stable ID	Transcript stable ID	Gene name
            WBGene00000001	Y110A7A.10.1	aap-1
            """
            header = next(f)

            for li in f:
                s = li.rstrip('\n').split('\t')

                if s[2] != '':
                    s[2], s[1] = (cls.simplify_locus_id(s[2]), cls.simplify_locus_id(s[1]))
                    worm_translator.add_mapping([s[2]], [s[1]], multiple_homologs_in_native_language_possible=False)
        
        return worm_translator

    @staticmethod
    def simplify_locus_id(locus_id):
        txpt = locus_id.split('.')
        if len(txpt) > 1:
            locus_id = '.'.join(txpt[:2])
        return re.sub('[a-z]$', '', locus_id)