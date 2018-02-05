import language
import targetSet
import complexTargetSet
import human_gene_ensmbl_map

class pum2Helper():
    
    @staticmethod
    def pum2_targets_as_list_of_ensmbl_sets():

        human_gene_symbol_to_human_ensmbl, human_ensmbl_to_human_gene_symbol = \
                                       human_gene_ensmbl_map.human_gene_ensmbl_map()

        pum2_targets = set([x.rstrip('\n') for x in open('pum2_targets.txt', 'r')])
        pum2_targets -= set(['NULL', '', '-'])

        symb_to_ens_tr = language.translation(
            'Gene symbol', 'Ensmbl', human_gene_symbol_to_human_ensmbl)

        return [symb_to_ens_tr.translate(gene_symbol) for gene_symbol in pum2_targets]

    @classmethod
    def pum2_targets_as_ensembl(cls, orthos):

        # PUM2 targets.
        pum2_in_ensmbl = cls.pum2_targets_as_list_of_ensmbl_sets()

        # Create PUM2 target set in Ensmbl language.
        pum2_ensembl_ts = complexTargetSet.complexTargetSet(pum2_in_ensmbl, 'Ensmbl', name='PUM2')
        pum2_ensembl_ts.stats()
        
        # Get the translator object ortho_cmap (worm locus ID -> Ensmbl) to remove non-translatable.
        ortho_locus_to_ensmbl, human_ensembl_to_locus_id, ortho_cmap = orthos

        # Human ensmbl <-> worm locus
        print("Total number of loci in ortholist: {v}. With translation to ESNGxx: {z}.".format(
            v=len(ortho_locus_to_ensmbl),
            z=len([x for x in ortho_locus_to_ensmbl if len(ortho_locus_to_ensmbl[x])])))
        
        pum2_ensembl_ts.remove_nontranslatable_cgenes(
            ortho_cmap.get_reverse_translator())

        return pum2_ensembl_ts