import collections, scipy.stats, pandas, re

from language import *
from human_gene_ensmbl_map import human_gene_ensmbl_map
import complexTargetSet
import targetSet


class fbfTargetSet(complexTargetSet.complexTargetSet):

    def __init__(self, input_filename=None, input_df=None):
        
        if input_filename is not None:
            fbf_ts = self.fbf_from_peaks(input_filename=input_filename)

            self.fbf_cts = complexTargetSet.complexTargetSet(
                [frozenset([x]) for x in fbf_ts.genes], 'Loci', name='FBF')

            self.fbf_cts.stats()
            
        if input_df is not None:
            if 'transcript_id' in input_df.columns:
                fbf_targs = set([self.simplify_locus_id(x) for x in input_df.transcript_id])
                self.fbf_cts = complexTargetSet.complexTargetSet(
                    [frozenset([x]) for x in fbf_targs], 'Loci', name='FBF')
                self.fbf_cts.stats()

    @staticmethod
    def simplify_locus_id(locus_id):
        txpt = locus_id.split('.')
        if len(txpt) > 1:
            locus_id = '.'.join(txpt[:2])
        return re.sub('[a-z]$', '', locus_id)
    
    def fbf_from_peaks(self, input_filename='lib/combined_fbf.txt'):
        """Returns targetSet object."""

        _targets = pandas.read_csv(input_filename, sep='\t')

        fbf_targs = set([simplify_locus_id(x) for x in _targets.transcript_id])        
                 
        return targetSet.targetSet(fbf_targs, 'Loci')

    def annotate_fbf_targs_df(self, pum2_ensembl_ts, orthos, use_this_df=None):
        
        ortho_locus_to_ensmbl, human_ensembl_to_locus_id, ortho_cmap = orthos
        
        human_gene_symbol_to_human_ensmbl, human_ensmbl_to_human_gene_symbol = \
                                           human_gene_ensmbl_map()
        
        if use_this_df is not None:
            df = use_this_df
        else:
            # .to_dataframe() uses the fbf_cts.annotated_cgenes information populated
            # by a fbf_cts.overlap_with_complexTargetSet() call. Specifically, this is
            # a_target, translations, and shared values. These denote whether the gene
            # is an FBF target, ENSG homologs, and whether it is also a PUM2 target,
            # respectively.
            df = self.fbf_cts.to_dataframe()  
            df['Locus ID'] = df['a_names']
        
        fbf_targs = df[df['a_target']==1].copy()  # FBF targets only.
        
        #fbf_targs = ann.annotate_a_dataframe(fbf_targs, key='Locus ID')        
        pum2_ensg, pum2_symbol = ([], [])
        
        pum2_targets_symbol = set(
            [x.rstrip('\n') for x in open('pum2_targets.txt', 'r')]) - set(['NULL', '', '-'])
        
        pum2_targ_set = pum2_ensembl_ts.set_of_all_targets(language='Ensmbl')
        
        for ensg_ids_of_fbf_target_homologs in fbf_targs['translations']:
            
            ensg_ids_of_fbf_target_homologs = set(ensg_ids_of_fbf_target_homologs.split(', '))
            
            # PUM2 target symbol -> ENSGs, overlaped with FBF locus ID -> ENSG homolog group.
            # This includes ENSG groups that all both homologs of the FBF locus, but only
            # one of which may be a PUM2 target.
            pum2_ensg.append(ensg_ids_of_fbf_target_homologs & pum2_targ_set)  
           
            pum2_symbol.append(set())
            for ensg_id in pum2_ensg[-1]:
                if ensg_id in human_ensmbl_to_human_gene_symbol:
                    pum2_symbol[-1].add(human_ensmbl_to_human_gene_symbol[ensg_id])  # All homologs.
                    
            pum2_symbol[-1] = pum2_symbol[-1] & pum2_targets_symbol  # Only the homologs actually targetted.

        fbf_targs['Human ortholog Ensembl IDs'] = [x for x in fbf_targs['translations'].tolist()]
        
        fbf_targs['PUM2 target HGNC symbol'] = [", ".join(x) for x in pum2_symbol]

        gtf = pandas.read_csv('/opt/lib/gtf_with_names_column.txt', sep='\t')

        gtf.transcript = [self.simplify_locus_id(x) for x in gtf['transcript_id'].tolist()]

        loci_to_gene_name = dict(list(zip(gtf.transcript, gtf.gene_name)))
            
        fbf_targs['Gene name'] = [loci_to_gene_name[x] for x in fbf_targs['Locus ID'].tolist()]
        fbf_targs['pum2 shared'] = fbf_targs['shared']
        fbf_targs.drop(['shared'], axis=1, inplace=True)
        fbf_targs['Worm locus ID'] = fbf_targs['Locus ID']
        fbf_targs['Worm gene name'] = fbf_targs['Gene name']
        fbf_targs['with_human'] = 1
        return fbf_targs

