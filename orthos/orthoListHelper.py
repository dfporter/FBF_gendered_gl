import pandas, os, re, sys, importlib, random, functools
import numpy as np

import simple.ManyToManyTranslator
importlib.reload(simple.ManyToManyTranslator)

class orthoListHelper():
    
    @staticmethod
    def read_ortholist_file(fname, max_orthologs=2):
        """Read csv into a dict of str(locus_id)<->set(other ids of orthologs).
        """

        # Read file into dict of locus ID -> row of data.
        df = pandas.read_csv(fname, sep=',')
        name_to_row = dict([(re.sub(' ', '', row['Common Name']), row) for row in df.to_dict('records')])
        
        # All columns in cols list are in ENSG, or 'Ensmbl' language.
        ensg_cols_of_orthologs = ['Ensembl Compara', 'InParanoid', 'Homologene', 'OrthoMCL']
        
        input_pairs = []
        # Populate the complexTranslation object, and a couple dicts for translation.
        for i, (worm_name, row) in enumerate(name_to_row.items()):

            if not(i % 1e3):
                print('Reading ortholist file, line: ', i)

            other_ids = set()

            ensg_ids = functools.reduce(
                lambda x, y: ','.join([x, y]), [str(row[col]) for col in ensg_cols_of_orthologs])
            ensg_ids =  [x.replace(' ', '').replace('nan', '') for x in re.split('[ ,]', ensg_ids)]

            other_ids |= set(ensg_ids) - set(['NULL', '', '-'])

            input_pairs.append([set([worm_name]), set(other_ids)])

        translator = simple.ManyToManyTranslator.ManyToManyTranslator(
            'Common Name', 'Ensmbl')
        pairs = translator.collapse_list_of_paired_sets(input_pairs)
        
        #worm_names = translator.expand_list_of_sets_to_include_all_in_group(
        #    [x[0] for x in pairs], [x[1] for x in pairs])
        #ensg = translator.expand_list_of_sets_to_include_all_in_group(
        #    [x[1] for x in pairs], [x[0] for x in pairs])
#        translator.define_mappings_from_list_of_paired_sets(pairs, max_orthologs=max_orthologs)

#        pairs = zip(worm_names, ensg)
        reverse_pairs = [[x[1], x[0]] for x in pairs]
        reverse_pairs = translator.collapse_list_of_paired_sets(reverse_pairs)
        pairs = [[x[1], x[0]] for x in reverse_pairs]

        translator.define_mappings_from_list_of_paired_sets(pairs, max_orthologs=max_orthologs)
        
        print("Finished making ortholist translator (worm names <-> human ENSG).")
        return translator
