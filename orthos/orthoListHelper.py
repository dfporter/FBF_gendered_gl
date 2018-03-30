import pandas, os, re, sys, importlib, random, functools
import numpy as np

import language
importlib.reload(language)


class orthoListHelper():
    
    @staticmethod
    def read_ortholist_file(fname, max_orthologs=2):
        """Read csv into a dict of str(locus_id)<->set(other ids of orthologs).
        """

        # Read file into dict of locus ID -> row of data.
        df = pandas.read_csv(fname, sep=',')
        locus_id_to_row = dict([(re.sub(' ', '', row['Common Name']), row) for row in df.to_dict('records')])
        
        # All columns in cols list are in ENSG, or 'Ensmbl' language.
        ensg_cols_of_orthologs = ['Ensembl Compara', 'InParanoid', 'Homologene', 'OrthoMCL']
        
        human_ensembl_to_locus_id = {}

        # Initialize a complexTranslation object.
        complexTransl = language.complexTranslation('Common Name', 'Ensmbl', {}, max_orthologs=max_orthologs)

        # Populate the complexTranslation object, and a couple dicts for translation.
        for i, (locus_id, row) in enumerate(locus_id_to_row.items()):

            if not(i % 1e3):
                print('Reading ortholist file, line: ', i)

            other_ids = set()

            ensg_ids = functools.reduce(
                lambda x, y: ','.join([x, y]), [str(row[col]) for col in ensg_cols_of_orthologs])
            ensg_ids =  [x.replace(' ', '').replace('nan', '') for x in re.split('[ ,]', ensg_ids)]

            for ensg in ensg_ids:
                human_ensembl_to_locus_id[ensg] = locus_id

            other_ids |= set(ensg_ids) - set(['NULL', '', '-'])

            if len(other_ids) > 0:
                complexTransl.add_mapping(frozenset([locus_id]), frozenset(other_ids), verbose=False)

            locus_id_to_row[locus_id] = other_ids

        print(complexTransl.info())
        return locus_id_to_row, human_ensembl_to_locus_id, complexTransl
