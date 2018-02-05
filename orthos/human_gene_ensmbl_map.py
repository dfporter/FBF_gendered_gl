import pandas


def human_gene_ensmbl_map(filename='lib/human_gene_symbol_to_ensmbl.txt'):
    """Symbol to Ensmbl.
    """
    
    records = pandas.read_csv(filename, sep='\t').to_dict('records')
    
    human_gene_symbol_to_human_ensmbl = dict(
        (row['HGNC symbol'], row['Ensembl Gene ID']) for row in records)
    
    human_ensmbl_to_human_gene_symbol = dict(
        (row['Ensembl Gene ID'], row['HGNC symbol']) for row in records)
    
    return human_gene_symbol_to_human_ensmbl, human_ensmbl_to_human_gene_symbol
