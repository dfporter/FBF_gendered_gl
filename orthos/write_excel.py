import pandas
import re
import os


def fix_ensembl_list(x):
    
    if type(x) == type(''):
        return re.sub('\s*', '', x)
    
    else: return ''

def clean(fbf_targs):
    
    if not os.path.exists('output/'):
        os.system('mkdir output/')

    fbf_targs['Human ortholog Ensembl IDs'] = [
        fix_ensembl_list(x) for x in \
        fbf_targs['Human ortholog Ensembl IDs'].tolist()]
    
    fbf_targs.replace(to_replace='', value='-', inplace=True)
    fbf_targs.sort_values(by=['Gene name'], inplace=True)
    
def write_tab3(fbf_targs, sheet_name=None, writer=None):

    clean(fbf_targs)
    
    if writer is None:
        writer = pandas.ExcelWriter('output/Table S5.xls')
        
    tab3_columns = [
        'Worm locus ID', 'Worm gene name','Human ortholog Ensembl IDs',
        'PUM2 target HGNC symbol']

    fbf_targs['Worm locus ID'] = fbf_targs['Locus ID']
    fbf_targs['Worm gene name'] = fbf_targs['Gene name']

    with_pum2_targ = fbf_targs[fbf_targs['pum2 shared']==1]

    try:
        with_pum2_targ.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
            columns=tab3_columns
            )
    except:
        with_pum2_targ.to_excel(
            writer,
            sheet_name=sheet_name[:25],
            index=False,
            columns=tab3_columns
            )
def write_excel(fbf_targs, tabs=[1, 2, 3, 4, 5], writer=None):
    
    if not os.path.exists('output/'):
        os.system('mkdir output/')
    
    tabs = [str(x) for x in tabs]
    
    if writer is None:
        writer = pandas.ExcelWriter('output/Table S5.xls')
    
    clean(fbf_targs)
    
    # Tab 1: sex determination genes.
    # Tab 2:
    if '2' in tabs:
        tab2_columns = [
            'Worm locus ID', 'Yeast Ensembl ID',
            'Yeast gene name',
            'Human ortholog Ensembl IDs']
        with_yeast_targ = fbf_targs[fbf_targs['yeast shared']==1]
        with_yeast_targ.to_excel(
            writer,
            sheet_name='tab2 FBF and yeast PUF shared',
            index=False,
            columns=tab2_columns)
        
    # Tab 3:
    if '3' in tabs:
        tab3_columns = [
            'Worm locus ID', 'Worm gene name','Human ortholog Ensembl IDs',
            'PUM2 target HGNC symbol']
        
        fbf_targs['Worm locus ID'] = fbf_targs['Locus ID']
        fbf_targs['Worm gene name'] = fbf_targs['Gene name']
        
        with_pum2_targ = fbf_targs[fbf_targs['pum2 shared']==1]
        
        with_pum2_targ.to_excel(
            writer,
            sheet_name='tab3 FBF and PUM2 shared',
            index=False,
            columns=tab3_columns
            )
        
    # Tab 4:
    if '4' in tabs:
        tab4_columns = [
            'Worm locus ID', 'Worm gene name',
            'Yeast Ensembl ID', 'Yeast gene name',
            'Human ortholog Ensembl IDs']
        with_both = fbf_targs[fbf_targs['pum2 shared']==1]
        with_both = with_both[with_both['yeast shared']==1]
        with_both.to_excel(
            writer,
            sheet_name='tab4 FBF, yeast, human shared',
            index=False,
            columns=tab4_columns
        )
    # Tab 5:
    if '5' in tabs:
        tab5_columns = [
            'Gene name', 'Locus ID',
            'Human ortholog Ensembl IDs',
            'Ensembl Compara',
            'InParanoid', 'Homologene', 'OrthoMCL']
        with_human = fbf_targs[fbf_targs['with_human']==1]
        with_human.to_excel(
            writer,
            sheet_name='tab5 FBF trgts with hum ortho.',
            index=False,
            columns=tab5_columns
            )
    writer.save()
