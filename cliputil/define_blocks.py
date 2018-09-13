import pandas, sys


def define_blocks(
	fname='../v2/tables/File S5 Complex frequencies.xls',
	outfname='../v2/tables/File S6 Blocks.xls'):
    
    df = pandas.read_excel(fname)#, sep='\t')

    if 'Gene name' in df.columns:
        df.index = df['Gene name']
    else:
        df.index = df['gene']
    #del df['gene']
    print(df.head())
    
    (blocki, blockii, blockiii) = ([], [], [])
    start_of_table = True
    in_ii, in_iii = (False, False)
    
    for gene in df.index:
        if start_of_table:
            blocki.append(gene)
        
        if gene == 'mib-1':
            start_of_table = False
            
        if gene == 'perm-1':
            in_ii = True
        
        if gene == 'gld-1':
            in_iii = True
        
        if in_ii:
            blockii.append(gene)
        
        if in_iii:
            blockiii.append(gene)
        
        if gene == 'gyg-2':
            in_ii = False
        
        if gene == 'htp-1':
            in_iii = False
    
    dfi = pandas.DataFrame([{'Gene name': x} for x in blocki])
    dfii = pandas.DataFrame([{'Gene name': x} for x in blockii])
    dfiii = pandas.DataFrame([{'Gene name': x} for x in blockiii])

    writer = pandas.ExcelWriter(outfname)
    
    dfi.to_excel(
        writer, sheet_name='Block I', index=False)
    dfii.to_excel(
        writer, sheet_name='Block II', index=False)
    dfiii.to_excel(
        writer, sheet_name='Block III', index=False)
    
    writer.save()


if __name__ == '__main__':
    print("Usage:")
    print("define_blocks.py <input read counts xls file> <output xls filename>")

    print("Loading heatmap data in {} and writing to {}.".format(
        sys.argv[1], sys.argv[2]))

    define_blocks(fname=sys.argv[1], outfname=sys.argv[2])

    print("Loading peaks")