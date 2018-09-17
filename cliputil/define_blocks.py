import pandas, sys, utils, excelHelper, peaksList

def get_wbid(dfs=[]):

    transl = utils.translator()
    transl.transl()

    for df in dfs:
        df['Wormbase ID'] = [transl.name_to_wbid.get(x, '') for x in df['Gene name']]

    return dfs

def subset_block_i(blocki_df, counts_df, peaks_xls='tables/Table S1 Peaks.xls'):
    #sp_df = pandas.read_excel(peaks_xls, sheet_name='SP FBF (25°)')

    print("Subsetting block I...")

    block_i_names = set(blocki_df['Gene name'])

    ##############
    # Counts filtering.
    counts_i = counts_df[[(x in block_i_names) for x in counts_df.index]].copy()
    counts_i['mini'] = [min(tup) for tup in zip(
        counts_i['25°C SP FBF rep 1'], counts_i['25°C SP FBF rep 2'], counts_i['25°C SP FBF rep 3'])]

    counts_cutoff = 3
    counts_i = counts_i[counts_i['mini']>counts_cutoff]

    print("{} RNAs in Block I have a reads-per-gene count >{} in all three replicates.".format(
        counts_i.shape[0], counts_cutoff))

    ##############
    # Program filtering.
    sp_pk = peaksList.peaksList()
    sp_pk.read_csv('combined_filtered/sp_both.txt')
    sp_pk.read_sp_vs_oo_as_programs()

    sp_pk.add_programs_to_df(blocki_df)


    print("Programs of block I RNAs:\n{}".format(blocki_df['Program'].value_counts()))

    print("Oogenic only: {}".format(blocki_df[blocki_df['Program']=='Oogenic only']))
    print("No program: {}".format(blocki_df[blocki_df['Program']=='']))

    blocki_df = blocki_df[[(x in set(counts_i.index)) for x in blocki_df['Gene name']]]
    sp_prog = blocki_df[blocki_df['Program']!='Oogenic only'].copy()

    print("{} RNAs in Block I belonged to the spermatogenic program.".format(sp_prog.shape[0]))
    print(sp_prog['Gene name'].tolist())

    print(blocki_df)

    ##############
    # Peaks filtering.
    sp_i = sp_pk.df[[(x in block_i_names) for x in sp_pk.df[sp_pk.gene_name_col]]].copy()

    oo_pk = peaksList.peaksList()
    oo_pk.read_csv('combined_filtered/oo_both.txt')
    oo_i = oo_pk.df[[(x in block_i_names) for x in oo_pk.df[oo_pk.gene_name_col]]].copy()

    print("{} SP peaks are in Block I.".format(sp_i.shape[0]))
    print("{} OO peaks are in Block I.".format(oo_i.shape[0]))

    sp_peak_only = set(sp_i.gene_name) - set(oo_i.gene_name)

    print("{} Block I RNAs have peaks in SP and not OO.".format(len(sp_peak_only)))
    cutoff = 25

    ##############
    # FBE filtering.
    print(sp_i['has_fbe'].value_counts())
    print(sp_i['minus_one_or_two_c'].value_counts())
    #fbe = sp_i[sp_i['Has UGUNNNAU (FBE)?']].copy()
    #w_fbe = set(sp_i[sp_i['minus_one_or_two_c']>0][sp_pk.gene_name_col])
    w_fbe = set(sp_i[sp_i['has_fbe']>0][sp_pk.gene_name_col])
    
    w_fbe_df = blocki_df[[ (x in w_fbe) for x in blocki_df['Gene name'] ]].copy()
    
    w_fbe_df.to_excel('tables/subset_of_block_I_with_fbe.xlsx')
    ##############
    # Peak height filtering.
    sp_i_robust = sp_i[[(x>cutoff) for x in sp_i['exp_reads']]].copy()

    print("{} of those peaks have >{} exp_reads.".format(sp_i_robust.shape[0], cutoff))

    robust_sp_plus = set(sp_i_robust[sp_pk.gene_name_col]) & w_fbe

    print("{} Block I RNAs have _, and have >{} reads.".format(
        len(robust_sp_plus), cutoff))
    #pk = peaksList.peaksList(dataframe=blocki_df, gene_name_col='Gene name')

    robust_sp_plus_df = sp_prog[[(x in robust_sp_plus) for x in sp_prog['Gene name']]].copy()

    robust_sp_plus_df.to_excel('tables/subset_of_block_I_with_robust_peaks_and_only_in_SP_program.xlsx')
    
    #pk.read_sp_vs_oo_as_programs()


def define_blocks(
	fname='../v2/tables/File S5 Complex frequencies.xls',
	outfname='../v2/tables/File S6 Blocks.xls'):
    
    counts_df = pandas.read_excel(fname)#, sep='\t')

    if 'Gene name' in counts_df.columns:
        counts_df.index = counts_df['Gene name']
    else:
        counts_df.index = counts_df['gene']

    print(counts_df.head())
    
    (blocki, blockii, blockiii, blockiv) = ([], [], [], [])
    start_of_table = True
    in_ii, in_iii, in_iv = (False, False, False)
    
    print("Identifying blocks in {}...".format(fname))

    for gene in counts_df.index:
        if start_of_table:
            blocki.append(gene)
        
        if gene == 'mib-1':
            start_of_table = False
            
        if gene == 'retr-1': #'perm-1':
            in_ii = True
        
        if gene == 'egg-6': #'gld-1':
            in_iii = True

        if gene == 'eef-1A.1':
            in_iv = True
        
        if in_ii:
            blockii.append(gene)
        
        if in_iii:
            blockiii.append(gene)
        
        if in_iv:
            blockiv.append(gene)

        if gene == 'fog-1': #'gyg-2':
            in_ii = False
        
        if gene == 'larp-1': #'htp-1':
            in_iii = False
    
        if gene == 'F43G6.3':
            in_iv = False

    # GLD-1 goes in block III based on PCA, despite the heatmap.
    blockiii.append('gld-1')

    dfi = pandas.DataFrame([{'Gene name': x} for x in blocki])
    dfii = pandas.DataFrame([{'Gene name': x} for x in blockii])
    dfiii = pandas.DataFrame([{'Gene name': x} for x in blockiii])
    dfiv = pandas.DataFrame([{'Gene name': x} for x in blockiv])

    (dfi, dfii, dfiii, dfiv) = get_wbid(dfs=[dfi, dfii, dfiii, dfiv])

    
    
    print("Adding gene descriptions...")

    helper = excelHelper.excelHelper()
    
    for df in [dfi, dfii, dfiii, dfiv]:
        #helper.add_gene_description(df)

        if 'transcript_id' in df.columns:
            del df['transcript_id']

    subset_block_i(dfi.copy(), counts_df)
    
    print("Writing excel file {}.".format(outfname))
    writer = pandas.ExcelWriter(outfname)

    dfi.to_excel(
        writer, sheet_name='Block I', index=False)
    dfii.to_excel(
        writer, sheet_name='Block II', index=False)
    dfiii.to_excel(
        writer, sheet_name='Block III', index=False)
    dfiv.to_excel(
        writer, sheet_name='Block IV', index=False)

    writer.save()


if __name__ == '__main__':
#python ../cliputil/define_blocks.py tables/File\ S5\ Complex\ frequencies.xls tables/File\ S6\ Blocks.xls

    print("Usage:")
    print("define_blocks.py <input read counts xls file> <output xls filename>")

    if(len(sys.argv) < 3):
        print("Input/output filenames were not given.")
        print("Using defaults {} and {}".format(
            'tables/File S5 Complex frequencies.xls',
            'tables/File S6 Blocks.xls'))
        sys.argv = [sys.argv[0]]
        sys.argv.append('tables/File S5 Complex frequencies.xls')
        sys.argv.append('tables/File S6 Blocks.xls')

    print("Loading heatmap data in {} and writing to {}.".format(
        sys.argv[1], sys.argv[2]))

    define_blocks(fname=sys.argv[1], outfname=sys.argv[2])

    print("Loading peaks")