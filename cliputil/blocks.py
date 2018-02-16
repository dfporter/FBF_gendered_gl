import pandas

in_file = '/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/tables/File S6 Blocks.xlsx'

i_df = pandas.read_excel(in_file, sheetname='Block I', index=False)
blocki = set(i_df['Gene name'].tolist()) - set([''])
ii_df = pandas.read_excel(in_file, sheetname='Block II', index=False)
blockii = set(ii_df['Gene name'].tolist()) - set([''])
iii_df = pandas.read_excel(in_file, sheetname='Block III', index=False)
blockiii = set(iii_df['Gene name'].tolist()) - set([''])