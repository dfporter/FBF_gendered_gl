config_text = """
ortholist=lib/OrthoList_20130531.csv

#-----------

test1=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/tables/test1.txt
    language=gene_name
    species=worm
    read_csv=True
    read_every_line_as_gene=True
    lowercase=True
    >]

test2=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/tables/test2.txt
    language=Gene symbol
    species=human
    read_csv=True
    read_every_line_as_gene=True
    capitalize=True
    >]

#-----------

FbfOO25C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=OO FBF (25°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    >]

Top500FbfOO25C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=OO FBF (25°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    subset=True
    subset_function=lambda x: x['Rank']<501
    >]

FbfSP25C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=SP FBF (25°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    >]

Top500FbfSP25C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=SP FBF (25°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    subset=True
    subset_function=lambda x: x['Rank']<501
    >]
    
Fbf1_20C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=OO FBF-1 (20°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    >]

Top500Fbf1_20C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=OO FBF-1 (20°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    subset=True
    subset_function=lambda x: x['Rank']<501
    >]

Fbf2_20C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=OO FBF-2 (20°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    >]

Top500Fbf2_20C=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/Table S1 Peaks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=OO FBF-2 (20°C)#('OO FBF (25°C)', 'SP FBF (25°C)', 'OO FBF-1 (20°C)', 'OO FBF-2 (20°C)')
    column_name=Gene name
    subset=True
    subset_function=lambda x: x['Rank']<501
    >]
    
BlockI=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/File S6 Blocks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=Block I
    column_name=Gene name
    >]
    
BlockII=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/File S6 Blocks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=Block II
    column_name=Gene name
    >]
    
BlockIII=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/File S6 Blocks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=Block III
    column_name=Gene name
    >]

BlockIV=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/tables/File S6 Blocks.xls
    language=gene_name
    species=worm
    load_sheets_of_excel=True
    sheetnames=Block IV
    column_name=Gene name
    >]
#-----------

Hafner=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/orthos/lib/pum2_targets.txt
    language=Gene symbol
    species=human
    read_csv=True
    read_every_line_as_gene=True
    >]

#-----------

gld1_raj=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/orthos/lib/gld1_raj.xlsx
    language=Gene symbol
    species=worm
    load_sheets_of_excel=True
    column_name=Gene name
    sheetnames=NIHMS339375-supplement-02
    >]
    
#-----------

Lin_pum1=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/orthos/lib/supp_gad.298752.117_SUPPLEMENTAL_Table2sites.xlsx
    load_sheets_of_excel=True
    language=Gene symbol
    species=human
    sheetnames=WT Pum1
    column_name=overlap gene symbol
    capitalize=True
    >]

Lin_pum2=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/orthos/lib/supp_gad.298752.117_SUPPLEMENTAL_Table2sites.xlsx
    load_sheets_of_excel=True
    language=Gene symbol
    species=human
    sheetnames=WT Pum2
    column_name=overlap gene symbol
    capitalize=True
    >]
    
Lin_pum=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/orthos/lib/supp_gad.298752.117_SUPPLEMENTAL_Table2sites.xlsx
    load_sheets_of_excel=True
    language=Gene symbol
    species=human
    sheetnames=Pum1 or Pum2
    column_name=overlap gene symbol
    capitalize=True
    >]
    
Lin_pum_with_motif=[<
    filename=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/orthos/lib/supp_gad.298752.117_SUPPLEMENTAL_Table2sites.xlsx
    load_sheets_of_excel=True
    language=Gene symbol
    species=human
    sheetnames=Pum1 or Pum2
    column_name=overlap gene symbol
    capitalize=True
    subset=True
    subset_function=lambda x: (x['motif'] != 'unidentified')
    >]
"""