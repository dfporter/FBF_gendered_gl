import functools, pandas, collections

import simple.overlapDictInfo

def simple_overlap(linguist, name='', col1='', col2='', verbose=False):
    stats = {
        'Worm public name': set(),
        'Human gene symbol': set(),
        'Human gene symbols with homologs': set(),
        'Worm public names with homologs': set(),
        'FBF targets (Worm public name)': set(),
        'Orthologs of FBF targets (Ensembl)': set(),
        'Orthologs of FBF targets (Human gene symbol)': set(),
        'PUM targets (Human gene symbols)': set(),
        'Orthologs of PUM targets (Ensembl)': set(),
        'Orthologs of PUM targets (Worm public name)': set(),
        'PUM targets shared (Human gene symbols)': set(),
        'Orthologs of shared targets (Ensembl)': set(),
        'FBF targets shared (Worm public name)': set()
    }
    stats_homologs_collapsed = {
        'Worm public name': 0,
        'Human gene symbol': 0,
        'Human gene symbols with homologs': 0,
        'Worm public names with homologs': 0,
        'FBF targets (Worm public name)': 0,
        'Orthologs of FBF targets (Ensembl)': 0,
        'Orthologs of FBF targets (Human gene symbol)': 0,
        'PUM targets (Human gene symbols)': 0,
        'Orthologs of PUM targets (Ensembl)': 0,
        'Orthologs of PUM targets (Worm public name)': 0,
        'PUM targets shared (Human gene symbols)': 0,
        'Orthologs of shared targets (Ensembl)': 0,
        'FBF targets shared (Worm public name)': 0
    }
    
    stats_homologs_uncollapsed = {
        'Worm public name': 0,
        'Human gene symbol': 0,
        'Human gene symbols with homologs': 0,
        'Worm public names with homologs': 0,
        'FBF targets (Worm public name)': 0,
        'Orthologs of FBF targets (Ensembl)': collections.defaultdict(int),
        'Orthologs of FBF targets (Human gene symbol)': 0,
        'PUM targets (Human gene symbols)': collections.defaultdict(int),
        'Orthologs of PUM targets (Ensembl)': collections.defaultdict(int),
        'Orthologs of PUM targets (Worm public name)': 0,
        'PUM targets shared (Human gene symbols)': 0,
        'Orthologs of shared targets (Ensembl)': 0,
        'FBF targets shared (Worm public name)': 0
    }
    
    linguist_rows = linguist.df.to_dict('records')
    
    for n, row in enumerate(linguist_rows):
        #print(row)
        
        stats['Worm public name'] |= row['Worm public name']
        stats_homologs_collapsed['Worm public name'] += 1
        stats_homologs_uncollapsed['Worm public name'] += len(row['Worm public name'])
        
        stats['Human gene symbol'] |= row['Human gene symbol']
        stats_homologs_collapsed['Human gene symbol'] += 1
        stats_homologs_uncollapsed['Human gene symbol'] += len(row['Human gene symbol'])

        if len(row['Human gene symbol']) and len(row['Worm public name']):
            stats['Human gene symbols with homologs'] |= row['Human gene symbol']
            stats_homologs_collapsed['Human gene symbols with homologs'] += 1
            stats_homologs_uncollapsed['Human gene symbols with homologs'] += len(row['Human gene symbol'])
            
            stats['Worm public names with homologs'] |= row['Worm public name']
            stats_homologs_collapsed['Worm public names with homologs'] += 1
            stats_homologs_uncollapsed['Worm public names with homologs'] += len(row['Worm public name'])
        
        if len(row[col1]):
            stats['FBF targets (Worm public name)'] |= row[col1]
            stats['Orthologs of FBF targets (Ensembl)'] |= row['Human ENSG']
            stats_homologs_collapsed['FBF targets (Worm public name)'] += 1
            stats_homologs_uncollapsed['FBF targets (Worm public name)'] += len(row[col1])
            stats_homologs_collapsed['Orthologs of FBF targets (Ensembl)'] += 1
            for item in list(row['Human ENSG']):
                stats_homologs_uncollapsed['Orthologs of FBF targets (Ensembl)'][item] += 1
            
            stats['Orthologs of FBF targets (Human gene symbol)'] |= row['Human gene symbol']
            stats_homologs_collapsed['Orthologs of FBF targets (Human gene symbol)'] += 1
            stats_homologs_uncollapsed['Orthologs of FBF targets (Human gene symbol)'] += len(
                row['Human gene symbol'])
            
        if len(row[col2]):
            stats['PUM targets (Human gene symbols)'] |= row[col2]
            stats_homologs_collapsed['PUM targets (Human gene symbols)'] += 1
            for item in list(row[col2]):
                stats_homologs_uncollapsed['PUM targets (Human gene symbols)'][item] += 1

            stats['Orthologs of PUM targets (Ensembl)'] |= row['Human ENSG']
            stats_homologs_collapsed['Orthologs of PUM targets (Ensembl)'] += 1
            for item in list(row['Human ENSG']):
                stats_homologs_uncollapsed['Orthologs of PUM targets (Ensembl)'][item] += 1
            #stats_homologs_uncollapsed['Orthologs of PUM targets (Ensembl)'] += len(row['Human ENSG'])

            stats['Orthologs of PUM targets (Worm public name)'] |= row['Worm public name']
            stats_homologs_collapsed['Orthologs of PUM targets (Worm public name)'] += 1
            stats_homologs_uncollapsed['Orthologs of PUM targets (Worm public name)'] += len(
                row['Worm public name'])

        if len(row[col1]) and len(row[col2]):
            stats['PUM targets shared (Human gene symbols)'] |= row[col2]
            stats['Orthologs of shared targets (Ensembl)'] |= row['Human ENSG']
            stats['FBF targets shared (Worm public name)'] |= row[col1]        
            stats_homologs_collapsed['FBF targets shared (Worm public name)'] += 1
            stats_homologs_uncollapsed['FBF targets shared (Worm public name)'] += len(row[col1])
            stats_homologs_collapsed['Orthologs of shared targets (Ensembl)'] += 1
            stats_homologs_uncollapsed['Orthologs of shared targets (Ensembl)'] += len(row['Human ENSG'])
            stats_homologs_collapsed['PUM targets shared (Human gene symbols)'] += 1
            stats_homologs_uncollapsed['PUM targets shared (Human gene symbols)'] += len(row[col2])
            
        #if not(n % 1000):
        #    print("Calc overlap for {0}:".format(n))
            #print(stats_homologs_uncollapsed)
            #break

    #for k,v in stats_homologs_collapsed.items():
    #    print(k, v)
#    df = pandas.DataFrame(stats)
    #ov = simple.overlapDictInfo.overlapDictInfo()
    #df_of_sets, stats_dict = ov.add_info(name, df, targs1, targs2)
    
    return stats, stats_homologs_collapsed,  stats_homologs_uncollapsed#df_of_sets, stats_dict