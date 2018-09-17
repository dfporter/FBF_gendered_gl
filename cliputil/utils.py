import pandas, collections, re

class translator():
    
    def transl(self, gtfname='/opt/lib/gtf_with_names_column.txt', use_gene_biotype=True):
        """Link WB IDs to public names, and public names to biotypes."""
        
        gtf = pandas.read_csv(gtfname, sep='\t')
        
        self.wbid_to_name = collections.defaultdict(str)
        self.wbid_to_name.update(dict(zip(gtf['gene_id'].tolist(), gtf['gene_name'].tolist())))
        
        self.name_to_wbid = collections.defaultdict(str)
        self.name_to_wbid.update(dict(zip(gtf['gene_name'].tolist(), gtf['gene_id'].tolist())))
        
        if use_gene_biotype:
            biotypes = [re.search('gene_biotype "([^"]+)"', x).group(1) for x in gtf['8']]
            self.name_to_biotype = dict(zip(gtf['gene_name'].tolist(), biotypes))
            
        else:  # Using transctrip biotype.
            self.name_to_biotype = dict(zip(gtf['gene_name'].tolist(), gtf['biotype'].tolist()))


def overlaps_in_one_list(_list):

    def get_overlaps_if_exist(_list):
        overlapping = []
        a_only = []
        if len(_list) < 2:
            return _list
        for x, a in enumerate(_list[:-1]):
            overlaps = False
            c = list(a)
            for y, b in enumerate(_list[x + 1:]):
                if overlap_of_tups(a, b):
                    c[1] = min([c[1], b[1]])
                    c[2] = max([c[2], b[2]])
                    overlaps = True

            if overlaps:
                overlapping.append(tuple(c))
            else:
                a_only.append(a)

        return (a_only, overlapping)

    overlapping = ['init']
    merged = _list
    while len(overlapping) > 0:
        a_only, overlapping = get_overlaps_if_exist(merged)
        merged = a_only + overlapping
        print('Merging {0} uinique, {1} overlapping ){2}...'.format(len(a_only), len(overlapping), str(overlapping)[:100]))

    return merged


def overlaps_in_lists(list_a, list_b):
    overlapping = []
    a_only = []
    b_only = []
    for a in list_a:
        c = list(a)
        overlaps = False
        for b in list_b:
            if overlap_of_tups(a, b):
                c[1] = min([c[1], b[1]])
                c[2] = max([c[2], b[2]])
                overlaps = True

        if overlaps:
            overlapping.append(tuple(c))
        else:
            a_only.append(a)

    for b in list_b:
        overlaps = False
        for a in list_a:
            if overlap_of_tups(a, b):
                overlaps = True

        if not overlaps:
            b_only.append(b)

    return (a_only, overlapping, b_only)


def overlap_of_tups(a, b):
    if not (a[1] <= b[1] <= a[2] or a[1] <= b[2] <= a[2] or b[1]) <= a[1] <= b[2]:
        return b[1] <= a[2] <= b[2] and True
    else:
        return False


def seq_from_iv(chrm, start, end, strand, sequences):
    seq = sequences[chrm][start:end]
    if strand == '-':
        return rc(seq)
    else:
        return seq


def complement(s):
    basecomplement = {'A': 'T',
     'C': 'G',
     'G': 'C',
     'T': 'A',
     'N': 'N',
     'a': 't',
     'c': 'g',
     'g': 'c',
     't': 'a',
     'n': 'n'}
    return ''.join([ basecomplement[base] for base in list(s) ])


def rc(s):
    return complement(s[::-1])