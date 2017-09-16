# Python bytecode 2.7 (62211) disassembled from Python 2.7
# Embedded file name: /opt/rh/python27/root/usr/lib/python2.7/site-packages/cliputil/peaksList.py
# Compiled at: 2016-08-25 17:30:14
import pandas
import numpy as np
import re
import sys
import collections
import HTSeq
from locatedPeak import locatedPeak

class peaksList(object):

    def __init__(self, dataframe=pandas.DataFrame(), name='Unnamed', gene_name_col='gene_name'):
        self.name = name
        self.df = dataframe
        self.gene_name_col = gene_name_col

    def __str__(self):
        li = '\n{a}:\n{b} peaks\n{c} genes\n'.format(a=self.name, b=len(self.df.index), c=len(set(self.df[self.gene_name_col].tolist())))
        if 'Classification OO/SP GL' in self.df.columns:
            li += '\nClassifications  (from Ortiz Deseq file):\n{d}'.format(d=str(self.df['Classification OO/SP GL'].value_counts()))
        if 'Program' in self.df.columns:
            li += '\nPrograms:\n{d}'.format(d=str(self.df['Program'].value_counts()))
        return li

    def read_csv(self, fname):
        self.df = pandas.read_csv(fname, sep='\t', index_col=False)
        if 'Wormbase ID' not in self.df.columns and 'WB ID' in self.df.columns:
            self.df['Wormbase ID'] = self.df['WB ID']
        if 'WB ID' not in self.df.columns and 'Wormbase ID' in self.df.columns:
            self.df['WB ID'] = self.df['Wormbase ID']
        if 'Wormbase ID' not in self.df.columns and 'WB Gene ID' in self.df.columns:
            self.df['Wormbase ID'] = self.df['WB Gene ID']
        if 'chrom' in self.df.columns and 'chrm' not in self.df.columns:
            self.df['chrm'] = self.df['chrom']
        if 'exp_reads' in self.df.columns:
            self.df.sort(columns='exp_reads', ascending=False, inplace=True)

    def to_csv(self, fname):
        self.df.to_csv(fname, sep='\t', index=False)

    def dict_by_first_instance(self, list_of_tups):
        if len(list_of_tups) < 1:
            return collections.defaultdict(str)
        if len(list_of_tups[0]) < 2:
            return collections.defaultdict(str)
        to_y = collections.defaultdict(type(list_of_tups[1]))
        for x in list_of_tups:
            if x[0] in to_y:
                continue
            to_y[x[0]] = x[1]

        return to_y

    def add_rip(self, fname='/opt/lib/fog3rip/TableS4_fog3_rip_targets.txt'):
        rip = pandas.read_csv(fname, sep='\t', index_col=False)
        rip['SAM rank'] = [ int(x) for x in rip['FOG-3 Rank'] ]
        to_sam = self.dict_by_first_instance(zip(rip['WB ID'].tolist(), rip['SAM rank'].tolist()))
        rip_wb = set(rip['WB ID'].tolist())
        self.df['Is RIP target?'] = [ x in rip_wb for x in self.df['Wormbase ID'].tolist() ]
        self.df['SAM rank'] = [ to_sam[x] for x in self.df['Wormbase ID'].tolist() ]

    def add_permutation_peaks(self, fname):
        peaks = pandas.read_csv(fname, sep='\t', index_col=False)
        to_peaks = self.dict_by_first_instance(zip(peaks['gene_id'].tolist(), peaks['height'].tolist()))
        clip_wb = set(peaks['gene_id'].tolist())
        print('clip wb IDs: {a}'.format(a=len(clip_wb)))
        self.df['Is CLIP target?'] = [ x in clip_wb for x in self.df['Wormbase ID'] ]
        self.df['CLIP height'] = [ to_peaks[x] for x in self.df['Wormbase ID'] ]

    def read_sp_vs_oo_as_programs(self, fname_oo='/opt/lib/ortiz/TableS1_oogenic.txt', fname_sp='/opt/lib/ortiz/TableS2_spermatogenic.txt', gtf_fname='/opt/lib/gtf_with_names_column.txt'):
        #print 'self.gene_name_col='
        #print self.gene_name_col
        if not hasattr(self, 'wbid_to_name'):
            self.transl(gtf_fname)
        oo_df = pandas.read_csv(open(fname_oo), sep='\t', index_col=False)
        sp_df = pandas.read_csv(fname_sp, sep='\t', index_col=False)
        sp_df['Program'] = [ re.sub('Spermatogenic and Oogenic', 'Oogenic and Spermatogenic', x) for x in sp_df['Program'].tolist() ]
        wb_program_oo = dict(zip(oo_df['Wormbase ID'].tolist(), oo_df['Program'].tolist()))
        wb_program_sp = dict(zip(sp_df['Wormbase ID'].tolist(), sp_df['Program'].tolist()))
        for k in set(wb_program_oo) & set(wb_program_sp):
            if wb_program_oo[k] != wb_program_sp[k]:
                print('Contradiction between %s and %s: %s (%s vs %s).' % (fname_oo,
                 fname_sp,
                 k,
                 wb_program_oo[k],
                 wb_program_sp[k]))

        wb_program_oo.update(wb_program_sp)
        self.program = wb_program_oo
        self.add_programs_to_df(self.df)
        print("""From Nobel et al., loaded {0} SP genes, {1} OO genes,
and {2} total GL genes.""".format(len(sp_df['Wormbase ID'].tolist()),
    len(oo_df['Wormbase ID'].tolist()), len(self.program)))
        #print self.df
        #sys.exit()
        if hasattr(self, 'clipdf'):
            self.add_programs_to_df(self.clipdf)
            print("CLIP peaks in a GL program: {0}".format(
            self.clipdf['Program'].value_counts()))

    def add_programs_to_df(self, _df):
        if 'Wormbase ID' not in _df.columns and self.gene_name_col in _df.columns:
            _df['Wormbase ID'] = [ self.name_to_wbid[x] for x in _df[self.gene_name_col].tolist() ]
        if 'Wormbase ID' in _df.columns:
            _df['Program'] = [ self.what_program(x) for x in _df['Wormbase ID'].tolist() ]
        
    def what_program(self, x):
        if x in self.program:
            return self.program[x]
        else:
            return ''

    def transl(self, gtfname='/opt/lib/gtf_with_names_column.txt'):
        gtf = pandas.read_csv(gtfname, sep='\t')
        self.wbid_to_name = collections.defaultdict(str)
        self.wbid_to_name.update(dict(zip(gtf['gene_id'].tolist(), gtf['gene_name'].tolist())))
        self.name_to_wbid = collections.defaultdict(str)
        self.name_to_wbid.update(dict(zip(gtf['gene_name'].tolist(), gtf['gene_id'].tolist())))
        self.name_to_biotype = dict(zip(gtf['gene_name'].tolist(), gtf['biotype'].tolist()))

    def programs_as_public_names(self):
        self.programs_as_pub_name = {}
        for wbid, prog in self.program.items():
            if wbid in self.wbid_to_name:
                self.programs_as_pub_name[self.wbid_to_name[wbid]] = \
                    self.program[wbid]
        return self.programs_as_pub_name
        
    def add_gene_name(self, gtfname):
        if not hasattr(self, 'wbid_to_name'):
            self.transl(gtfname)
        self.df['gene_name'] = [ self.wbid_to_name[x] for x in self.df.gene_id ]

    def add_wbid(self, gtfname='/opt/lib/gtf_with_names_column.txt'):
        if not hasattr(self, 'name_to_wbid'):
            self.transl(gtfname=gtfname)
        self.df['Wormbase ID'] = [ self.name_to_wbid[x] for x in self.df[self.gene_name_col].tolist() ]

    def add_biotype(self, gtfname):
        if not hasattr(self, 'name_to_biotype'):
            self.transl(gtfname)
        self.df['biotype'] = [ self.name_to_biotype[x] for x in self.df.gene_name ]

    def read_sp_vs_oo(self, fname='/opt/lib/ortiz/DESeq_genes_in_gonad.txt'):
        """/opt/lib/ortiz/DESeq_genes_in_gonad.txt
        """
        self.ortiz = pandas.read_csv(fname, sep='\t')
        fc = self.ortiz['Fold change (of normalized reads)']
        self.loci_to_fc = dict(zip(self.ortiz['Gene ID'].tolist(), self.ortiz['log2 fold change (of normalized reads)'].tolist()))
        self.gene_name_to_fc = dict(zip(self.ortiz['Gene name'].tolist(), self.ortiz['log2 fold change (of normalized reads)'].tolist()))
        self.gene_name_to_class = dict(zip(self.ortiz['Gene name'].tolist(), self.ortiz['Gene expression'].tolist()))
        self.loci_to_class = dict(zip(self.ortiz['Gene ID'].tolist(), self.ortiz['Gene expression'].tolist()))
        self.gene_name_to_oo_rpkm = dict(zip(self.ortiz['Gene name'].tolist(), self.ortiz['Expression in  fog-2 gonads (RPKM)'].tolist()))
        self.loci_to_oo_rpkm = dict(zip(self.ortiz['Gene ID'].tolist(), self.ortiz['Expression in  fog-2 gonads (RPKM)'].tolist()))
        self.gene_name_to_sp_rpkm = dict(zip(self.ortiz['Gene name'].tolist(), self.ortiz['Expression in fem-3 gonads (RPKM)'].tolist()))
        self.loci_to_sp_rpkm = dict(zip(self.ortiz['Gene ID'].tolist(), self.ortiz['Expression in fem-3 gonads (RPKM)'].tolist()))
        self.wb_id_to_oo_rpkm = dict(zip(self.ortiz['WormBase ID (WS240)'].tolist(), self.ortiz['Expression in  fog-2 gonads (RPKM)'].tolist()))
        self.wb_id_to_sp_rpkm = dict(zip(self.ortiz['WormBase ID (WS240)'].tolist(), self.ortiz['Expression in fem-3 gonads (RPKM)'].tolist()))

    def name_to_fc(self, name):
        if name in self.gene_name_to_fc:
            return self.gene_name_to_fc[name]
        elif name in self.loci_to_fc:
            return self.loci_to_fc[name]
        else:
            return np.nan

    def name_to_sp_oo_class(self, name):
        if name in self.gene_name_to_class:
            return self.gene_name_to_class[name]
        elif name in self.loci_to_class:
            return self.loci_to_class[name]
        else:
            return np.nan

    def name_to_oo_rpkm(self, name):
        if name in self.gene_name_to_oo_rpkm:
            return self.gene_name_to_oo_rpkm[name]
        elif name in self.loci_to_oo_rpkm:
            return self.loci_to_oo_rpkm[name]
        elif name in self.wb_id_to_oo_rpkm:
            return self.wb_id_to_oo_rpkm[name]
        else:
            return np.nan

    def name_to_sp_rpkm(self, name):
        if name in self.gene_name_to_sp_rpkm:
            return self.gene_name_to_sp_rpkm[name]
        elif name in self.loci_to_sp_rpkm:
            return self.loci_to_sp_rpkm[name]
        elif name in self.wb_id_to_sp_rpkm:
            return self.wb_id_to_sp_rpkm[name]
        else:
            return np.nan

    def annotate_sp_vs_oo(self):
        self.read_sp_vs_oo()
        print('%s annotate_sp_vs_oo, df.columns:' % self.name)
        print(self.df.columns)
        self.df['Fold change SP/OO GL'] = [ self.name_to_fc(x) for x in self.df[self.gene_name_col] ]
        self.df['Classification OO/SP GL'] = [ self.name_to_sp_oo_class(x) for x in self.df[self.gene_name_col] ]
        self.df['Oo. expression'] = [ self.name_to_oo_rpkm(x) for x in self.df[self.gene_name_col] ]
        self.df['Sp. expression'] = [ self.name_to_sp_rpkm(x) for x in self.df[self.gene_name_col] ]
        if 'fog_TGGC' in self.df.columns:
            t = zip(self.df['fog_TGGC'].tolist(), self.df['fog_GGTT'].tolist(), self.df['fog_GGCA'].tolist(), self.df['fog_CGGA'].tolist())
            self.df['height'] = [ np.mean(x) for x in t ]
            print(self.df['height'])

    def add_reads(self, ga=HTSeq.GenomicArray('auto'), name='unnamed array'):
        if self.df is None:
            print('Asked to set reads on an empty peaks object!')
            return False
        else:
            peaks = self.df
            ivs = zip(peaks['chrm'].tolist(), peaks['left'].tolist(), peaks['right'].tolist(), peaks['strand'].tolist())
            try:
                ivs = zip(peaks['chrm'].tolist(), peaks['left'].tolist(), peaks['right'].tolist(), peaks['strand'].tolist())
            except:
                print('Could not create ivs from %s.' % self.name)
                return False

            self.df[name] = [ self.get_val(ga, HTSeq.GenomicInterval(*iv)) for iv in ivs ]
            return

    def get_val(self, ga, iv):
        return np.max(np.fromiter(ga[iv], dtype=np.float))

    def add_location_from_integral(self, gtf, ga, use_name=None):
        if use_name is None:
            use_name = self.gene_name_col
        if 'chrm' not in self.df.columns:
            self.df['chrm'] = self.df['chrom']
        ivs = zip(self.df.chrm, self.df.left, self.df.right, self.df.strand, self.df[use_name].tolist())
        located_peaks = [ locatedPeak(*iv) for iv in ivs ]
        locations = [ x.add_location_from_integral(gtf, ga) for x in located_peaks ]
        self.df['location'] = locations
        return

    def add_seqs(self, sequences):
        seq = [ seq_from_iv(tup[0], tup[1], tup[2], tup[3], sequences) for tup in zip(self.df.chrm, self.df.left, self.df.right, self.df.strand) ]
        self.df['seq'] = seq

    def merge_ranges_with_other_peaksList(self, _p, use_name='gene_name'):
        assert type(_p) == type(self)
        ivs_a_tups = zip(self.df.chrm, self.df.left, self.df.right, self.df.strand, self.df[use_name].tolist())
        ivs_b_tups = zip(_p.df.chrm, _p.df.left, _p.df.right, _p.df.strand, _p.df[use_name].tolist())
        ivs_a = collections.defaultdict(list)
        ivs_b = collections.defaultdict(list)
        for tup in ivs_a_tups:
            ivs_a[tup[0], tup[3]].append(tup)

        for tup in ivs_b_tups:
            ivs_b[tup[0], tup[3]].append(tup)

        merged = collections.defaultdict(list)
        for k in ivs_a:
            if k in ivs_b:
                a, b, c = overlaps_in_lists(ivs_a[k], ivs_b[k])
                merged[k] = a + b + c
            else:
                merged[k] = ivs_a[k]

        for k in [ x for x in ivs_b.keys() if x not in ivs_a.keys() ]:
            merged[k] = ivs_b[k]

        merged_list = []
        for k in merged:
            merged_list = merged_list + overlaps_in_one_list(merged[k])

        return merged_list


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
