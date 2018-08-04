#!/usr/bin/python
"""
Annotate the */*/null_hyp_4.txt files with reads,
 and write copies with fewer columns.
 
 Combine peak calling info and add reads-in-peaks/annotations.
 ---
 
 This will create combined_filtered/ and combined_unfiltered/ directories:
 
 ```bash
 python cliputil/filter.py -i individual_clip
 # Outputs combined_unfiltered/.
 # Annotate these peaks with seq, FBE and biotype.
 python cliputil/annotate_peaks.py -i combined_unfiltered/
 ```
 
 This adds reads-in-peak numbers for both normalized and unnormalized numbers.
 
 Filtering
 ---
 
 Filter by the hardcoded ratios:
 ```bash
 python cliputil/filter.py -f -i combined_unfiltered/
 ```
 
 The -f, -r, and -p parameters are mutually exclusive.
"""
import argparse
from __init__ import *
import add_reads_to_peaks
from peaks import *


class peakSet(object):

    def __init__(self):
        self._peaks = []
        pass

    def add_peak(self, peak_obj):
        self._peaks.append(peak_obj)

    def cut_columns(self, remove_those_with_patterns=[
                '_cor$', '_rej$', '_poisson',
                '_gene_norm$', '_mu$', 'pvalues',
                '\Aexons$', '_nb$', '_std$', '\Alocal$',
                'clip_local_norm', 'rna_seq_local_norm']):
     
        for _p in self._peaks:
            _p.cut_columns(
                remove_those_with_patterns=remove_those_with_patterns)

    def write_tables(self, dirname):
        pass

def filename_and_label_mappings():
    
    unfiltered = './combined_unfiltered/'
    filtered = './combined_filtered/'

    for dirname in [unfiltered, filtered]:
        if os.path.exists(dirname):
            _in = input("Output directory {} already exists. Overwrite? [Y/N]: ".format(dirname))
            if _in[0].upper() == 'Y':
                pass
            else:
                raise IOError("Aborting because output folder exists.")
    
    fname_to_label = {
# Uncommented before initial submission
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/old_fbf2.txt': 'old_fbf2',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/old_fbf1.txt': 'old_fbf1',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/oo_both.txt': 'oo_both',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/oo_fbf1.txt': 'oo_fbf1',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/oo_fbf2.txt': 'oo_fbf2',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/sp_both.txt': 'sp_both',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/sp_fbf1.txt': 'sp_fbf1',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/sp_fbf2.txt': 'sp_fbf2',
#'/groups/Kimble/Common/fbf_celltype/combined_unfiltered/old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/old_fbf1_to_fbf2_n2/09_39_PM/peaks/combined_09_39_PM/null_hyp_4.txt': 'old_fbf1_to_fbf2_n2',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/old_fbf1/11_01_PM/peaks/combined_11_01_PM/null_hyp_4.txt': 'old_fbf1',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/old_fbf2/old_fbf2/peaks/combined_old_fbf2/null_hyp_4.txt': 'old_fbf2',

# Commented out before initial submission.
#'/groups/Kimble/Common/fbf_celltype/individual_clip/oo_both/08_51_PM/peaks/combined_08_51_PM/null_hyp_4.txt': 'oo_both',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/oo_fbf1/08_17_PM/peaks/combined_08_17_PM/null_hyp_4.txt': 'oo_fbf1',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/oo_fbf2/08_18_PM/peaks/combined_08_18_PM/null_hyp_4.txt': 'oo_fbf2',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/sp_both/09_08_PM/peaks/combined_09_08_PM/null_hyp_4.txt': 'sp_both',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/sp_fbf1/08_19_PM/peaks/combined_08_19_PM/null_hyp_4.txt': 'sp_fbf1',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/sp_fbf2/08_20_PM/peaks/combined_08_20_PM/null_hyp_4.txt': 'sp_fbf2'}

# Uncommented before initial submission.
#'/groups/Kimble/Common/fbf_celltype/recall_using_old_sam_w_multimappers_removed/oo_both/08_51_PM/peaks/combined_08_51_PM/null_hyp_4.txt': 'oo_both',
#'/groups/Kimble/Common/fbf_celltype/recall_using_old_sam_w_multimappers_removed/oo_fbf1/12_15_AM/peaks/combined_12_15_AM/null_hyp_4.txt': 'oo_fbf1',
#'/groups/Kimble/Common/fbf_celltype/recall_using_old_sam_w_multimappers_removed/oo_fbf2/12_15_AM/peaks/combined_12_15_AM/null_hyp_4.txt': 'oo_fbf2',
#'/groups/Kimble/Common/fbf_celltype/recall_using_old_sam_w_multimappers_removed/sp_both/09_08_PM/peaks/combined_09_08_PM/null_hyp_4.txt': 'sp_both',
#'/groups/Kimble/Common/fbf_celltype/recall_using_old_sam_w_multimappers_removed/sp_fbf1/12_14_AM/peaks/combined_12_14_AM/null_hyp_4.txt': 'sp_fbf1',
#'/groups/Kimble/Common/fbf_celltype/recall_using_old_sam_w_multimappers_removed/sp_fbf2/12_14_AM/peaks/combined_12_14_AM/null_hyp_4.txt': 'sp_fbf2',

# Commented out before initial submission.
#'/groups/Kimble/Common/fbf_celltype/individual_clip/spermatogenic/bak_fbf_sp/peaks/combined_fbf_sp/null_hyp_4.txt': 'bak_spermatogenic',
#'/groups/Kimble/Common/fbf_celltype/individual_clip/spermatogenic/fbf_sp/peaks/combined_fbf_sp/null_hyp_4.txt': 'spermatogenic',
#'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/oo_fbf1/05_40_PM/peaks/combined_05_40_PM/null_hyp_4.txt': 'oo_fbf1',
#'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/sp_fbf1/05_16_PM/peaks/combined_05_16_PM/null_hyp_4.txt': 'sp_fbf1',
#'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/sp_fbf/05_05_PM/peaks/combined_05_05_PM/null_hyp_4.txt': 'sp_fbf',
#'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/oo_fbf/05_06_PM/peaks/combined_05_06_PM/null_hyp_4.txt': 'oo_fbf',
#'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/oo_fbf2/07_04_PM/peaks/combined_07_04_PM/null_hyp_4.txt': 'oo_fbf2',
#'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/sp_fbf2/04_13_PM/peaks/combined_04_13_PM/null_hyp_4.txt': 'sp_fbf2',

'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/old_fbf2.txt': 'old_fbf2',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/old_fbf1.txt': 'old_fbf1',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/oo_both.txt': 'oo_both',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/oo_fbf1.txt': 'oo_fbf1',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/oo_fbf2.txt': 'oo_fbf2',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/sp_both.txt': 'sp_both',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/sp_fbf1.txt': 'sp_fbf1',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/sp_fbf2.txt': 'sp_fbf2',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/sp_fbf.txt': 'sp_fbf',
'/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/individual_clip/combined_unfiltered/oo_fbf.txt': 'oo_fbf',
}
    
    label_to_unfiltered_output_file = dict([
        (x, '{0}/{1}.txt'.format(unfiltered, x)) for x in list(fname_to_label.values())])
        
    label_to_filtered_output_file = dict([
        (x, '{0}/{1}.txt'.format(filtered, x)) for x in list(fname_to_label.values())])
        
    for x in [unfiltered, filtered]:
        if not os.path.exists(x):
            os.system('mkdir {0}'.format(x))
    
    return unfiltered, filtered, fname_to_label, label_to_unfiltered_output_file, label_to_filtered_output_file


def run(args):
    
    (unfiltered, filtered, fname_to_label, label_to_unfiltered_output_file,
     label_to_filtered_output_file) = filename_and_label_mappings()
     
    # For the null_hyp_4.txt files found in args.input, make peaks objects, add bedgraph data,
    # and return a list of the peaks objects.
    peak_objs = add_reads_to_peaks.run(args)

    print("Loaded all peaks files.")
    
    # Clean.
    for _p in peak_objs:

        print('~ * ' * 10)
        print('{} (filepath {}):\ncutting columns, writing to a table {}'.format(
            fname_to_label[os.path.realpath(_p.file)], 
            _p.file,
            label_to_unfiltered_output_file[fname_to_label[os.path.realpath(_p.file)]]))
        
        _p.cut_columns(
            remove_those_with_patterns=[
                '_cor$', '_rej$', '_poisson',
                '_gene_norm$', '_mu$', 'pvalues',
                '\Aexons$', '_nb$', '_std$', '\Alocal$',
                'clip_local_norm', 'rna_seq_local_norm'])
                
        _p.write_table('test.txt')
        
        if os.path.realpath(_p.file) not in fname_to_label:
            print("{0} not in fname_to_label {1}.".format(_p.file, fname_to_label))
            continue
        
        _p.label = fname_to_label[os.path.realpath(_p.file)]
        
        print(os.path.realpath(_p.file))
        print(" .____. ")
        print(_p.label)
        
        _p.write_table(label_to_unfiltered_output_file[_p.label])

    for _p in peak_objs:
        
        if not hasattr(_p, 'label'):
            print("Skipping the generation of a filtered {0}, as it has no label attr".format(str(_p)))
            continue
        
        filt_p = _p.get_filtered_obj(col='unnorm_ratio', cutoff=10)
        filt_p = filt_p.get_filtered_obj(col='ratio', cutoff=2)
        filt_p.write_table(label_to_filtered_output_file[_p.label])


def report(dirname):

    for f in glob.glob(dirname + '/*txt'):
        f = os.path.realpath(f)
        df = pandas.read_csv(f, sep='\t')

        print('{}:'.format(f))

        for col in ['unnorm_ratio', 'ratio']:
            if col not in df.columns:
                print("Missing expected column {} in {}".format(
                    col, os.path.basename(f)))
                continue

            _min = np.min(df[col].tolist())
            _max = np.max(df[col].tolist())

            print("{}, {}: minimum {} maximum {}".format(
                os.path.basename(f), col, _min, _max))

def read_dir(dirname):
    
    print(dirname)
    
    (unfiltered, filtered, fname_to_label, label_to_unfiltered_output_file,
     label_to_filtered_output_file) = filename_and_label_mappings()
     
    for f in glob.glob(dirname + '/*txt'):
        f = os.path.realpath(f)
        
        if re.search('old_fbf1.txt', f):
            cutoff = 4#4
            unnorm_cutoff = 4#6
        
        elif re.search('old_fbf2.txt', f):
            cutoff = 2#30
            unnorm_cutoff = 30#30
        
        elif re.search('old_fbf1_to_fbf2_n2.txt', f):
            cutoff = 3#40
            unnorm_cutoff = 40#60
        
        #elif re.search('sp_both.txt', f):
        #    cutoff = 5#2#5
        #    unnorm_cutoff = 60#10#50

        elif re.search('sp_fbf.txt', f):
            cutoff = 5#5
            unnorm_cutoff = 20#50     

#        elif re.search('oo_both.txt', f):
#            cutoff = 2#3
#            unnorm_cutoff = 5#10

        elif re.search('oo_fbf.txt', f):
            cutoff = 2.5#3
            unnorm_cutoff = 20#5#10

        elif re.search('oo_fbf\d.txt', f):
            cutoff = 2
            unnorm_cutoff = 10
        
        else:
            cutoff = 2.5#3
            unnorm_cutoff = 20#20
        
        print(f)
        print(cutoff, unnorm_cutoff)
        
        _p = peaks(file=f)
        _p.label = fname_to_label[os.path.realpath(_p.file)]
        
        filtp = _p.get_filtered_obj(col='ratio', cutoff=cutoff)
        
        if 'unnorm_ratio' in _p.data.columns:
            filtp = filtp.get_filtered_obj(
                col='unnorm_ratio', cutoff=unnorm_cutoff)
        
        filtp.write_table(label_to_filtered_output_file[_p.label])

def find_lowest_cutoffs_of_positive_controls(dirname):

    positives = "gld-1,htp-1,htp-2,mpk-1,him-3,fbf-1,lip-1,syp-2,fbf-2,fog-1,fem-3,syp-3,gld-3,fog-3,egl-4".split(',')

    positives.remove('egl-4')
    positives.remove('fog-3')

    for f in glob.glob(dirname + '/*txt'):
        f = os.path.realpath(f)

        df = pandas.read_csv(f, sep='\t')
        positive_rows = df[[(x in positives) for x in df.gene_name]]
        positive_rows.sort_values(by=['ratio'], inplace=True, ascending=False)
        positive_rows.drop_duplicates(subset='gene_name', inplace=True)
        print('>{}'.format(os.path.basename(f)))
        print(positive_rows[['gene_name', 'ratio', 'unnorm_ratio']])
        
if __name__ == '__main__':          
    
    p = argparse.ArgumentParser(__doc__)

    p.add_argument(
        '-i', '--input', help='Directory of peaks files, which are searched \
recursively.')

    #p.add_argument('-c', '--config', help='config.ini file',
    #               default='auto.ini')
    p.add_argument('-o', '--output', help='(Optional) A directory to output \
peaks with coverages added and columns simplified.',
                   default=None)

    p.add_argument('-f', '--filter',
                   help='Just filter the input dir, not recursively.',
                   action='store_true', default=False)

    p.add_argument('-r', '--report',
                   help='Report on the cutoffs observed in a directory. Do not filter.',
                   action='store_true', default=False)

    p.add_argument('-p', '--positives',
                   help='Report on the cutoffs on postibe control genes. Do not filter.',
                   action='store_true', default=False)


    args = p.parse_args()
    args.input = os.path.realpath(args.input)
    #lib = config(args.config)

    if args.report:
        report(args.input)

    elif args.positives:
        find_lowest_cutoffs_of_positive_controls(args.input)
    
    elif args.filter:
        read_dir(args.input)
    
    else:
        run(args)#, lib)
    