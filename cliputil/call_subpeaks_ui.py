import sys
sys.path.insert(
    0, 'fbf_clip/analysis/src/feature_locations/')
import subpeaks_ui
import subpeaks
import config
import determine_feature_locations

args = subpeaks_ui.read_args()
lib = config.config(filepath=args.config_ini)
#if ('coverage_wigs' not in lib) and ('bedgraphs_folder' in lib):
#    lib['coverage_wigs'] = lib['bedgraphs_folder']
lib['coverage_wigs'] = '/groups/Kimble/Common/fbf_celltype/header_wigs/'
lib['gtf_one_txpt_per_gene'] = '/opt/lib/gtf_one_txtp_per_gene.txt'
lib['gtf_pickle'] = '/opt/lib/gtf_as_dict.p'
lib['feat_loc_dir'] = 'data/'
lib['txpt_obj_path'] = 'data/feat.p'
lib['chr_sizes'] = '/opt/lib/chr_sizes.txt'
pk_list = [
    'old_fbf1_to_fbf2_n2',
    'old_fbf2',
    'oo_both',
    #'sp_fbf1',
    #'sp_fbf2',
    'sp_both']
to_wig = {
    'old_fbf1_to_fbf2_n2': ('exp_fbf1_+.wig', 'exp_fbf1_-.wig'),
    'old_fbf2': ('exp_fbf2_+.wig', 'exp_fbf2_-.wig'),
    'oo_both': ('exp_fbf_oo_+.wig', 'exp_fbf_oo_-.wig'),
    'sp_both': ('exp_fbf_sp_+.wig', 'exp_fbf_sp_-.wig')
    }
for label in pk_list:
    args.input = 'combined_filtered/{0}.txt'.format(label)
    lib['bedgraph_exp_plus'] = lib['coverage_wigs'] + \
                               '/' + to_wig[label][0]
    lib['bedgraph_exp_minus'] = lib['coverage_wigs'] + \
                                '/' + to_wig[label][1]
    ga = subpeaks.get_bedgraph(
        bedgraphs_folder=lib['coverage_wigs'],
        args=args, lib=lib)
    (peaks, txpts, chr_lens) = determine_feature_locations.get_data(args, lib=lib)
    (sequences, rc_sequences, chr_lens) = subpeaks.get_sequences(lib)
    subpeaks.run(
            args.input, ga, lib, txpts, sequences, rc_sequences,
            chr_lens, args=args)
