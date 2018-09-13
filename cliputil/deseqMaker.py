
import pandas, os, sys, collections, glob, re

top_dir = '/groups/Kimble/Common/fbf_celltype/'  # Server.
top_dir = '/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/v2/combined_filtered/'

fname_to_label = {
top_dir + '/old_fbf2.txt': 'old_fbf2',
top_dir + '/old_fbf1.txt': 'old_fbf1',
top_dir + '/oo_both.txt': 'oo_both',
top_dir + '/oo_fbf1.txt': 'oo_fbf1',
top_dir + '/oo_fbf2.txt': 'oo_fbf2',
top_dir + '/sp_both.txt': 'sp_both',
top_dir + '/sp_fbf1.txt': 'sp_fbf1',
top_dir + '/sp_fbf2.txt': 'sp_fbf2',
top_dir + '/old_fbf1_to_fbf2_n2.txt': 'old_fbf1_to_fbf2_n2',
}


def _mk(d):
    if not os.path.exists(d): os.system('mkdir {}'.format(d))


class deseqMaker(object):

    def yield_sp_oo(self, indir):
        for f in glob.glob(indir + '/*'):
            f = os.path.basename(f)
            # Only sp/oo.
            if not (
                re.search('exp_fbf.*_oo', f) or re.search('exp_fbf.*_sp', f)):
                continue
            yield f
            
    def remove_underscores(self, _f):
            m = re.match('.*lane(\d+)_rt(\d+)_counts.*', _f)
            if m is not None:
                new_seq = "lane{}rt{}".format(m.group(1), m.group(2))
                _f = re.sub('lane\d+_rt\d+_counts', new_seq, _f)
    
    def fix_count_filenames_for_deseq_of_sp_vs_oo_clip_11_reps(
        self, indir='counts/', outdir='counts_alt_fnames/'):
        _mk(outdir)
        to_move = []
        print("----")
        for f in self.yield_sp_oo(indir):
            # Only individual replicates.
            if re.search('and', f): continue
            # Skip the tiny dataset
            if re.search('fbf1_oo_lane2_rt1', f): continue
            print(f,) 
            _f = re.sub('exp_fbf1_oo', 'oo-fbf1', f)
            _f = re.sub('exp_fbf2_oo', 'oo-fbf2', _f)
            _f = re.sub('exp_fbf1_sp', 'sp-fbf1', _f)
            _f = re.sub('exp_fbf2_sp', 'sp-fbf2', _f)
            print(">" + _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)
        
    def cp_files(self, to_move, indir, outdir):
        print("Copying {} files from {} to {}".format(
            len(to_move), indir, outdir))
        for f, _f in to_move:
            in_f = os.path.join(indir, f)
            out_f = os.path.join(outdir, _f)
            print(('cp {} {}'.format(in_f, out_f)))
            os.system('cp {} {}'.format(in_f, out_f))

    def mk_count_dir_for_deseq_of_sp_vs_oo_clip_6_reps(
        self, indir='counts/', outdir='counts_6_reps/'):
        _mk(outdir)
        to_move = []
        print("----")
        for f in self.yield_sp_oo(indir):
            # Only combined replicates.
            if not re.search('and', f): continue
            print(f)
            _f = re.sub('exp_fbf_oo_', 'oo-', f)
            _f = re.sub('exp_fbf_sp_', 'sp-', _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)
        
    def mk_count_dir_for_deseq_of_sp_vs_oo_clip_6_reps_from_combined_counts_file(
        self, infile='combined_counts.txt', outdir='counts_6_reps/'):
        _mk(outdir)
        to_move = []
        
        print("----")
        df = pandas.read_csv(infile, sep='\t', index_col=False)
        
        out_cols = [x for x in df.columns if re.search('exp_fbf_sp.*', x) or re.search('exp_fbf_oo.*', x)]
        df.index = df['gene']
        #df = df[['gene'] + out_cols].copy()
        print(df)
        
        if '_ambiguous' in df.index:
            print('amb')
            df.drop('_ambiguous', inplace=True)
            df.drop('_no_feature', inplace=True)
            
        for col in out_cols:
            to_out_df = df[['gene', col]]
            
            _f = re.sub('exp_fbf_oo_', 'oo-', col)
            _f = re.sub('exp_fbf_sp_', 'sp-', _f)
            
            to_out_df[col] = [int(x) for x in to_out_df[col]]
            to_out_df.to_csv(outdir + '/' + _f, sep='\t', index=False, header=False)
        
    def mk_count_dir_for_deseq_of_lt_fbf1_vs_fbf2(
        self, indir='counts/', outdir='counts_lt_fbf/'):
        _mk(outdir)
        to_move = []
        for f in [
            'exp_fbf1_CGGA_counts.txt',
            'exp_fbf1_GGTT_counts.txt',
            'exp_fbf1_TGGC_counts.txt',
            'exp_fbf2_CGGA_counts.txt',
            'exp_fbf2_GGTT_counts.txt',
            'exp_fbf2_TGGC_counts.txt']:
            print(f)
            _f = re.sub('exp_fbf1_', 'fbf1-', f)
            _f = re.sub('exp_fbf2_', 'fbf2-', _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)

    def mk_count_dir_for_deseq_of_lt_fbf1_vs_ht_fbf1(
        self, indir='counts/', outdir='counts/counts_lt_fbf1_vs_ht_fbf1/'):
        _mk(outdir)
        to_move = []
        for f in [
            'exp_fbf1_CGGA_counts.txt',
            'exp_fbf1_GGTT_counts.txt',
            'exp_fbf1_TGGC_counts.txt',
            'exp_fbf1_oo_lane2_rt1_counts.txt',
            'exp_fbf1_oo_lane2_rt6_counts.txt',
            'exp_fbf1_oo_lane2_rt9_counts.txt']:
            print(f)
            _f = re.sub('exp_fbf1_oo_', 'fbf1_ht-', f)
            _f = re.sub('exp_fbf1_', 'fbf1_lt-', _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)

    def mk_count_dir_for_deseq_of_lt_fbf1_vs_ht_fbf(
        self, indir='counts/', outdir='counts/counts_lt_fbf1_vs_ht_fbf/'):
        _mk(outdir)
        to_move = []
        for f in [
            'exp_fbf1_CGGA_counts.txt',
            'exp_fbf1_GGTT_counts.txt',
            'exp_fbf1_TGGC_counts.txt',
            'exp_fbf_oo_rt11_and_6_counts.txt',
            'exp_fbf_oo_rt2_and_13_counts.txt',
            'exp_fbf_oo_rt9_and1_counts.txt']:
            print(f)
            _f = re.sub('exp_fbf_oo_', 'fbf1_ht-', f)
            _f = re.sub('exp_fbf1_', 'fbf1_lt-', _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)

    def mk_count_dir_for_deseq_of_lt_fbf2_vs_ht_fbf(
        self, indir='counts/', outdir='counts/counts_lt_fbf2_vs_ht_fbf/'):
        _mk(outdir)
        to_move = []
        for f in [
            'exp_fbf2_CGGA_counts.txt',
            'exp_fbf2_GGTT_counts.txt',
            'exp_fbf2_TGGC_counts.txt',
            'exp_fbf_oo_rt11_and_6_counts.txt',
            'exp_fbf_oo_rt2_and_13_counts.txt',
            'exp_fbf_oo_rt9_and1_counts.txt']:
            print(f)
            _f = re.sub('exp_fbf_oo_', 'fbf1_ht-', f)
            _f = re.sub('exp_fbf2_', 'fbf2_lt-', _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)

    def mk_count_dir_for_deseq_of_lt_fbf1_and_2_vs_ht_fbf(
        self, indir='counts/', outdir='counts/counts_lt_fbf1_and_2_vs_ht_fbf/'):
        _mk(outdir)
        to_move = []
        for f in [
            'exp_fbf1_CGGA_counts.txt',
            'exp_fbf1_GGTT_counts.txt',
            'exp_fbf1_TGGC_counts.txt',
            'exp_fbf2_CGGA_counts.txt',
            'exp_fbf2_GGTT_counts.txt',
            'exp_fbf2_TGGC_counts.txt',
            'exp_fbf_oo_1_counts.txt',
            'exp_fbf_oo_2_counts.txt',
            'exp_fbf_oo_3_counts.txt']:
            print(f)
            _f = re.sub('exp_fbf_oo_', 'fbf_ht-', f)
            _f = re.sub('exp_fbf1_', 'fbf_lt-1', _f)
            _f = re.sub('exp_fbf2_', 'fbf_lt-2', _f)
            to_move.append((f, _f))
        self.cp_files(to_move, indir, outdir)
        
    def mk_count_dir_for_deseq_of_lt_fbf1_and_2_vs_ht_fbf_from_combined_counts_file(
        self, infile='combined_counts.txt', outdir='counts_lt_fbf1_and_2_vs_ht_fbf/'):
        _mk(outdir)
        to_move = []
        
        print("----")
        df = pandas.read_csv(infile, sep='\t', index_col=False)
        
        out_cols = [
            #'old_fbf1',
            #'old_fbf2',
            #'old_fbf1_to_fbf2_n2',
            #'exp_fbf'
            'exp_fbf1_CGGA_counts.txt',
            'exp_fbf1_GGTT_counts.txt',
            'exp_fbf1_TGGC_counts.txt',
            'exp_fbf2_CGGA_counts.txt',
            'exp_fbf2_GGTT_counts.txt',
            'exp_fbf2_TGGC_counts.txt',
            'exp_fbf_oo_1_counts.txt',
            'exp_fbf_oo_2_counts.txt',
            'exp_fbf_oo_3_counts.txt',
            #'exp_fbf_oo_rt11_and_6_counts.txt',
            #'exp_fbf_oo_rt2_and_13_counts.txt',
            #'exp_fbf_oo_rt9_and1_counts.txt'
            ]
        df.index = df['gene']
        #df = df[['gene'] + out_cols].copy()
        print(df)
        
        if '_ambiguous' in df.index:
            print('amb')
            df.drop('_ambiguous', inplace=True)
            df.drop('_no_feature', inplace=True)
            
        for col in out_cols:
            to_out_df = df[['gene', col]]
            
            _f = re.sub('exp_fbf_oo_', 'fbf_ht-', col)
            _f = re.sub('exp_fbf1_', 'fbf_lt-1', _f)
            _f = re.sub('exp_fbf2_', 'fbf_lt-2', _f)
            
            to_out_df[col] = [int(x) for x in to_out_df[col]]
            to_out_df.to_csv(outdir + '/' + _f, sep='\t', index=False, header=False)


if __name__ == '__main__':
    a = deseqMaker()
    #a.mk_count_dir_for_deseq_of_lt_fbf1_vs_ht_fbf1()
    #a.mk_count_dir_for_deseq_of_lt_fbf1_vs_ht_fbf()
    #a.mk_count_dir_for_deseq_of_lt_fbf2_vs_ht_fbf()
    a.mk_count_dir_for_deseq_of_lt_fbf1_and_2_vs_ht_fbf_from_combined_counts_file()
    #a.mk_count_dir_for_deseq_of_lt_fbf1_vs_fbf2()
    #a.fix_count_filenames_for_deseq_of_sp_vs_oo_clip_11_reps(
    #    indir='counts/', outdir='counts_11_reps/')
    #a.mk_count_dir_for_deseq_of_sp_vs_oo_clip_6_reps()
    a.mk_count_dir_for_deseq_of_sp_vs_oo_clip_6_reps_from_combined_counts_file()

    
rscript = """

library(DESeq2)
library("xlsx")

# DESeq expects tab separated input files with no header, and genes in the left column.

get_res = function(directory) {
        sf = grep('txt', list.files(directory), value=TRUE)
        sampleCondition = sub("(.*)-.*", "\\1", sf)
        sampleTable = data.frame(sampleName=sf, fileName=sf, condition=sampleCondition)
        ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,  directory=directory, design= ~ condition)
        ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 50, ]
        ddsHTSeq = DESeq(ddsHTSeq)
        res = results(ddsHTSeq)
        resO = res[order(res$padj),]
        return(resO)
}


res_wt = get_res('counts_lt_fbf1_and_2_vs_ht_fbf/')
write.table(res_wt, file='tables/lt_fbf1_and_2_vs_ht_fbf_deseq.txt', quote=FALSE, sep='\t')
write.xlsx(res_wt, file='tables/lt_fbf1_and_2_vs_ht_fbf_deseq.xlsx', sheetName = "Sheet1", 
  col.names = TRUE, row.names = TRUE, append = FALSE)

res_wt = get_res('counts_6_reps/')
write.table(res_wt, file='tables/6_reps_sp_vs_oo.txt', quote=FALSE, sep='\t')
write.table(res_wt, file='tables/File S4 unedited Deseq2 SP vs OO.txt', quote=FALSE, sep='\t')
write.xlsx(res_wt, "tables/File S4 raw DESeq2 output.xlsx", sheetName = "Sheet1", 
  col.names = TRUE, row.names = TRUE, append = FALSE)
  


"""
