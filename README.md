Comparing FBF iCLIP in differnt cell types
======

Code for a paper on FBF in the gendered transcriptome. 

Data download
======

Data was downloaded from archive SRP077050 using fastq-dump in the NCBI toolkit.

Fastq files are 51 nt long reads split by barcode (the four bases read_str[4:7])
by the script clip-preprocess/split_by_barcode.py.

From CLIP-preprocess readme:
```bash
python src/split_by_barcode.py <fastq>

Several steps are done to get to uncollapsed bed files. The script split_fastq_to_bed.py will move barcodes to read names, remove linkers, map and take the best map for each read, and finally create bams and beds from the sam output by STAR. The sequence of folders is:

input dir -> temp_adapter_in_name/ -> temp_fastq/ -> temp_clipped/ -> sams/ -> temp_sam_collapse_multi_map/ -> sams/ -> bed_uncollapsed/

bed_collapsed is the only output needed for everything else.

This is run by:

$ python split_fastq_to_bed.py -i input_dir

That script does the following:

The barcodes were moved to the read header:

$ python move_barcode_to_name_in_fastq.py <fastq directory/>

The three prime linker was removed:

$ python clip_adapter.py -t -i <in_dir> -o <out_dir>

The RT primer was removed:

$ python clip_adapter.py -r -i <in_dir> -o <out_dir>

Reads were mapped to the genome:

$ python map_with_star.py -i <in_dir>
```

For CSEQ parameters, the sam headers are (for example):

```bash
@HD	VN:1.4
@SQ	SN:I	LN:15072434
@SQ	SN:II	LN:15279421
@SQ	SN:III	LN:13783801
@SQ	SN:IV	LN:17493829
@SQ	SN:V	LN:20924180
@SQ	SN:X	LN:17718942
@SQ	SN:MtDNA	LN:13794
@PG	ID:STAR	PN:STAR	VN:STAR_2.4.2a	CL:/home/users/dfporter/dfp_scratch/very_temp/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR   --runThreadN 16   --genomeDir /home/users/dfporter/dfp_scratch/very_temp/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/   --readFilesIn ../temp_clipped/SRR3715268.fastq      --outFileNamePrefix SRR3715268_   --outSAMunmapped Within   --outFilterMultimapNmax 3   --outFilterMismatchNmax 2   --winAnchorMultimapNmax 10000   --seedSearchStartLmax 6   --alignIntronMax 1   --alignTranscriptsPerReadNmax 50000   --alignEndsType Local   --sjdbGTFfile /home/users/dfporter/dfp_scratch/very_temp/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/sjdb.txt   --sjdbGTFtagExonParentTranscript transcript_id
@CO	user command line: /home/users/dfporter/dfp_scratch/very_temp/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR --alignIntronMax 1 --sjdbGTFfile /home/users/dfporter/dfp_scratch/very_temp/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/sjdb.txt --genomeDir /home/users/dfporter/dfp_scratch/very_temp/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/ --readFilesIn ../temp_clipped/SRR3715268.fastq --outSAMunmapped Within --outFilterMultimapNmax 3 --outFilterMismatchNmax 2 --seedSearchStartLmax 6 --winAnchorMultimapNmax 10000 --alignEndsType Local --sjdbGTFtagExonParentTranscript transcript_id --runThreadN 16 --outFileNamePrefix SRR3715268_ --alignTranscriptsPerReadNmax 50000
```

These are converted to bed files of format:

```
I	928	959	SRR3715251.4196778#AATGA	255	+	0	31M	*	0	0	TGTCTGCCCAAATATAGGCAACCAAAAATAA	IIIIIIIIIIIIIIIIIIHIIIIIIIIHIIINH:i:1	HI:i:1	AS:i:30	nM:i:0
```

map_with_star.py -c is used to subset the bed columns.
It has some parameters for specific usages:

```bash
$ python ../../../clip-preprocess/map_with_star.py -h
usage: map_with_star.py [-h] [-i INPUT_DIR] [-s SAMS_INPUT]
                        [-c SUBSET_BED_COLUMNS]

Always writes to a ./sam/ directory.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
  -s SAMS_INPUT, --sams_input SAMS_INPUT
                        Don't map, just convert this sams directory to bed.
  -c SUBSET_BED_COLUMNS, --subset_bed_columns SUBSET_BED_COLUMNS
                        Just subset the column of the bed files in this
                        directory to a format for CIMS.
```

The bed files now have the format:
```bash
$ head ../../../subset_bed_columns/SRR3715251.bed 
I	928	959	SRR3715251.4196778#AATGA	0	+
```

This is the correct format for duplicate collapsing.
The commands for duplicate collapsing are generated:

```bash
$ python clip-preprocess/write_commands_for_collapsing.py
```

The collapsed bed files then need their filenames translated from SRR* numbers:
```bash
$ python clip-preprocess/translate_filenames.py -i <input folder> -o <output folder>
```

At this point, bedgraph files can be created:

```bash
$ python ../../../clip-preprocess/bed_to_wig.py -h
usage: bed_to_wig.py [-h] [-b INPUT_BED] [-o OUTPUT_BEDGRAPH_NORM]
                     [-u OUTPUT_BEDGRAPH_UNNORM]

optional arguments:
  -h, --help            show this help message and exit
  -b INPUT_BED, --input_bed INPUT_BED
                        Folder of bed files.
  -o OUTPUT_BEDGRAPH_NORM, --output_bedgraph_norm OUTPUT_BEDGRAPH_NORM
                        Output folder for normalized bedgraphs.
  -u OUTPUT_BEDGRAPH_UNNORM, --output_bedgraph_unnorm OUTPUT_BEDGRAPH_UNNORM
                        Output folder for unnormalized bedgraphs.
```

Individual CLIP analysis folders for peak calling are then created.

```bash
$ python ../../../clip-preprocess/create_individual_clip_folders.py -h
../../../clip-preprocess/create_individual_clip_folders.py
usage: create_individual_clip_folders.py [-h] [-i INPUT] [-o OUTPUT] [-r]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input folder.
  -o OUTPUT, --output OUTPUT
                        Output folder.
  -r, --run_commands    Run commands, rather than just print them.
```

Peaks are then called using cwt-callpeaks for each set individually.

Analysis
======

This analysis begins with the filtered, separately called oo/sp/20deg peaks
 files and associated wig files.

 It also includes library files of DESeq2 output - quite a few library files,
 in fact.

Essentially, CLIP data constitutes:

1. Peaks files in combined_filtered/
2. bedgraphs, norm (per million reads)
3. bedgraphs, unnorm (raw read counts)
4. bed files

The second type of input into this analysis is a variety of library files,
 placed under /opt/lib/.

CLIP datasets used:

1. 20 deg FBF-1
2. 20 deg FBF-2
3. Oo. FBF-1
4. Oo. FBF-2
5. Sp. FBF-1
6. Sp. FBF-2
7. Oo. FBF
8. Sp. FBF
9. 20 deg FBF



Peak calls
---

Changed wd into each directory under individual_clip/ and generared auto.ini files.

```bash
python ../spermatogenic/cwt-peakcaller/generate_config_ini.py .
```

Change min_rep_number to 3 (except, oo_fbf1).

Call peaks.

```bash
nohup python /groups/Kimble/Common/fbf_celltype/individual_clip/spermatogenic/cwt-peakcaller/main.py -c auto.ini -m -n -x > t &
# -m: no NB. -n: no UI. -x: clobber.
```

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

Final results:
```bash
python cliputil/score_metrics.py -c auto.ini -p combined_filtered/
```

Assigning to genes
---
This creates the combined_counts.txt file and the counts/ directories that are used for 
reads per gene analysis, including DESeq2.

```bash
python cliputil/assign_to_genes.py --gtf lib/gtf_with_names_column --bed all_bed_collapsed/
```


Supplementary tables
---

```bash
# File S1 Peaks
python cliputil/create_excel.py

# File S2 Stats.
# (When run with bed_files accessible.)
python cliputil/create_excel.py

# File S3 GO terms:
# ipynb for formating:
# cliputil/Convert IDs to WB ID for GO analysis on DAVID.ipynb

# File S4
python cliputil/deseqMaker.py
# Then run the R script given at the bottom of the .py file.

# File S5 Complex frequencies
python cliputil/heatmap.py

# File S6 Blocks
# Selected from the appropriate range in File S5 by hand.

# File S7 iPython notebook for conservation.
# orthos/pum2_vs_fbf.ipynb
```



