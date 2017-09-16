Comparing FBF iCLIP in differnt cell types
======

Code for a paper on FBF in the gendered transcriptome

This analysis begins with the filtered, separately called oo/sp/20deg peaks
 files and associated wig files.
 It also includes library files of DESeq2 output - quite a few library files,
 in fact.i

Essentially, CLIP data constitutes:
    1. Peaks files in combined_filtered/
    2. bedgraphs, norm
    3. bedgraphs, unnorm
    4. bed files

The second type of input into this analysis is a variety of library files,
 ALL OF WHICH will be placed under /opt/lib/.
 The analysis files are:
     1. DESeq2 of sp/oo germlines
     2. GL mitotic SAM.

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



