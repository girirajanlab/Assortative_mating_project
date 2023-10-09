#!/bin/bash

#Merge BED files of all batches together, to create one fileset with all samples
mkdir -p plink_files/1_merge

/data5/software/plink --merge-list list_files/plink_filenames.list --out plink_files/1_merge/16p12_initial_merge

# Perform strand flipping on Hudson Alpha and Yale and re-merge
/data5/software/plink --bfile 16p12_hudson --flip plink_files/1_merge/16p12_initial_merge.missnp --make-bed --out plink_files/1_merge/ha_flip
/data5/software/plink --bfile 16p12_yale --flip plink_files/1_merge/16p12_initial_merge.missnp --make-bed --out plink_files/1_merge/yale_flip
/data5/software/plink --merge-list list_files/plink_filenames_flipped.list --out plink_files/1_merge/16p12_flip_merge

# Remove excluded samples
/data5/software/plink --bfile plink_files/1_merge/16p12_flip_merge --keep 16p12_microarray_samples.keep.txt --make-bed --out plink_files/1_merge/16p12_keep
/data5/software/plink --bfile plink_files/1_merge/16p12_keep --update-ids 16p12_microarray_samples.update.txt --make-bed --out plink_files/1_merge/16p12_merge
