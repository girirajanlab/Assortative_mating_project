#!/bin/bash

# Merge individual chromosome files
find plink_files/2_subset_bed/*.bed | cut -f 1 -d . > list_files/chr_files.list

mkdir plink_files/3_merge

/data5/software/plink --merge-list list_files/chr_files.list --make-bed --out plink_files/3_merge/ukb
