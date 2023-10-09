#!/bin/bash

mkdir -p plink_files/6_final_filter

/data5/software/plink --bfile plink_files/4_filter/UKB_QC_final --make-bed --keep list_files/5_european_pairs.fam --out plink_files/6_final_filter/UKB_assort_mating_final --extract plink_files/4_filter/UKB_QC_final.snplist
/data5/software/plink --bfile plink_files/6_final_filter/UKB_assort_mating_final --make-bed --update-ids list_files/5_update_fid.list --out plink_files/6_final_filter/UKB_king_input
