#!/bin/bash

mkdir -p plink_files/2_filter_samples

/data5/software/plink --bfile ../1_Kinship/plink_files/4_filter/UKB_QC_final --make-bed --keep list_files/1_european_samples.fam --out plink_files/2_filter_samples/UKB_european
