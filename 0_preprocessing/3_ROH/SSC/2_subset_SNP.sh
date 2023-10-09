#!/bin/bash

# Subset the SNP data to only include the proband files we want to run
mkdir Analysis_files/2_subset_probands

# Take PLINK files from assortative mating analysis
/data5/software/plink --bfile ../../2_kinship_PCs/SSC/plink_files/3_ld_het_filter/SSC_QC_final --keep Analysis_files/1_roh_samples.fam --make-bed --out Analysis_files/2_subset_probands/2_roh_samples
