#!/bin/bash

mkdir -p plink_files/7_run_king

# Filter for European spouse pairs and run KING
/data5/software/plink --bfile plink_files/3_ld_het_filter/SSC_QC_final --make-bed --keep list_files/6_european_spouses.fam --out plink_files/7_run_king/SSC_european

# Update sample IDs
/data5/software/plink --bfile plink_files/7_run_king/SSC_european --update-ids SSC_microarray_samples.update.txt --make-bed --out plink_files/7_run_king/SSC_KING


#Step 2: Calculate kinship coefficients between pairs of parents
#Use KING software to calculate pairwise kinship coefficients for all parents (http://people.virginia.edu/~wc9c/KING/manual.html)
/data5/software/king -b plink_files/7_run_king/SSC_KING.bed --kinship --degree 1
