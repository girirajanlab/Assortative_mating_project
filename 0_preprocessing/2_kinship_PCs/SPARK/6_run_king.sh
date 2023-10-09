#!/bin/bash

mkdir -p plink_files/6_run_king

# Filter for European spouse pairs and run KING
/data5/software/plink --bfile plink_files/3_ld_het_filter/SPARK_QC_final --make-bed --keep 5_FASTSTRUCTURE/Results/2_european_spouses.txt --out plink_files/6_run_king/SPARK_KING

#Step 2: Calculate kinship coefficients between pairs of parents
#Use KING software to calculate pairwise kinship coefficients for all parents (http://people.virginia.edu/~wc9c/KING/manual.html)
/data5/software/king -b plink_files/6_run_king/SPARK_KING.bed --kinship --degree 1
