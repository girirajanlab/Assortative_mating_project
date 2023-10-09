#!/bin/bash

mkdir -p plink_files/2_qc_filter

#PLINK QC filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-5 (under selection), geno<0.05 (SNPs missing in >5% of subjects)
#Sample-level filtering: mind<0.01 (Individuals with >1% missing genotypes), me>0.05 (Mendelian error rate for both variants and samples)
/data5/software/plink --bfile plink_files/1_merge/16p12_merge --maf 0.05 --hwe 1e-5 --geno 0.05 --chr 1-22 --make-bed --out plink_files/2_qc_filter/16p12_SNP_QC
/data5/software/plink --bfile plink_files/2_qc_filter/16p12_SNP_QC --mind 0.01 --me 0.05 0.05 --write-snplist --make-bed --out plink_files/2_qc_filter/16p12_FAM_QC

#Check for duplicate SNPs from PLINK files
cut -f 2 plink_files/2_qc_filter/16p12_FAM_QC.bim | sort | uniq -d > list_files/2_dup_snps.list
#No duplicate SNPs
