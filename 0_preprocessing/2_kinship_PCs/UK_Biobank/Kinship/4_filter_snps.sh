#!/bin/bash

#PLINK QC filters: Cohort MAF <0.05, Hardy-Weinberg equilibrium <1e-5 (under selection), geno<0.05 (SNPs missing in >5% of subjects)
#Sample-level filtering: mind<0.01 (Individuals with >1% missing genotypes), me>0.05 (Mendelian error rate for both variants and samples)
#Note: parameters are slightly different from PRS calculations
/data5/software/plink --bfile plink_files/3_merge/ukb --maf 0.05 --hwe 1e-5 --geno 0.05 --chr 1-22 --make-bed --out plink_files/4_filter/UKB_SNP_QC
/data5/software/plink --bfile plink_files/4_filter/UKB_SNP_QC --mind 0.01 --me 0.05 0.05 --write-snplist --make-bed --out plink_files/4_filter/UKB_FAM_QC

#Remove SNPs within long-range LD regions (Price et al, AJHG 2007)--224,283 SNPs remain
#File format: CHR START END ID
/data5/software/plink --bfile plink_files/4_filter/UKB_FAM_QC --exclude 'range' Price_2007_blacklist_hg19.txt --make-bed  --write-snplist --out plink_files/4_filter/UKB_QC_final

#Remove samples with extremely high heterozygosity
#Select SNPs in 200 SNP/50 SNP sliding windows that have LD r^2 scores of >0.25
/data5/software/plink --bfile plink_files/4_filter/UKB_QC_final --keep plink_files/4_filter/UKB_QC_final.fam --extract plink_files/4_filter/UKB_QC_final.snplist --indep-pairwise 200 50 0.25 --out plink_files/4_filter/UKB_QC_final
#Calculate F coefficient samples to estimate heterozygosity rate in each sample
/data5/software/plink --bfile plink_files/4_filter/UKB_QC_final --extract plink_files/4_filter/UKB_QC_final.prune.in --keep plink_files/4_filter/UKB_QC_final.fam --het --out plink_files/4_filter/UKB_QC_final
#Remove samples with >2SD heterozygosity rate:
Rscript scripts/heterzygosity_filter.R plink_files/4_filter/UKB_QC_final.het plink_files/4_filter/UKB_QC_het.valid.sample plink_files/4_filter/UKB_QC_het.invalid.sample
