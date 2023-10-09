#!/bin/bash

mkdir -p plink_files/3_ld_het_filter

#Remove SNPs within long-range LD regions (Price et al, AJHG 2007)
#File format: CHR START END ID
/data5/software/plink --bfile plink_files/2_qc_filter/SPARK_FAM_QC --exclude 'range' Price_2007_blacklist_hg19.txt --make-bed  --write-snplist --out plink_files/3_ld_het_filter/SPARK_QC_final

#Remove samples with extremely high heterozygosity
#Select SNPs in 200 SNP/50 SNP sliding windows that have LD r^2 scores of >0.25
/data5/software/plink --bfile plink_files/3_ld_het_filter/SPARK_QC_final --keep plink_files/3_ld_het_filter/SPARK_QC_final.fam --extract plink_files/3_ld_het_filter/SPARK_QC_final.snplist --indep-pairwise 200 50 0.25 --out plink_files/3_ld_het_filter/SPARK_QC_final
#Calculate F coefficient samples to estimate heterozygosity rate in each sample
/data5/software/plink --bfile plink_files/3_ld_het_filter/SPARK_QC_final --extract plink_files/3_ld_het_filter/SPARK_QC_final.prune.in --keep plink_files/3_ld_het_filter/SPARK_QC_final.fam --het --out plink_files/3_ld_het_filter/SPARK_QC_final
#Remove samples with >2SD heterozygosity rate:
python scripts/3_heterozygosity_filter.py plink_files/3_ld_het_filter/SPARK_QC_final.het list_files/3.valid.sample list_files/3.invalid.sample
