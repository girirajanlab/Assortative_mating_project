#!/bin/bash

#Determine basic ethnicity of samples, to later remove families with spouses who are from different populations
#HapMap3 data source (hg19 liftOver): https://www.dropbox.com/sh/ycvtoner1d2bujs/AACtKzZwsL9ln0vI5Fnxuof0a
#Data filtering and merging pipeline adapted from: https://cran.r-project.org/web/packages/plinkQC/vignettes/AncestryCheck.pdf
#Pipeline source: http://rajanil.github.io/fastStructure/

#Find overlapping SNPs shared with HapMap3 array data
cut -f2 plink_files/3_ld_het_filter/SSC_QC_final.bim > list_files/4_SSC_SNP.list
cut -f2 hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim > list_files/4_hapmap3_SNP.list
grep -w -f list_files/4_hapmap3_SNP.list list_files/4_SSC_SNP.list  | sort | uniq > list_files/4_shared_snps.list

mkdir -p plink_files/4_merge_hapmap

#Use PLINK to generate BED/BIM files with common SNPs
/data5/software/plink --bfile plink_files/3_ld_het_filter/SSC_QC_final  --extract list_files/4_shared_snps.list --make-bed --out plink_files/4_merge_hapmap/SSC
/data5/software/plink --bfile hapmap3/hapmap3_r1_b37_fwd_consensus.qc.poly.recode --extract list_files/4_shared_snps.list --make-bed --out plink_files/4_merge_hapmap/hapmap3

#Remove ambiguous SNPs (A/C or G/T).
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' plink_files/4_merge_hapmap/SSC.bim > list_files/4_ambig_snps.list
/data5/software/plink --bfile plink_files/4_merge_hapmap/SSC --exclude list_files/4_ambig_snps.list --make-bed --out plink_files/4_merge_hapmap/SSC_noambig
/data5/software/plink --bfile plink_files/4_merge_hapmap/hapmap3 --exclude list_files/4_ambig_snps.list --make-bed --out plink_files/4_merge_hapmap/hapmap3_noambig

#Prune SNPs in linkage disequilibrium (R^2 >0.2 in 50x5-SNP sliding windows)
/data5/software/plink --bfile plink_files/4_merge_hapmap/SSC_noambig --indep-pairwise 50 5 0.2 --out plink_files/4_merge_hapmap/ld_prune
/data5/software/plink --bfile plink_files/4_merge_hapmap/hapmap3_noambig --extract plink_files/4_merge_hapmap/ld_prune.prune.in --make-bed --out plink_files/4_merge_hapmap/hapmap3_prune
/data5/software/plink --bfile plink_files/4_merge_hapmap/SSC_noambig --extract plink_files/4_merge_hapmap/ld_prune.prune.in --make-bed --out plink_files/4_merge_hapmap/SSC_prune

#Correct allele flips and SNP chromosome and position mismatches in HapMap file
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} ($2 in a && a[$2] != $1) {print a[$2],$2}' plink_files/4_merge_hapmap/SSC_prune.bim plink_files/4_merge_hapmap/hapmap3_prune.bim  > list_files/update_chr.list
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} ($2 in a && a[$2] != $4) {print a[$2],$2}' plink_files/4_merge_hapmap/SSC_prune.bim plink_files/4_merge_hapmap/hapmap3_prune.bim  > list_files/update_pos.list
#No chr mismatches identified
/data5/software/plink --bfile plink_files/4_merge_hapmap/hapmap3_prune --update-map list_files/update_pos.list 1 2 --make-bed --out plink_files/4_merge_hapmap/hapmap3_update_pos

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' plink_files/4_merge_hapmap/SSC_prune.bim plink_files/4_merge_hapmap/hapmap3_update_pos.bim > list_files/update_flip.list
/data5/software/plink --bfile plink_files/4_merge_hapmap/hapmap3_update_pos --flip list_files/update_flip.list --make-bed --out plink_files/4_merge_hapmap/hapmap3_flip

#Remove SNPs with alleles that don't match from both datasets
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' plink_files/4_merge_hapmap/SSC_prune.bim plink_files/4_merge_hapmap/hapmap3_flip.bim > list_files/4_mismatch.list
/data5/software/plink --bfile plink_files/4_merge_hapmap/hapmap3_flip --exclude list_files/4_mismatch.list --make-bed --out plink_files/4_merge_hapmap/hapmap3_final
/data5/software/plink --bfile plink_files/4_merge_hapmap/SSC_prune --exclude list_files/4_mismatch.list --make-bed --out plink_files/4_merge_hapmap/SSC_final

#Finally, merge HapMap3 and SSC datasets together
/data5/software/plink --bfile plink_files/4_merge_hapmap/SSC_final --bmerge plink_files/4_merge_hapmap/hapmap3_final --make-bed --out plink_files/4_merge_hapmap/SSC_hapmap_combined
