#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=2_subset_bed
#SBATCH -o logs/2_subset_bed.log
#SBATCH -e logs/2_subset_bed.log
#SBATCH --ntasks=1
#SBATCH --time=720:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --array 1-25

chrom=`head -n $SLURM_ARRAY_TASK_ID list_files/chr.list | tail -n1`

/data5/software/plink --bed download/ukb_cal_chr$chrom'_v2.bed' \
  --bim download/ukb_snp_chr$chrom'_v2.bim' \
 --fam download/ukb45023_s488288.fam \
 --keep-fam list_files/spouses.fam \
 --make-bed \
 --out plink_files/2_subset_bed/spouse_chr$chrom
