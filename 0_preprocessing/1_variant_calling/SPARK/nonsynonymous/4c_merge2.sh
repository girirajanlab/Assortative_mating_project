#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=merge2
#SBATCH -o logs/4c_merge2.log
#SBATCH -e logs/4c_merge2.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G

echo `date` starting job on $HOSTNAME

bcftools merge -m none vcfs/4b_merge1/*.vcf.gz | bgzip > vcfs/4c_merge2.vcf.gz

tabix -p vcf vcfs/4c_merge2.vcf.gz

echo `date` finished
