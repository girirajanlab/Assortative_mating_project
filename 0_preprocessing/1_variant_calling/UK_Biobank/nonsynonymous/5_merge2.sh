#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=5_merge2
#SBATCH -o logs/5_merge2.log
#SBATCH -e logs/5_merge2.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G

echo `date` starting job on $HOSTNAME

bcftools merge -m none vcfs/4_merge1/*.vcf.gz | bgzip > vcfs/all.vcf.gz
tabix -p vcf vcfs/all.vcf.gz

echo `date` finished
