#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_gnomad
#SBATCH -o logs/10_gnomad.log
#SBATCH -e logs/10_gnomad.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G

echo `date` starting job on $HOSTNAME

invcf=vcfs/9_all_annovar_exonic_splicing_cadd.vcf.gz
outvcf=vcfs/10_gnomad.vcf.gz

vcfanno_linux64.1 -p 15 gnomad.toml $invcf | bgzip > $outvcf
tabix -p vcf $outvcf

echo `date` finished
