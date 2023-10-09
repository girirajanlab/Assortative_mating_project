#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=9_cadd
#SBATCH -o logs/9_cadd.log
#SBATCH -e logs/9_cadd.log
#SBATCH --ntasks=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=6G

echo `date` starting job on $HOSTNAME

# Annotate with CADD scores

invcf=vcfs/8_all_annovar_exonic_splicing.vcf
outvcf=vcfs/9_all_annovar_exonic_splicing_cadd.vcf.gz

vcfanno_linux64.1 -p 15 cadd.toml $invcf | bgzip > $outvcf
tabix -p vcf $outvcf

echo `date` finished
