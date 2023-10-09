#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=clinvar
#SBATCH -o logs/8_annotate_clinvar/batch_%a.log
#SBATCH -e logs/8_annotate_clinvar/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=5G
#SBATCH --array 2-346%10

echo `date` starting job on $HOSTNAME

SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=exonic_vcfs/cadd/$SAMPLE.vcf.gz
OUTPUT_VCF=vcfs/8_annotate_clinvar/$SAMPLE.vcf.gz

echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_VCF

vcfanno_linux64.1 ClinVar.toml $INPUT_VCF | bgzip > $OUTPUT_VCF

echo `date` finished
