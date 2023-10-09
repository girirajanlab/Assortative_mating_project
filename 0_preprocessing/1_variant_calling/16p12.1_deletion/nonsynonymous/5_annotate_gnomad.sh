#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gnomad
#SBATCH -o logs/5_annotate_gnomad/batch_%a.log
#SBATCH -e logs/5_annotate_gnomad/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 2-346

echo `date` starting job on $HOSTNAME

SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=exonic_vcfs/filter1/$SAMPLE.vcf.gz
OUTPUT_VCF=exonic_vcfs/gnomad/$SAMPLE.vcf.gz

echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_VCF

vcfanno_linux64.1 gnomad.toml $INPUT_VCF | bgzip > $OUTPUT_VCF

echo `date` finished
