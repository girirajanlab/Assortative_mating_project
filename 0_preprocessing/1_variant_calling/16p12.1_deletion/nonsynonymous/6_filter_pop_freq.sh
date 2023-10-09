#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gnomad
#SBATCH -o logs/6_filter_pop_freq/batch_%a.log
#SBATCH -e logs/6_filter_pop_freq/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 2-346%100

echo `date` starting job on $HOSTNAME

SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=exonic_vcfs/gnomad/$SAMPLE.vcf.gz
OUTPUT_VCF=exonic_vcfs/filter2/$SAMPLE.vcf.gz

echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_VCF

bcftools view -i 'gnomad_exome_AF<=0.001 | gnomad_exome_AF="."' $INPUT_VCF | bcftools view -i 'gnomad_genome_AF<=0.001 | gnomad_genome_AF="."' | bgzip > $OUTPUT_VCF

echo `date` finished
