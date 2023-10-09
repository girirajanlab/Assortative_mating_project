#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=filter1
#SBATCH -o logs/4_filter_exonic/batch_%a.log
#SBATCH -e logs/4_filter_exonic/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=2G
#SBATCH --array 2-346%100

echo `date` starting job on $HOSTNAME

SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=vcfs/annovar_genes/$SAMPLE.hg19_multianno.vcf
OUTPUT_FILE=exonic_vcfs/filter1/$SAMPLE.vcf.gz

echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_FILE

# note this will also include variants with Func.refGene=ncRNA_exonic and Func.refGene=ncRNA_splcing
# I will need to filter these out at a later step.. in 10_filter_variant_types.py
bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV19~"exonic" | INFO/Func.wgEncodeGencodeBasicV19~"splicing"' $INPUT_VCF | bgzip > $OUTPUT_FILE

echo `date` finished
