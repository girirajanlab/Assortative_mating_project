#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=annovar
#SBATCH -o logs/3_annovar_genes/batch_%a.log
#SBATCH -e logs/3_annovar_genes/batch_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=30G
#SBATCH --array 2-346%100

echo `date` starting job on $HOSTNAME

SAMPLE=`cat samples.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n1 | cut -f2 -d,`
INPUT_VCF=vcfs/filter1/$SAMPLE.vcf.gz
OUTPUT_PREFIX=vcfs/annovar_genes/$SAMPLE

echo $SAMPLE
echo input file: $INPUT_VCF
echo output file: $OUTPUT_PREFIX

perl annovar/table_annovar.pl $INPUT_VCF hg19/annotations/annovar \
 -buildver hg19 \
 -out $OUTPUT_PREFIX \
 -remove \
 -protocol refGene,wgEncodeGencodeBasicV19 \
 -operation g,g \
 -nastring . \
 -vcfinput \
 -arg '-hgvs','-hgvs'

echo `date` finished
