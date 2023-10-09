#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=merge1
#SBATCH -o logs/4b_merge1/%a.log
#SBATCH -e logs/4b_merge1/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 1-100

echo `date` $HOSTNAME

INDIR=vcfs/2_quality_filtered/
OUTDIR=vcfs/4b_merge1/

i=`cat task_arrays.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f1 -d' '`
j=`cat task_arrays.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f2 -d' '`
echo $i $j
bcftools merge -m none $INDIR/$j/*${i}${j}.vcf.gz | bcftools view -G | bgzip > $OUTDIR/${i}${j}.vcf.gz

tabix -p vcf $OUTDIR/${i}${j}.vcf.gz


echo `date` $HOSTNAME
