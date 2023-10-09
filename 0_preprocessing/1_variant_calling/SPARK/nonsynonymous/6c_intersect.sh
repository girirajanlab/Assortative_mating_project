#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=intersect
#SBATCH -o logs/6_intersect/%a.log
#SBATCH -e logs/6_intersect/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 1-200%100

echo `date` $HOSTNAME

ANNOTATED_FILE=vcfs/5_annotations/all.hg38_multianno.exonic.gnomad.cadd.vcf.gz
IN_DIR=vcfs/2_quality_filtered
OUT_DIR=vcfs/6_intersect

#Set the number of runs that each SLURM task should do
PER_TASK=137

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  #Do your stuff here
  LINE=$(sed -n "$run"p samples.csv)
  echo $LINE
  sample=`echo $LINE | cut -d, -f1`
  last_digit=`echo $LINE | cut -d, -f2`
  in=$IN_DIR/$last_digit/$sample.vcf.gz
  out1=$OUT_DIR/$last_digit/$sample.vcf.gz
  out2=$OUT_DIR/$last_digit/$sample.tsv
  echo $sample $in $out

  # annotate the variants in a sample
  bcftools annotate -a $ANNOTATED_FILE -c INFO $in | bcftools view -i "ANNOVAR_DATE!='.'" | bgzip > $out1
  # to table
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV38\t%Gene.wgEncodeGencodeBasicV38\t%GeneDetail.wgEncodeGencodeBasicV38\t%ExonicFunc.wgEncodeGencodeBasicV38\t%AAChange.wgEncodeGencodeBasicV38\t%gnomad_exome_AF\t%gnomad_genome_AF\t%CADD_PHRED\t%CADD_RawScore[\t%SAMPLE\t%GT\t%DP\t%AD\t%GQ\t%PL]\n' $out1 > $out2

done

echo `date` $HOSTNAME
