#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=filter2
#SBATCH -o logs/2_filter_vcfs/%a.log
#SBATCH -e logs/2_filter_vcfs/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 1-136

echo `date` $HOSTNAME

FASTA=Resources/genome.fa
OUT_DIR=vcfs/2_quality_filtered/

#Set the number of runs that each SLURM task should do
PER_TASK=201

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
  in=`echo $LINE | cut -d, -f3`
  out=$OUT_DIR/$last_digit/$sample.vcf.gz
  echo $sample $in $out
 
  bcftools norm -f ${FASTA} -m-both $in | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' -o $out -Oz
  bcftools index $out

done


echo `date` $HOSTNAME
