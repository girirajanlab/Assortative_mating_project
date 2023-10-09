#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_tables
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=4G
#SBATCH -o logs/12_anno_samples/%a.log
#SBATCH -e logs/12_anno_samples/%a.log
#SBATCH --array 1-1000%150
#SBATCH --exclude qingyu,laila

# Annotate variants within each sample and format as table
# Create table out of VCFs

echo `date` starting job on $HOSTNAME

#Set the number of runs that each SLURM task should do
# 999 jobs x 201 samples per job = 200,709 > 200604 total samples
PER_TASK=201

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

if [ $END_NUM -gt 200605 ];
  then
    END_NUM=200605
fi

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  # Get sample info
  SAMPLE=$(head -n $run samples.csv | tail -n 1 | cut -f 2 -d ,)
  LAST_TWO=$(head -n $run samples.csv | tail -n 1 | cut -f 3 -d ,)
  SAMPLE_VCF=vcfs/by_sample/$LAST_TWO/$SAMPLE.vcf.gz
  TMP_VCF=vcfs/$SAMPLE.tmp.vcf.gz
  OUTPUT_TABLE=tables/by_sample/$LAST_TWO/$SAMPLE.tsv
  echo $SAMPLE $SAMPLE_VCF $OUTPUT_TABLE
  # Annotate VCF and convert to table
  bcftools annotate -a ../../2023_02_24/vcfs/11_filter3.vcf.gz -c INFO -O z -o $TMP_VCF $SAMPLE_VCF
  echo Finished annotation
  # ExonicFunc.wgEncodeGencodeBasicV38
  bcftools view -i "ANNOVAR_DATE!='.'" $TMP_VCF | bcftools view -i "INFO/ExonicFunc.wgEncodeGencodeBasicV38='synonymous_SNV'" | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%Func.wgEncodeGencodeBasicV38\t%Gene.wgEncodeGencodeBasicV38\t%GeneDetail.wgEncodeGencodeBasicV38\t%ExonicFunc.wgEncodeGencodeBasicV38\t%AAChange.wgEncodeGencodeBasicV38\t%CADD_PHRED[\t${sample_name}\t%GT\t%DP\t%AD]\n" -o $OUTPUT_TABLE
  rm $TMP_VCF
done

echo `date` ending job
