#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=intersect
#SBATCH -o logs/7_filter/%a.log
#SBATCH -e logs/7_filter/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 1-200%100

echo `date` $HOSTNAME

# Apply the following filters to all files
#       gnomAD_genome <= 0.001
#       gnomAD_exome <= 0.001
#	Func: `exonic` or  `splicing`

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
  input=vcfs/6_intersect/$last_digit/$sample.vcf.gz
  output1=vcfs/7_filter_variants/$last_digit/$sample.vcf.gz
  output2=tables/7_filter_variants/$sample.tsv
  echo $sample $input $output

  # annotate the variants in a sample
  bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' $input | bcftools view -i 'gnomad_exome_AF<=0.001 | gnomad_exome_AF="."' | bcftools view -i 'gnomad_genome_AF<=0.001 | gnomad_genome_AF="."' | bgzip > $output1
  # to table
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV38\t%Gene.wgEncodeGencodeBasicV38\t%GeneDetail.wgEncodeGencodeBasicV38\t%ExonicFunc.wgEncodeGencodeBasicV38\t%AAChange.wgEncodeGencodeBasicV38\t%gnomad_exome_AF\t%gnomad_genome_AF\t%CADD_PHRED\t%CADD_RawScore[\t%SAMPLE\t%GT\t%DP\t%AD\t%GQ\t%PL]\n' $output1 > $output2
done

echo `date` $HOSTNAME
