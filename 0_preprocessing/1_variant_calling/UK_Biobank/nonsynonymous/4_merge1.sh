#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=4_merge1
#SBATCH -o logs/4_merge1/%a.log
#SBATCH -e logs/4_merge1/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --array 1-100

echo `date` starting job on $HOSTNAME

num_minus1=$(($SLURM_ARRAY_TASK_ID-1))
formatted_num=$(printf "%02d\n" $num_minus1)

echo $num_minus1
echo $formatted_num

# Can't merge all files at once - too many files for bcftools to handle
# Split into smaller groups based on first digit
for i in {1..6};
do
echo $i
bcftools merge -m none vcfs/by_sample/$formatted_num/$i*.vcf.gz | bcftools view -G | bgzip > vcfs/4_merge1/$formatted_num.$i.vcf.gz
tabix -p vcf vcfs/4_merge1/$formatted_num.$i.vcf.gz
done

echo Merging all
bcftools merge -m none vcfs/4_merge1/$formatted_num.*.vcf.gz | bcftools view -G | bgzip > vcfs/4_merge1/$formatted_num.vcf.gz
tabix -p vcf vcfs/4_merge1/$formatted_num.vcf.gz

rm vcfs/4_merge1/$formatted_num.*.vcf.gz*

echo `date` finished
