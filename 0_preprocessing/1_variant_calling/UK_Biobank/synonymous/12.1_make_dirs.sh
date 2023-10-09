#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=12.1_make_dirs
#SBATCH -o logs/12.1_make_dirs.log
#SBATCH -e logs/12.1_make_dirs.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=1G
#SBATCH --chdir /data5/UK_Biobank/annotations/annovar/2023_08_16/synonymous
#SBATCH --array 1-100

num_minus1=$(($SLURM_ARRAY_TASK_ID-1))
formatted_num=$(printf "%02d\n" $num_minus1)

echo $num_minus1
echo $formatted_num

mkdir tables/by_sample/$formatted_num
