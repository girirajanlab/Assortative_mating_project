#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=6_merge_calls
#SBATCH -o logs/6_merge.log
#SBATCH -e logs/6_merge.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G

python 6_merge_calls.py Analysis_files/1Mv1_snp_pos.txt raw_out/1Mv1_autosome.txt raw_out/1Mv1_chrX.txt tables/1Mv1_merged.txt
python 6_merge_calls.py Analysis_files/1Mv3_snp_pos.txt raw_out/1Mv3_autosome.txt raw_out/1Mv3_chrX.txt tables/1Mv3_merged.txt
python 6_merge_calls.py Analysis_files/Omni2.5_snp_pos.txt raw_out/Omni2.5_autosome.txt raw_out/Omni2.5_chrX.txt tables/Omni2.5_merged.txt
