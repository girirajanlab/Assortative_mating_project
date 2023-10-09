#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=6_merge_calls
#SBATCH -o logs/6_merge/%a.log
#SBATCH -e logs/6_merge/%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --array 1-888

echo `date` $HOSTNAME

#Set the number of runs that each SLURM task should do
PER_TASK=4

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

for (( run=$START_NUM; run<=END_NUM; run++ )); do
	samp=`head -n $run Analysis_files/1_samples.csv | tail -n1 | cut -f 2 -d ,`
	grep $samp raw_out/iwes_autosome_adj.txt > tables/tmp/$samp.autosome.txt
	grep $samp raw_out/iwes_chrX_adj.txt >> tables/tmp/$samp.chrX.txt

	python 6_merge_calls.py Analysis_files/iWES_snp_pos.txt tables/tmp/$samp.autosome.txt tables/tmp/$samp.chrX.txt tables/6_merge_individual/$samp.txt

	rm tables/tmp/$samp.autosome.txt
	rm tables/tmp/$samp.chrX.txt
done

echo `date` $HOSTNAME
