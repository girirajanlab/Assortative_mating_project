#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=spark_faststructure
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH -o logs/1_faststructure_subsets/%a.log
#SBATCH -e logs/1_faststructure_subsets/%a.log
#SBATCH --array 1-1000%200
#SBATCH --nodelist ramona

# To save time, run faststructure on just a few samples at once

echo `date` starting job on $HOSTNAME

PER_TASK=57

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

if [ $END_NUM -gt 56259 ];
  then
    END_NUM=56259
fi

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Subset the desired samples
head -n $END_NUM ../plink_files/4_merge_hapmap/SPARK_final.fam | tail -n $PER_TASK > plink_files/$SLURM_ARRAY_TASK_ID.tmp.fam
cat ../plink_files/4_merge_hapmap/hapmap3_final.fam >> plink_files/$SLURM_ARRAY_TASK_ID.tmp.fam
/data5/software/plink --bfile ../plink_files/4_merge_hapmap/SPARK_hapmap_combined --keep plink_files/$SLURM_ARRAY_TASK_ID.tmp.fam --make-bed --out plink_files/$SLURM_ARRAY_TASK_ID.tmp2

# Run faststructure
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/bin/gsl/lib
export CFLAGS="-I/usr/bin/gsl/include"
export LDFLAGS="-L/usr/bin/gsl/lib"

python2 fastStructure/structure.py -K 6 --input=plink_files/$SLURM_ARRAY_TASK_ID.tmp2 --output=plink_files/faststructure_output/$SLURM_ARRAY_TASK_ID
paste plink_files/$SLURM_ARRAY_TASK_ID.tmp2.fam plink_files/faststructure_output/$SLURM_ARRAY_TASK_ID.6.meanQ > Results/faststructure_output/$SLURM_ARRAY_TASK_ID'_output.txt'

rm plink_files/$SLURM_ARRAY_TASK_ID.tmp*
