#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=am_simulation
#SBATCH -o logs/1_%a.log
#SBATCH -e logs/1_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:0:0
#SBATCH --mem-per-cpu=150G
#SBATCH --chdir /data6/corrine/10_Generational_Model/NEW_MODEL/2023_09_28
#SBATCH --array 1-28
#SBATCH --exclude qingyu,sarah,laila

echo Started on $hostname at `date`

mkdir -p Simulation_outputs/$SLURM_ARRAY_TASK_ID

parameters=$(head -n $SLURM_ARRAY_TASK_ID 0_model_parameters.txt | tail -1)
echo $SLURM_ARRAY_TASK_ID $parameters

python 1_varying_parameters.py $SLURM_ARRAY_TASK_ID $parameters

echo Finished at `date`
