#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=add_inheritance
#SBATCH -o logs/2_add_inheritance.log
#SBATCH -e logs/2_add_inheritance.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/SPARK_WES/annotations/2022_03_21/inheritance

echo `date` starting job on $HOSTNAME

python 2_add_inheritance.py

echo `date` finished
