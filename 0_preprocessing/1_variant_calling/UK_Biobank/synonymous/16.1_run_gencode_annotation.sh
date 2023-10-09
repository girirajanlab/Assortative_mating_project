#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_gencode
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=100G
#SBATCH -o logs/16_anno_gencode.log
#SBATCH -e logs/16_anno_gencode.log

echo `date` starting job on $HOSTNAME

# Annotate GENCODE genes in the UK Biobank calls
python 16_annotate_gencode_genes.py

echo `date` ending job
