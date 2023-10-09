#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ssc_pipeline_final
#SBATCH -o logs/ssc_pipeline_final.log
#SBATCH -e logs/ssc_pipeline_final.err.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G

mkdir intermediate_files/final

for FILE in `ls intermediate_files/*/*.11_variants_merged.txt`; do cat $FILE >> intermediate_files/final/12_variants_merged_all.txt; done
python 13_filter_variant_types.py
python 14_intracohort_filter.py
python 15_gencode_gene_ids.py
python 16_loeuf_scores.py
