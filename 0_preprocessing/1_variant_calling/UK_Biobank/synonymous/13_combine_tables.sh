#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_combine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=4G
#SBATCH -o logs/13_combine_tables.log
#SBATCH -e logs/13_combine_tables.log

# Cat variants into a single file

echo `date` starting job on $HOSTNAME

for (( run=2; run<=200605; run++ )); do
  # Get sample information
  SAMPLE=$(head -n $run samples.csv | tail -n 1 | cut -f 2 -d ,)
  LAST_TWO=$(head -n $run samples.csv | tail -n 1 | cut -f 3 -d ,)
  OUTPUT_DIR=tables/by_sample/$LAST_TWO
  echo $SAMPLE $OUTPUT_DIR
  # Combine variants into one table
  sed "s/$/\t$SAMPLE/" "$OUTPUT_DIR"/"$SAMPLE".tsv > tmp/$SAMPLE.tsv
  cat tmp/$SAMPLE.tsv >> tables/13_filtered_combined.tsv
  rm tmp/$SAMPLE.tsv
done

echo `date` ending job
