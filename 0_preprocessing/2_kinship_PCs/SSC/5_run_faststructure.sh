#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=faststructure
#SBATCH -o logs/5_faststructure.log
#SBATCH -e logs/5_faststructure.log
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --nodelist ramona

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/bin/gsl/lib
export CFLAGS="-I/usr/bin/gsl/include"
export LDFLAGS="-L/usr/bin/gsl/lib"

mkdir -p plink_files/5_faststructure

python2 fastStructure/structure.py -K 6 --input=plink_files/4_merge_hapmap/SSC_hapmap_combined --output=plink_files/5_faststructure/SSC_hapmap_combined
paste plink_files/4_merge_hapmap/SSC_hapmap_combined.fam plink_files/5_faststructure/SSC_hapmap_combined.6.meanQ > list_files/5_SSC_hapmap_combined_output.txt
