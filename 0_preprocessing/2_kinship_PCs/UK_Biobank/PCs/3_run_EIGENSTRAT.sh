#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=eigenstrat
#SBATCH -o logs/3_EIGENSTRAT.log
#SBATCH -e logs/3_EIGENSTRAT.log
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=300G

echo Started on $HOSTNAME at `date`

# Perform PCA using EIGENSTRAT

PATH=/data5/software/EIG/bin:$PATH

# Rename files for EIGENSTRAT
cp plink_files/2_filter_samples/UKB_european.bim Analysis_files/3_eigen_input.pedsnp
cp plink_files/2_filter_samples/UKB_european.fam Analysis_files/3_eigen_input.pedind

# All samples have -9 for phenotype - replace with 3: https://github.com/DReichLab/AdmixTools/issues/56
sed -i 's/-9/3/g' Analysis_files/3_eigen_input.pedind

perl /data5/software/EIG/bin/smartpca.perl -i plink_files/2_filter_samples/UKB_european.bed -a Analysis_files/3_eigen_input.pedsnp -b Analysis_files/3_eigen_input.pedind -k 20 -o Results/3_eigenstrat.pca -p Results/3_eigenstrat.plot -e Results/3_eigenstrat.eval -m 0 -l logs/3_eigenstrat.log

echo Ended at `date`
