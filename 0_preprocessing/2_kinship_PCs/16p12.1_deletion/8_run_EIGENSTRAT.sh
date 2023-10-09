#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=eigenstrat
#SBATCH -o logs/8_EIGENSTRAT.log
#SBATCH -e logs/8_EIGENSTRAT.log
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=5G

echo Started on $HOSTNAME at `date`

# Perform PCA using EIGENSTRAT

PATH=/data5/software/EIG/bin:$PATH

# Rename files for EIGENSTRAT
cp plink_files/7_run_king/16p12_KING.bim analysis_files/8_eigen_input.pedsnp
cp plink_files/7_run_king/16p12_KING.fam analysis_files/8_eigen_input.pedind

# All samples have -9 for phenotype - replace with 3: https://github.com/DReichLab/AdmixTools/issues/56
sed -i 's/-9/3/g' analysis_files/8_eigen_input.pedind

perl /data5/software/EIG/bin/smartpca.perl -i plink_files/7_run_king/16p12_KING.bed -a analysis_files/8_eigen_input.pedsnp -b analysis_files/8_eigen_input.pedind -k 20 -o Results/8_eigenstrat.pca -p Results/8_eigenstrat.plot -e Results/8_eigenstrat.eval -m 0 -l logs/8_eigenstrat.log

echo Ended at `date`
