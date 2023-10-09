#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=2_make_pfb
#SBATCH -o logs/2_make_pfb.log
#SBATCH -e logs/2_make_pfb.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G

# Compile a Population frequency of B allele file
export PENNCNV_PATH=/data5/software/PennCNV-1.0.5/

$PENNCNV_PATH/compile_pfb.pl --listfile Analysis_files/1Mv1_files.list --snpposfile Analysis_files/1Mv1_snp_pos.txt  --output Analysis_files/1Mv1.pfb
$PENNCNV_PATH/compile_pfb.pl --listfile Analysis_files/1Mv3_files.list --snpposfile Analysis_files/1Mv3_snp_pos.txt  --output Analysis_files/1Mv3.pfb
$PENNCNV_PATH/compile_pfb.pl --listfile Analysis_files/Omni2.5_files.list --snpposfile Analysis_files/Omni2.5_snp_pos.txt  --output Analysis_files/Omni2.5.pfb
