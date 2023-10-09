#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=4_run_penncnv
#SBATCH -o logs/4_chrX_cnvs.log
#SBATCH -e logs/4_chrX_cnvs.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G

export PENNCNV_PATH=/data5/software/PennCNV/PennCNV-1.0.5

echo Started on $HOSTNAME at `date`

# Inputs needed:
# 1. HMM file
#	- Supplies the HMM model to PennCNV
#	- PennCNV provides a model, which is used here
# 2. PFB file
#	- Provides the population frequency of B allele information for each marker
#	- This file would need to be re-generated for every new array used
#	- This was generated from 1000 random UKB samples (see script 6_make_pfb.sh)
# 3. Input signal intensity files
#	- File shows Log R ratio and B allele frequency at each SNP for each sample
#	- Input files can be provided as a list (as here)
perl $PENNCNV_PATH/detect_cnv.pl -test --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Analysis_files/1Mv1.pfb --list Analysis_files/1Mv1_files.list -log logs/4_1mv1_chrX.log -out raw_out/1Mv1_chrX.txt --sexfile Analysis_files/1Mv1_sexfile.list
perl $PENNCNV_PATH/detect_cnv.pl -test --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Analysis_files/1Mv3.pfb --list Analysis_files/1Mv3_files.list -log logs/4_1mv3_chrX.log -out raw_out/1Mv3_chrX.txt --sexfile Analysis_files/1Mv3_sexfile.list
perl $PENNCNV_PATH/detect_cnv.pl -test --chrx -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Analysis_files/Omni2.5.pfb --list Analysis_files/Omni2.5_files.list -log logs/4_omni2.5_chrX.log -out raw_out/Omni2.5_chrX.txt --sexfile Analysis_files/Omni2.5_sexfile.list

echo Ended at `date`
