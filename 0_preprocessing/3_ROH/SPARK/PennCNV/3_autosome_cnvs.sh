#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=3_run_penncnv
#SBATCH -o logs/3_autosome_cnvs.log
#SBATCH -e logs/3_autosome_cnvs.log
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
# (OPTIONAL) 4. GC model file
#	- File has the GC content in a 1Mb region around each SNP
#	- Only needed if "a significant fraction" of samples have a WF > 0.04 or < -0.04
perl $PENNCNV_PATH/detect_cnv.pl -test -hmm $PENNCNV_PATH/lib/hhall.hmm -pfb Analysis_files/iWES.pfb --list Analysis_files/iwes_command.list -log logs/3_iwes_autosome.log -out raw_out/iwes_autosome_adj.txt -gcmodel Analysis_files/iWES_gcmod.txt
