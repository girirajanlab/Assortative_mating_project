#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=2.1_gc_mod_adj
#SBATCH -o logs/2.1_gc_mod_adj.log
#SBATCH -e logs/2.1_gc_mod_adj.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G

# Compile a Population frequency of B allele file
export PENNCNV_PATH=/data5/software/PennCNV-1.0.5/

gunzip -c /data5/software/PennCNV/PennCNV-1.0.5/gc_file/hg19.gc5Base.txt.gz > Analysis_files/hg19.gc5Base.txt

GCfile='Analysis_files/hg19.gc5Base.txt'

$PENNCNV_PATH/cal_gc_snp.pl $GCfile Analysis_files/iWES_snp_pos.txt -output Analysis_files/iWES_gcmod.txt

