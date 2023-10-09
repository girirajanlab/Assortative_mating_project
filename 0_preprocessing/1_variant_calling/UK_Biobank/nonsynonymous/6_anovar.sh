#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=6_annovar
#SBATCH -o logs/6_annovar.log
#SBATCH -e logs/6_annovar.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=40G

echo `date` starting job on $HOSTNAME

invcf=vcfs/all.vcf.gz
outprefix=vcfs/all_annovar

perl annovar/table_annovar.pl $invcf \
	annovar/humandb/ \
	-buildver hg38 -out $outprefix -remove \
	-protocol wgEncodeGencodeBasicV38 -operation g \
	-nastring . \
	-vcfinput \
	-arg '-hgvs'

echo `date` finished
