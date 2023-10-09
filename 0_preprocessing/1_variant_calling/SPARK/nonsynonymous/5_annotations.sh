#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=annotations
#SBATCH -o logs/5_annotations.log
#SBATCH -e logs/5_annotations.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G

echo `date` starting job on $HOSTNAME

in=vcfs/4c_merge2.vcf.gz

# annotate with annovar
perl annovar/table_annovar.pl $in annovar/humandb/ \
-buildver hg38 \
-out vcfs/5_annotations/all \
-remove \
-protocol wgEncodeGencodeBasicV38 \
-operation g \
-nastring . \
-vcfinput \
-arg '-hgvs'

# filter for exonic regions
bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' -o vcfs/5_annotations/all.hg38_multianno.exonic.vcf.gz -Oz vcfs/5_annotations/all.hg38_multianno.vcf

# annotate gnomad
vcfanno_linux64.1 gnomad.toml vcfs/5_annotations/all.hg38_multianno.exonic.vcf.gz | bgzip > vcfs/5_annotations/all.hg38_multianno.exonic.gnomad.vcf.gz

# annotate cadd
vcfanno_linux64.1 cadd.toml vcfs/5_annotations/all.hg38_multianno.exonic.gnomad.vcf.gz | bgzip > vcfs/5_annotations/all.hg38_multianno.exonic.gnomad.cadd.vcf.gz

echo `date` finished
