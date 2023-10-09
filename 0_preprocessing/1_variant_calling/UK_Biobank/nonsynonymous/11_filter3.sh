#!/bin/bash

# Filter for:
# Exonic or splicing mutations
# gnomAD genome and exome frequency <=0.001
bcftools view -i "ANNOVAR_DATE!='.'" vcfs/10_gnomad.vcf.gz  | bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' | bcftools view -i 'gnomad_exome_AF<=0.001 | gnomad_exome_AF="."' | bcftools view -i 'gnomad_genome_AF<=0.001 | gnomad_genome_AF="."' | bgzip > vcfs/11_filter3.vcf.gz
tabix -p vcf vcfs/11_filter3.vcf.gz
