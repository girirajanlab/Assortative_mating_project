#!/bin/bash
invcf=vcfs/7_remove_dot.vcf
outvcf=vcfs/8_all_annovar_exonic_splicing.vcf

bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' $invcf > $outvcf
