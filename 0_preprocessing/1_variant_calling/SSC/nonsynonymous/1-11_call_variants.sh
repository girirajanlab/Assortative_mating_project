#!/bin/bash

a=$1
in=Project_SSC_9209samples.JGvariants.2019-06-21/CCDG_9000JG_B01_GRM_WGS_2019-03-21_chr${a}.recalibrated_variants.flagged.vcf.gz

# intermediate files will be stored here
intermediate_dir=./intermediate_files/${a}
mkdir -p $intermediate_dir

# first bcftools split multi allelic records and left-normalize variants; then pipe to filter for quality
FASTA=hg38/Homo_sapiens_assembly38.fasta
bcftools norm -f ${FASTA} -m-both $in | bcftools view -i 'GT="alt" & QUAL>=50 & FORMAT/DP>=8 & FORMAT/AD[:1]>0 & (FORMAT/AD[:1])/(FORMAT/DP)>=0.25 & QUAL/(FORMAT/AD[:1])>=1.5' | bcftools view -i '(FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9' -o $intermediate_dir/${a}.1_quality_filter.vcf.gz -Oz
bcftools index $intermediate_dir/${a}.1_quality_filter.vcf.gz

# strip samples from vcf
bcftools view -G -o $intermediate_dir/${a}.2_samples_stripped.vcf.gz -Oz $intermediate_dir/${a}.1_quality_filter.vcf.gz

# annotate with annovar
perl annovar/table_annovar.pl $intermediate_dir/${a}.2_samples_stripped.vcf.gz annovar/humandb/ \
 -buildver hg38 \
 -out $intermediate_dir/${a}_3_annovar \
 -remove \
 -protocol wgEncodeGencodeBasicV38 \
 -operation g \
 -nastring . \
 -vcfinput \
 -arg '-hgvs'

# filter for exonic regions
bcftools view -i 'INFO/Func.wgEncodeGencodeBasicV38~"exonic" | INFO/Func.wgEncodeGencodeBasicV38~"splicing"' -o $intermediate_dir/${a}.4_filter_exon.vcf.gz -Oz $intermediate_dir/${a}_3_annovar.hg38_multianno.vcf

# annotate gnomad
vcfanno_linux64.1 gnomad.toml $intermediate_dir/${a}.4_filter_exon.vcf.gz | bgzip > $intermediate_dir/${a}.5_gnomad.vcf.gz

# filter for rarity
bcftools view -i 'gnomad_exome_AF<=0.001 | gnomad_exome_AF="."' $intermediate_dir/${a}.5_gnomad.vcf.gz | bcftools view -i 'gnomad_genome_AF<=0.001 | gnomad_genome_AF="."' -o $intermediate_dir/${a}.6_population_frequency.vcf.gz -Oz

# annotate cadd
/data5/software/vcfanno_linux64.1 cadd.toml $intermediate_dir/${a}.6_population_frequency.vcf.gz | bgzip > $intermediate_dir/${a}.7_cadd.vcf.gz

# Extract and write records from A shared by both A and B using exact allele match
bcftools index $intermediate_dir/${a}.7_cadd.vcf.gz
bcftools isec $intermediate_dir/${a}.1_quality_filter.vcf.gz $intermediate_dir/${a}.7_cadd.vcf.gz -p $intermediate_dir/${a}.8_intersect -n =2 -w 1
bgzip $intermediate_dir/${a}.8_intersect/0000.vcf
bcftools index $intermediate_dir/${a}.8_intersect/0000.vcf.gz

# merge the information from both vcfs
bcftools merge $intermediate_dir/${a}.8_intersect/0000.vcf.gz $intermediate_dir/${a}.7_cadd.vcf.gz > $intermediate_dir/${a}.9_merged.vcf

# vcf to table
bgzip $intermediate_dir/${a}.9_merged.vcf
mkdir ./intermediate_files/${a}/10.by_sample

java -Xmx50G -jar jvarkit/dist/biostar130456.jar -x -z -p "$intermediate_dir/10.by_sample/__SAMPLE__.vcf" $intermediate_dir/${a}.9_merged.vcf.gz

while read sample; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Func.wgEncodeGencodeBasicV38\t%Gene.wgEncodeGencodeBasicV38\t%GeneDetail.wgEncodeGencodeBasicV38\t%ExonicFunc.wgEncodeGencodeBasicV38\t%AAChange.wgEncodeGencodeBasicV38\t%gnomad_exome_AF\t%gnomad_genome_AF\t%CADD_PHRED\t%CADD_RawScore[\t%SAMPLE\t%GT\t%DP\t%AD\t%GQ\t%PL\t%RGQ\t%SB]\n' $intermediate_dir/10.by_sample/${sample}.vcf > $intermediate_dir/10.by_sample/${sample}.${a}.variants.txt; done < ssc_sample_list.txt
cat $intermediate_dir/10.by_sample/*.variants.txt > $intermediate_dir/${a}.11_variants_merged.txt
