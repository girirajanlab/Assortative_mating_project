#!/bin/python3

import pandas as pd

df = pd.read_csv('tables/1_synonymous_variants.csv')

samples_df = pd.read_csv('16p12_All_Participants_v5.csv')
samples_df = samples_df[samples_df['WGS'] == 'X']

df = df[df.Sample.isin(samples_df.Sample)]

# Reorder columns
columns = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
		'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19',
		'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'ClinVar_CLNDN', 'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG', 'ClinVar_ALLELEID',
		'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB']
df = df[columns]

# Filter for intracohort frequency
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
variant_counts = df['variant_id'].value_counts()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

df.to_csv('tables/2_intracohort.csv', index=False)
