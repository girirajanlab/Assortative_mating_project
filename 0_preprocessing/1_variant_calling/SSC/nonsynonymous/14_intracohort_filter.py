#!/bin/python3

import pandas as pd

df = pd.read_csv('intermediate_files/final/13_rare_deleterious_exonic_variants.csv')

# Remove samples that do not have consent forms or we are excluding for WGS
# Skip this step for SSC--keep all 9209 samples

#samples_df = pd.read_csv('16p12_All_Participants_v5.csv')
#samples_df = samples_df[samples_df['No_consent_forms'] != 'X']
#samples_df = samples_df[samples_df['WGS'] == 'X']

#df = df[df.Sample.isin(samples_df.Sample)]

# Reorder columns
columns = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Mut_type', 'variant_filter', 'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19',
		'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB']
df = df[columns]

# Filter for intracohort frequency
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
variant_counts = df['variant_id'].value_counts()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

df.to_csv('intermediate_files/final/14_rare_deleterious_exonic_variants_cohort_count.csv', index=False)
