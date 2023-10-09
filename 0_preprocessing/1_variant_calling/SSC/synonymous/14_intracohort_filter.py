#!/bin/python3

import pandas as pd

df = pd.read_csv('tables/13_synonymous.csv')

# Reorder columns
columns = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Mut_type', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38', 'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38',
		'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB']
df = df[columns]

# Filter for intracohort frequency
df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
variant_counts = df['variant_id'].value_counts()
df['cohort_count'] = df['variant_id'].map(variant_counts)

df = df[df['cohort_count'] <= 10]

df.to_csv('tables/14_intracohort_filter.csv', index=False)
