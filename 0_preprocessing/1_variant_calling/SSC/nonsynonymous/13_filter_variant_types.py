#!/bin/python3

import pandas as pd
import numpy as np

df = pd.read_csv('intermediate_files/final/12_variants_merged_all.txt', sep = '\t',
		names = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19',
		'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore','Sample', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB'])
print(df)

def dot2na(s):
	if s == '.':
		return np.nan
	return float(s)

df['CADD_PHRED'] = df['CADD_PHRED'].apply(dot2na)

df['Func.wgEncodeGencodeBasicV19'] = df['Func.wgEncodeGencodeBasicV19'].apply(lambda s: s.replace('\\x3b', ';'))
df['Gene.wgEncodeGencodeBasicV19'] = df['Gene.wgEncodeGencodeBasicV19'].apply(lambda s: s.replace('\\x3b', ';'))
df['ExonicFunc.wgEncodeGencodeBasicV19'] = df['ExonicFunc.wgEncodeGencodeBasicV19'].apply(lambda s: s.replace('\\x3b', ';'))

keep_locations = ['exonic', 'splicing']
def should_keep_locations(locations):
	locations = locations.split(';')
	for location in locations:
		if location in keep_locations:
			return True
	return False

df = df[df['Func.wgEncodeGencodeBasicV19'].apply(lambda s: should_keep_locations(s))]

# here I check to see if there are any cases where a variant has
# more than one ExonicFunc. There aren't any so I don't have to worry about it
print(df['ExonicFunc.wgEncodeGencodeBasicV19'].value_counts())

missense_variant_classes = ['nonsynonymous_SNV', 'nonframeshift_deletion', 'nonframeshift_insertion', 'nonframeshift_substitution']
lof_variant_functions = ['stopgain', 'stoploss', 'frameshift_deletion', 'frameshift_insertion', 'frameshift_substitution']
def get_mutation_type(functions, exonic_function):
	functions = functions.split(';')
	if 'exonic' in functions:
		if exonic_function in missense_variant_classes:
			return 'missense'
		if exonic_function in lof_variant_functions:
			return 'lof'
	if 'splicing' in functions:
		return 'splice'
	return exonic_function

df['Mut_type'] = df[['Func.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19']].apply(lambda s: get_mutation_type(s[0], s[1]), axis=1)

df = df[df['Mut_type'].isin(['splice', 'missense', 'lof'])]

def filter_variants(mut_type, cadd):
	if mut_type=='lof':
		return True
	if cadd>=25:
		return True
	return False

df['variant_filter'] = df[['Mut_type', 'CADD_PHRED']].apply(lambda s: filter_variants(s[0], s[1]), axis=1)

df = df[df['variant_filter']]

# in some cases a variant has Func="exonic;splicing" region or
# Func="ncRNA_exonic;splicing"
# (see below)
pd.set_option('display.max_rows', 400)
print(df[df['Func.wgEncodeGencodeBasicV19'].apply(lambda s: ';' in s)][['Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19']])
print(df[df['Gene.wgEncodeGencodeBasicV19'].apply(lambda s: ';' in s)][['Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19']])

# in which case, only keep the gene whose function we're more interested in
# exonic > splicing > ncRNA_exonic
for i, row in df.iterrows():
	locations = row['Func.wgEncodeGencodeBasicV19'].split(';')
	genes = row['Gene.wgEncodeGencodeBasicV19'].split(';')

	if len(locations) == 1:
		continue

	if 'exonic' in locations:
		index = locations.index('exonic')
		df.at[i, 'Func.wgEncodeGencodeBasicV19'] = locations[index]
		df.at[i, 'Gene.wgEncodeGencodeBasicV19'] = genes[index]

	elif 'splicing' in locations:
		index = locations.index('splicing')
		df.at[i, 'Func.wgEncodeGencodeBasicV19'] = locations[index]
		df.at[i, 'Gene.wgEncodeGencodeBasicV19'] = genes[index]

# at this point there are no more Func, Gene, or ExonicFunc with ;
df.to_csv('intermediate_files/final/13_rare_deleterious_exonic_variants.csv', index=False)
