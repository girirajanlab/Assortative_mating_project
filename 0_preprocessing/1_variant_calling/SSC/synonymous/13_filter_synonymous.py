#!/bin/python3

import pandas as pd
import numpy as np

df = pd.read_csv('nonsynonymous/intermediate_files/final/12_variants_merged_all.txt', sep = '\t',
		names = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38', 'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38',
		'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore','Sample', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB'])
# GT="alt"
df=df[df.GT.str.contains('1')]

# QUAL>=50
df=df[df.Qual>=50]

# FORMAT/DP>=8
df=df[df.DP!='.']
df.DP=df.DP.astype(int)
df=df[df.DP>=8]

# FORMAT/AD[:1]>0
df['AD_count']=df.AD.str.count(',')
print(df[df.AD_count>2])

df['alt_depth']=df.AD.str.split(',', expand=True)[1].astype(int)
df=df[df.alt_depth > 0]

# FORMAT/AD[:1])/(FORMAT/DP)>=0.25
df['alt_ratio']=df.alt_depth/df.DP
df=df[df.alt_ratio>=0.25]

# (FORMAT/AD[:1])/(FORMAT/DP)<=0.75 | (FORMAT/AD[:1])/(FORMAT/DP)>=0.9
df=df[(df.alt_ratio<=0.75) | (df.alt_ratio>=0.9)]

# QUAL/(FORMAT/AD[:1])>=1.5
df['qual_ratio']=df.Qual/df.alt_depth
df=df[df.qual_ratio>=1.5]

df['Func.wgEncodeGencodeBasicV38'] = df['Func.wgEncodeGencodeBasicV38'].apply(lambda s: s.replace('\\x3b', ';'))
df['Gene.wgEncodeGencodeBasicV38'] = df['Gene.wgEncodeGencodeBasicV38'].apply(lambda s: s.replace('\\x3b', ';'))
df['ExonicFunc.wgEncodeGencodeBasicV38'] = df['ExonicFunc.wgEncodeGencodeBasicV38'].apply(lambda s: s.replace('\\x3b', ';'))

print(df['ExonicFunc.wgEncodeGencodeBasicV38'].value_counts())
df=df[(df['ExonicFunc.wgEncodeGencodeBasicV38'].str.contains('synonymous_SNV')) & (df['ExonicFunc.wgEncodeGencodeBasicV38']!='nonsynonymous_SNV')]
print(df['ExonicFunc.wgEncodeGencodeBasicV38'].value_counts())

df['Mut_type'] = 'synonymous'

# at this point there are no more Func, Gene, or ExonicFunc with ;
df.to_csv('tables/13_synonymous.csv', index=False)
