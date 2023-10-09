import pandas as pd
import numpy as np

# Filter variants for mutation type and frequency

# Read in data
cols=['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38',
	'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38',
	'gnomAD_exome', 'gnomAD_genome', 'CADD_PHRED', 'CADD_Raw',
	'Sample', 'GT', 'DP', 'AD', 'GQ', 'PL']

df=pd.read_csv('nonsynonymous/tables/8_combined.tsv', sep='\t', header=None, names=cols)
print(df.shape)

# Reformat chromosome
df['Chr']='chr'+df['Chr'].astype(str)

# Assign a mutation type to each variant
df=df[(df['ExonicFunc.wgEncodeGencodeBasicV38'].str.contains('synonymous_SNV')) & (df['ExonicFunc.wgEncodeGencodeBasicV38']!='nonsynonymous_SNV')]
print(df['ExonicFunc.wgEncodeGencodeBasicV38'].value_counts())
df['Mut_type']='synonymous'

# Intracohort frequency filter
df['variant_id']=df.Chr.astype(str)+'_'+df.Pos.astype(str)+'_'+df.Ref.astype(str)+'_'+df.Alt.astype(str)
counts=df.variant_id.value_counts().to_dict()

sample_count=len(list(df.Sample.unique()))
print(sample_count)
df['frequency']=df.variant_id.map(counts)/sample_count

# Filter frequency
df=df[df.frequency<=0.001]
print(df.shape)

# Save to file
df.to_csv('tables/1_synonymous.csv', index=False)
