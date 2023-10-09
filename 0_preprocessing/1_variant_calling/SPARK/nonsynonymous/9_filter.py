import pandas as pd
import numpy as np

# Filter variants for mutation type and frequency

# Read in data
cols=['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38',
	'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38',
	'gnomAD_exome', 'gnomAD_genome', 'CADD_PHRED', 'CADD_Raw',
	'Sample', 'GT', 'DP', 'AD', 'GQ', 'PL']

df=pd.read_csv('tables/8_combined.tsv', sep='\t', header=None, names=cols)
print(df.shape)

# Reformat chromosome
df['Chr']='chr'+df['Chr'].astype(str)

# Assign a mutation type to each variant
missense_variant_classes = ['nonsynonymous_SNV', 'nonframeshift_deletion', 'nonframeshift_insertion', 'nonframeshift_substitution']
lof_variant_classes = ['stopgain', 'stoploss', 'frameshift_deletion', 'frameshift_insertion', 'frameshift_substitution']
df['Mut_type'] = '.'
df.loc[df['Func.wgEncodeGencodeBasicV38'].str.contains('splicing'), 'Mut_type']='splice'
df.loc[(df['Func.wgEncodeGencodeBasicV38'].str.contains('exonic')) & (df['ExonicFunc.wgEncodeGencodeBasicV38'].isin(missense_variant_classes)), 'Mut_type']='missense'
df.loc[(df['Func.wgEncodeGencodeBasicV38'].str.contains('exonic')) & (df['ExonicFunc.wgEncodeGencodeBasicV38'].isin(lof_variant_classes)), 'Mut_type']='lof'

# Filter missense and splice for only those with CADD>=25
# Replace '.' in CADD with NA
df.loc[df.CADD_PHRED=='.', 'CADD_PHRED']=np.nan
df.CADD_PHRED=df.CADD_PHRED.astype(float)
df['Mutation_filter']=False
df.loc[df.Mut_type=='lof', 'Mutation_filter']=True
df.loc[(df.Mut_type.isin(['missense', 'splice'])) & (df.CADD_PHRED>=25), 'Mutation_filter']=True
df=df[df.Mutation_filter]
print(df.shape)

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
df.to_csv('tables/9_filtered.csv', index=False)
