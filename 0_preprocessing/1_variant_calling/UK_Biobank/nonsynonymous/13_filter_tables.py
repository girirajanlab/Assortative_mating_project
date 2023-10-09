import pandas as pd
import numpy as np
import sys

# Filter tables to include only
# LOF, missense, and splice mutations
# missense and splice must have a CADD score >= 25

sample=sys.argv[1]
prefix=sys.argv[2]

cols=['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38',
	'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38', 'CADD_PHRED',
	'Sample', 'GT', 'DP', 'AD']
# Read data
df=pd.read_csv(prefix+'/'+sample+'.tsv', sep='\t', header=None, names=cols)

# Sample field did not output correctly - adjust here
df.Sample=sample

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

# Save to file
df.to_csv(prefix+'/'+sample+'_filtered.csv', index=False, header=False)
