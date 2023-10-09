import pandas as pd

# Filter variants for those found in <= 0.1% of the cohort

# Read in data
cols=['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38',
	'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38', 'CADD_PHRED',
	'BLANK', 'GT', 'DP', 'AD', 'Sample']
df=pd.read_csv('tables/13_filtered_combined.tsv', header=None, names=cols, sep='\t')
print(df)

df['variant_id']=df.Chr.astype(str)+'_'+df.Pos.astype(str)+'_'+df.Ref.astype(str)+'_'+df.Alt.astype(str)
counts=df.variant_id.value_counts().to_dict()

# Number of samples is total number of samples we called variants on - in samples.csv
samples=pd.read_csv('samples.csv')
sample_count=samples.shape[0]

df['frequency']=df.variant_id.map(counts)/sample_count

# Filter frequency
df=df[df.frequency<=0.001]

# Save to file
df[['Chr', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38',
        'GeneDetail.wgEncodeGencodeBasicV38', 'ExonicFunc.wgEncodeGencodeBasicV38', 'AAChange.wgEncodeGencodeBasicV38',
        'Sample', 'GT', 'DP', 'AD', 'frequency']].to_csv('tables/14_intracohort_filter.csv', index=False)
