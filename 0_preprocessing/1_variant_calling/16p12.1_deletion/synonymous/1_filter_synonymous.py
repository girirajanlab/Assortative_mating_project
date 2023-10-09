import pandas as pd

# Get only the synonymous variants for the cohort
df=pd.read_csv('../nonsynonymous/tables/9_exonic_variants_gnomad_cadd_clinvar.txt', sep='\t', header=None,
		names = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
                'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19',
                'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'ClinVar_CLNDN', 'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG', 'ClinVar_ALLELEID',
                'Sample', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB'])

# Get exonic/splicing
print(df['Func.wgEncodeGencodeBasicV19'].value_counts())
df=df[((df['Func.wgEncodeGencodeBasicV19'].str.contains('exonic')) | (df['Func.wgEncodeGencodeBasicV19'].str.contains('splicing'))) & ~(df['Func.wgEncodeGencodeBasicV19'].isin(['ncRNA_exonic', 'ncRNA_splicing']))]
print(df['Func.wgEncodeGencodeBasicV19'].value_counts())

# Get synonymous variants
print(df['ExonicFunc.wgEncodeGencodeBasicV19'].value_counts())
df=df[(df['ExonicFunc.wgEncodeGencodeBasicV19'].str.contains('synonymous_SNV')) & (df['ExonicFunc.wgEncodeGencodeBasicV19']!='nonsynonymous_SNV')]
print(df['ExonicFunc.wgEncodeGencodeBasicV19'].value_counts())

# Save to file
df.to_csv('tables/1_synonymous_variants.csv', index=False)


