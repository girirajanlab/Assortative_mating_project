import pandas as pd

# Annotate inherited variants with disease gene sets
df=pd.read_csv('tables/2_inherited_variants_slim.csv')
df.sort_values(by=['Sample', 'inheritance'], inplace=True)

# Get the number of genes in each disease gene set for each sample
gene_annos=pd.read_csv('Gene_annotations.csv')
genesets=[i for i in gene_annos.columns.to_list()[2:15] if i not in ['BrainExpressed_Kang2011', 'Constrained_LOEUF']]
print(genesets)

# Annotate variants with geneset
for gs in genesets:
	df[gs]=0
	df.loc[df.gene_id.isin(gene_annos[gene_annos[gs]!='.']['gene_id'].to_list()), gs]=1
	print(df[gs].value_counts())

# Check the number of sets each gene is annotated with
df['set_num']=df[genesets].sum(axis=1)
df=df[df.set_num>0]

# Restrict to probands who inherited at least one variant from each parent
def keep_pros(df):
	fpros=df[df.inheritance=='F']['Sample'].to_list()
	mpros=df[df.inheritance=='M']['Sample'].to_list()
	keep_pros=list(set(fpros).intersection(set(mpros)))
	return keep_pros

df=df[df.Sample.isin(keep_pros(df))]

for gs in genesets:
	print(df[gs].value_counts())

# Save to file
df.to_csv('tables/3_inherited_geneset_annotated.csv', index=False)
