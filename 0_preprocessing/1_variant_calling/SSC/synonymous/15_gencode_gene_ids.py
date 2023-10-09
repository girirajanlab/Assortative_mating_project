#!/bin/python3

import pandas as pd

df = pd.read_csv('tables/14_intracohort_filter.csv')

print(df[df['Gene.wgEncodeGencodeBasicV38'].apply(lambda s: ';' in s)][['Func.wgEncodeGencodeBasicV38', 'Gene.wgEncodeGencodeBasicV38']])

gencode_df = pd.read_csv('../../gencode.v38.parsed.genes.csv')

df['Chr_Pos']=df.Chrom+'_'+df.Pos.astype(str)
uniq_locations=list(df.Chr_Pos.unique())

# Annotate unique locations with Gencode information
gencode_df = pd.read_csv('/data5/SSC_WGS/SSC_annotations_MJ/gencode.v38.parsed.genes.csv')

def get_gene_attribute(chr_pos, attributes=['gene_name', 'gene_id', 'gene_type'], offset=0):
	chrom=chr_pos.split('_')[0]
	pos=int(chr_pos.split('_')[1])
	# Get overlapping lines
	sub_gencode=gencode_df[(gencode_df.Chrom==chrom) & (gencode_df.Start<=pos+offset) & (gencode_df.End>=pos-offset)]

	# increase search space if none found
	if sub_gencode.shape[0]==0:
		print(chr_pos)
		off=offset+10
		if off>100:
			print('Exceeded recursion limit')
			return ['.' for i in attributes]
		return get_gene_attribute(chr_pos, attributes, offset=offset+10)

	out=[]
	for att in attributes:
		out.append(';'.join(sub_gencode[att].to_list()))
	return out

symbol_dict={}
id_dict={}
biotype_dict={}
for loc in uniq_locations:
	out=get_gene_attribute(loc)
	symbol_dict[loc]=out[0]
	id_dict[loc]=out[1]
	biotype_dict[loc]=out[2]

df['Gene_symbol'] = df.Chr_Pos.map(symbol_dict)
df['Gene_id'] = df.Chr_Pos.map(id_dict)
df['Gene_biotype'] = df.Chr_Pos.map(biotype_dict)

df[['Gene.wgEncodeGencodeBasicV38','Gene_symbol']]

# remove gene symbols that annovar did not annotate
for i, row in df.iterrows():
	annovar_genes = row['Gene.wgEncodeGencodeBasicV38'].split(';')
	gencode_genes = row['Gene_symbol'].split(';')
	gencode_ids = row['Gene_id'].split(';')
	gencode_biotypes = row['Gene_biotype'].split(';')

	# create dict of gene ids and biotypes
	symbol2id = {}
	symbol2biotype = {}
	for j in range(len(gencode_genes)):
		gencode_gene = gencode_genes[j]
		gene_id = gencode_ids[j]
		gene_biotype = gencode_biotypes[j]

		symbol2id[gencode_gene] = gene_id
		symbol2biotype[gencode_gene] = gene_biotype

	# remove gencode genes that aren't annotated by annovar
	genes_to_remove = []
	for gene in gencode_genes:
		if gene not in annovar_genes:
			genes_to_remove.append(gene)

	for gene in genes_to_remove:
		gencode_genes.remove(gene)

	# sort gencode and annovar genes
	gencode_genes = sorted(list(set(gencode_genes)))
	annovar_genes = sorted(list(set(annovar_genes)))

	# make sure gene ids and biotypes match order of genes
	gencode_ids = [symbol2id[s] for s in gencode_genes]
	biotype_ids = [symbol2biotype[s] for s in gencode_genes]

	df.at[i, 'Gene_symbol'] = ';'.join(gencode_genes)
	df.at[i, 'Gene.wgEncodeGencodeBasicV38'] = ';'.join(annovar_genes)
	df.at[i, 'Gene_id'] = ';'.join(gencode_ids)
	df.at[i, 'Gene_biotype'] = ';'.join(biotype_ids)

# assert that all genes found by annovar match genes found by this script
assert((df['Gene.wgEncodeGencodeBasicV38'] != df['Gene_symbol']).sum() == 0)

# Keep only variants in protein-coding genes
df = df[df['Gene_biotype'].apply(lambda s: 'protein_coding' in s)]

df.to_csv('tables/15_gencode.csv', index=False)
