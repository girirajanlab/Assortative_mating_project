#!/bin/python3

# gnomad.v2.1.1.lof_metrics.by_gene.tsv from Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/gnomad_annotations

import pandas as pd

df = pd.read_csv('intermediate_files/final/15_rare_deleterious_exonic_variants_cohort_count_gene_ids.csv')

gnomad_df = pd.read_csv('gnomad.v2.1.1.lof_metrics.by_gene.tsv', sep='\t')
gnomad_df = gnomad_df.set_index('gene_id')

def remove_dot(gene_ids):
	gene_ids = gene_ids.split(';')
	new_gene_ids = []
	for gene_id in gene_ids:
		new_gene_ids.append(gene_id.split('.')[0])
	return ';'.join(new_gene_ids)


# remove dot (ENSG00000215910.3 > ENSG00000215910)
df['Gene_id_'] = df['Gene_id'].apply(remove_dot)

df[df['Gene_id'].apply(lambda s: ';' in s)][['Gene_id', 'Gene_id_']]

# get the lower loeuf score
df['LOEUF'] = ''
for i, row in df.iterrows():
	gene_ids = row['Gene_id_'].split(';')
	gene_ids_in_gnomad = [s for s in gene_ids if s in gnomad_df.index]

	if len(gene_ids_in_gnomad) == 0:
		# none of the gene ids are in gnomad
		print(gene_ids)
		continue

	lowest_loeuf = 1000
	for gene_id in gene_ids_in_gnomad:
		loeuf = gnomad_df.at[gene_id, 'oe_lof_upper']
		if loeuf < lowest_loeuf:
			lowest_loeuf = loeuf


	df.at[i, 'LOEUF'] = lowest_loeuf

df.to_csv('intermediate_files/final/16_rare_deleterious_exonic_variants_cohort_count_gene_ids_loeuf.csv', index=False)

