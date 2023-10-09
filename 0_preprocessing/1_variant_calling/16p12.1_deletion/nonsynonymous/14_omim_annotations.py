#!/bin/python3

import pandas as pd

df = pd.read_csv('tables/13_rare_deleterious_exonic_variants_cohort_count_gene_ids_loeuf.csv')

omim_df = pd.read_csv('genemap2.txt', sep='\t', comment='#', header=None)
omim_columns = ['Chromosome', 'Genomic Position Start', 'Genomic Position End', 'Cyto Location', 'Computed Cyto Location', 'MIM Number',
		'Gene Symbols', 'Gene Name', 'Approved Gene Symbol', 'Entrez Gene ID', 'Ensembl Gene ID', 'Comments', 'Phenotypes', 'Mouse Gene Symbol/ID']
omim_df.columns = omim_columns

ensembl_genes_in_omim = omim_df['Ensembl Gene ID'].to_list()

df['MIM_number'] = ''
for i, row in df.iterrows():
	gene_ids = row['Gene_id_'].split(';')
	gene_ids_in_omim = [s for s in gene_ids if s in ensembl_genes_in_omim]

	if len(gene_ids_in_omim) == 0:
		# none of the gene ids are in omim
		continue

	subomim_df = omim_df[omim_df['Ensembl Gene ID'].isin(gene_ids)]

	# a gene id may have more than one omim number
	# in which case annotate all seperated by semicolons
	# (600179;607206)
	omim_numbers = subomim_df['MIM Number'].to_list()
	omim_numbers = [str(s) for s in omim_numbers]
	df.at[i, 'MIM_number'] = ';'.join(omim_numbers)

df.to_csv('tables/14_rare_deleterious_exonic_variants_cohort_count_gene_ids_loeuf_omim.csv', index=False)










