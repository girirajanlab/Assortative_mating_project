#!/bin/python3
import pandas as pd

filename = 'intermediate_files/final/16_rare_deleterious_exonic_variants_cohort_count_gene_ids_loeuf.csv'

df = pd.read_csv(filename)

def get_depth_of_alternative_allele(ad):
	sad = ad.split(',')
	if len(sad) != 2:
		print('problem ', ad)
	return int(sad[1])


df['alternative_depth'] = df.AD.apply(get_depth_of_alternative_allele)

df = df[df.DP >= 8]
df = df[df.alternative_depth > 0]
df = df[(df.alternative_depth / df.DP) >= 0.25]
df = df[(df.Qual / df.alternative_depth) >= 1.5]
df = df[((df.alternative_depth / df.DP) <= 0.75) | ((df.alternative_depth / df.DP) >= 0.9)]


# drop the column that I created
df = df.drop('alternative_depth', axis=1)


filename = 'intermediate_files/final/16_rare_deleterious_exonic_variants_cohort_count_gene_ids_loeuf_qc_filtered.csv'
df.to_csv(filename, index=False)
