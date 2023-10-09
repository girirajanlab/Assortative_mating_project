import pandas as pd

# Combine the nonsynonymous and synonymous burden counts into a single file
syn=pd.read_csv('tables/15_gencode.csv')
sampleid_map=pd.read_csv('documentation/nygc_sfari_id_map.csv')
sampleid_map.index=sampleid_map['Repository Id']

syn['Sample']=syn.Sample.map(sampleid_map['SFARI ID'].to_dict())
syn_burden=pd.DataFrame(syn.Sample.value_counts())
syn_burden.columns=['burden']

nonsyn=pd.read_csv('../nonsynonymous/intermediate_files/final/16_rare_deleterious_exonic_variants_cohort_count_gene_ids_loeuf_qc_filtered_sample_names.csv')
nonsyn_burden=pd.DataFrame(nonsyn.Sample.value_counts())
nonsyn_burden.columns=['burden']

df=pd.merge(syn_burden, nonsyn_burden, left_index=True, right_index=True, how='outer', suffixes=['_synonymous', '_nonsynonymous'])
df['Sample']=df.index.to_list()
df.fillna(0, inplace=True)
df.sort_values(by=['burden_nonsynonymous', 'burden_synonymous'], inplace=True)
print(df)


# Save to file
df.to_csv('tables/16_burden_table.csv', index=False)
