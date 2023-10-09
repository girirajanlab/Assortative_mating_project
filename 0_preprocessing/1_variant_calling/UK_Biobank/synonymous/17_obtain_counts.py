import pandas as pd

# Get variant burden counts for each sample
syn=pd.read_csv('tables/16_gene_ids.csv')

syn_burden=pd.DataFrame(syn.Sample.value_counts())
syn_burden.columns=['burden']

# Nonsynonymous
nonsyn=pd.read_csv('../nonsynonymous/tables/17_gene_ids.csv')

nonsyn_burden=pd.DataFrame(nonsyn.Sample.value_counts())
nonsyn_burden.columns=['burden']

df=pd.merge(nonsyn_burden, syn_burden, left_index=True, right_index=True, how='outer', suffixes=['_nonsynonymous', '_synonymous'])
df['Sample']=df.index.to_list()
df.fillna(0, inplace=True)

df=df.astype(int)
df.sort_values(by=['burden_nonsynonymous', 'burden_synonymous', 'Sample'], inplace=True)
print(df)

# Save to file
df.to_csv('tables/17_burden_table.csv', index=False)
