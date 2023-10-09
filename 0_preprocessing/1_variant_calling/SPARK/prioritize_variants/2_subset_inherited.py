import pandas as pd
import numpy as np

# Get the inherited variants in probands
df=pd.read_csv('SPARK_Variants_wInheritance.txt', sep='\t')
print(df.inheritance.value_counts())
df=df[df.inheritance.isin(['M', 'F'])][['Sample', 'Chr', 'Pos', 'Ref', 'Alt', 'Gene_symbol', 'Gene_id', 'AAChange.wgEncodeGencodeBasicV38', 'Mut_type', 'inheritance', 'batch']]
print(df)

# Restrict to samples with variants from both parents
m_samps=set(df[df.inheritance=='M']['Sample'].to_list())
f_samps=set(df[df.inheritance=='F']['Sample'].to_list())
keep_samps=list(m_samps.intersection(f_samps))
df=df[df.Sample.isin(keep_samps)]
print(df)

# Split up variants that affect multiple genes
df['vid']=df.Chr+'_'+df.Pos.astype(str)+'_'+df.Ref+'_'+df.Alt
df['gene_id']=df.Gene_id.str.split('.', expand=True)[0]

multi_gene=df[df.Gene_id.str.contains(';')].copy()
multi_gene['gene1']=df.Gene_id.str.split(';', expand=True)[0]
multi_gene['gene2']=df.Gene_id.str.split(';', expand=True)[1]

g1=multi_gene.copy()
g1['gene_id']=g1.gene1.str.split('.', expand=True)[0]
g2=multi_gene.copy()
g2['gene_id']=g2.gene2.str.split('.', expand=True)[0]

df=df[~df.Gene_id.str.contains(';')].copy()
df=pd.concat([df, g1, g2])

df=df[['Sample', 'Chr', 'Pos', 'Ref', 'Alt', 'vid', 'Gene_symbol', 'Gene_id', 'gene_id', 'AAChange.wgEncodeGencodeBasicV38', 'Mut_type', 'inheritance', 'batch']]

print(df)

# Save variants
df.to_csv('tables/2_inherited_variants_slim.csv', index=False)
