import pandas as pd

# Merge SNV (synonymous and other) and EIGENSTRAT data into one table
snv=pd.read_csv('../../1_variant_calling/16p12.1_deletion/synonymous/tables/5_burden_table.csv', index_col=0)
snv['Sample']=snv.index.to_list()
print(snv)

eig=pd.read_csv('Results/8_eigenstrat.pca.evec', delim_whitespace=True)
eig.columns=['PC'+str(i) for i in range(1, 21)]+['Pheno']
eig['Sample']=eig.index.to_list()
eig.Sample=eig.Sample.str.split(':', expand=True)[1]
print(eig)

df=pd.merge(snv, eig[[i for i in eig.columns.to_list() if 'PC' in i or i=='Sample']], on='Sample', how='inner')

# Save to file
df.to_csv('analysis_files/9_snv_eig.csv', index=False)
