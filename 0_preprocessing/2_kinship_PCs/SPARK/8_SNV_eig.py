import pandas as pd

# Merge SNV burden and EIGENSTRAT PCs into one table
snv=pd.read_csv('../../1_variant_calling/SPARK/synonymous/analysis_files/4_burden_table.csv')
snv['variant_burden']=snv.burden_nonsynonymous
snv['synonymous_burden']=snv.burden_synonymous

eig=pd.read_csv('Results/7_eigenstrat.pca.evec', delim_whitespace=True)
eig.columns=['PC'+str(i) for i in range(1, 21)]+['Phenotype']
eig['Sample']=eig.index.to_list()
eig.Sample=eig.Sample.str.split(':', expand=True)[1]

tab=pd.merge(snv[['Sample', 'variant_burden', 'synonymous_burden']], eig[['Sample']+['PC'+str(i) for i in range(1, 21)]], on='Sample', how='inner')
print(tab)

# Save to file
tab.to_csv('SPARK_snv_table.csv', index=False)
