import pandas as pd

# Merge SNV burden and EIGENSTRAT PCs into one table
snv=pd.read_csv('../../0_preprocessing/1_variant_calling/UK_Biobank/synonymous/tables/17_burden_table.csv')
snv['variant_burden']=snv.burden_nonsynonymous
snv['synonymous_burden']=snv.burden_synonymous

eig=pd.read_csv('../../0_preprocessing/2_kinship_PCs/UK_Biobank/PCs/Results/3_eigenstrat.pca.evec', delim_whitespace=True)
eig.columns=['PC'+str(i) for i in range(1, 21)]+['Phenotype']
eig['Sample']=eig.index.to_list()
eig.Sample=eig.Sample.str.split(':', expand=True)[0]
eig.Sample=eig.Sample.astype(int)

tab=pd.merge(snv[['Sample', 'variant_burden', 'synonymous_burden']], eig[['Sample']+['PC'+str(i) for i in range(1, 21)]], on='Sample', how='inner')
print(tab)

# Save to file
tab.to_csv('UKB_snv_table.csv', index=False)
