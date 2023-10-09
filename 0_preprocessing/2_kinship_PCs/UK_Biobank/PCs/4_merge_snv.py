import pandas as pd

# Merge the PC results with the SNV burden
pc=pd.read_csv('Results/3_eigenstrat.pca.evec', delim_whitespace=True, index_col=0)
pc.columns=['PC'+str(i) for i in range(1, 21)]+['Pheno']
pc['Sample']=pc.index.to_list()
pc.Sample=pc.Sample.str.split(':', expand=True)[1]
print(pc)

snv=pd.read_csv('../../../1_variant_calling/UK_Biobank/nonsynonymous/tables/18_burden_table.csv', header = None, names = ['IID', 'SNV_burden'])
