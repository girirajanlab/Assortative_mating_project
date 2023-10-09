import pandas as pd

# Filter spouses for only those of European ancestry
# Do this using Data Field 22006 from the ukb30075.csv file
# This field is "Genetic ethnic group" and
# "Indicates samples who self-identified as 'White British' according to Field 21000 and have very similar genetic ancestry based on a principal components analysis of the genotypes."

samps=pd.read_csv('../Kinship/plink_files/4_filter/UKB_QC_het.valid.sample', sep='\t')

pheno_df = pd.read_csv('download/ukb30075.csv', encoding='unicode_escape', usecols = ['eid', '22006-0.0'])
pheno_df.columns=['IID', 'European']

df=pd.merge(samps, pheno_df, on='IID', how='left')
df=df[df.European==1]
print(df.shape)

# Also filter spouses for those with SNV data available
snv=pd.read_csv('../../../1_variant_calling/UK_Biobank/nonsynonymous/tables/18_burden_table.csv', header = None, names = ['IID', 'SNV_burden'])
df=df[df.IID.isin(snv.IID.to_list())]
print(df.shape)

# Filter for complete pairs
pairs=pd.read_csv('../../../../4_UK\ Biobank/2_spouse_table/UKB_family_table.csv')
couple_dict=dict(zip(pairs.Male.to_list()+pairs.Female.to_list(), pairs.Couple.to_list()+pairs.Couple.to_list()))
df['Couple']=df.IID.map(couple_dict)
comp_pairs=df.Couple.value_counts()
comp_pairs=comp_pairs[comp_pairs==2]
df=df[df.Couple.isin(comp_pairs.index.to_list())]
print(df.shape)

# Save needed samples as FAM
df[['FID', 'IID']].to_csv('list_files/1_european_samples.fam', sep=' ', index=False)
